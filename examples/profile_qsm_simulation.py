"""
Profile a single QSM aero-structural simulation with detailed timing breakdown.

This script:
1. Runs a complete aero-structural simulation using the QSM solver
2. Instruments key functions to track time spent in each component
3. Generates a detailed performance report

Usage:
    python examples/profile_qsm_simulation.py

Optional arguments:
    --steering STR      Override steering_tape_final_extension [m]
    --profile-type PT   'cprofile', 'line_profiler', or 'instrumented' (default)
                        - cprofile: expensive but comprehensive
                        - instrumented: lightweight, focused on main loops
    --output DIR        Save timing report to DIR (default: ./profile_results)
"""

import numpy as np
from pathlib import Path
import argparse
import copy
import cProfile
import pstats
import io
import sys
import json
from datetime import datetime
from collections import defaultdict
from time import perf_counter

from kitesim.logging_config import *
from kitesim.utils import (
    load_and_save_config_files,
    load_yaml,
    save_results,
    rotate_geometry,
)
from kitesim import (
    aero2struc_level_1,
    aerodynamic_vsm,
    aerostructural_coupled_solver_qsm,
    read_struc_geometry_yaml_level_1,
    structural_pss,
)
from awetrim.system.system_model import SystemModel
from awetrim.system.tether import RigidLumpedTether


class TimingTracker:
    """Track timing for different components during simulation."""

    def __init__(self):
        self.timings = defaultdict(list)
        self.markers = {}

    def start(self, name: str) -> None:
        """Mark start of a timing block."""
        self.markers[name] = perf_counter()

    def end(self, name: str) -> None:
        """Mark end of a timing block and record elapsed time."""
        if name not in self.markers:
            return
        elapsed = perf_counter() - self.markers[name]
        self.timings[name].append(elapsed)
        del self.markers[name]

    def report(self) -> dict:
        """Generate a timing report."""
        stats = {}
        for name, times in self.timings.items():
            stats[name] = {
                "count": len(times),
                "total_s": sum(times),
                "mean_s": np.mean(times),
                "min_s": np.min(times),
                "max_s": np.max(times),
                "std_s": np.std(times) if len(times) > 1 else 0.0,
            }
        return stats

    def print_report(self) -> None:
        """Print timing report to console."""
        stats = self.report()
        print("\n" + "=" * 70)
        print("TIMING REPORT")
        print("=" * 70)
        total_measured = sum(s["total_s"] for s in stats.values())
        for name in sorted(stats.keys()):
            s = stats[name]
            pct = 100 * s["total_s"] / total_measured if total_measured > 0 else 0
            if s["count"] == 1:
                print(f"{name:40s}: {s['total_s']:8.2f}s")
            else:
                print(
                    f"{name:40s}: {s['total_s']:8.2f}s "
                    f"({s['count']:3d}x, {s['mean_s']:6.2f}s avg, "
                    f"min={s['min_s']:.2f}s, max={s['max_s']:.2f}s)"
                )
            print(f"  {pct:5.1f}% of total measured time")
        print("=" * 70)


def wrap_coupled_solver_with_timing(solver_fn, tracker: TimingTracker):
    """Wrap the coupled solver to track timing at the iteration level."""
    original_fn = solver_fn

    def wrapper(*args, **kwargs):
        tracking_data, meta = original_fn(*args, **kwargs)

        # Extract per-iteration timing from tracking_data if available
        if tracking_data is not None and hasattr(tracking_data, "get"):
            n_iter = meta.get("n_iter", 0)
            if n_iter > 0:
                # Log that solver completed
                tracker.timings["solver_iterations"].append(n_iter)

        return tracking_data, meta

    return wrapper


def run_profiled_simulation(
    config_path,
    struc_geometry_path,
    aero_geometry_path,
    kite_name="TUDELFT_V3_KITE",
    steering_override=None,
    output_dir=None,
    profile_type="instrumented",
):
    """Run a single QSM simulation with detailed profiling."""

    if output_dir is None:
        output_dir = Path("./profile_results")
    else:
        output_dir = Path(output_dir)

    output_dir.mkdir(parents=True, exist_ok=True)

    # Load configuration
    config = load_yaml(config_path)
    struc_geometry = load_yaml(struc_geometry_path)

    # Override steering if provided
    if steering_override is not None:
        config["steering_tape_final_extension"] = float(steering_override)

    PROJECT_DIR = Path(config_path).parents[1]

    print(f"\n{'='*70}")
    print(f"QSM Simulation Profiler")
    print(f"{'='*70}")
    print(f"Config: {config_path}")
    print(f"Steering: {config.get('steering_tape_final_extension', 0.0)} m")
    print(
        f"Aero update interval: {config['aero_structural_solver'].get('aero_update_interval', 1)}"
    )
    print(f"VSM max iterations: {config['aerodynamic'].get('max_iterations', 100)}")
    print(f"{'='*70}\n")

    # Track timing
    tracker = TimingTracker()

    # Load structural geometry
    tracker.start("load_structural_geometry")
    (
        struc_nodes_base,
        m_arr,
        struc_node_le_indices,
        struc_node_te_indices,
        power_tape_index,
        steering_tape_indices,
        pulley_node_indices,
        kite_connectivity_arr,
        bridle_connectivity_arr,
        bridle_diameter_arr,
        l0_arr,
        k_arr,
        c_arr,
        linktype_arr,
        pulley_line_indices,
        pulley_line_to_other_node_pair_dict,
    ) = read_struc_geometry_yaml_level_1.main(struc_geometry, config=config)
    tracker.end("load_structural_geometry")

    # Apply initial rotation
    tracker.start("rotate_initial_geometry")
    struc_nodes_base = rotate_geometry(struc_nodes_base)
    tracker.end("rotate_initial_geometry")

    # Initialize aerodynamic solver
    tracker.start("initialize_aerodynamic_solver")
    n_wing_struc_nodes = len(struc_geometry["wing_particles"]["data"])
    n_struc_ribs = n_wing_struc_nodes / 2
    n_panels_aero = (n_struc_ribs - 1) * config["aerodynamic"][
        "n_aero_panels_per_struc_section"
    ]
    bridle_path = (
        struc_geometry_path if config.get("is_with_aero_bridle", False) else None
    )
    body_aero_init, vsm_solver_init, vel_app, initial_polar_data_init = (
        aerodynamic_vsm.initialize(
            aero_geometry_path, config, n_panels_aero, bridle_path=bridle_path
        )
    )
    tracker.end("initialize_aerodynamic_solver")

    # Initialize aero-to-structural mapping
    tracker.start("initialize_aero2struc_mapping")
    aero2struc_mapping = aero2struc_level_1.initialize_mapping(
        body_aero_init.panels,
        struc_nodes_base,
        struc_node_le_indices,
        struc_node_te_indices,
    )
    tracker.end("initialize_aero2struc_mapping")

    # Initialize structural solver
    tracker.start("instantiate_structural_solver")
    struc_nodes = struc_nodes_base.copy()
    (psystem, pss_initial_conditions, pss_params, struc_nodes_initial) = (
        structural_pss.instantiate(
            config,
            struc_nodes,
            m_arr,
            kite_connectivity_arr,
            l0_arr,
            k_arr,
            c_arr,
            linktype_arr,
            pulley_line_to_other_node_pair_dict,
        )
    )
    tracker.end("instantiate_structural_solver")

    # Set up system model
    tracker.start("setup_system_model")
    tether = RigidLumpedTether(diameter=config["tether"]["diameter"])
    system_model = SystemModel(tether=tether)
    system_model.mass_wing = float(np.sum(m_arr))
    system_model.angle_elevation = np.deg2rad(
        float(config.get("angle_elevation_deg", 30.0))
    )
    system_model.angle_azimuth = np.deg2rad(
        float(config.get("angle_azimuth_deg", 20.0))
    )
    system_model.angle_course = np.deg2rad(float(config.get("angle_course_deg", 90.0)))
    system_model.speed_radial = float(config.get("speed_radial", 0.0))
    system_model.distance_radial = float(config.get("distance_radial", 200.0))
    system_model.wind.speed_wind_ref = float(config.get("wind_speed_wind_ref", 4.0))
    system_model.timeder_speed_tangential = float(
        config.get("timeder_speed_tangential", 0.0)
    )
    system_model.timeder_speed_radial = float(config.get("timeder_speed_radial", 0.0))
    tracker.end("setup_system_model")

    # Run coupled solver with profiling
    print("Running coupled aero-structural solver...")
    tracker.start("coupled_solver_main")

    if profile_type == "cprofile":
        # Use cProfile for detailed profiling
        pr = cProfile.Profile()
        pr.enable()

    tracking_data, meta = aerostructural_coupled_solver_qsm.main(
        m_arr=m_arr,
        struc_nodes=struc_nodes,
        struc_nodes_initial=struc_nodes_initial,
        system_model=system_model,
        config=config,
        initial_length_power_tape=l0_arr[power_tape_index],
        n_power_tape_steps=0,
        power_tape_final_extension=0.0,
        power_tape_extension_step=0.0,
        kite_connectivity_arr=kite_connectivity_arr,
        bridle_connectivity_arr=bridle_connectivity_arr,
        pulley_line_indices=pulley_line_indices,
        pulley_line_to_other_node_pair_dict=pulley_line_to_other_node_pair_dict,
        struc_node_le_indices=struc_node_le_indices,
        struc_node_te_indices=struc_node_te_indices,
        body_aero=copy.deepcopy(body_aero_init),
        vsm_solver=copy.deepcopy(vsm_solver_init),
        vel_app=vel_app,
        initial_polar_data=copy.deepcopy(initial_polar_data_init),
        bridle_diameter_arr=bridle_diameter_arr,
        aero2struc_mapping=aero2struc_mapping,
        power_tape_index=power_tape_index,
        psystem=psystem,
        kite_fem_structure=None,
        initial_qs_guess=None,
    )

    if profile_type == "cprofile":
        pr.disable()

    tracker.end("coupled_solver_main")

    # Save results
    tracker.start("save_results")
    case_folder = (
        f"steering_{config.get('steering_tape_final_extension', 0.0)*1000:.0f}mm"
    )
    case_dir = output_dir / case_folder
    case_dir.mkdir(parents=True, exist_ok=True)
    h5_path = case_dir / "sim_output.h5"
    save_results(tracking_data, meta, h5_path)
    tracker.end("save_results")

    # Generate reports
    print("\n" + "=" * 70)
    print("SIMULATION COMPLETED")
    print("=" * 70)
    print(f"Total iterations: {meta.get('n_iter', 0)}")
    print(f"Converged: {meta.get('converged', False)}")
    print(f"Total time (meta): {meta.get('total_time_s', 0):.2f}s")
    print("=" * 70)

    tracker.print_report()

    # Save detailed timing report
    timing_report = {
        "timestamp": datetime.now().isoformat(),
        "config": {
            "steering_m": float(config.get("steering_tape_final_extension", 0.0)),
            "wind_speed_ms": float(config.get("wind_speed_wind_ref", 4.0)),
            "aero_update_interval": int(
                config["aero_structural_solver"].get("aero_update_interval", 1)
            ),
            "vsm_max_iterations": int(config["aerodynamic"].get("max_iterations", 100)),
        },
        "results": {
            "n_iter": int(meta.get("n_iter", 0)),
            "converged": bool(meta.get("converged", False)),
            "total_time_s": float(meta.get("total_time_s", 0)),
        },
        "timings": tracker.report(),
    }

    report_path = case_dir / "timing_report.json"
    with open(report_path, "w") as f:
        json.dump(timing_report, f, indent=2)

    print(f"\nTiming report saved to: {report_path}")
    print(f"Simulation results saved to: {h5_path}")

    # If using cProfile, save detailed profile
    if profile_type == "cprofile":
        prof_path = case_dir / "cprofile_stats.txt"
        with open(prof_path, "w") as f:
            ps = pstats.Stats(pr, stream=f)
            ps.sort_stats("cumulative")
            ps.print_stats(50)  # Top 50 functions
        print(f"cProfile statistics saved to: {prof_path}")

    return tracking_data, meta, timing_report


def main():
    parser = argparse.ArgumentParser(
        description="Profile a QSM aero-structural simulation"
    )
    parser.add_argument(
        "--steering",
        type=float,
        default=None,
        help="Override steering_tape_final_extension [m]",
    )
    parser.add_argument(
        "--profile-type",
        choices=["cprofile", "instrumented"],
        default="instrumented",
        help="Type of profiling to use",
    )
    parser.add_argument(
        "--output",
        type=str,
        default=None,
        help="Output directory for profiling results",
    )

    args = parser.parse_args()

    PROJECT_DIR = Path(__file__).resolve().parents[1]
    kite_name = "TUDELFT_V3_KITE"

    config_path = PROJECT_DIR / "data" / kite_name / "config.yaml"
    struc_geometry_path = (
        PROJECT_DIR / "data" / kite_name / "struc_geometry_level_1_manual_JULIA.yaml"
    )
    aero_geometry_path = PROJECT_DIR / "data" / kite_name / "aero_geometry.yaml"

    try:
        tracking_data, meta, timing_report = run_profiled_simulation(
            config_path=config_path,
            struc_geometry_path=struc_geometry_path,
            aero_geometry_path=aero_geometry_path,
            kite_name=kite_name,
            steering_override=args.steering,
            output_dir=args.output,
            profile_type=args.profile_type,
        )
    except Exception as e:
        print(f"\nERROR: {e}", file=sys.stderr)
        import traceback

        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
