#!/usr/bin/env python3
from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import yaml


SCRIPT_DIR = Path(__file__).resolve().parent
KITE_DIR = SCRIPT_DIR.parent
RESULTS_DIR = SCRIPT_DIR / "results"

AERO_GEOMETRY_PATH = (
    KITE_DIR
    / "2D_airfoils_polars_plots_BEST"
    / "aero_geometry_CFD_CAD_derived.yaml"
)
FEM_GEOMETRY_PATH = RESULTS_DIR / "struc_geometry_FEM_full_stretched.yaml"
PSM_GEOMETRY_PATH = (
    RESULTS_DIR / "struc_geometry_PSM_reduced_stretched_mapped_to_FEM_stretched.yaml"
)
OUTPUT_MD = RESULTS_DIR / "area_consistency_report.md"


@dataclass
class Station:
    label: str
    forward: np.ndarray
    rearward: np.ndarray


@dataclass
class PanelArea:
    label: str
    projected_area_m2: float
    flat_area_m2: float


@dataclass
class AreaSummary:
    name: str
    path: Path
    stations: list[Station]
    panels: list[PanelArea]

    @property
    def projected_area_m2(self) -> float:
        return sum(panel.projected_area_m2 for panel in self.panels)

    @property
    def flat_area_m2(self) -> float:
        return sum(panel.flat_area_m2 for panel in self.panels)


def load_yaml(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as stream:
        data = yaml.safe_load(stream)
    if not isinstance(data, dict):
        raise ValueError(f"{path} did not parse as a YAML mapping.")
    return data


def header_map(table: dict[str, Any]) -> dict[str, int]:
    return {str(name): idx for idx, name in enumerate(table["headers"])}


def as_point(values: Iterable[Any]) -> np.ndarray:
    return np.array([float(value) for value in values], dtype=float)


def ordered_forward_rearward(
    label: str,
    point_a: np.ndarray,
    point_b: np.ndarray,
) -> Station:
    if point_a[0] <= point_b[0]:
        return Station(label=label, forward=point_a, rearward=point_b)
    return Station(label=label, forward=point_b, rearward=point_a)


def aero_stations(path: Path) -> list[Station]:
    geometry = load_yaml(path)
    wing_sections = geometry["wing_sections"]
    headers = header_map(wing_sections)
    stations: list[Station] = []

    for section_i, row in enumerate(wing_sections["data"], start=1):
        airfoil_id = int(row[headers["airfoil_id"]])
        le = as_point(
            [
                row[headers["LE_x"]],
                row[headers["LE_y"]],
                row[headers["LE_z"]],
            ]
        )
        te = as_point(
            [
                row[headers["TE_x"]],
                row[headers["TE_y"]],
                row[headers["TE_z"]],
            ]
        )
        stations.append(
            ordered_forward_rearward(
                f"section_{section_i:02d}_airfoil_{airfoil_id}", le, te
            )
        )

    return stations


def structural_wing_stations(path: Path) -> list[Station]:
    geometry = load_yaml(path)
    wing_particles = geometry["wing_particles"]
    headers = header_map(wing_particles)
    rows = sorted(wing_particles["data"], key=lambda row: int(row[headers["id"]]))

    if len(rows) % 2 != 0:
        raise ValueError(
            f"{path} has {len(rows)} wing_particles; expected an even number."
        )

    stations: list[Station] = []
    for row_a, row_b in zip(rows[0::2], rows[1::2]):
        node_a = int(row_a[headers["id"]])
        node_b = int(row_b[headers["id"]])
        point_a = as_point([row_a[headers["x"]], row_a[headers["y"]], row_a[headers["z"]]])
        point_b = as_point([row_b[headers["x"]], row_b[headers["y"]], row_b[headers["z"]]])
        stations.append(
            ordered_forward_rearward(f"nodes_{node_a}_{node_b}", point_a, point_b)
        )

    return stations


def project_points(points: np.ndarray, projection_plane: str) -> np.ndarray:
    if projection_plane == "xy":
        return points[:, [0, 1]]
    if projection_plane == "xz":
        return points[:, [0, 2]]
    if projection_plane == "yz":
        return points[:, [1, 2]]
    raise ValueError(f"Unsupported projection plane: {projection_plane}")


def polygon_area_2d(points: np.ndarray) -> float:
    x = points[:, 0]
    y = points[:, 1]
    return 0.5 * abs(float(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1))))


def triangle_area_3d(a: np.ndarray, b: np.ndarray, c: np.ndarray) -> float:
    return 0.5 * float(np.linalg.norm(np.cross(b - a, c - a)))


def quad_flat_area_3d(points: np.ndarray) -> float:
    return triangle_area_3d(points[0], points[1], points[2]) + triangle_area_3d(
        points[0], points[2], points[3]
    )


def panel_points(station_a: Station, station_b: Station) -> np.ndarray:
    return np.array(
        [
            station_a.forward,
            station_a.rearward,
            station_b.rearward,
            station_b.forward,
        ],
        dtype=float,
    )


def panel_areas(stations: list[Station], projection_plane: str) -> list[PanelArea]:
    panels: list[PanelArea] = []
    for panel_i, (station_a, station_b) in enumerate(
        zip(stations[:-1], stations[1:]),
        start=1,
    ):
        points_3d = panel_points(station_a, station_b)
        points_2d = project_points(points_3d, projection_plane)
        panels.append(
            PanelArea(
                label=f"panel_{panel_i:02d}: {station_a.label} -> {station_b.label}",
                projected_area_m2=polygon_area_2d(points_2d),
                flat_area_m2=quad_flat_area_3d(points_3d),
            )
        )
    return panels


def build_summary(
    name: str,
    path: Path,
    stations: list[Station],
    projection_plane: str,
) -> AreaSummary:
    if len(stations) < 2:
        raise ValueError(f"{name} has fewer than two wing stations.")
    return AreaSummary(
        name=name,
        path=path,
        stations=stations,
        panels=panel_areas(stations, projection_plane),
    )


def pct_delta(value: float, reference: float) -> float:
    if reference == 0.0:
        return float("nan")
    return 100.0 * (value - reference) / reference


def format_delta(value: float, reference: float) -> str:
    return f"{value - reference:+.6f} ({pct_delta(value, reference):+.2f}%)"


def markdown_report(summaries: list[AreaSummary], projection_plane: str) -> str:
    truth = summaries[0]
    lines = [
        "# Area Consistency Report",
        "",
        f"Projection plane: `{projection_plane}`. For the default `xy` plane, "
        "`z` is dropped and the result is the planform area.",
        "",
        "Flat area is computed as the sum of two 3D triangle areas per "
        "LE/TE quadrilateral. The panel vertex order is "
        "`forward_i -> rearward_i -> rearward_j -> forward_j`.",
        "",
        "The aero geometry is treated as the truth. The PSM area is a "
        "reduced-panel approximation because the reduced PSM wing has fewer "
        "stations than the CAD/FEM geometry.",
        "",
        "## Summary",
        "",
        "| Model | Source | Stations | Panels | Projected area [m2] | Projected delta vs aero | Flat area [m2] | Flat delta vs aero | Flat/projected |",
        "| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |",
    ]

    for summary in summaries:
        projected_delta = (
            "-"
            if summary is truth
            else format_delta(summary.projected_area_m2, truth.projected_area_m2)
        )
        flat_delta = (
            "-"
            if summary is truth
            else format_delta(summary.flat_area_m2, truth.flat_area_m2)
        )
        ratio = summary.flat_area_m2 / summary.projected_area_m2
        lines.append(
            "| "
            f"{summary.name} | "
            f"`{summary.path}` | "
            f"{len(summary.stations)} | "
            f"{len(summary.panels)} | "
            f"{summary.projected_area_m2:.6f} | "
            f"{projected_delta} | "
            f"{summary.flat_area_m2:.6f} | "
            f"{flat_delta} | "
            f"{ratio:.6f} |"
        )

    for summary in summaries:
        lines.extend(
            [
                "",
                f"## {summary.name} Panels",
                "",
                "| Panel | Projected area [m2] | Flat area [m2] | Flat/projected |",
                "| --- | ---: | ---: | ---: |",
            ]
        )
        for panel in summary.panels:
            ratio = panel.flat_area_m2 / panel.projected_area_m2
            lines.append(
                "| "
                f"{panel.label} | "
                f"{panel.projected_area_m2:.6f} | "
                f"{panel.flat_area_m2:.6f} | "
                f"{ratio:.6f} |"
            )

    lines.append("")
    return "\n".join(lines)


def print_console_summary(summaries: list[AreaSummary]) -> None:
    truth = summaries[0]
    print("Area consistency summary:")
    print(
        f"- {truth.name}: projected={truth.projected_area_m2:.6f} m2, "
        f"flat={truth.flat_area_m2:.6f} m2"
    )
    for summary in summaries[1:]:
        print(
            f"- {summary.name}: projected={summary.projected_area_m2:.6f} m2 "
            f"({format_delta(summary.projected_area_m2, truth.projected_area_m2)}), "
            f"flat={summary.flat_area_m2:.6f} m2 "
            f"({format_delta(summary.flat_area_m2, truth.flat_area_m2)})"
        )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Compare projected and flat wing panel area between the CAD-derived "
            "aero geometry and stretched FEM/PSM structural geometries."
        )
    )
    parser.add_argument("--aero", type=Path, default=AERO_GEOMETRY_PATH)
    parser.add_argument("--fem", type=Path, default=FEM_GEOMETRY_PATH)
    parser.add_argument("--psm", type=Path, default=PSM_GEOMETRY_PATH)
    parser.add_argument("--output", type=Path, default=OUTPUT_MD)
    parser.add_argument(
        "--projection-plane",
        choices=("xy", "xz", "yz"),
        default="xy",
        help="Plane used for projected area. Default drops z and uses x-y.",
    )
    parser.add_argument(
        "--no-write",
        action="store_true",
        help="Print the summary only; do not write the markdown report.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    summaries = [
        build_summary(
            "Aero CAD-derived truth",
            args.aero,
            aero_stations(args.aero),
            args.projection_plane,
        ),
        build_summary(
            "FEM stretched",
            args.fem,
            structural_wing_stations(args.fem),
            args.projection_plane,
        ),
        build_summary(
            "PSM reduced mapped stretched",
            args.psm,
            structural_wing_stations(args.psm),
            args.projection_plane,
        ),
    ]

    print_console_summary(summaries)

    if not args.no_write:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(
            markdown_report(summaries, args.projection_plane),
            encoding="utf-8",
        )
        print(f"[compare_aero_structural_area] wrote {args.output}")


if __name__ == "__main__":
    main()
