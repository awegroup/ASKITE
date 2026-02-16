import numpy as np
import matplotlib.pyplot as plt
from kite_fem.FEMStructure import FEM_structure
from kite_fem.Functions import adapt_stiffnesses
from kite_fem.Plotting import plot_structure


def fix_nodes(kite: FEM_structure, indices):
    """Fix all 6 DOFs for a list of node indices."""
    for node_id in indices:
        kite.bc[6 * node_id : 6 * node_id + 6] = False
        kite.fixed[6 * node_id : 6 * node_id + 6] = True
    return kite


# TODO: should go back to kite_fem
def relaxbridles(
    kite: FEM_structure,
    canopy_nodes,
    origin,
    pull_down_force_z=-1500.0,
    settle_force_z=-1500.0,
):
    """
    Relax bridles by fixing canopy nodes and applying a vertical force on origin nodes.

    Args:
        kite: FEM structure instance.
        canopy_nodes: Node ids fixed during relaxation.
        origin: Node ids that get the vertical pull force.
        pull_down_force_z: Main Z-force used to pull origin node(s).
        settle_force_z: Smaller Z-force used in a second settle solve.
    """
    kite = fix_nodes(kite, canopy_nodes)
    initial_conditions = kite.initial_conditions
    pulley_matrix = kite.pulley_matrix
    spring_matrix = kite.spring_matrix
    beam_matrix = kite.beam_matrix
    fe = np.zeros(kite.N)

    for node_id in origin:
        kite.bc[6 * node_id + 2] = True
        kite.fixed[6 * node_id + 2] = False
        fe[6 * node_id + 2] = float(pull_down_force_z)

    kite.solve(fe, max_iterations=300, tolerance=0.01, print_info=False)
    for node_id in origin:
        fe[6 * node_id + 2] = float(settle_force_z)
    kite.solve(fe, max_iterations=300, tolerance=0.01, print_info=False)

    newcoords = np.reshape(kite.coords_current, (-1, 3))

    initial_conditions_new = []
    for node_id, (_, vel, mass, fixed) in enumerate(initial_conditions):
        posnew = np.array(newcoords[node_id], copy=True)
        initial_conditions_new.append([posnew, vel, mass, fixed])

    kite_relaxed = FEM_structure(
        initial_conditions_new, spring_matrix, pulley_matrix, beam_matrix
    )

    # Ensure node 0 (or provided origin node) ends exactly at [0, 0, 0].
    return _recenter_structure_node_to_origin(kite_relaxed, node_idx=origin[0])


def _recenter_structure_node_to_origin(kite_fem_structure, node_idx=0):
    """
    Translate the full structure so node_idx is exactly at [0, 0, 0].
    Rebuilds FEM_structure so all internal arrays remain consistent.
    """
    coords = kite_fem_structure.coords_init.reshape(-1, 3).copy()
    translation = -coords[node_idx]
    coords += translation

    initial_conditions_new = []
    for i, (_, vel, mass, fixed) in enumerate(kite_fem_structure.initial_conditions):
        initial_conditions_new.append(
            [coords[i], np.array(vel, copy=True), mass, fixed]
        )

    return FEM_structure(
        initial_conditions=initial_conditions_new,
        spring_matrix=kite_fem_structure.spring_matrix,
        pulley_matrix=kite_fem_structure.pulley_matrix,
        beam_matrix=kite_fem_structure.beam_matrix,
    )


def instantiate(
    config,
    struc_geometry,
    struc_nodes,
    kite_connectivity_arr,
    l0_arr,
    k_arr,
    c_arr,
    m_arr,
    linktype_arr,
    pulley_line_to_other_node_pair_dict,
    canopy_sections,
    strut_sections,
):

    # --- initial conditions ---
    initial_conditions = []
    if config.get("is_with_initial_point_velocity"):
        raise ValueError("Error: initial point velocity has never been defined")
    vel_ini = np.zeros((len(struc_nodes), 3))

    fixed_set = set(int(i) for i in struc_geometry.get("fixed_point_indices", []))
    for i in range(len(struc_nodes)):
        fixed = i in fixed_set
        initial_conditions.append([struc_nodes[i], vel_ini[i], m_arr[i], fixed])

    pulley_matrix = []
    spring_matrix = []
    beam_matrix = []
    # Deduplicate pulleys: remember which (ci,cj,ck) we’ve emitted already
    seen_pulley_triplets = set()

    for idx, (cicj, k, c, l0, linktype) in enumerate(
        zip(kite_connectivity_arr, k_arr, c_arr, l0_arr, linktype_arr)
    ):
        ci, cj = int(cicj[0]), int(cicj[1])
        lt = str(linktype).lower()

        if lt == "pulley":
            # Expect mapping: { str(idx) : [cj, ck, l0_other, ci, cj, ck] }
            map_val = pulley_line_to_other_node_pair_dict.get(str(idx))
            if map_val is None:
                # No mapping (indexing mismatch) → skip or raise
                # Here we skip gracefully and treat as spring to avoid crashing:
                spring_matrix.append([ci, cj, float(k), float(c), float(l0), lt])
                continue

            # Extract the full triplet from the mapping to avoid ambiguity
            # (the last three elements are [ci, cj, ck] in your initializer)
            try:
                ci_map, cj_map, ck = int(map_val[3]), int(map_val[4]), int(map_val[5])
            except Exception:
                # Fallback to using [cj, ck] plus current ci
                cj_map, ck = int(map_val[0]), int(map_val[1])
                ci_map = ci

            triplet = (ci_map, cj_map, ck)
            if triplet in seen_pulley_triplets:
                continue
            seen_pulley_triplets.add(triplet)

            # Recover EA and consistent effective properties
            l0_total = float(l0)  # total rest length of the *whole* line
            k1 = float(k)  # equals EA / l0_total
            EA = k1 * l0_total
            k_eff = 0.0 if l0_total == 0.0 else EA / l0_total  # == k1
            alpha = (float(c) / k1) if k1 != 0.0 else 0.0  # damping per stiffness
            c_eff = alpha * k_eff

            # pyfe3d pulley: [ci, cj, ck, k_eff, c_eff, l0_total]
            ##TODO: fix not so clean solution
            if ci_map != cj_map:
                pulley_matrix.append([ci_map, cj_map, ck, k_eff, c_eff, l0_total])
        elif lt == "inflatable_beam":
            diameter = k
            pressure = c
            beam_matrix.append([ci, cj, float(diameter), float(pressure), float(l0)])
        else:
            # Regular spring: [ci, cj, k, c, l0, springtype]
            spring_matrix.append([ci, cj, float(k), float(c), float(l0), lt])

    # initial_conditions = initial_conditions  # [[x,y,z,vel_x,vel_y,vel_z,m,fixed]]
    # pulley_matrix = pulley_matrix  # [[ci, cj, ck, k_eff, c_eff, l0_total], ...]
    # spring_matrix = spring_matrix  # [[ci, cj, k, c, l0, springtype], ...]

    kite_fem_structure = FEM_structure(
        initial_conditions=initial_conditions,
        spring_matrix=spring_matrix,
        pulley_matrix=pulley_matrix,
        beam_matrix=beam_matrix,
    )

    # Relax the bridle lines
    canopy_nodes = list(
        set([node for section in canopy_sections + strut_sections for node in section])
    )
    pull_down_force_z = config.get("structural_kite_fem", {}).get(
        "relaxbridles_pull_down_force_z", -100.0
    )
    settle_force_z = config.get("structural_kite_fem", {}).get(
        "relaxbridles_settle_force_z", -1.0
    )
    kite_fem_structure = relaxbridles(
        kite_fem_structure,
        canopy_nodes,
        [0],
        pull_down_force_z=pull_down_force_z,
        settle_force_z=settle_force_z,
    )
    struc_nodes_initial = kite_fem_structure.coords_init.reshape(-1, 3)

    if config.get("is_with_initial_structure_plot", False):
        ax, fig = plot_structure(kite_fem_structure, plot_node_numbers=True)
        ax.set_title("Initial structure")
        plt.show()

    return (
        kite_fem_structure,
        initial_conditions,
        pulley_matrix,
        spring_matrix,
        struc_nodes_initial,
    )


def get_rest_lengths(kite_fem_structure, kite_connectivity_arr):
    # Build lookup of ASKITE connectivity key -> rest length.
    # Use both spring and beam elements so inflatable_beam entries are covered.
    l0_map = {}

    # Spring elements (including pulley elements)
    spring_l0s = kite_fem_structure.modify_get_spring_rest_length()
    for spring_element, l0 in zip(kite_fem_structure.spring_elements, spring_l0s):
        n1 = int(spring_element.spring.n1)
        n2 = int(spring_element.spring.n2)
        key = (min(n1, n2), max(n1, n2))
        if spring_element.springtype == "pulley":
            l0 = l0 / 2
        l0_map[key] = float(l0)

    # Beam elements (inflatable beams)
    for beam_element in kite_fem_structure.beam_elements:
        n1 = int(beam_element.beam.n1)
        n2 = int(beam_element.beam.n2)
        key = (min(n1, n2), max(n1, n2))
        l0_map[key] = float(beam_element.L)

    # Map to ASKITE connectivity ordering. Use NaN for unmatched keys to keep
    # the array numeric and HDF5-compatible.
    mapped_l0s = []
    for connectivity in kite_connectivity_arr:
        n1c, n2c = int(connectivity[0]), int(connectivity[1])
        key = (min(n1c, n2c), max(n1c, n2c))
        mapped_l0s.append(l0_map.get(key, np.nan))

    return np.asarray(mapped_l0s, dtype=np.float64)


def run_kite_fem(
    kite_fem_structure,
    f_ext_flat,
    config_structural_kite_fem,
):

    # [fx, fy, fz, mx, my, mz] for each node
    f_ext_reshaped = f_ext_flat.reshape(-1, 3)
    fe_6d = [[fe[0], fe[1], fe[2], 0, 0, 0] for fe in f_ext_reshaped]
    fe_6d = np.array(fe_6d).flatten()

    is_structural_converged, residual = kite_fem_structure.solve(
        fe=fe_6d,
        max_iterations=config_structural_kite_fem["max_iterations"],
        tolerance=config_structural_kite_fem["tolerance"],
        step_limit=config_structural_kite_fem["step_limit"],
        relax_init=config_structural_kite_fem["relax_init"],
        relax_update=config_structural_kite_fem["relax_update"],
        k_update=config_structural_kite_fem["k_update"],
        I_stiffness=config_structural_kite_fem["I_stiffness"],
        pseudo_dt=config_structural_kite_fem.get("pseudo_dt", None),
        k_reg_min=config_structural_kite_fem.get("k_reg_min", 0.0),
        print_info=config_structural_kite_fem["print_info"],
    )

    if config_structural_kite_fem["update_stiffness"]:
        adapt_stiffnesses(kite_fem_structure,max_stiffness=config_structural_kite_fem["max_stiffness"])  #increases k every iter for >1% strain springs
        
    struc_nodes = kite_fem_structure.coords_current
    # reshape from flat to (n_nodes, 3)
    struc_nodes = struc_nodes.reshape(-1, 3)
    f_int = -kite_fem_structure.fi
    # set fixed nodes to the values of -fe_6d
    f_int = np.where(kite_fem_structure.bc == True, f_int, -fe_6d)
    # remove moments
    f_int = f_int.reshape(-1, 6)[:, :3].flatten()

    if config_structural_kite_fem.get("is_with_diagnostics", False):
        node_idx = config_structural_kite_fem.get("num_node_diagnostics", 34)
        _diagnose_node_force_balance_fem(
            kite_fem_structure, f_ext_flat, node_idx=node_idx
        )

    return kite_fem_structure, is_structural_converged, struc_nodes, f_int


def _diagnose_node_force_balance_fem(kite_fem_structure, f_ext_flat, node_idx=34):
    """Print per-spring force contributions at a given node for kite_fem diagnostics."""
    n3 = node_idx * 3
    f_ext_node = f_ext_flat[n3 : n3 + 3]
    coords = kite_fem_structure.coords_current  # flat, length = 3 * num_nodes

    f_int_node = np.zeros(3)
    print(f"\n{'='*80}")
    print(f"FORCE BALANCE DIAGNOSTIC (kite_fem) \u2014 Node {node_idx}")
    print(f"  Position: {coords[n3:n3+3]}")
    print(f"  f_ext:    {f_ext_node}  (|f_ext| = {np.linalg.norm(f_ext_node):.4f} N)")
    print(f"{'='*80}")

    for se_idx, se in enumerate(kite_fem_structure.spring_elements):
        n1 = se.spring.n1
        n2 = se.spring.n2
        if n1 != node_idx and n2 != node_idx:
            continue

        # Compute current length and unit vector
        unit, l_current = se.unit_vector(coords)

        # Handle pulley coupling
        l_other_pulley = 0.0
        if se.springtype == "pulley":
            other = kite_fem_structure.spring_elements[se.i_other_pulley]
            _, l_other = other.unit_vector(coords)
            l_other_pulley = l_other

        # Get spring force via the element's own method
        fi_elem = se.spring_internal_forces(coords, l_other_pulley)
        fi_3d = fi_elem[:3]  # [fx, fy, fz] \u2014 force from n1 toward n2

        # Convention: fi is subtracted at n1, added at n2
        if n1 == node_idx:
            contrib = -fi_3d
        else:
            contrib = fi_3d

        f_int_node += contrib

        l0 = se.l0
        k = se.k
        strain = (l_current - l0) / l0 if l0 > 0 else 0
        sign_str = "+" if n2 == node_idx else "-"
        slack_str = " [SLACK]" if se.slack else ""
        pulley_str = (
            f", l_other={l_other_pulley:.6f}" if se.springtype == "pulley" else ""
        )
        f_mag = np.linalg.norm(fi_3d)

        print(
            f"  SE {se_idx:3d} [{n1:2d}\u2192{n2:2d}] {se.springtype:16s} "
            f"l={l_current:.6f} l0={l0:.6f} strain={strain:+.6f} k={k:.1f} "
            f"|F|={f_mag:.4f}N ({sign_str}){pulley_str}{slack_str}"
        )
        print(
            f"         F_contrib = [{contrib[0]:+.4f}, {contrib[1]:+.4f}, {contrib[2]:+.4f}]"
        )

    f_residual_node = f_int_node + f_ext_node
    print(f"{'\u2500'*80}")
    print(
        f"  SUM f_int:      [{f_int_node[0]:+.4f}, {f_int_node[1]:+.4f}, {f_int_node[2]:+.4f}]  |f_int| = {np.linalg.norm(f_int_node):.4f} N"
    )
    print(
        f"  f_ext:          [{f_ext_node[0]:+.4f}, {f_ext_node[1]:+.4f}, {f_ext_node[2]:+.4f}]  |f_ext| = {np.linalg.norm(f_ext_node):.4f} N"
    )
    print(
        f"  f_residual:     [{f_residual_node[0]:+.4f}, {f_residual_node[1]:+.4f}, {f_residual_node[2]:+.4f}]  |f_res| = {np.linalg.norm(f_residual_node):.4f} N"
    )
    print(f"{'='*80}\n")
