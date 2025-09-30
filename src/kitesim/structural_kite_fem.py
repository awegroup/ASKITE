import numpy as np
from kite_fem.FEMStructure import FEM_structure
from matplotlib import pyplot as plt


def instantiate(
    config,
    struc_geometry,
    struc_nodes,
    conn_arr,
    l0_arr,
    k_arr,
    c_arr,
    m_arr,
    linktype_arr,
    pulley_line_to_other_node_pair_dict,
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

    # Deduplicate pulleys: remember which (ci,cj,ck) we’ve emitted already
    seen_pulley_triplets = set()

    for idx, (cicj, k, c, l0, linktype) in enumerate(
        zip(conn_arr, k_arr, c_arr, l0_arr, linktype_arr)
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
            ##TODO: fix this beun oplossing
            if ci_map != cj_map:
                pulley_matrix.append([ci_map, cj_map, ck, k_eff, c_eff, l0_total])

        else:
            # Regular spring: [ci, cj, k, c, l0, springtype]
            spring_matrix.append([ci, cj, float(k), float(c), float(l0), lt])

    ##TODO: use the defined inputs to initialize the pyfe3d
    initial_conditions = initial_conditions  # [[x,y,z,vel_x,vel_y,vel_z,m,fixed]]
    pulley_matrix = pulley_matrix  # [[ci, cj, ck, k_eff, c_eff, l0_total], ...]
    spring_matrix = spring_matrix  # [[ci, cj, k, c, l0, springtype], ...]

    kite_fem_structure = FEM_structure(
        initial_conditions=initial_conditions,
        spring_matrix=spring_matrix,
        pulley_matrix=pulley_matrix,
    )

    # print(f"initial_conditions: {initial_conditions[0:3]}")
    # breakpoint()
    # print(f"\n pulley_matrix: {pulley_matrix[0:3]}")
    # for pulley in pulley_matrix:
    # print(f"  pulley: {pulley}")

    # print(f"\n spring_matrix: {spring_matrix[0:10]}")

    # breakpoint()

    # if is_plot is True:
    #     kite_fem_structure.plot_3D(color="blue")
    #     plt.show()
    #     plt.close()

    return kite_fem_structure, initial_conditions, pulley_matrix, spring_matrix


def extract_rest_length(kite_fem_structure):
    """
    Extracts the rest lengths of the spring elements in the FEM structure.

    Args:
        kite_fem_structure (FEM_structure): The FEM structure object containing spring elements.
    Returns:

    """
    return np.array([link.l0 for link in kite_fem_structure.spring_elements])
