import numpy as np
from kite_fem.FEMStructure import FEM_structure
from matplotlib import pyplot as plt


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
    )
    struc_nodes_initial = np.array([node_data[0] for node_data in initial_conditions])
    kite_fem_structure.update_internal_forces()
    
    return (
        kite_fem_structure,
        initial_conditions,
        pulley_matrix,
        spring_matrix,
        struc_nodes_initial,
    )

def get_rest_lengths(kite_fem_structure,kite_connectivity_arr):
    #get connectivity data kite_fem and kite_connectivity
    n1s = []
    n2s = []
    n1_conn = []
    n2_conn = []
    springtypes = []
    for spring_element in kite_fem_structure.spring_elements:
        n1s.append(spring_element.spring.n1)
        n2s.append(spring_element.spring.n2)
        springtypes.append(spring_element.springtype)
    for connectivity in kite_connectivity_arr:
        n1_conn.append(connectivity[0])
        n2_conn.append(connectivity[1])
    
    #get rest lengths
    l0s = kite_fem_structure.modify_get_spring_rest_length()
    
    for i,springtype in enumerate(springtypes):
        if springtype == "pulley":
            l0s[i] /= 2
    
    #Map l0s from kitefem output to askite input
    l0_map = {(min(n1, n2), max(n1, n2)): l0 for n1, n2, l0 in zip(n1s, n2s, l0s)}
    mapped_l0s = []
    for n1c, n2c in zip(n1_conn, n2_conn):
        key = (min(n1c, n2c), max(n1c, n2c))
        l0_val = l0_map.get(key)
        mapped_l0s.append(l0_val)
    mapped_l0s = np.array(mapped_l0s)
    return mapped_l0s

def run_kite_fem(
    kite_fem_structure,
    f_ext_flat,
    config_structural_kite_fem,
):
    # reset to this iteration (sets coords_current to initial coords)
    kite_fem_structure.reset()

    # [fx, fy, fz, mx, my, mz] for each node
    f_ext_reshaped = f_ext_flat.reshape(-1, 3)
    fe_6d = [[fe[0], fe[1], fe[2], 0, 0, 0] for fe in f_ext_reshaped]
    fe_6d = np.array(fe_6d).flatten()

    is_structural_converged = kite_fem_structure.solve(
        fe=fe_6d,
        max_iterations=config_structural_kite_fem["max_iterations"],
        tolerance=config_structural_kite_fem["tolerance"],
        step_limit=config_structural_kite_fem["step_limit"],
        relax_init=config_structural_kite_fem["relax_init"],
        relax_update=config_structural_kite_fem["relax_update"],
        k_update=config_structural_kite_fem["k_update"],
        I_stiffness=config_structural_kite_fem["I_stiffness"],
        print_info=config_structural_kite_fem["print_info"]
    )
    struc_nodes = kite_fem_structure.coords_current
    # reshape from flat to (n_nodes, 3)
    struc_nodes = struc_nodes.reshape(-1, 3)
    f_int = -kite_fem_structure.fi
    #set fixed nodes to the values of -fe_6d 
    f_int = np.where(kite_fem_structure.bu==True,f_int,-fe_6d)
    # remove moments
    f_int = f_int.reshape(-1, 6)[:, :3].flatten()
    return kite_fem_structure, is_structural_converged, struc_nodes, f_int
