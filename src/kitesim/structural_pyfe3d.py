import numpy as np


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
    initial_conditions = []
    if config["is_with_initial_point_velocity"]:
        raise ValueError("Error: initial point velocity has never been defined")
    else:
        vel_ini = np.zeros((len(struc_nodes), 3))

    for i in range(len(struc_nodes)):
        fixed = i in struc_geometry["fixed_point_indices"]
        initial_conditions.append([struc_nodes[i], vel_ini[i], m_arr[i], fixed])

    pulley_matrix = []
    spring_matrix = []

    for idx, (cicj, k, c, l0, linktype) in enumerate(
        zip(conn_arr, k_arr, c_arr, l0_arr, linktype_arr)
    ):
        ci, cj = int(cicj[0]), int(cicj[1])
        lt = str(linktype).lower()

        if lt == "pulley":
            # mapping stores: { str(idx) -> [pulley_node(=cj), ck, l0_other] }
            pulley_node, ck, l0_other = pulley_line_to_other_node_pair_dict[str(idx)]
            ck = int(ck)

            # This row already encodes a single pulley with two halves:
            # current half: ci -- cj (rest length l0_1 = l0)
            # other  half: cj -- ck (rest length l0_2 = l0_other)
            l0_1 = float(l0)
            l0_2 = float(l0_other)
            k1 = float(k)
            c1 = float(c)

            # Effective series properties (assuming same EA on both halves):
            # EA = k1 * l0_1  ->  k_eff = EA / (l0_1 + l0_2)
            l0_tot = l0_1 + l0_2
            k_eff = (k1 * l0_1) / l0_tot if l0_tot != 0.0 else 0.0
            alpha = (c1 / k1) if k1 != 0.0 else 0.0  # damping per stiffness
            c_eff = alpha * k_eff

            # pyfe3d pulley format: [n1, pulley_node(=cj), n3, k_eff, c_eff, l0_total]
            pulley_matrix.append([ci, cj, ck, k_eff, c_eff, l0_tot])

        else:
            # spring format: [n1, n2, k, c, l0, springtype]
            spring_matrix.append([ci, cj, float(k), float(c), float(l0), lt])
