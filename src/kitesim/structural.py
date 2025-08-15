import logging
import numpy as np
import matplotlib.pyplot as plt
from PSS.particleSystem.SpringDamper import SpringDamperType
from PSS.particleSystem import ParticleSystem


def plot_3d_kite_structure(
    struc_nodes, connectivity, power_tape_index, fixed_nodes=None, pulley_nodes=None
):
    """
    Plot the 3D structure of a kite with enhanced visualization features.

    Args:
        struc_nodes (np.ndarray): Array of 3D coordinates for each node (n_nodes, 3).
        connectivity (list): List of [i, j, k, c, type] for each connection.
        fixed_nodes (iterable, optional): Indices of fixed nodes.
        pulley_nodes (iterable, optional): Indices of pulley nodes.

    Returns:
        None. Displays a 3D plot.
    """
    fig = plt.figure(figsize=(12, 10))
    ax = fig.add_subplot(111, projection="3d")

    # Create sets for fixed and pulley nodes if not provided
    if fixed_nodes is None:
        fixed_nodes = set()
    else:
        fixed_nodes = set(np.atleast_1d(fixed_nodes))

    if pulley_nodes is None:
        pulley_nodes = set()
    else:
        pulley_nodes = set(np.atleast_1d(pulley_nodes))

    # Extract node masses from connectivity (using the m_array would be better if available)
    node_masses = {}
    for conn in connectivity:
        i, j = int(conn[0]), int(conn[1])
        if hasattr(conn[4], "value"):
            link_type = conn[4].value
        else:
            link_type = conn[4]

        # Initialize masses if not already in dictionary
        if i not in node_masses:
            node_masses[i] = 0
        if j not in node_masses:
            node_masses[j] = 0

    # Create sets to track which elements are tubular frame, te_lines, or other noncompressive
    tubular_frame_nodes = set()
    te_line_nodes = set()
    pulley_line_nodes = set()

    # Line style mapping
    line_styles = {
        "default": {
            "color": "black",
            "linestyle": "-",
            "linewidth": 2.5,
            "label": "Tubular Frame",
        },
        "noncompressive": {"color": "green", "linestyle": "-", "linewidth": 1.5},
        "pulley": {
            "color": "purple",
            "linestyle": "-",
            "linewidth": 1.5,
            "label": "Pulley Lines",
        },
    }

    # Track which labels have been used
    used_labels = set()

    # First pass to identify TE lines and bridle lines
    # This is necessary because we need to know which noncompressive lines are TE lines before plotting
    for conn in connectivity:
        i, j = int(conn[0]), int(conn[1])

        if hasattr(conn[4], "value"):
            link_type = conn[4].value
        else:
            link_type = conn[4]

        # Mark te_line_idx_list nodes (this is a placeholder - in actual code,
        # we would use the te_line_idx_list parameter to identify TE lines)
        # For now, we're just propagating the te_line_nodes set
        if link_type.lower() == "noncompressive":
            if i in te_line_nodes or j in te_line_nodes:
                te_line_nodes.add(i)
                te_line_nodes.add(j)

    # Plot connections with appropriate styling
    for conn in connectivity:
        i, j = int(conn[0]), int(conn[1])
        k, c = float(conn[2]), float(conn[3])

        if hasattr(conn[4], "value"):
            link_type = conn[4].value
        else:
            link_type = conn[4]

        x_vals = [struc_nodes[i][0], struc_nodes[j][0]]
        y_vals = [struc_nodes[i][1], struc_nodes[j][1]]
        z_vals = [struc_nodes[i][2], struc_nodes[j][2]]

        # Default styling
        style = line_styles.get(
            link_type.lower(), {"color": "gray", "linestyle": "-", "linewidth": 1}
        )

        # Separate noncompressive elements into TE lines and bridle lines
        if link_type.lower() == "noncompressive":
            if i in te_line_nodes or j in te_line_nodes:
                style["color"] = "orange"
                if "Canopy TE" not in used_labels:
                    style["label"] = "Canopy TE"
                    used_labels.add("Canopy TE")
                else:
                    style.pop("label", None)

            else:
                style["color"] = "blue"
                if "Bridle Lines" not in used_labels:
                    style["label"] = "Bridle Lines"
                    used_labels.add("Bridle Lines")
                else:
                    style.pop("label", None)

        # Track nodes for tubular frame and pulley lines
        if link_type.lower() == "default":
            tubular_frame_nodes.add(i)
            tubular_frame_nodes.add(j)
            if "Tubular Frame" not in used_labels:
                used_labels.add("Tubular Frame")
            else:
                style.pop("label", None)

        if link_type.lower() == "pulley":
            pulley_line_nodes.add(i)
            pulley_line_nodes.add(j)
            if "Pulley Lines" not in used_labels:
                used_labels.add("Pulley Lines")
            else:
                style.pop("label", None)

        # Include damping in the label if requested
        if "label" in style and "damping" not in style["label"]:
            style["label"] += f" (k={k:.1f}, c={c:.2f})"

        # Plot the line
        ax.plot(x_vals, y_vals, z_vals, **style)

    # Create legend labels for nodes
    node_handles = []
    node_labels = []

    # Plot nodes - separate loop to ensure nodes are drawn on top of lines
    for i, point in enumerate(struc_nodes):
        # Plot the index of the node
        ax.text(
            point[0] + 0.02,
            point[1] + 0.02,
            point[2] + 0.02,
            str(i),
            color="black",
            fontsize=6,
        )
        if i in fixed_nodes:
            marker = ax.scatter(
                point[0],
                point[1],
                point[2],
                color="red",
                s=5,
                label="",  # We'll add to legend separately
            )
            if "Fixed Node" not in used_labels:
                node_handles.append(marker)
                node_labels.append("Fixed Node")
                used_labels.add("Fixed Node")
        elif i in pulley_nodes:
            marker = ax.scatter(
                point[0],
                point[1],
                point[2],
                color="purple",
                s=8,
                label="",  # We'll add to legend separately
            )
            if "Pulley Node" not in used_labels:
                node_handles.append(marker)
                node_labels.append("Pulley Node")
                used_labels.add("Pulley Node")
        else:
            marker = ax.scatter(
                point[0],
                point[1],
                point[2],
                color="black",
                s=8,
                label="",  # We'll add to legend separately
            )
            if "Free Node" not in used_labels:
                node_handles.append(marker)
                node_labels.append("Free Node")
                used_labels.add("Free Node")

    for idx, (i, j, k, _, line_type) in enumerate(connectivity):
        # Get coordinates
        p1 = np.array(struc_nodes[i])
        p2 = np.array(struc_nodes[j])

        # Midpoint for label
        midpoint = (p1 + p2) / 2
        # label = f"{line_type.name}\nk={k:.1e}"
        label = f"{idx}"
        # compute distance between p1 and p2
        distance = np.linalg.norm(p2 - p1)
        label = f"{1e3*distance:.1f} mm"

        # Add label slightly offset from midpoint
        offset = 0.02 * np.linalg.norm(p2 - p1)
        ax.text(
            midpoint[0] + offset,
            midpoint[1] + offset,
            midpoint[2] + offset,
            label,
            fontsize=6,
            color="blue",
        )
        if power_tape_index == idx:
            # Highlight the power tape line
            ax.plot(
                [p1[0], p2[0]],
                [p1[1], p2[1]],
                [p1[2], p2[2]],
                color="red",
                linestyle="-",
                linewidth=3,
                label="Power Tape",
            )
            if "Power Tape" not in used_labels:
                used_labels.add("Power Tape")

    # Set labels and title
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_zlabel("Z")
    ax.set_title("3D Kite Structure")

    # Equal aspect ratio
    if hasattr(struc_nodes, "max") and hasattr(struc_nodes, "min"):
        bb = struc_nodes.max(axis=0) - struc_nodes.min(axis=0)
        ax.set_box_aspect(bb)
    else:
        # If points is not a numpy array with max/min methods
        struc_nodes_arr = np.array(struc_nodes)
        bb = struc_nodes_arr.max(axis=0) - struc_nodes_arr.min(axis=0)
        ax.set_box_aspect(bb)

    # Add legend - use a separate legend for nodes
    # Get existing handles and labels from the lines
    handles, labels = ax.get_legend_handles_labels()

    # Combine with node handles and labels
    all_handles = handles + node_handles
    all_labels = labels + node_labels

    # Create legend outside the plot area to ensure visibility
    plt.legend(
        all_handles,
        all_labels,
        loc="upper left",
        bbox_to_anchor=(1.05, 1),
        borderaxespad=0,
    )

    # Adjust layout to make room for the legend
    plt.tight_layout(rect=[0, 0, 0.85, 1])  # Leave space on the right for the legend

    plt.show()


def instantiate_psystem(
    config_dict,
    geometry_dict,
    struc_nodes,
    wing_connectivity,
    kite_connectivity,
    rest_lengths,
    m_array,
    tubular_frame_line_idx_list,
    te_line_idx_list,
    pulley_point_indices,
    pulley_line_indices,
    pulley_line_to_other_node_pair_dict,
    power_tape_index,
):
    """
    Instantiate the particle system for the structural solver.

    Args:
        config_dict (dict): Main configuration dictionary.
        geometry_dict (dict): Kite-specific configuration dictionary.
        struc_nodes (np.ndarray): Initial node positions (n_nodes, 3).
        wing_connectivity (np.ndarray): Wing connectivity matrix.
        kite_connectivity (np.ndarray): Full kite connectivity matrix.
        rest_lengths (np.ndarray): Rest lengths for all springs.
        m_array (np.ndarray): Mass array for all nodes.
        tubular_frame_line_idx_list (list): Indices for tubular frame lines.
        te_line_idx_list (list): Indices for trailing edge lines.
        pulley_point_indices (list): Indices for pulley points.
        pulley_data (dict): Pulley connectivity and mapping data.

    Returns:
        psystem (ParticleSystem): Instantiated particle system.
        pss_param_dict (dict): Parameters for the particle system.
        pss_kite_connectivity (list): Connectivity matrix for the particle system.
    """

    pss_param_dict = {
        "pulley_other_line_pair": pulley_line_to_other_node_pair_dict,
        "l0": rest_lengths,
        "c": config_dict["structural"]["damping_constant"],
        "is_with_visc_damping": config_dict["structural"]["is_with_visc_damping"],
        "dt": config_dict["structural"]["dt"],
        "t_steps": config_dict["structural"]["n_internal_time_steps"],
        "abs_tol": config_dict["structural"]["abs_tol"],
        "rel_tol": config_dict["structural"]["rel_tol"],
        "max_iter": config_dict["structural"]["max_iter"],
        "n": len(struc_nodes),
        "g": -config_dict["grav_constant"][2],
    }

    # creating PSS style connectivity matrix
    pss_kite_connectivity = []
    for idx, _ in enumerate(kite_connectivity):
        if idx in tubular_frame_line_idx_list:
            # compression and tension
            k = geometry_dict["tubular_frame_stiffness"]
            c = geometry_dict["tubular_frame_damping"]
            linktype = "default"
        elif idx in te_line_idx_list:
            # only tension
            k = geometry_dict["trailing_edge_stiffness"]
            c = geometry_dict["trailing_edge_damping"]
            linktype = "noncompressive"
        elif idx in pulley_line_indices:
            # pulley
            k = geometry_dict["bridle_line_stiffness"]
            c = geometry_dict["bridle_line_damping"]
            linktype = "pulley"
        else:
            # only compression
            k = geometry_dict["bridle_line_stiffness"]
            c = geometry_dict["canopy_damping"]
            linktype = "noncompressive"

        pss_kite_connectivity.append(
            [
                int(kite_connectivity[idx][0]),
                int(kite_connectivity[idx][1]),
                float(k),
                float(c),
                SpringDamperType(linktype.lower()),
            ]
        )

    ## INITIAL CONDITIONS
    if config_dict["is_with_initial_point_velocity"]:
        raise ValueError("Error: initial point velocity has never been defined")
    else:
        vel_ini = np.zeros((len(struc_nodes), 3))

    fixed_nodes = np.array(geometry_dict["fixed_node_indices"])
    pss_initial_conditions = []
    n = len(struc_nodes)
    for i in range(n):
        if i in fixed_nodes:
            pss_initial_conditions.append(
                [struc_nodes[i], vel_ini[i], m_array[i], True]
            )
        else:
            pss_initial_conditions.append(
                [struc_nodes[i], vel_ini[i], m_array[i], False]
            )

    psystem = ParticleSystem(
        pss_kite_connectivity,
        pss_initial_conditions,
        pss_param_dict,
    )

    # TODO: deal with the rest length problem
    # updating all the rest lengths to the user input, instead of based on initial node-to-node distance
    set_rest_lengths = psystem.extract_rest_length

    # update rest lengths to original wing rest lengths
    # wing_rest_lengths = np.array(
    #     [
    #         0.34867407,
    #         1.51367733,
    #         0.88927689,
    #         1.72113939,
    #         1.32426864,
    #         1.32890591,
    #         2.24038103,
    #         1.32801374,
    #         1.31377741,
    #         2.49418968,
    #         1.32335236,
    #         1.31128284,
    #         2.61551116,
    #         1.328504,
    #         1.321442,
    #         2.61551116,
    #         1.32335236,
    #         1.31128284,
    #         2.49418968,
    #         1.32801374,
    #         1.31377741,
    #         2.24038103,
    #         1.32426864,
    #         1.32890591,
    #         1.72113939,
    #         1.51367733,
    #         0.88927689,
    #         0.34867407,
    #     ]
    # )
    # rest_lengths[1:29] = wing_rest_lengths

    for idx, curr_set_rest_length in enumerate(psystem.extract_rest_length):
        #     # delta = rest_lengths[idx] - set_res_len
        #     # if np.abs(delta) > 0.25:
        #     #     print(f"\nci,cj: {kite_connectivity[idx][0]}, {kite_connectivity[idx][1]}")
        #     #     print(f"set res len: {set_res_len}")
        #     #     print(f"rest length: {rest_lengths[idx]}")
        #     #     print(f"Delta l0: {rest_lengths[idx] - set_res_len}")
        psystem.update_rest_length(idx, rest_lengths[idx] - curr_set_rest_length)
    # print(
    #     f"ci,cj: {kite_connectivity[idx][0]}, {kite_connectivity[idx][1]},curr rest length: {set_res_len:.3f}, set rest length: {rest_lengths[idx]:.3f}"
    # )

    # # 3% of TE canopy billowing induced extra rest length
    # if idx in te_line_idx_list:
    #     psystem.update_rest_length(idx, 0.03 * curr_set_rest_length)

    if config_dict["is_with_initial_structure_plot"]:
        plot_3d_kite_structure(
            struc_nodes,
            pss_kite_connectivity,
            power_tape_index,
            fixed_nodes=fixed_nodes,
            pulley_nodes=pulley_point_indices,
        )

    return psystem, pss_param_dict, pss_kite_connectivity


def run_pss(psystem, params, f_ext):
    """
    Run the particle system simulation with kinetic damping until convergence.

    Args:
        psystem (ParticleSystem): The particle system to simulate.
        params (dict): Simulation parameters.
        f_ext (np.ndarray): Flattened external force vector (n_nodes*3,).

    Returns:
        psystem (ParticleSystem): The updated particle system after simulation.
    """
    t_vector_internal = np.linspace(
        params["dt"], params["t_steps"] * params["dt"], params["t_steps"]
    )
    E_kin = []
    f_int = []
    E_kin_tol = 1e-3  # 1e-29

    logging.debug(f"Running PS simulation, f_int: {psystem.f_int}")

    # And run the simulation
    for step_internal in t_vector_internal:
        psystem.kin_damp_sim(f_ext)

        E_kin.append(np.linalg.norm(psystem.x_v_current[1] ** 2))
        f_int.append(np.linalg.norm(psystem.f_int))

        converged = False
        if step_internal > 10:
            if np.max(E_kin[-10:-1]) <= E_kin_tol:
                converged = True
        if converged and step_internal > 1:
            # print("Kinetic damping PS converged", step_internal)
            break
    return psystem
