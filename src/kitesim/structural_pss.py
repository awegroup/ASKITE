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
                s=25,
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


def instantiate(
    config,
    struc_geometry,
    struc_nodes,
    m_arr,
    conn_arr,
    l0_arr,
    k_arr,
    c_arr,
    linktype_arr,
    pulley_line_to_other_node_pair_dict,
):
    # TODO: add l0 to the instantiate method and change ParticleSystem accordingly
    pss_connectivity = []
    for cicj, k, c, l0, linktype in zip(conn_arr, k_arr, c_arr, l0_arr, linktype_arr):
        pss_connectivity.append(
            [
                int(cicj[0]),
                int(cicj[1]),
                float(k),
                float(c),
                # float(l0),
                SpringDamperType(linktype.lower()),
            ]
        )

    pss_initial_conditions = []
    if config["is_with_initial_point_velocity"]:
        raise ValueError("Error: initial point velocity has never been defined")
    else:
        vel_ini = np.zeros((len(struc_nodes), 3))

    for i in range(len(struc_nodes)):
        if i in struc_geometry["fixed_point_indices"]:
            pss_initial_conditions.append([struc_nodes[i], vel_ini[i], m_arr[i], True])
        else:
            pss_initial_conditions.append([struc_nodes[i], vel_ini[i], m_arr[i], False])

    pss_params = {
        "pulley_other_line_pair": pulley_line_to_other_node_pair_dict,
        # "l0": l0_arr,
        # "c": config["structural"]["damping_constant"],
        # "is_with_visc_damping": config["structural"]["is_with_visc_damping"],
        "dt": config["structural_pss"]["dt"],
        "t_steps": config["structural_pss"]["n_internal_time_steps"],
        "abs_tol": config["structural_pss"]["abs_tol"],
        "rel_tol": config["structural_pss"]["rel_tol"],
        "max_iter": config["structural_pss"]["max_iter"],
        # "n": len(struc_nodes),
        # "g": -config["grav_constant"][2],
    }

    psystem = ParticleSystem(
        pss_connectivity,
        pss_initial_conditions,
        pss_params,
    )
    return (
        psystem,
        pss_connectivity,
        pss_initial_conditions,
        pss_params,
    )


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
    return psystem, converged
