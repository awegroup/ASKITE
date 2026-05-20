import numpy as np


DEFAULT_KCU = {
    "node_index": 0,
    "length": 1.0,
    "diameter": 0.48,
    "cd_perpendicular": 0.69,
    "cd_parallel": 0.83,
}


def _unit(vec, fallback=None):
    vec = np.asarray(vec, dtype=float).reshape(3)
    norm = np.linalg.norm(vec)
    if norm <= 1e-12:
        if fallback is None:
            fallback = [1.0, 0.0, 0.0]
        return np.asarray(fallback, dtype=float).reshape(3)
    return vec / norm


def infer_kcu_to_kite_axis(struc_nodes, bridle_connectivity_arr, kcu_node_index=0):
    """Infer the KCU-to-kite direction from bridle lines connected to the KCU node."""
    struc_nodes = np.asarray(struc_nodes, dtype=float)
    kcu_node_index = int(kcu_node_index)
    kcu_position = struc_nodes[kcu_node_index]
    connected = []

    if bridle_connectivity_arr is not None:
        for ci, cj in np.asarray(bridle_connectivity_arr, dtype=int).reshape(-1, 2):
            if ci == kcu_node_index:
                connected.append(int(cj))
            elif cj == kcu_node_index:
                connected.append(int(ci))

    if connected:
        target = np.mean(struc_nodes[np.asarray(connected, dtype=int)], axis=0)
    else:
        target = np.mean(struc_nodes, axis=0)
    return _unit(target - kcu_position, fallback=[0.0, 0.0, 1.0])


def compute_kcu_aerodynamic_force(vel_app, kcu_axis, rho, kcu_config=None):
    """
    Compute KCU parasitic drag as a separate finite-cylinder force vector.

    This follows the EKF-AWE two-point kite model convention: the KCU is pitched
    90 degrees relative to the tether/bridle axis. The apparent-wind component
    along the tether therefore uses the cylinder perpendicular coefficient and
    side area; the remaining component uses the tangential coefficient and
    circular frontal area.
    """
    cfg = dict(DEFAULT_KCU)
    if kcu_config:
        cfg.update(kcu_config)

    vel_app = np.asarray(vel_app, dtype=float).reshape(3)
    axis = _unit(kcu_axis, fallback=[0.0, 0.0, 1.0])
    rho = float(rho)
    length = float(cfg["length"])
    diameter = float(cfg["diameter"])
    cd_perpendicular = float(cfg["cd_perpendicular"])
    cd_parallel = float(cfg["cd_parallel"])

    area_perpendicular = diameter * length
    area_parallel = np.pi * (diameter / 2.0) ** 2

    vel_along_tether = np.dot(vel_app, axis) * axis
    vel_across_tether = vel_app - vel_along_tether

    force_perpendicular_to_kcu = (
        0.5
        * rho
        * np.linalg.norm(vel_along_tether)
        * vel_along_tether
        * cd_perpendicular
        * area_perpendicular
    )
    force_parallel_to_kcu = (
        0.5
        * rho
        * np.linalg.norm(vel_across_tether)
        * vel_across_tether
        * cd_parallel
        * area_parallel
    )
    return force_perpendicular_to_kcu + force_parallel_to_kcu


def main(struc_nodes, bridle_connectivity_arr, vel_app, rho, kcu_config=None):
    cfg = dict(DEFAULT_KCU)
    if kcu_config:
        cfg.update(kcu_config)
    node_index = int(cfg["node_index"])
    force = compute_kcu_aerodynamic_force(
        vel_app,
        infer_kcu_to_kite_axis(struc_nodes, bridle_connectivity_arr, node_index),
        rho,
        cfg,
    )
    f_aero_kcu = np.zeros_like(np.asarray(struc_nodes, dtype=float))
    f_aero_kcu[node_index] = force
    return f_aero_kcu
