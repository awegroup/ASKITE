import numpy as np


EPS = 1e-12


def _as_index_array(indices):
    return np.asarray(indices, dtype=int).reshape(-1)


def _safe_norm(vec):
    return float(np.linalg.norm(np.asarray(vec, dtype=float)))


def _unit(vec, fallback):
    vec = np.asarray(vec, dtype=float).reshape(3)
    norm = np.linalg.norm(vec)
    if norm <= EPS:
        return np.asarray(fallback, dtype=float).reshape(3)
    return vec / norm


def wind_axis_basis(vel_app):
    """
    Build a simple wind-axis basis from the apparent-wind vector.

    Drag is positive along the apparent-wind vector. Lift is perpendicular to
    drag in the x-z plane where possible and is oriented toward positive z.
    """
    drag_axis = _unit(vel_app, [1.0, 0.0, 0.0])
    side_ref = np.array([0.0, 1.0, 0.0])
    side_axis = side_ref - np.dot(side_ref, drag_axis) * drag_axis
    if np.linalg.norm(side_axis) <= EPS:
        side_axis = np.cross([0.0, 0.0, 1.0], drag_axis)
    side_axis = _unit(side_axis, [0.0, 1.0, 0.0])
    lift_axis = np.cross(drag_axis, side_axis)
    lift_axis = _unit(lift_axis, [0.0, 0.0, 1.0])
    if lift_axis[2] < 0:
        lift_axis *= -1.0
        side_axis *= -1.0
    return drag_axis, lift_axis, side_axis


def compute_global_coefficients(
    f_wing_total,
    f_bridle_total=None,
    f_kcu_total=None,
    vel_app=None,
    rho=1.225,
    s_ref_m2=1.0,
):
    """
    Compute integrated force coefficients from wing plus optional bridle/KCU forces.
    """
    f_wing_total = np.asarray(f_wing_total, dtype=float).reshape(3)
    if f_bridle_total is None:
        f_bridle_total = np.zeros(3)
    f_bridle_total = np.asarray(f_bridle_total, dtype=float).reshape(3)
    if f_kcu_total is None:
        f_kcu_total = np.zeros(3)
    f_kcu_total = np.asarray(f_kcu_total, dtype=float).reshape(3)
    f_total = f_wing_total + f_bridle_total + f_kcu_total
    vel_app = np.asarray(vel_app, dtype=float).reshape(3)
    speed = _safe_norm(vel_app)
    s_ref_m2 = float(s_ref_m2)
    rho = float(rho)

    q_s = 0.5 * rho * speed**2 * s_ref_m2
    drag_axis, lift_axis, side_axis = wind_axis_basis(vel_app)

    def _coeffs(force):
        lift = float(np.dot(force, lift_axis))
        drag = float(np.dot(force, drag_axis))
        side = float(np.dot(force, side_axis))
        resultant = _safe_norm(force)
        if q_s <= EPS:
            return lift, drag, side, np.nan, np.nan, np.nan, np.nan
        c_l = lift / q_s
        c_d = drag / q_s
        c_r = resultant / q_s
        glide = c_l / c_d if abs(c_d) > EPS else np.nan
        return lift, drag, side, c_l, c_d, c_r, glide

    lift, drag, side, c_l, c_d, c_r, glide = _coeffs(f_wing_total)
    (
        lift_total,
        drag_total,
        side_total,
        c_l_total,
        c_d_total,
        c_r_total,
        glide_total,
    ) = _coeffs(f_total)

    return {
        "aero_force_wing_total": f_wing_total,
        "aero_force_bridle_total": f_bridle_total,
        "aero_force_kcu_total": f_kcu_total,
        "aero_force_total": f_total,
        "lift_wing": lift,
        "drag_wing": drag,
        "side_force_wing": side,
        "C_L_wing": c_l,
        "C_D_wing": c_d,
        "C_R_wing": c_r,
        "glide_ratio_wing": glide,
        "lift_total_aero": lift_total,
        "drag_total_aero": drag_total,
        "side_force_total_aero": side_total,
        "C_L_total_aero": c_l_total,
        "C_D_total_aero": c_d_total,
        "C_R_total_aero": c_r_total,
        "glide_ratio_total_aero": glide_total,
        "V_a": speed,
    }


def _paired_le_te(nodes, le_indices, te_indices):
    nodes = np.asarray(nodes, dtype=float)
    le_indices = _as_index_array(le_indices)
    te_indices = _as_index_array(te_indices)
    n = min(le_indices.size, te_indices.size)
    le_indices = le_indices[:n]
    te_indices = te_indices[:n]

    order = np.argsort(nodes[le_indices, 1])
    return le_indices[order], te_indices[order]


def compute_projected_span(nodes, le_indices, te_indices):
    nodes = np.asarray(nodes, dtype=float)
    wing_indices = np.unique(
        np.concatenate([_as_index_array(le_indices), _as_index_array(te_indices)])
    )
    if wing_indices.size == 0:
        return np.nan
    y = nodes[wing_indices, 1]
    return float(np.nanmax(y) - np.nanmin(y))


def compute_projected_area_from_le_te(nodes, le_indices, te_indices, axes=(0, 1)):
    nodes = np.asarray(nodes, dtype=float)
    le_idx, te_idx = _paired_le_te(nodes, le_indices, te_indices)
    if le_idx.size < 2 or te_idx.size < 2:
        return np.nan

    polygon = np.vstack([nodes[le_idx][:, axes], nodes[te_idx[::-1]][:, axes]])
    x = polygon[:, 0]
    y = polygon[:, 1]
    return float(0.5 * abs(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1))))


def compute_chord_pitch_angles_deg(nodes, le_indices, te_indices):
    nodes = np.asarray(nodes, dtype=float)
    le_idx, te_idx = _paired_le_te(nodes, le_indices, te_indices)
    if le_idx.size == 0:
        return np.array([], dtype=float)
    chords = nodes[te_idx] - nodes[le_idx]
    return np.rad2deg(np.arctan2(chords[:, 2], chords[:, 0]))


def compute_geometry_metrics(nodes, le_indices, te_indices):
    angles = compute_chord_pitch_angles_deg(nodes, le_indices, te_indices)
    if angles.size == 0:
        midspan_pitch = np.nan
        mean_twist = np.nan
        max_abs_twist = np.nan
    else:
        mid_idx = int(np.argmin(np.abs(np.arange(angles.size) - (angles.size - 1) / 2)))
        midspan_pitch = float(angles[mid_idx])
        twist = angles - midspan_pitch
        mean_twist = float(np.nanmean(twist))
        max_abs_twist = float(np.nanmax(np.abs(twist)))

    return {
        "projected_span_m": compute_projected_span(nodes, le_indices, te_indices),
        "projected_area_m2": compute_projected_area_from_le_te(
            nodes, le_indices, te_indices
        ),
        "midspan_pitch_deg": midspan_pitch,
        "mean_twist_deg": mean_twist,
        "max_abs_twist_deg": max_abs_twist,
    }


def resolve_reference_area(config, nodes, le_indices, te_indices):
    """
    Resolve coefficient reference area from config, falling back to initial projection.
    """
    candidates = []
    if isinstance(config, dict):
        candidates.extend(
            [
                ("S_ref_m2", config.get("S_ref_m2")),
                ("reference_area_m2", config.get("reference_area_m2")),
                ("area_ref_m2", config.get("area_ref_m2")),
            ]
        )
        aero_cfg = config.get("aerodynamic", {})
        if isinstance(aero_cfg, dict):
            candidates.extend(
                [
                    ("aerodynamic.S_ref_m2", aero_cfg.get("S_ref_m2")),
                    ("aerodynamic.reference_area_m2", aero_cfg.get("reference_area_m2")),
                    ("aerodynamic.area_ref_m2", aero_cfg.get("area_ref_m2")),
                ]
            )

    for source, value in candidates:
        try:
            value = float(value)
        except (TypeError, ValueError):
            continue
        if np.isfinite(value) and value > EPS:
            return value, source

    projected_area = compute_projected_area_from_le_te(nodes, le_indices, te_indices)
    if np.isfinite(projected_area) and projected_area > EPS:
        return projected_area, "initial_projected_area_xy_from_le_te"

    return 1.0, "fallback_1_m2"
