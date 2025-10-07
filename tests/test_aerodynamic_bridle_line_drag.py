import numpy as np

from kitesim.aerodynamic_bridle_line_drag import compute_line_aerodynamic_force


def _expected_force(va, p1, p2, diameter, cd_cable, cf_cable, rho):
    if p1[2] > p2[2]:
        p1, p2 = p2, p1

    va = np.asarray(va, dtype=float)
    p1 = np.asarray(p1, dtype=float)
    p2 = np.asarray(p2, dtype=float)

    segment = p2 - p1
    length = np.linalg.norm(segment)
    ej = segment / length
    speed = np.linalg.norm(va)
    dir_d = va / speed

    cos_theta = np.clip(np.dot(va, ej) / speed, -1.0, 1.0)
    theta = np.arccos(cos_theta)
    sin_theta = np.sin(theta)

    cd_t = cd_cable * sin_theta**3 + np.pi * cf_cable * cos_theta**3
    cl_t = cd_cable * sin_theta**2 * cos_theta - np.pi * cf_cable * sin_theta * cos_theta**2

    dir_l = np.cross(dir_d, np.cross(dir_d, ej))
    dynamic_pressure_area = 0.5 * rho * speed**2 * length * diameter

    return dynamic_pressure_area * (cl_t * dir_l + cd_t * dir_d)


def test_compute_line_aerodynamic_force_parallel_flow():
    cd_cable = 1.1
    cf_cable = 0.01
    rho = 1.225
    diameter = 0.01

    va = np.array([0.0, 0.0, 10.0])
    p1 = np.array([0.0, 0.0, 0.0])
    p2 = np.array([0.0, 0.0, 10.0])

    result = compute_line_aerodynamic_force(va, (p1, p2, diameter), cd_cable, cf_cable, rho)

    length = np.linalg.norm(p2 - p1)
    speed = np.linalg.norm(va)
    dynamic_pressure_area = 0.5 * rho * speed**2 * length * diameter
    expected_drag_coeff = np.pi * cf_cable
    expected = dynamic_pressure_area * expected_drag_coeff * (va / speed)

    np.testing.assert_allclose(result, expected, rtol=1e-12, atol=1e-12)


def test_compute_line_aerodynamic_force_general_case_matches_expected():
    cd_cable = 1.1
    cf_cable = 0.01
    rho = 1.225
    diameter = 0.02

    p1 = np.array([0.0, 0.0, 0.0])
    p2 = np.array([0.0, 4.0, 3.0])
    va = np.array([3.0, 4.0, 0.0])

    expected = _expected_force(va, p1, p2, diameter, cd_cable, cf_cable, rho)
    result = compute_line_aerodynamic_force(va, (p1, p2, diameter), cd_cable, cf_cable, rho)
    swapped = compute_line_aerodynamic_force(va, (p2, p1, diameter), cd_cable, cf_cable, rho)

    np.testing.assert_allclose(result, expected, rtol=1e-12, atol=1e-12)
    np.testing.assert_allclose(swapped, expected, rtol=1e-12, atol=1e-12)
