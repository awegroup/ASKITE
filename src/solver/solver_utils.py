import numpy as np
import scipy.optimize
from src.coupling import coupling_aero2struc, coupling_struc2aero
from src.aerodynamic import VSM, bridle_line_system_aero, tether_aero


def calculate_fx(
    vel_kite_x,
    vel_app,
    points,
    config,
    input_VSM,
    input_bridle_aero,
):
    vel_app[0] = config.vel_wind[0] - float(vel_kite_x)

    # Struc --> aero
    points_left_to_right = coupling_struc2aero.order_struc_nodes_right_to_left(
        points, config.kite.connectivity.plate_point_indices
    )
    # Wing Aerodynamic
    (
        force_aero_wing_VSM,
        moment_aero_wing_VSM,
        F_rel,
        ringvec,
        controlpoints,
        wingpanels,
        rings,
        coord_L,
        coord_refined,
    ) = VSM.calculate_force_aero_wing_VSM(points_left_to_right, vel_app, input_VSM)
    # Aero --> struc
    force_aero_wing = coupling_aero2struc.aero2struc(
        points,
        config.kite.connectivity.wing_ci,
        config.kite.connectivity.wing_cj,
        config.kite.connectivity.plate_point_indices,
        force_aero_wing_VSM,
        moment_aero_wing_VSM,
        ringvec,
        controlpoints,
    )  # Put in the new positions of the points
    # Bridle Aerodynamics
    if config.is_with_aero_bridle:
        force_aero_bridle = (
            bridle_line_system_aero.calculate_force_aero_bridle_thedens2022(
                points, vel_app, input_bridle_aero
            )
        )
    else:
        force_aero_bridle = [0]

    force_aero = force_aero_wing + force_aero_bridle

    force_aero_x = np.sum(force_aero[:, 0])

    drag_tether = tether_aero.calculate_tether_drag(
        config.tether.diameter,
        config.tether.length,
        config.rho,
        config.aero.cd_cylinder,
        vel_app,
    )
    # print(f"force_aero_x: {force_aero_x:.1f}N and drag_tether: {drag_tether:.1f}N")

    force_resultant_x = drag_tether + force_aero_x
    # print(f"force_resultant_x: {force_resultant_x:.1f}N")
    return force_resultant_x


def optimalisation_of_vk_for_fx_0(
    vk_x_initial_guess,
    vel_app,
    points,
    config,
    input_VSM,
    input_bridle_aero,
    tol_vk_optimization,
):
    vel_wind = config.vel_wind
    n_segments = config.kite.n_segments
    plate_point_indices = config.kite.connectivity.plate_point_indices

    # TODO: remove hardcoding
    # initialising
    # Trying again with conventional functions
    # define the bounds
    lower_bound = -6 * vel_wind[2]
    upper_bound = -vel_wind[2]
    boundData = scipy.optimize.Bounds(lower_bound, upper_bound)

    # #define the linear constraints
    # lowerBoundLinearConstraints=0
    # upperBoundLinearConstraints=12
    # linearConstraints = LinearConstraint(lowerBoundLinearConstraints, upperBoundLinearConstraints)

    # define the nonlinear constraints
    lowerBoundNonLinearConstraints = [0, 0]
    # TODO: remove hardcoding, make this an input parameter
    upperBoundNonLinearConstraints = [12, 0]

    # defining angle of mid-chord
    points_LE_segment = points[plate_point_indices[int(n_segments / 2)]][:2]
    points_TE_segment = points[plate_point_indices[int(n_segments / 2)]][2:]
    LE_point_diff = points_LE_segment[1] - points_LE_segment[0]
    TE_point_diff = points_TE_segment[1] - points_TE_segment[0]
    mid_point_LE = points_LE_segment[0] + 0.5 * LE_point_diff
    mid_point_TE = points_TE_segment[0] + 0.5 * TE_point_diff
    # trim_angle   = np.arctan(mid_point_LE[0]/mid_point_LE[2])
    delta_z_chord = mid_point_LE[2] - mid_point_TE[2]
    delta_x_chord = np.abs(mid_point_LE[0] - mid_point_TE[0])
    angle_of_mid_chord = np.rad2deg(np.arctan(delta_z_chord / delta_x_chord))

    def nonlinearConstraintMiddle(vel_kite_x):
        # unpacking args
        vel_app[0] = vel_wind[0] - vel_kite_x

        # getting the aoa
        angle_va_with_xy = np.rad2deg(np.arctan(vel_app[2] / vel_app[0]))
        aoa = angle_of_mid_chord + angle_va_with_xy

        fx = calculate_fx(
            vel_kite_x,
            vel_app,
            points,
            config,
            input_VSM,
            input_bridle_aero,
        )

        return [aoa, fx]

    nonlinearConstraints = scipy.optimize.NonlinearConstraint(
        nonlinearConstraintMiddle,
        lowerBoundNonLinearConstraints,
        upperBoundNonLinearConstraints,
    )
    # nonlinearconstraint_function_with_args = lambda x: nonlinearConstraintMiddle(x, vel_app, points, arguments)

    # see all the options for the solver
    # scipy.optimize.show_options(solver='minimize',method='trust-constr')

    result = scipy.optimize.minimize(
        calculate_fx,
        x0=vk_x_initial_guess,
        method="trust-constr",
        constraints=[nonlinearConstraints],
        tol=tol_vk_optimization,
        # method="SLSQP",
        args=(vel_app, points, config, input_VSM, input_bridle_aero),
        # options={'verbose': 1},
        bounds=boundData,
    )

    # here is the solution
    # print(f' ')
    print(
        f"--- success: {result.success}, result.x[0]: {result.x[0]:.5f}, fx: {result.fun:.5f}N, aoa: {nonlinearConstraintMiddle(result.x[0])[0]:.5f}deg"
    )
    # print(f' fx check: {calculate_fx(result.x[0],vel_app,points,arguments):.5f}N')

    return result.x[0]


def optimize_fc_to_get_mx_to_zero(moment_aero, points, force_centrifugal):
    args_mc_x = [np.sum(moment_aero[:, 0]), points, force_centrifugal]
    fc_scaling = 1.0
    boundData = scipy.optimize.Bounds(0.95, 1.05)

    def calculate_mx(fc_scaling, args):
        ma_x_sum, points, force_centrifugal = args
        force_centrifugal = fc_scaling * force_centrifugal
        mc_x = np.cross(points, force_centrifugal)[:, 0]
        mc_x_sum = np.sum(mc_x)
        print(f"ma_x_sum: {ma_x_sum}, mc_x_sum:{mc_x_sum}")
        return ma_x_sum + mc_x_sum

    result = scipy.optimize.minimize(
        calculate_mx,
        x0=fc_scaling,
        method="trust-constr",
        options={"disp": True},
        tol=1e-1,
        args=args_mc_x,
        bounds=boundData,
    )
    print(
        f"--- success: {result.success}, result.x: {result.x}, result.f: {result.fun}"
    )
    return result.x * force_centrifugal


# defining several functions
def calculate_center_of_gravity(points, mass_points):
    """calculates the center of gravity of the kite, given the points and the mass of each point"""
    numerator_x = np.sum([mass_points[i] * point[0] for i, point in enumerate(points)])
    numerator_y = np.sum([mass_points[i] * point[1] for i, point in enumerate(points)])
    numerator_z = np.sum([mass_points[i] * point[2] for i, point in enumerate(points)])
    denominator = np.sum(mass_points)
    return np.array([numerator_x, numerator_y, numerator_z]) / denominator


def calculate_force_centrifugal_distribution(force_aero, mass_points, vel_app, points):
    """calculates the centrifugal force distribution, i.e. the force on each node Fc,i"""
    # calculate center of gravity
    center_of_gravity = calculate_center_of_gravity(points, mass_points)

    # initialize some variables
    force_aero_side = np.sum(force_aero[:, 1])
    sum_of_masses = np.sum(mass_points)
    v_k = vel_app[0]

    # calculate turning radius of the whole kite r_k
    r_k = (sum_of_masses * v_k**2) / -force_aero_side

    # calculate distance from turning-axis to x,y plane r_0
    r_0 = r_k - center_of_gravity[1]

    # Calculate kite velocity of each node v_k,i
    r_ratio_array = np.array([((r_0 + point[1]) / r_k) for point in points])
    v_k_i_array = v_k * r_ratio_array

    # calculate centrifugal force distribution, i.e. the force on each node Fc,i
    force_centrifugal_i_array = [
        ((mass_points[i] * (v_k_i_array[i] ** 2)) / (r_0 + point[1]))
        for i, point in enumerate(points)
    ]

    return force_centrifugal_i_array, r_0, v_k_i_array, center_of_gravity


def calculate_force_centrifugal_distribution_moment_based(
    moment_aero, mass_points, vel_app, points
):
    """calculates the centrifugal force distribution, i.e. the force on each node Fc,i"""
    # calculate center of gravity
    center_of_gravity = calculate_center_of_gravity(points, mass_points)

    # initialize some variables
    sum_of_masses = np.sum(mass_points)
    v_k = vel_app[0]

    # calculat the roll moment
    moment_aero_roll = np.sum(moment_aero[:, 0])

    # calculate Fc that would create the same moment as the moment_aero
    force_centrifugal = moment_aero_roll / np.linalg.norm(center_of_gravity)

    # calculate turning radius of the whole kite r_k
    r_k = (sum_of_masses * v_k**2) / force_centrifugal

    # calculate distance from turning-axis to x,y plane r_0
    r_0 = r_k - center_of_gravity[1]

    # Calculate kite velocity of each node v_k,i
    r_ratio_array = np.array([((r_0 + point[1]) / r_k) for point in points])
    v_k_i_array = v_k * r_ratio_array

    # calculate centrifugal force distribution, i.e. the force on each node Fc,i
    force_centrifugal_i_array = [
        ((mass_points[i] * (v_k_i_array[i] ** 2)) / (r_0 + point[1]))
        for i, point in enumerate(points)
    ]

    return force_centrifugal_i_array, r_0, v_k_i_array, center_of_gravity


def calculate_force_centrifugal_distribution_with_tether_tension(
    r_0_initial_guess,
    force_aero,
    mass_points,
    vel_app,
    points,
    length_tether,
    is_print_intermediate_results=True,
):
    """calculates the centrifugal force distribution, i.e. the force on each node Fc,i"""
    # calculate center of gravity
    center_of_gravity = calculate_center_of_gravity(points, mass_points)

    # initialize some variables
    force_aero_sum = np.sum(force_aero, axis=0)
    v_k = -vel_app[0]
    sum_of_masses = np.sum(mass_points)

    def return_zero_for_correct_r_0(
        r_0, sum_of_masses, v_k, length_tether, force_aero_sum
    ):
        # Should be negative (oriented in negative-side/y direction)
        f_aero_side = force_aero_sum[1]
        f_tether_side = -force_aero_sum[2] * np.tan(np.arcsin(r_0 / length_tether))
        # Should be positive (oriented in positive-side/y direction)
        f_centrifugal_side = sum_of_masses * (v_k**2 / r_0)

        ## now r_0 will be chosen such that it balances these two contributions
        return f_centrifugal_side + f_tether_side + f_aero_side

    # Find the root numerically and get full output
    result = scipy.optimize.root(
        return_zero_for_correct_r_0,
        r_0_initial_guess,
        args=(sum_of_masses, v_k, length_tether, force_aero_sum),
        method="hybr",
    )

    if result.x > 0:  # if the root is positive, use it
        r_0 = result.x[0]
    else:  # if the root is negative, use the previous value
        r_0 = r_0_initial_guess
        print("!! WARNING !! ; a negative r_0 was found, previous value is used")

    # calculate turning radius of the whole kite r_k
    r_k = r_0 + center_of_gravity[1]

    # Calculate kite velocity of each node v_k,i
    r_ratio_array = np.array([((r_0 + point[1]) / r_k) for point in points])
    v_k_i_array = v_k * r_ratio_array

    # calculate centrifugal force distribution, i.e. the force on each node Fc,i
    force_centrifugal_side = [
        (mass_points[i] * (v_k_i_array[i] ** 2)) / (r_0 + point[1])
        for i, point in enumerate(points)
    ]
    force_centrifugal = np.array(
        [np.array([0, 1, 0]) * force for i, force in enumerate(force_centrifugal_side)]
    )

    # ## printing for debugging
    # if is_print_intermediate_results:
    #     print(
    #         f"r_0: {r_0:.2f}m (r_k: {r_k:.2f}m, cg_side: {center_of_gravity[1]:.2f}m)"
    #     )
    #     # Should be negative (oriented in negative-side/y direction)
    #     f_aero_side = force_aero_sum[1]
    #     f_tether_side = -force_aero_sum[2] * np.tan(np.arcsin(r_0 / length_tether))
    #     # Should be positive (oriented in positive-side/y direction)
    #     f_centrifugal_side = sum_of_masses * (v_k**2 / r_0)
    #     print(
    #         f"Fa_side: {f_aero_side:.3f}N, Ft_side: {f_tether_side:.3f}N, Fc_side: {f_centrifugal_side:.3f}N"
    #     )
    #     # print(
    #     #     f"SUM (Fc+Fa+Ft): {f_aero_side+f_tether_side+f_centrifugal_side:.3f}N (should be close to zero)"
    #     # )
    #     Ma_side = np.sum(np.cross(points, force_aero), axis=0)[0]
    #     Mc_side = np.sum(np.cross(points, force_centrifugal), axis=0)[0]
    #     print(
    #         f"Mx: {Ma_side+Mc_side:.1f}Nm (Ma_side: {Ma_side:.1f}Nm, Mc_side: {Mc_side:.1f}Nm)"
    #     )

    return force_centrifugal, r_0, v_k_i_array, center_of_gravity, r_k


def calculate_moment(distance_array, force_array):
    """calculates the moment given the distance and force arrays"""
    moment_array = [
        np.cross(point, force)
        for i, (point, force) in enumerate(zip(distance_array, force_array))
    ]
    moment_0, moment_1, moment_2 = 0, 0, 0
    for moment in moment_array:
        moment_0 += moment[0]
        moment_1 += moment[1]
        moment_2 += moment[2]
    return (
        np.array([moment_0, moment_1, moment_2]),
        np.sum(moment_array, axis=0),
        np.sum(np.cross(distance_array, force_array), axis=0),
    )


def calculate_vel_app_distribution(force_aero, mass_points, vel_app, points, vel_wind):
    """calculates the apparant wind speed at each node, given the velocity of the kite and the turning radius"""
    # calculate center of gravity
    center_of_gravity = calculate_center_of_gravity(points, mass_points)

    # initialize some variables
    force_aero_side = np.sum(force_aero[:, 1])
    sum_of_masses = np.sum(mass_points)
    v_k = vel_app[0]

    # calculate turning radius of the whole kite r_k
    r_k = (sum_of_masses * v_k**2) / -force_aero_side

    # calculate distance from turning-axis to x,y plane r_0
    r_0 = r_k - center_of_gravity[1]

    # Calculate kite velocity of each node v_k,i
    r_ratio_array = np.array([((r_0 + point[1]) / r_k) for point in points])
    v_k_i_array = v_k * r_ratio_array

    # calculate new apparant wind speed
    vel_app_array = [(vel_wind - v_k_i) for v_k_i in v_k_i_array]

    return vel_app_array, r_0, v_k_i_array, center_of_gravity
