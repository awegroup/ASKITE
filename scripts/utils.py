#%%
# making things autorelad - needed for Jupyter Kernel
%load_ext autoreload
%autoreload 2

# %% Checking the path
# checking the path and updating the path, accordingly (it was one folder to high)
import sys
import os

%matplotlib widget

# TODO: remove this hardcoding
folder_path = '/home/jellepoland/surfdrive/phd/code/phd_Jelle_Poland/Simulations'
os.chdir(folder_path)  # This should not be needed
sys.path.append(os.getcwd())

# %%
from src.coupling import coupling_aero2struc, coupling_struc2aero
from src.structural import structural_model, structural_mesher
from src.solver import solver_utils as solver_utils
from initialisation import load_surfplan, pulley_connectivity
from post_processing import functions_print, functions_plot
from post_processing import post_processing_utils as post_processing_utils
from test import test_main
from src import mass_distribution, actuation_relations
from src.solver import solver_main
from src.aerodynamic import VSM, breukels_2D, plate_aero, bridle_line_system_aero
import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import time
import scipy.optimize
import yaml
import importlib
import pytest
# {} gives units. {{}} is normal {} in Latex
from IPython.display import display, Latex

from src.initialisation.yaml_loader import config

#TODO:remove dicts
AERO_CONFIG = {'RING_GEOMETRY': config.aero.ring_geometry,
 'MODEL': config.aero.model,
 'N_ITER': config.aero.n_iter,
 'ERROR_LIMIT': config.aero.error_limit,
 'RELAX_FACTOR': config.aero.relax_factor,
 'N_SPLITS': config.aero.n_splits,
 'BRIDLE_CD': config.aero.cd_cylinder,
 'BRIDLE_CD_SHEAR': config.aero.cd_shear_cylinder,
 'CD_MULTIPLIER': config.aero.cd_multiplier,}

# Structural
SOLVER_CONFIG = {'METHOD': config.solver.method,
 'TOL': config.solver.tol,
 'MAXFEV': config.solver.max_fev,
 'INTEGRATION_METHOD': config.solver.integration_method,
 'RTOL_newton': config.solver.newton_rtol,
 'MAXITER': config.solver.newton_max_iter,
 'DISP': config.solver.newton_disp}

# AeroStructural  Solver
AEROSTRUCTURAL_SETTINGS = {'MAX_ITER': config.aero_structural.max_iter,
 'IT_CHECK': config.aero_structural.it_check,
 'TOL': config.aero_structural.tol,
 'MAX_ITER_VK': config.aero_structural.crosswind_max_iter,
 'TOL_VK': config.aero_structural.crosswind_tol,
 'RELAX_FACTOR_VK': config.aero_structural.crosswind_relax_factor}

# Stiffness Information
STIFFNESS_DATA = {'BRIDLE': config.kite.stiffness.bridle,
 'TUBE': config.kite.stiffness.tube,
 'TE': config.kite.stiffness.trailing_edge,
 'CANOPY': config.kite.stiffness.canopy,
 'IS_WITH_ELONGATION_LIMIT': config.kite.is_with_elongation_limit,
 'ELONGATION_LIMIT': config.kite.elongation_limit,
 'IS_WITH_COMPRESSION_LIMIT': config.kite.is_with_compression_limit,
 'COMPRESSION_LIMIT': config.kite.compression_limit,
 'LIMIT_STIFFNESS_FACTOR': config.kite.limit_stiffness_factor}

if config.kite_name == 'V3_25':  # TODO: Should be generated/imported from the surfplan file instead

    BILLOWING_ANGLES = config.kite.billowing_angles
    TUBE_DIAMETERS = config.kite.airfoil.tube_diameters
    CANOPY_MAX_HEIGHTS = config.kite.airfoil.canopy_max_heights

    BRIDLE_DATA = {'DIAMETER': config.kite.bridle.diameter, 
                   'RHO': config.kite.bridle.density, 
                   'BRIDLE_POINT_INDEX': config.kite.bridle.bridle_point_index}
    PULLEY_DATA =   {'POINT_INDICES': config.kite.pulley.point_indices,
                     'MASS': config.kite.pulley.mass,
                    'NUMBER_OF_PULLEYS_IN_BACK_LINES': config.kite.pulley.number_of_pulleys_in_back_lines}
    KCU_DATA = {'CD': config.kite.kcu.drag_coefficient,
                'DIAMETER': config.kite.kcu.diameter,
                'INDEX': config.kite.kcu.index,
                'MASS': config.kite.kcu.mass}

if config.kite_name == 'V9_60C':  # TODO: Should be generated/imported directly from the surfplan file instead
    surfplan_file = None
    simulation_directory = os.getcwd()

    extract_rib_db_whole = importlib.import_module(
        f'initialisation.{config.kite_name}.analyze_surfplan').extract_rib_db_whole
    RIB_DB_WHOLE_STRUTS = extract_rib_db_whole(
        simulation_directory, surfplan_file)
    extract_airfoil_geometry_data_from_ribs = importlib.import_module(
        f'initialisation.{config.kite_name}.analyze_surfplan').extract_airfoil_geometry_data_from_ribs
    TUBE_DIAMETERS, CANOPY_MAX_HEIGHTS, BILLOWING_ANGLES = extract_airfoil_geometry_data_from_ribs(
        RIB_DB_WHOLE_STRUTS, simulation_directory, config.kite_name)
    extract_bridle_line_system_data = importlib.import_module(
        f'initialisation.{config.kite_name}.analyze_surfplan').extract_bridle_line_system_data
    BRIDLE_DATA, PULLEY_DATA, KCU_DATA = extract_bridle_line_system_data(
        simulation_directory, surfplan_file)

#TODO: changing these caused convergence issues, be wary when getting to it.
## Adding the Connectivity to the Pulley_data dict
PULLEY_DATA = pulley_connectivity.extract_pulley_connectivity(
    config.kite.points_ini, config.kite.connectivity.bridle_ci, config.kite.connectivity.bridle_cj, PULLEY_DATA)

# Aerodynamic Model Setup
vel_app = config.vel_wind - config.vel_kite  # defining vel_app (vector of vel_app_norm)

# TODO: should be imported instead of defined here
airfoil_data = {'TUBE_DIAMETERS': TUBE_DIAMETERS,
                'IS_TUBE_DIAMETER_DIMENSIONLESS': config.kite.airfoil.is_tube_diameter_dimensionless,
                'CANOPY_MAX_HEIGHTS': CANOPY_MAX_HEIGHTS,
                'IS_MAX_CANOPY_HEIGHT_DIMENSIONLESS': config.kite.airfoil.is_canopy_max_height_dimensionless,}

# TODO: write such that also works for other kites
# Plotting geometry to check
if config.is_with_initial_plot and config.kite_name == 'V9_60C':
    plot_points_initial = [2, config.kite.points_ini, 'black', 2]
    plot_points_pulley = [
        2, config.kite.points_ini[PULLEY_DATA['POINT_INDICES']], 'red', 3]

    # acquiring the attachment points that lie in between
    points_attachment_inbetween = []
    points_index_attachment_inbetween = []
    for idx in config.kite.connectivity.wing_ci[config.kite.connectivity.tube_line_indices]:
        # if it is not a TE point
        if idx not in config.kite.connectivity.wing_ci[config.kite.connectivity.te_line_indices] \
                and idx not in config.kite.connectivity.wing_cj[config.kite.connectivity.te_line_indices] \
                and idx not in config.kite.connectivity.plate_point_indices:  # and not a plate point
            points_attachment_inbetween.append(config.kite.points_ini[idx])
            points_index_attachment_inbetween.append(idx)

    # need to sort these points per strut
    points_between_dict = {}
    tol_linepoint = .2
    for plate_index in config.kite.connectivity.plate_point_indices:
        left_le = config.kite.points_ini[plate_index[0]]
        left_te = config.kite.points_ini[plate_index[3]]
        for idx, point in enumerate(points_attachment_inbetween):
            if abs(structural_model.distance_point_to_line(point, [left_le, left_te])[0]) < tol_linepoint:
                points_between_dict[str(points_index_attachment_inbetween[idx])] = [
                    plate_index[0], plate_index[3]]

    # treating last strut separately
    right_le = config.kite.points_ini[config.kite.connectivity.plate_point_indices[-1][1]]
    right_te = config.kite.points_ini[config.kite.connectivity.plate_point_indices[-1][2]]
    for idx, point in enumerate(points_attachment_inbetween):
        if abs(structural_model.distance_point_to_line(point, [right_le, right_te])[0]) < tol_linepoint:
            points_between_dict[str(points_index_attachment_inbetween[idx])] = [
                config.kite.connectivity.plate_point_indices[-1][1], config.kite.connectivity.plate_point_indices[-1][2]]

    points_between_IN_dict = [
        config.kite.points_ini[int(key)] for key in points_between_dict.keys()]
    points_edges_dict_values_0 = [
        config.kite.points_ini[int(value[0])] for value in points_between_dict.values()]
    points_edges_dict_values_1 = [
        config.kite.points_ini[int(value[1])] for value in points_between_dict.values()]
    plot_point_edges_0 = [2, points_edges_dict_values_0, 'purple', 10]
    plot_point_edges_1 = [2, points_edges_dict_values_1, 'yellow', 10]
    plot_points_inbetween = [2, points_between_IN_dict, 'blue', 8]

    # plot_point_edges_0 = [2,points_edges_dict_values_0[0:4], 'purple',10]
    # plot_point_edges_1 = [2,points_edges_dict_values_1[0:4], 'yellow',10]
    # plot_points_inbetween =  [2,points_between_IN_dict[0:4], 'blue',8]

    plot_lines_pulley = [config.kite.points_ini, PULLEY_DATA['ci'],
                         PULLEY_DATA['cj'], 'red', 4, []]
    plot_lines_bridle = [config.kite.points_ini, config.kite.connectivity.bridle_ci, config.kite.connectivity.bridle_cj, 'black', 1.5, []]
    plot_lines_wing = [config.kite.points_ini, config.kite.connectivity.wing_ci, config.kite.connectivity.wing_cj, 'blue', 2, []]
    plot_lines_tube = [config.kite.points_ini, config.kite.connectivity.wing_ci[config.kite.connectivity.tube_line_indices],
                       config.kite.connectivity.wing_cj[config.kite.connectivity.tube_line_indices], 'grey', 10, []]
    functions_plot.plot_kite([plot_points_initial, plot_points_pulley, plot_points_inbetween, plot_point_edges_0, plot_point_edges_1], [
                             plot_lines_pulley, plot_lines_bridle, plot_lines_wing, plot_lines_tube])

# TODO: left-off here
# Strut bending resistance
# Realize that the initial meshing is not such that attachment points align
# Therefore we will displace the initial meshing; such that it aligns.
# Otherwise enforcing the bending resistance is hard

# but before doing so we will try with an initialisation
# does not seem to work, so we should try to get them to align.

#TODO: work on getting points_ini, differently.
points = config.kite.points_ini

## Alinging the points
points_between_dict_ini = {}
if config.kite_name == 'V9_60C':
    force_strut = np.zeros(points.shape)
    for key in points_between_dict.keys():
        point = points[int(key)]
        line_end1 = points[points_between_dict[key][0]]
        line_end2 = points[points_between_dict[key][1]]
        distance, distance_vector, angle = structural_model.distance_point_to_line(
            point, [line_end1, line_end2])
        points_between_dict_ini[key] = [
            points_between_dict[key][0], points_between_dict[key][1], angle, distance_vector]

        # putting points back
        points[int(key)] += distance_vector

# initialising solver
# Spring rest-lengths
bridle_rest_lengths = config.kite.bridle_rest_lengths_initial 
wing_rest_lengths = config.kite.wing_rest_lengths_initial

if config.kite_name == 'V3_25':  # TODO: Should be generated/imported from the surfplan file instead

    # ACTUATION
    depower_tape_extension = actuation_relations.up_to_ld(
        config.u_p, config.depower_tape_extension_percentage)  # [mm] depower-tape extension
    # bridle_rest_lengths[1] corresponds to the depower tape
    bridle_rest_lengths[1] += depower_tape_extension

# start with slightly shorter springs, for numerical reasons ##TODO: figure out if this is needed
bridle_rest_lengths = np.array(
    bridle_rest_lengths)*config.bridle_initial_compression_factor
wing_rest_lengths = np.array(wing_rest_lengths) * \
    config.bridle_initial_compression_factor

# config.kite.mass_points = mass_distribution.calculate_mass_distribution(
#     config.kite.points_ini, config.kite.connectivity.bridle_ci, config.kite.connectivity.bridle_cj, config.kite.connectivity.wing_ci, config.kite.connectivity.wing_cj, config.kite.wing_mass, BRIDLE_DATA, KCU_DATA, PULLEY_DATA)

# config.kite.mass_points = config.kite.config.kite.mass_points


if config.is_with_gravity:
    gravitational_constant = np.array([0, 0, -config.grav_constant])
else:
    gravitational_constant = np.array([0, 0, 0])
aerostructural_start_time = time.time()
calculate_force_spring_arguments = [bridle_rest_lengths, wing_rest_lengths, config.kite.connectivity.bridle_ci,
                                    config.kite.connectivity.bridle_cj, config.kite.connectivity.wing_ci, config.kite.connectivity.wing_cj, config.kite.connectivity.te_line_indices, config.kite.connectivity.tube_line_indices, STIFFNESS_DATA, PULLEY_DATA]
solver_arguments = [calculate_force_spring_arguments,
                    gravitational_constant, config.acc_kite, config.kite.mass_points, BRIDLE_DATA]
args_bridle = [config.is_with_aero_bridle, config.kite.connectivity.bridle_ci, config.kite.connectivity.bridle_cj, config.rho,
               BRIDLE_DATA, KCU_DATA]  # needed for optimilisation of Vk

# #%% STATIC AERO
# static aerodynamic analysis
points_left_to_right = coupling_struc2aero.order_struc_nodes_right_to_left(
    config.kite.points_ini, config.kite.connectivity.plate_point_indices)
force_aero_wing_VSM, moment_aero_wing_VSM, F_rel, ringvec, controlpoints, wingpanels, rings, coord_L, coord_refined = VSM.calculate_force_aero_wing_VSM(
    points_left_to_right, BILLOWING_ANGLES, vel_app, AERO_CONFIG, airfoil_data, config.kite.n_segments)
force_aero_wing = coupling_aero2struc.aero2struc(config.kite.points_ini, config.kite.connectivity.wing_ci, config.kite.connectivity.wing_cj, config.kite.connectivity.plate_point_indices,
                                                 force_aero_wing_VSM, moment_aero_wing_VSM, ringvec, controlpoints)  # Put in the new positions of the points

# Bridle Aerodynamics
if config.is_with_aero_bridle:
    force_aero_bridle = bridle_line_system_aero.calculate_force_aero_bridle_thedens2022(
        config.vel_wind, config.kite.connectivity.bridle_ci, config.kite.connectivity.bridle_cj, config.kite.points_ini, vel_app, config.rho, AERO_CONFIG, BRIDLE_DATA, KCU_DATA)
else:
    force_aero_bridle = [0]
force_aero = force_aero_wing+force_aero_bridle
force_total = solver_main.calculate_total_force(
    config.kite.points_ini, solver_arguments, force_aero)
print(
    f'force total: Fx: {np.sum(force_total[:,0]):.2f}N, Fy: {np.sum(force_total[:,1]):.2f}N, Fz: {np.sum(force_total[:,2]):.2f}N')
print(f' ')
functions_print.print_aero(config.kite.points_ini, config.kite.connectivity.plate_point_indices, config.kite.n_segments, vel_app, config.vel_wind,
                           force_aero, config.kite.area_projected, AERO_CONFIG, [0], KCU_DATA['INDEX'], config.rho, config.kite.ref_chord, config.mu)
# functions_plot.plot_aero(config.vel_wind, config.vel_kite, vel_app, wingpanels,
#                          controlpoints, rings, coord_L, config.kite.n_segments, F_rel, elev=10, azim=230)

#TODO: could write this into a test?
print(f' ')
print(f'-------- SHOULD BE ----------------')
print(f'force total: Fx: 188.36N, Fy: 0.00N, Fz: 4148.40N')
print(f'Aero Force  --> Fx 204.93 N, Fy 0.00 N, Fz 4152.67 N')
print(f'Angles      --> AoA: 8.70deg , trim-angle: -2.40deg')


# TODO: add back test
# TODO: write this as one liner, should take in all items and ONLY print those results that are none zero.
# if actually all zero, should return test is passed, kite is symmetrical and no forces are present
# main_test.print_test_results(config.kite.points_ini,vel_app,force_total,force_aero,force_aero_wing,calculate_force_spring_arguments,AERO_CONFIG,config.bridle_initial_compression_factor,
#                        config.kite.connectivity.wing_ci,config.kite.connectivity.wing_cj,config.kite.connectivity.bridle_ci,config.kite.connectivity.bridle_cj,wing_rest_lengths,bridle_rest_lengths,config.kite.connectivity.tube_line_indices,config.kite.area_projected)

# Plate-aero static aero
# force_aero_wing = plate_aero.calculate_force_aero_plate(config.kite.connectivity.plate_point_indices,config.kite.points_ini,vel_app_norm,config.kite.area_projected,config.rho,equal_boolean=False)
# force_aero_bridle = bridle_line_system_aero.calculate_force_aero_bridle_thedens2022(config.vel_wind,config.kite.connectivity.bridle_ci,config.kite.connectivity.bridle_cj,config.kite.points_ini,vel_app,config.rho,AERO_CONFIG,BRIDLE_DATA,KCU_DATA)
# force_aero = force_aero_wing + force_aero_bridle
# print(f' ')
# print(f'    plate-aero')
# coord_refined, controlpoints, rings, wingpanels, ringvec, coord_L, n_chordwise_elements_refined, config.kite.n_segments_refined = coupling_struc2aero.struc2aero(
# config.kite.points_ini,vel_app,n_chordwise_elements,config.kite.connectivity.plate_point_indices,BILLOWING_ANGLES,AERO_CONFIG)
# functions_print.print_aero(coord_refined,config.vel_wind,vel_app,vel_app_norm,force_aero,force_aero_bridle,
#    force_aero_wing,AERO_CONFIG,config.kite.ref_chord,MU,config.rho,config.kite.area_projected,KCU_DATA['INDEX'],
#    config.kite.points_ini[config.kite.connectivity.plate_point_indices[int(config.kite.n_segments/2)]][:2])


# %% AEROSTRUCTURAL SOLVER --- non-maneuvering, stationairy near zenith (symmetric) aerostructural solver
# AEROSTRUCTURAL SOLVER --- non-maneuvering, stationairy near zenith (symmetric) aerostructural solver

# initialisation
points = config.kite.points_ini
functions_print.print_settings(
    config.u_p, vel_app, AEROSTRUCTURAL_SETTINGS, AERO_CONFIG, config.kite.mass_points, SOLVER_CONFIG)
it_check = 1  # number of iterations between printing out results

for it in range(AEROSTRUCTURAL_SETTINGS['MAX_ITER']+1):

    time_iter = time.time()  # to print out each iteration time
    # Wing Aerodynamics
    points_left_to_right = coupling_struc2aero.order_struc_nodes_right_to_left(
        points, config.kite.connectivity.plate_point_indices)
    force_aero_wing_VSM, moment_aero_wing_VSM, F_rel, ringvec, controlpoints, wingpanels, rings, coord_L, coord_refined = VSM.calculate_force_aero_wing_VSM(
        points_left_to_right, BILLOWING_ANGLES, vel_app, AERO_CONFIG, airfoil_data, config.kite.n_segments)
    force_aero_wing = coupling_aero2struc.aero2struc(points, config.kite.connectivity.wing_ci, config.kite.connectivity.wing_cj, config.kite.connectivity.plate_point_indices, force_aero_wing_VSM,
                                                     moment_aero_wing_VSM, ringvec, controlpoints)  # Put in the new positions of the points
    # Bridle Aerodynamics
    if config.is_with_aero_bridle:
        force_aero_bridle = bridle_line_system_aero.calculate_force_aero_bridle_thedens2022(
            config.vel_wind, config.kite.connectivity.bridle_ci, config.kite.connectivity.bridle_cj, points, vel_app, config.rho, AERO_CONFIG, BRIDLE_DATA, KCU_DATA)
    else:
        force_aero_bridle = [0]

    force_aero = force_aero_wing+force_aero_bridle
    time_iter_aero = time.time()

    # Strut bending resistance
    if config.kite_name == 'V9_60C':
        force_strut = np.zeros(config.kite.points_ini.shape)
        for key in points_between_dict.keys():
            point = points[int(key)]
            line_end1 = points[points_between_dict[key][0]]
            line_end2 = points[points_between_dict[key][1]]
            distance, distance_vector, angle = structural_model.distance_point_to_line(
                point, [line_end1, line_end2])
            # adding a force from point to line
            # direction is verified and correct
            # adding force to force_aero, as an external force
            force_strut[int(key)] += STIFFNESS_DATA['STRUT']*distance_vector

        # adding strut force to aero force
        force_aero += force_strut

    # solver (here the structural solver is called)
    points_new, solver_max_error = solver_main.calculate_points_new(
        points, solver_arguments, force_aero, SOLVER_CONFIG)

    # Convergence
    # iter_difference with last iteration
    iter_diff = np.amax(np.linalg.norm(points_new-points, axis=1))

    # if converged
    if iter_diff < AEROSTRUCTURAL_SETTINGS['TOL'] and solver_max_error < SOLVER_CONFIG['TOL']:
        is_convergence = True
        break

    # if reached max-iterations
    elif it == AEROSTRUCTURAL_SETTINGS['MAX_ITER']:
        is_convergence = False
        break

    elif iter_diff > 2:  # when changes are bigger than 2meters
        print(f'iter_diff > 10, change is to large to be reasonable, stopping iterations')
        is_convergence = False
        break

    else:  # if not convergenced
        points = points_new
        if config.is_print_mid_results:
            if it % AEROSTRUCTURAL_SETTINGS['IT_CHECK'] == 0:
                trim_angle = functions_print.calculate_trim_angle(
                    coord_refined, AERO_CONFIG['N_SPLITS'], points_new[config.kite.connectivity.plate_point_indices[int(config.kite.n_segments/2)]][:2])
                print(f'i{it} - diff: {1e3*iter_diff:.3f} mm - error: {solver_max_error:.3f} - trim-angle: {trim_angle:.3f} deg - dt {(time.time() - time_iter):.2f} s (aero: {(time_iter_aero - time_iter):.2f} - struc: {(time.time() - time_iter_aero):.2f})')

# TODO: come up with a better name
np.save(
    f"results/V3_25/points_up_{str(int(100*config.u_p))}_symmetric.npy", points_new)

# post-processing
aerostructural_total_time = time.time()-aerostructural_start_time

# Computing trim-angle
if AERO_CONFIG['MODEL'] == ('VSM' or 'LLT'):
    trim_angle = trim_angle
elif AERO_CONFIG['MODEL'] == 'plate_aero':
    trim_angle = 10.  # TODO: write a function for this as well

# Printing out the solution nicely
functions_print.print_aerostructural(
    is_convergence, it, AEROSTRUCTURAL_SETTINGS['MAX_ITER'], aerostructural_total_time)
functions_print.print_structural(points_new, config.kite.points_ini, config.kite.connectivity.bridle_ci, config.kite.connectivity.bridle_cj, config.kite.connectivity.wing_ci,
                                 config.kite.connectivity.wing_cj, bridle_rest_lengths, wing_rest_lengths, config.u_p, STIFFNESS_DATA, config.kite.connectivity.tube_line_indices)

# computing the aero-again as its somehow not correctly stored ##TODO: fix-this
# Wing Aerodynamics
points_left_to_right = coupling_struc2aero.order_struc_nodes_right_to_left(
    points_new, config.kite.connectivity.plate_point_indices)
force_aero_wing_VSM, moment_aero_wing_VSM, F_rel, ringvec, controlpoints, wingpanels, rings, coord_L, coord_refined = VSM.calculate_force_aero_wing_VSM(
    points_left_to_right, BILLOWING_ANGLES, vel_app, AERO_CONFIG, airfoil_data, config.kite.n_segments)
force_aero_wing = coupling_aero2struc.aero2struc(points_new, config.kite.connectivity.wing_ci, config.kite.connectivity.wing_cj, config.kite.connectivity.plate_point_indices,
                                                 force_aero_wing_VSM, moment_aero_wing_VSM, ringvec, controlpoints)  # Put in the new positions of the points

# Bridle Aerodynamics
if config.is_with_aero_bridle:
    force_aero_bridle = bridle_line_system_aero.calculate_force_aero_bridle_thedens2022(
        config.vel_wind, config.kite.connectivity.bridle_ci, config.kite.connectivity.bridle_cj, points_new, vel_app, config.rho, AERO_CONFIG, BRIDLE_DATA, KCU_DATA)
else:
    force_aero_bridle = [0]
force_aero = force_aero_wing+force_aero_bridle
# Total force
force_total = solver_main.calculate_total_force(
    config.kite.points_ini, solver_arguments, force_aero)
print(
    f'force_total: Fx: {np.sum(force_total[:,0]):.2f}N, Fy: {np.sum(force_total[:,1]):.2f}N, Fz: {np.sum(force_total[:,2]):.2f}N')
# Spring force
force_spring = structural_model.calculate_force_spring(
    config.kite.points_ini, calculate_force_spring_arguments)
print(
    f'force_spring: Fx: {np.sum(force_spring[:,0]):.5f}N, Fy: {np.sum(force_spring[:,1]):.5f}N, Fz: {np.sum(force_spring[:,2]):.5f}N')
print(f' ---------------------------')
functions_print.print_aero(points_new, config.kite.connectivity.plate_point_indices, config.kite.n_segments, vel_app, config.vel_wind, force_aero,
                           config.kite.area_projected, AERO_CONFIG, force_aero_bridle, KCU_DATA['INDEX'], config.rho, config.kite.ref_chord, config.mu)
functions_plot.plot_aero(config.vel_wind, config.vel_kite, vel_app, wingpanels,
                         controlpoints, rings, coord_L, config.kite.n_segments, F_rel, elev=10, azim=230)
# functions_plot.plot_aero(config.vel_wind, config.vel_kite, vel_app, wingpanels, controlpoints, rings, coord_L, config.kite.n_segments, F_rel, elev=0, azim=-90)

# Plotting the solution - offline PLOTLY
plot_points_old = [2, config.kite.points_ini, 'black', 2]
plot_points_new = [2, points_new, 'green', 4]
elongation_values = post_processing_utils.calculate_elongation(
    points_new, config.kite.connectivity.wing_ci, config.kite.connectivity.wing_cj, config.kite.connectivity.bridle_ci, config.kite.connectivity.bridle_cj, wing_rest_lengths, bridle_rest_lengths, config.kite.connectivity.tube_line_indices)[2]
plot_lines_old = [config.kite.points_ini, np.concatenate(
    (config.kite.connectivity.bridle_ci, config.kite.connectivity.wing_ci)), np.concatenate((config.kite.connectivity.bridle_cj, config.kite.connectivity.wing_cj)), 'grey', 1.5, []]
plot_lines_new = [points_new, np.concatenate((config.kite.connectivity.bridle_ci, config.kite.connectivity.wing_ci)), np.concatenate(
    (config.kite.connectivity.bridle_cj, config.kite.connectivity.wing_cj)), 'purple', 4, elongation_values]
plot_surface_old = [True, config.kite.points_ini, config.kite.connectivity.plate_point_indices, 'purple', 0.0]
plot_surface_new = [True, points_new, config.kite.connectivity.plate_point_indices, 'lightgrey', 0.3]
functions_plot.plot_kite([plot_points_new, plot_points_old], [
                         plot_lines_new, plot_lines_old], f" ", [plot_surface_new, plot_surface_old])
# functions_plot.plot_kite([plot_points_new],[plot_lines_new],f"Light-grey: surfplan, Purple: simulation",[plot_surface_new])


#
#
#
#
#
#
#
#
#
#
#
#
#
#
#



# %% PERFECT CROSSWIND - AEROSTRUCTURAL SOLVER
# Loyd's perfect crosswind (symmetric) aerostructural solver

# PLOTTING fx vs vk_x
# Plotting the output of this function to see why the optimisation struggles...
# Define arrays to store data for plotting

# TODO: delete below, is now called from solver_utils
# def calculate_fx(vel_app, points, args):

#     # unpacking args
#     config.vel_wind, BILLOWING_ANGLES, AERO_CONFIG, airfoil_data, config.kite.n_segments, config.kite.connectivity.wing_ci, config.kite.connectivity.wing_cj, config.kite.connectivity.plate_point_indices = args

#     # Wing Aerodynamics
#     points_left_to_right = coupling_struc2aero.order_struc_nodes_right_to_left(
#         config.kite.points_ini, config.kite.connectivity.plate_point_indices)
#     force_aero_wing_VSM, moment_aero_wing_VSM, F_rel, ringvec, controlpoints, wingpanels, rings, coord_L, coord_refined = VSM.calculate_force_aero_wing_VSM(
#         points_left_to_right, BILLOWING_ANGLES, vel_app, AERO_CONFIG, airfoil_data, config.kite.n_segments)
#     force_aero_wing = coupling_aero2struc.aero2struc(config.kite.points_ini, config.kite.connectivity.wing_ci, config.kite.connectivity.wing_cj, config.kite.connectivity.plate_point_indices,
#                                                      force_aero_wing_VSM, moment_aero_wing_VSM, ringvec, controlpoints)  # Put in the new positions of the points

#     # Bridle Aerodynamics
#     if config.is_with_aero_bridle:
#         force_aero_bridle = bridle_line_system_aero.calculate_force_aero_bridle_thedens2022(
#             config.vel_wind, config.kite.connectivity.bridle_ci, config.kite.connectivity.bridle_cj, config.kite.points_ini, vel_app, config.rho, AERO_CONFIG, BRIDLE_DATA, KCU_DATA)
#     else:
#         force_aero_bridle = [0]
#     force_aero = force_aero_wing+force_aero_bridle
#     return np.sum(force_aero[:, 0])

vel_kite_x_values = []
fx_values = []
aoa_values = []

# getting the angle of mid-chord
points_LE_segment = points[config.kite.connectivity.plate_point_indices[int(config.kite.n_segments/2)]][:2]
points_TE_segment = points[config.kite.connectivity.plate_point_indices[int(config.kite.n_segments/2)]][2:]
LE_point_diff = points_LE_segment[1]-points_LE_segment[0]
TE_point_diff = points_TE_segment[1]-points_TE_segment[0]
mid_point_LE = points_LE_segment[0] + 0.5*LE_point_diff
mid_point_TE = points_TE_segment[0] + 0.5*TE_point_diff
trim_angle = np.arctan(mid_point_LE[0]/mid_point_LE[2])
delta_z_chord = mid_point_LE[2]-mid_point_TE[2]
delta_x_chord = np.abs(mid_point_LE[0]-mid_point_TE[0])
angle_of_mid_chord = np.rad2deg(np.arctan(delta_z_chord/delta_x_chord))

# TODO: make this a function in a separate file
for vel_kite_x_old in np.arange(15, 40, 1):
    vel_app[0] = config.vel_wind[0] + vel_kite_x_old
    args = (config.vel_wind, BILLOWING_ANGLES, AERO_CONFIG, airfoil_data,
            config.kite.n_segments, config.kite.connectivity.wing_ci, config.kite.connectivity.wing_cj, config.kite.connectivity.plate_point_indices)
    fx_sum = solver_utils.calculate_fx(-vel_kite_x_old,
                                       vel_app, points, args, args_bridle)

    # Append data to arrays for plotting
    vel_kite_x_values.append(-vel_kite_x_old)
    fx_values.append(fx_sum)

    # Also getting the aoa for plotting
    angle_va_with_xy = np.rad2deg(np.arctan(vel_app[2]/vel_app[0]))
    aoa = angle_of_mid_chord+angle_va_with_xy
    aoa_values.append(aoa)

    # print(f'vel_kite_x: {-vel_kite_x_old}, fx: {fx_sum}')

fig, (ax1, ax2, ax3) = plt.subplots(
    3, 1, figsize=(8, 15))  # Adjust the figsize as needed

# Plot the first data in the first subplot
ax1.plot(fx_values, vel_kite_x_values, marker='o', linestyle='--')
ax1.set_title('vel_kite_x vs fx')
ax1.set_xlabel('fx')
ax1.set_ylabel('vel_kite_x')
ax1.grid(True)
# ax1.set_xlim(-40, -10)

# create the second y-axis
ax1_twin = ax1.twinx()

# Plot the second data on the right y-axis
ax1_twin.plot(fx_values, aoa_values, marker='s',
              linestyle='--', color='r', label='aoa')
ax1_twin.set_ylabel('aoa', color='r')

# Add legend
lines1, labels1 = ax1.get_legend_handles_labels()
lines2, labels2 = ax1_twin.get_legend_handles_labels()
ax1.legend(lines1 + lines2, labels1 + labels2, loc='upper left')

# Plot the second data in the second subplot
ax2.plot(vel_kite_x_values, fx_values, marker='o', color='b',
         linestyle='--', label='aoa')  # Modify these with your data
ax2.set_title('vel_kite_x vs fx')
ax2.set_xlabel('vel_kite_x')
ax2.set_ylabel('fx')
ax2.grid(True)
# ax2.set_xlim(-40, -5)  # Modify the x-axis limits as needed

ax3.plot(aoa_values, fx_values, marker='o', color='r',
         linestyle='--', label='aoa')  # Modify these with your data
ax3.set_title('aoa vs fx')
ax3.set_xlabel('aoa')
ax3.set_ylabel('fx')
ax3.grid(True)

plt.tight_layout()  # Ensures proper spacing between subplots
plt.show()

# optimisation function
# special initialisation for perfect crosswind flight case

tol_vel_kite_optimization = 1e-1

# TODO: think about how I can further tune the Vk, such that the force_aero_x is zero
# looks like it doesn.t deal with after the point precision, should change that!

# from scipy.optimize import minimize

# %%

# initial guess
vk_x_initial_guess = -3*config.vel_wind[2]
# defining args
arguments = (config.vel_wind, BILLOWING_ANGLES, AERO_CONFIG, airfoil_data,
             config.kite.n_segments, config.kite.connectivity.wing_ci, config.kite.connectivity.wing_cj, config.kite.connectivity.plate_point_indices)

pre_simulation_of_vk = solver_utils.optimalisation_of_vk_for_fx_0(
    vk_x_initial_guess, vel_app, points, args, args_bridle)

# %% AEROSTRUCTURAL SOLVER --- PERFECT CROSSWIND
# AEROSTRUCTURAL SOLVER --- PERFECT CROSSWIND

# initialisation
vk_x_initial_guess = pre_simulation_of_vk
points = points_new  # use previous solution as initialisation
functions_print.print_settings(
    config.u_p, vel_app, AEROSTRUCTURAL_SETTINGS, AERO_CONFIG, config.kite.mass_points, SOLVER_CONFIG)
it_check = 1  # number of iterations between printing out results
# it_update_vk = 2 #numer of iterations at which Vk needs updating
tol_fx_ratio_to_fz = 0.005  # tolerance for fx, if fx>tol*fz, then update vk_x

is_convergence = False
for it in range(AEROSTRUCTURAL_SETTINGS['MAX_ITER']+1):

    # Updating the velocity of the kite in the x-direction, to make Fx go to zero
    # run if fx > tol_fx_ratio_to_fz*Fz
    force_resultant_x = solver_utils.calculate_fx(-vel_app[0], vel_app, points, args, args_bridle)
    if np.abs(force_resultant_x) > tol_fx_ratio_to_fz*np.abs(np.sum(force_aero[:, 2])):
        print(f'i{it}, updating vel_kite_x')
        # vel_app[0] = -finding_vk_which_makes_fx_zero(initial_guess_vk_x,points,config.vel_wind,BILLOWING_ANGLES, AERO_CONFIG, airfoil_data, config.kite.n_segments, config.kite.connectivity.wing_ci, config.kite.connectivity.wing_cj, config.kite.connectivity.plate_point_indices)
        args = (config.vel_wind, BILLOWING_ANGLES, AERO_CONFIG, airfoil_data,
                config.kite.n_segments, config.kite.connectivity.wing_ci, config.kite.connectivity.wing_cj, config.kite.connectivity.plate_point_indices)
        vel_app[0] = - \
            solver_utils.optimalisation_of_vk_for_fx_0(
                -vel_app[0], vel_app, points, args, args_bridle)

    time_iter = time.time()  # to print out each iteration time
    # Wing Aerodynamics
    points_left_to_right = coupling_struc2aero.order_struc_nodes_right_to_left(
        points, config.kite.connectivity.plate_point_indices)
    force_aero_wing_VSM, moment_aero_wing_VSM, F_rel, ringvec, controlpoints, wingpanels, rings, coord_L, coord_refined = VSM.calculate_force_aero_wing_VSM(
        points_left_to_right, BILLOWING_ANGLES, vel_app, AERO_CONFIG, airfoil_data, config.kite.n_segments)
    force_aero_wing = coupling_aero2struc.aero2struc(points, config.kite.connectivity.wing_ci, config.kite.connectivity.wing_cj, config.kite.connectivity.plate_point_indices, force_aero_wing_VSM,
                                                     moment_aero_wing_VSM, ringvec, controlpoints)  # Put in the new positions of the points
    # Bridle Aerodynamics
    if config.is_with_aero_bridle:
        force_aero_bridle = bridle_line_system_aero.calculate_force_aero_bridle_thedens2022(
            config.vel_wind, config.kite.connectivity.bridle_ci, config.kite.connectivity.bridle_cj, points, vel_app, config.rho, AERO_CONFIG, BRIDLE_DATA, KCU_DATA)
    else:
        force_aero_bridle = [0]

    force_aero = force_aero_wing+force_aero_bridle
    time_iter_aero = time.time()

    if config.kite_name == 'V9_60C':
        force_strut = np.zeros(config.kite.points_ini.shape)
        for key in points_between_dict.keys():
            point = points[int(key)]
            line_end1 = points[points_between_dict[key][0]]
            line_end2 = points[points_between_dict[key][1]]
            distance, distance_vector, angle = structural_model.distance_point_to_line(
                point, [line_end1, line_end2])
            # adding force to force_aero, as an external force
            force_strut[int(key)] += STIFFNESS_DATA['STRUT']*distance_vector
        # adding strut force to aero force
        force_aero += force_strut

    # solver (here the structural solver is called)
    points_new, solver_max_error = solver_main.calculate_points_new(
        points, solver_arguments, force_aero, SOLVER_CONFIG)

    # Convergence
    # iter_difference with last iteration
    iter_diff = np.amax(np.linalg.norm(points_new-points, axis=1))

    # if converged
    if iter_diff < AEROSTRUCTURAL_SETTINGS['TOL'] and solver_max_error < SOLVER_CONFIG['TOL']:
        print(f'Converged at: i{it}, final update of vel_kite_x')
        # vel_app[0] = -finding_vk_which_makes_fx_zero(initial_guess_vk_x,points,config.vel_wind,BILLOWING_ANGLES, AERO_CONFIG, airfoil_data, config.kite.n_segments, config.kite.connectivity.wing_ci, config.kite.connectivity.wing_cj, config.kite.connectivity.plate_point_indices)
        args = (config.vel_wind, BILLOWING_ANGLES, AERO_CONFIG, airfoil_data,
                config.kite.n_segments, config.kite.connectivity.wing_ci, config.kite.connectivity.wing_cj, config.kite.connectivity.plate_point_indices)
        vel_app[0] = - \
            solver_utils.optimalisation_of_vk_for_fx_0(
                -vel_app[0], vel_app, points, args, args_bridle)
        is_convergence = True
        break

    # elif it == AEROSTRUCTURAL_SETTINGS['MAX_ITER']: #if reached max-iterations
    #     is_convergence = False
    #     break

    elif iter_diff > 2:  # when changes are bigger than 2meters
        print(f'iter_diff > 10, change is to large to be reasonable, stopping iterations')
        is_convergence = False
        break

    else:  # if not convergenced
        points = points_new
        if config.is_print_mid_results:
            if it % AEROSTRUCTURAL_SETTINGS['IT_CHECK'] == 0:
                trim_angle = functions_print.calculate_trim_angle(
                    coord_refined, AERO_CONFIG['N_SPLITS'], points_new[config.kite.connectivity.plate_point_indices[int(config.kite.n_segments/2)]][:2])
                print(
                    f'i{it} - diff: {1e3*iter_diff:.3f} mm - fx: {force_resultant_x:.3f}N - error: {solver_max_error:.3f} - trim-angle: {trim_angle:.3f} deg - dt {(time.time() - time_iter):.2f} s (aero: {(time_iter_aero - time_iter):.2f} - struc: {(time.time() - time_iter_aero):.2f})')


np.save(
    f"results/V3_25/points_up_{str(int(100*config.u_p))}_perfectcrosswind.npy", points_new)

# post-processing
aerostructural_total_time = time.time()-aerostructural_start_time

# Computing trim-angle
if AERO_CONFIG['MODEL'] == ('VSM' or 'LLT'):
    trim_angle = trim_angle
elif AERO_CONFIG['MODEL'] == 'plate_aero':
    trim_angle = 10.  # TODO: write a function for this as well

# Printing out the solution nicely
functions_print.print_aerostructural(
    is_convergence, it, AEROSTRUCTURAL_SETTINGS['MAX_ITER'], aerostructural_total_time)
functions_print.print_structural(points_new, config.kite.points_ini, config.kite.connectivity.bridle_ci, config.kite.connectivity.bridle_cj, config.kite.connectivity.wing_ci,
                                 config.kite.connectivity.wing_cj, bridle_rest_lengths, wing_rest_lengths, config.u_p, STIFFNESS_DATA, config.kite.connectivity.tube_line_indices)

# computing the aero-again as its somehow not correctly stored ##TODO: fix-this
# Wing Aerodynamics
points_left_to_right = coupling_struc2aero.order_struc_nodes_right_to_left(
    points_new, config.kite.connectivity.plate_point_indices)
force_aero_wing_VSM, moment_aero_wing_VSM, F_rel, ringvec, controlpoints, wingpanels, rings, coord_L, coord_refined = VSM.calculate_force_aero_wing_VSM(
    points_left_to_right, BILLOWING_ANGLES, vel_app, AERO_CONFIG, airfoil_data, config.kite.n_segments)
force_aero_wing = coupling_aero2struc.aero2struc(points_new, config.kite.connectivity.wing_ci, config.kite.connectivity.wing_cj, config.kite.connectivity.plate_point_indices,
                                                 force_aero_wing_VSM, moment_aero_wing_VSM, ringvec, controlpoints)  # Put in the new positions of the points

# Bridle Aerodynamics
if config.is_with_aero_bridle:
    force_aero_bridle = bridle_line_system_aero.calculate_force_aero_bridle_thedens2022(
        config.vel_wind, config.kite.connectivity.bridle_ci, config.kite.connectivity.bridle_cj, points_new, vel_app, config.rho, AERO_CONFIG, BRIDLE_DATA, KCU_DATA)
else:
    force_aero_bridle = [0]
force_aero = force_aero_wing+force_aero_bridle
# Total force
force_total = solver_main.calculate_total_force(
    config.kite.points_ini, solver_arguments, force_aero)
print(
    f'force_total: Fx: {np.sum(force_total[:,0]):.2f}N, Fy: {np.sum(force_total[:,1]):.2f}N, Fz: {np.sum(force_total[:,2]):.2f}N')
# Spring force
force_spring = structural_model.calculate_force_spring(
    config.kite.points_ini, calculate_force_spring_arguments)
print(
    f'force_spring: Fx: {np.sum(force_spring[:,0]):.5f}N, Fy: {np.sum(force_spring[:,1]):.5f}N, Fz: {np.sum(force_spring[:,2]):.5f}N')
print(f' ---------------------------')
functions_print.print_aero(points_new, config.kite.connectivity.plate_point_indices, config.kite.n_segments, vel_app, config.vel_wind, force_aero,
                           config.kite.area_projected, AERO_CONFIG, force_aero_bridle, KCU_DATA['INDEX'], config.rho, config.kite.ref_chord, config.mu)
# functions_plot.plot_aero(config.vel_wind, config.vel_kite, vel_app, wingpanels, controlpoints, rings, coord_L, config.kite.n_segments, F_rel, elev=10, azim=230)
# functions_plot.plot_aero(config.vel_wind, config.vel_kite, vel_app, wingpanels, controlpoints, rings, coord_L, config.kite.n_segments, F_rel, elev=0, azim=-90)


# Plotting the solution - offline PLOTLY
plot_points_old = [2, config.kite.points_ini, 'black', 2]
plot_points_new = [2, points_new, 'green', 4]
elongation_values = post_processing_utils.calculate_elongation(
    points_new, config.kite.connectivity.wing_ci, config.kite.connectivity.wing_cj, config.kite.connectivity.bridle_ci, config.kite.connectivity.bridle_cj, wing_rest_lengths, bridle_rest_lengths, config.kite.connectivity.tube_line_indices)[2]
plot_lines_old = [config.kite.points_ini, np.concatenate(
    (config.kite.connectivity.bridle_ci, config.kite.connectivity.wing_ci)), np.concatenate((config.kite.connectivity.bridle_cj, config.kite.connectivity.wing_cj)), 'grey', 1.5, []]
plot_lines_new = [points_new, np.concatenate((config.kite.connectivity.bridle_ci, config.kite.connectivity.wing_ci)), np.concatenate(
    (config.kite.connectivity.bridle_cj, config.kite.connectivity.wing_cj)), 'purple', 4, elongation_values]
plot_surface_old = [True, config.kite.points_ini, config.kite.connectivity.plate_point_indices, 'purple', 0.0]
plot_surface_new = [True, points_new, config.kite.connectivity.plate_point_indices, 'lightgrey', 0.3]
functions_plot.plot_kite([plot_points_new, plot_points_old], [
                         plot_lines_new, plot_lines_old], f" ", [plot_surface_new, plot_surface_old])
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#


# %% Initialising asymmetrical state
# Initialising asymmetrical state, applying delta_ls and stop iteration before kite starts to move out of domain


points_new = np.load(f"results/V3_25/points_up_70_perfectcrosswind.npy")
#points_new = config.kite.points_ini

# applying a small delta to the bridle rest lengths
bridle_rest_lengths = structural_mesher.calculate_edge_lengths(
    config.kite.connectivity.bridle_ci, config.kite.connectivity.bridle_cj, config.kite.points_ini)  # [m]
wing_rest_lengths = structural_mesher.calculate_edge_lengths(
    config.kite.connectivity.wing_ci, config.kite.connectivity.wing_cj, config.kite.points_ini)  # [m]
# start with slightly shorter springs, for numerical reasons ##TODO: figure out if this is needed
bridle_rest_lengths = np.array(
    bridle_rest_lengths)*config.bridle_initial_compression_factor
wing_rest_lengths = np.array(wing_rest_lengths) * \
    config.bridle_initial_compression_factor
delta_ls = 0.01
bridle_rest_lengths[2] += delta_ls  # spring[1] is the left steering tape
bridle_rest_lengths[3] -= delta_ls  # spring[2] is the right steering tape

# updating the args
new_stiffness = 3e5
STIFFNESS_DATA['BRIDLE'] = new_stiffness
STIFFNESS_DATA['TUBE'] = new_stiffness
STIFFNESS_DATA['TE'] = new_stiffness

calculate_force_spring_arguments = [bridle_rest_lengths, wing_rest_lengths, config.kite.connectivity.bridle_ci,
                                    config.kite.connectivity.bridle_cj, config.kite.connectivity.wing_ci, config.kite.connectivity.wing_cj, config.kite.connectivity.te_line_indices, config.kite.connectivity.tube_line_indices, STIFFNESS_DATA, PULLEY_DATA]
solver_arguments = [calculate_force_spring_arguments,
                    gravitational_constant, config.acc_kite, config.kite.mass_points, BRIDLE_DATA]

# running the solver
points = points_new  # use previous solution as initialisation
max_iter_turning_initialisation = 1
it_check = 1
it = 0
SOLVER_CONFIG['method'] = 'newton'

is_convergence = False
for it in range(max_iter_turning_initialisation):

    time_iter = time.time()  # to print out each iteration time

    # Wing Aerodynamics
    points_left_to_right = coupling_struc2aero.order_struc_nodes_right_to_left(
        points, config.kite.connectivity.plate_point_indices)
    force_aero_wing_VSM, moment_aero_wing_VSM, F_rel, ringvec, controlpoints, wingpanels, rings, coord_L, coord_refined = VSM.calculate_force_aero_wing_VSM(
        points_left_to_right, BILLOWING_ANGLES, vel_app, AERO_CONFIG, airfoil_data, config.kite.n_segments)
    force_aero_wing = coupling_aero2struc.aero2struc(points, config.kite.connectivity.wing_ci, config.kite.connectivity.wing_cj, config.kite.connectivity.plate_point_indices, force_aero_wing_VSM,
                                                     moment_aero_wing_VSM, ringvec, controlpoints)  # Put in the new positions of the points
    # Bridle Aerodynamics
    if config.is_with_aero_bridle:
        force_aero_bridle = bridle_line_system_aero.calculate_force_aero_bridle_thedens2022(
            config.vel_wind, config.kite.connectivity.bridle_ci, config.kite.connectivity.bridle_cj, points, vel_app, config.rho, AERO_CONFIG, BRIDLE_DATA, KCU_DATA)
    else:
        force_aero_bridle = [0]

    force_aero = force_aero_wing+force_aero_bridle
    time_iter_aero = time.time()

    # solver (here the structural solver is called)
    points_new, solver_max_error = solver_main.calculate_points_new(
        points, solver_arguments, force_aero, SOLVER_CONFIG)

    # Convergence
    # iter_difference with last iteration
    iter_diff = np.amax(np.linalg.norm(points_new-points, axis=1))
    sideangle_leading_edge_midpoint = functions_print.calculate_side_angle_leading_edge_midpoint(
        points_new[config.kite.connectivity.plate_point_indices[int(config.kite.n_segments/2)]][:2])

    if iter_diff < AEROSTRUCTURAL_SETTINGS['TOL']:  # if convergenced
        # print(f'Converged at: i{it}, final update of vel_kite_x')
        # # vel_app[0] = -finding_vk_which_makes_fx_zero(initial_guess_vk_x,points,config.vel_wind,BILLOWING_ANGLES, AERO_CONFIG, airfoil_data, config.kite.n_segments, config.kite.connectivity.wing_ci, config.kite.connectivity.wing_cj, config.kite.connectivity.plate_point_indices)
        # args = (config.vel_wind,BILLOWING_ANGLES, AERO_CONFIG, airfoil_data, config.kite.n_segments, config.kite.connectivity.wing_ci, config.kite.connectivity.wing_cj, config.kite.connectivity.plate_point_indices)
        # vel_app[0] = -solver_utils.optimalisation_of_vk_for_fx_0(vk_x_initial_guess,vel_app,points,args,args_bridle))
        is_convergence = True
        break

    elif iter_diff > 2:  # when changes are bigger than 2meters
        print(f'iter_diff > 10, change is to large to be reasonable, stopping iterations')
        is_convergence = False
        break

    else:  # if not convergenced
        points = points_new
        if config.is_print_mid_results:
            if it % AEROSTRUCTURAL_SETTINGS['IT_CHECK'] == 0:
                trim_angle = functions_print.calculate_trim_angle(
                    coord_refined, AERO_CONFIG['N_SPLITS'], points_new[config.kite.connectivity.plate_point_indices[int(config.kite.n_segments/2)]][:2])
                print(f'i{it} - diff: {1e3*iter_diff:.3f}mm - error: {solver_max_error:.3f}N? - side: {sideangle_leading_edge_midpoint:.2f}deg - trim: {trim_angle:.3f}deg - dt: {(time.time() - time_iter):.2f}s (aero: {(time_iter_aero - time_iter):.2f} - struc: {(time.time() - time_iter_aero):.2f})')

# post-processing
aerostructural_total_time = time.time()-aerostructural_start_time

# Computing trim-angle
if AERO_CONFIG['MODEL'] == ('VSM' or 'LLT'):
    trim_angle = trim_angle
elif AERO_CONFIG['MODEL'] == 'plate_aero':
    trim_angle = 10.  # TODO: write a function for this as well

# Printing out the solution nicely
functions_print.print_aerostructural(
    is_convergence, it, AEROSTRUCTURAL_SETTINGS['MAX_ITER'], aerostructural_total_time)
functions_print.print_structural(points_new, config.kite.points_ini, config.kite.connectivity.bridle_ci, config.kite.connectivity.bridle_cj, config.kite.connectivity.wing_ci,
                                 config.kite.connectivity.wing_cj, bridle_rest_lengths, wing_rest_lengths, config.u_p, STIFFNESS_DATA, config.kite.connectivity.tube_line_indices)

# computing the aero-again as its somehow not correctly stored ##TODO: fix-this
# Wing Aerodynamics
points_left_to_right = coupling_struc2aero.order_struc_nodes_right_to_left(
    points_new, config.kite.connectivity.plate_point_indices)
force_aero_wing_VSM, moment_aero_wing_VSM, F_rel, ringvec, controlpoints, wingpanels, rings, coord_L, coord_refined = VSM.calculate_force_aero_wing_VSM(
    points_left_to_right, BILLOWING_ANGLES, vel_app, AERO_CONFIG, airfoil_data, config.kite.n_segments)
force_aero_wing = coupling_aero2struc.aero2struc(points_new, config.kite.connectivity.wing_ci, config.kite.connectivity.wing_cj, config.kite.connectivity.plate_point_indices,
                                                 force_aero_wing_VSM, moment_aero_wing_VSM, ringvec, controlpoints)  # Put in the new positions of the points

# Bridle Aerodynamics
if config.is_with_aero_bridle:
    force_aero_bridle = bridle_line_system_aero.calculate_force_aero_bridle_thedens2022(
        config.vel_wind, config.kite.connectivity.bridle_ci, config.kite.connectivity.bridle_cj, points_new, vel_app, config.rho, AERO_CONFIG, BRIDLE_DATA, KCU_DATA)
else:
    force_aero_bridle = [0]
force_aero = force_aero_wing+force_aero_bridle
# Total force
force_total = solver_main.calculate_total_force(
    config.kite.points_ini, solver_arguments, force_aero)
print(
    f'force_total: Fx: {np.sum(force_total[:,0]):.2f}N, Fy: {np.sum(force_total[:,1]):.2f}N, Fz: {np.sum(force_total[:,2]):.2f}N')
# Spring force
force_spring = structural_model.calculate_force_spring(
    config.kite.points_ini, calculate_force_spring_arguments)
print(
    f'force_spring: Fx: {np.sum(force_spring[:,0]):.5f}N, Fy: {np.sum(force_spring[:,1]):.5f}N, Fz: {np.sum(force_spring[:,2]):.5f}N')
print(f' ---------------------------')
functions_print.print_aero(points_new, config.kite.connectivity.plate_point_indices, config.kite.n_segments, vel_app, config.vel_wind, force_aero,
                           config.kite.area_projected, AERO_CONFIG, force_aero_bridle, KCU_DATA['INDEX'], config.rho, config.kite.ref_chord, config.mu)
# functions_plot.plot_aero(config.vel_wind, config.vel_kite, vel_app, wingpanels, controlpoints, rings, coord_L, config.kite.n_segments, F_rel, elev=10, azim=230)
# functions_plot.plot_aero(config.vel_wind, config.vel_kite, vel_app, wingpanels, controlpoints, rings, coord_L, config.kite.n_segments, F_rel, elev=0, azim=-90)

# Plotting the solution - offline PLOTLY
plot_points_old = [2, config.kite.points_ini, 'black', 2]
plot_points_new = [2, points_new, 'green', 4]
elongation_values = post_processing_utils.calculate_elongation(
    points_new, config.kite.connectivity.wing_ci, config.kite.connectivity.wing_cj, config.kite.connectivity.bridle_ci, config.kite.connectivity.bridle_cj, wing_rest_lengths, bridle_rest_lengths, config.kite.connectivity.tube_line_indices)[2]
plot_lines_old = [config.kite.points_ini, np.concatenate(
    (config.kite.connectivity.bridle_ci, config.kite.connectivity.wing_ci)), np.concatenate((config.kite.connectivity.bridle_cj, config.kite.connectivity.wing_cj)), 'grey', 1.5, []]
plot_lines_new = [points_new, np.concatenate((config.kite.connectivity.bridle_ci, config.kite.connectivity.wing_ci)), np.concatenate(
    (config.kite.connectivity.bridle_cj, config.kite.connectivity.wing_cj)), 'purple', 4, elongation_values]
plot_surface_old = [True, config.kite.points_ini, config.kite.connectivity.plate_point_indices, 'purple', 0.0]
plot_surface_new = [True, points_new, config.kite.connectivity.plate_point_indices, 'lightgrey', 0.3]
# functions_plot.plot_kite([plot_points_new,plot_points_old],[plot_lines_new,plot_lines_old],f" ",[plot_surface_new,plot_surface_old])
functions_plot.plot_kite(
    [plot_points_new], [plot_lines_new], f" ", [plot_surface_new])

# %%
# #%% Initializing Turning variables
# initializing turning variables


# updating vk
args = (config.vel_wind, BILLOWING_ANGLES, AERO_CONFIG, airfoil_data,
        config.kite.n_segments, config.kite.connectivity.wing_ci, config.kite.connectivity.wing_cj, config.kite.connectivity.plate_point_indices)
time_start_finding_vk = time.time()
vel_kite_x = solver_utils.optimalisation_of_vk_for_fx_0(
    -vel_app[0], vel_app, points, args, args_bridle)
print(f'time to find vk: {time.time()-time_start_finding_vk:.2f}s')

# defining several functions

def calculate_center_of_gravity(points, config.kite.mass_points):
    ''' calculates the center of gravity of the kite, given the points and the mass of each point'''
    numerator_x = np.sum([config.kite.mass_points[i]*point[0]
                         for i, point in enumerate(points)])
    numerator_y = np.sum([config.kite.mass_points[i]*point[1]
                         for i, point in enumerate(points)])
    numerator_z = np.sum([config.kite.mass_points[i]*point[2]
                         for i, point in enumerate(points)])
    denominator = np.sum(config.kite.mass_points)
    return np.array([numerator_x, numerator_y, numerator_z])/denominator


def calculate_force_centrifugal_distribution(force_aero, config.kite.mass_points, vel_app, points):
    ''' calculates the centrifugal force distribution, i.e. the force on each node Fc,i '''
    # calculate center of gravity
    center_of_gravity = calculate_center_of_gravity(points, config.kite.mass_points)

    # initialize some variables
    force_aero_side = np.sum(force_aero[:, 1])
    sum_of_masses = np.sum(config.kite.mass_points)
    v_k = vel_app[0]

    # calculate turning radius of the whole kite r_k
    r_k = (sum_of_masses*v_k**2) / -force_aero_side

    # calculate distance from turning-axis to x,y plane r_0
    r_0 = r_k - center_of_gravity[1]

    # Calculate kite velocity of each node v_k,i
    r_ratio_array = np.array([((r_0 + point[1])/r_k) for point in points])
    v_k_i_array = v_k * r_ratio_array

    # calculate centrifugal force distribution, i.e. the force on each node Fc,i
    force_centrifugal_i_array = [
        ((config.kite.mass_points[i] * (v_k_i_array[i]**2)) / (r_0 + point[1])) for i, point in enumerate(points)]

    return force_centrifugal_i_array, r_0, v_k_i_array, center_of_gravity

def calculate_force_centrifugal_distribution_moment_based(moment_aero, config.kite.mass_points, vel_app, points):
    ''' calculates the centrifugal force distribution, i.e. the force on each node Fc,i '''
    # calculate center of gravity
    center_of_gravity = calculate_center_of_gravity(points, config.kite.mass_points)

    # initialize some variables
    sum_of_masses = np.sum(config.kite.mass_points)
    v_k = vel_app[0]

    # calculat the roll moment
    moment_aero_roll = np.sum(moment_aero[:, 0])

    # calculate Fc that would create the same moment as the moment_aero
    force_centrifugal = moment_aero_roll / np.linalg.norm(center_of_gravity)

    # calculate turning radius of the whole kite r_k
    r_k = (sum_of_masses*v_k**2) / force_centrifugal

    # calculate distance from turning-axis to x,y plane r_0
    r_0 = r_k - center_of_gravity[1]

    # Calculate kite velocity of each node v_k,i
    r_ratio_array = np.array([((r_0 + point[1])/r_k) for point in points])
    v_k_i_array = v_k * r_ratio_array

    # calculate centrifugal force distribution, i.e. the force on each node Fc,i
    force_centrifugal_i_array = [
        ((config.kite.mass_points[i] * (v_k_i_array[i]**2)) / (r_0 + point[1])) for i, point in enumerate(points)]

    return force_centrifugal_i_array, r_0, v_k_i_array, center_of_gravity

def calculate_force_centrifugal_distribution_with_tether_tension(r_0_previous, force_aero, config.kite.mass_points, vel_app, points):
    ''' calculates the centrifugal force distribution, i.e. the force on each node Fc,i '''
    # calculate center of gravity
    center_of_gravity = calculate_center_of_gravity(points, config.kite.mass_points)

    # initialize some variables
    force_aero_sum = np.sum(force_aero,axis=0)
    sum_of_masses = np.sum(config.kite.mass_points)
    v_k = vel_app[0]
    length_tether = 200 # m

    ### CALCULATING r_0
    def equation_to_solve(r_0, sum_of_masses, v_k, length_tether, force_aero_sum):
        x_aero_induced_tether_component =  force_aero_sum[2] * np.tan(np.arcsin(r_0 / length_tether))
        z_aero_induced_component = force_aero_sum[1]
        inwards_directed_force = np.abs(x_aero_induced_tether_component) + np.abs(z_aero_induced_component)
        centrifugal_force = sum_of_masses * (v_k**2 / r_0)

        ## now r_0 will be chosen such that it balances these two contributions
        return centrifugal_force - inwards_directed_force

    # Find the root numerically and get full output
    r_0_initial_guess = r_0_previous
    result = scipy.optimize.root(equation_to_solve, r_0_initial_guess, args=(sum_of_masses, v_k, length_tether, force_aero_sum), method='hybr')

    if result.x > 0: # if the root is positive, use it
        r_0 = result.x[0]
    else: # if the root is negative, use the previous value
        r_0 = r_0_previous
        print('!! WARNING !! ; a negative r_0 was found, previous value is used')

    # calculate turning radius of the whole kite r_k
    r_k = r_0 + center_of_gravity[1]

    # Calculate kite velocity of each node v_k,i
    r_ratio_array = np.array([((r_0 + point[1])/r_k) for point in points])
    v_k_i_array = v_k * r_ratio_array

    # calculate centrifugal force distribution, i.e. the force on each node Fc,i
    force_centrifugal_i_array = [
        ((config.kite.mass_points[i] * (v_k_i_array[i]**2)) / (r_0 + point[1])) for i, point in enumerate(points)]

    return force_centrifugal_i_array, r_0, v_k_i_array, center_of_gravity

def calculate_moment(distance_array,force_array):
    ''' calculates the moment given the distance and force arrays'''
    moment_array = [np.cross(point,force) for i, (point,force) in enumerate(zip(distance_array,force_array))]
    moment_0,moment_1,moment_2 = 0,0,0
    for moment in moment_array:
        moment_0 += moment[0]
        moment_1 += moment[1]
        moment_2 += moment[2]
    return np.array([moment_0,moment_1,moment_2]), np.sum(moment_array,axis=0), np.sum(np.cross(distance_array,force_array),axis=0)

## Printing out results
force_aero_sum = np.sum(force_aero, axis=0)
print(f'sum(Fa):{force_aero_sum}N')
moment_aero = np.cross(points, force_aero)
moment_aero_sum = np.sum(moment_aero, axis=0)
center_of_fa =  calculate_center_of_gravity(points, force_aero)
print(f'Force center: {center_of_fa}m')
print(f'Ma at force-center:{np.cross(center_of_fa,force_aero_sum)[0]:.0f}Nm')

roll_moment = np.sum(np.cross(points, force_aero), axis=0)[0]  # Calculate the roll moment
print(f'Roll Moment: {roll_moment:.0f} Nm')
print(f'Ma_roll due to to Fa_side: {force_aero_sum[1]*center_of_fa[2]:.0f}Nm')
print(f'Ma[0] {np.sum(moment_aero[:,0]):.0f}Nm, Ma[1] {np.sum(moment_aero[:,1]):.0f}Nm, Ma[2] {np.sum(moment_aero[:,2]):.0f}Nm')


print(f'Ma individual: {calculate_moment(points,force_aero)}')

print(f'')
print(f'Equal forces')
force_centrifugal_i_array, r_0, v_k_i_array, center_of_gravity = calculate_force_centrifugal_distribution(
    force_aero, config.kite.mass_points, vel_app, points)
force_centrifugal = np.array(
    [np.array([0, 1, 0])*force for i, force in enumerate(force_centrifugal_i_array)])
moment_centrifugal = np.cross(points, force_centrifugal)
print(f'Fa_side: {np.sum(force_aero[:, 1]):.0f}N, Fc: {np.sum(force_centrifugal_i_array):.0f}N, diff: {np.abs(np.sum(force_centrifugal_i_array) + (np.sum(force_aero[:,1]))):.0f}N')
print(f'Ma_side: {np.sum(moment_aero[:, 0]):.0f}Nm, Mc: {np.sum(moment_centrifugal[:, 0]):.0f}Nm, diff: {np.abs(np.sum(moment_centrifugal[:, 0]) + (np.sum(moment_aero[:,0]))):.0f}Nm, %: {100* np.abs(np.sum(moment_centrifugal[:, 0]) + (np.sum(moment_aero[:,0])))/np.sum(moment_aero[:, 0]):.3f}')
print(f'r_0: {r_0:.2f}m')
print(f'Ma at ac:{np.sum(force_aero[:, 1])*center_of_fa[2]:.0f}Nm, {np.cross(center_of_fa,[0,np.sum(force_aero[:, 1]),0])[0]:.0f}Nm')
print(f'Mc at cg:{np.sum(force_centrifugal_i_array)*center_of_gravity[2]:.0f}Nm, {np.cross(center_of_gravity,[0,np.sum(force_centrifugal_i_array),0])[0]:.0f}Nm')
print(f'cg: {center_of_gravity}m')

print(f'')
print(f'Moment-based ')
force_centrifugal_i_array, r_0, v_k_i_array, center_of_gravity = calculate_force_centrifugal_distribution_moment_based(
    moment_aero, config.kite.mass_points, vel_app, points)
force_centrifugal = np.array(
    [np.array([0, 1, 0])*force for i, force in enumerate(force_centrifugal_i_array)])
moment_centrifugal = np.cross(points, force_centrifugal)
print(f'Fa_side: {np.sum(force_aero[:, 1]):.0f}N, Fc: {np.sum(force_centrifugal_i_array):.0f}N, diff: {np.abs(np.sum(force_centrifugal_i_array) + (np.sum(force_aero[:,1]))):.0f}N')
print(f'Ma_side: {np.sum(moment_aero[:, 0]):.0f}Nm, Mc: {np.sum(moment_centrifugal[:, 0]):.0f}Nm, diff: {np.abs(np.sum(moment_centrifugal[:, 0]) + (np.sum(moment_aero[:,0]))):.0f}Nm, %: {100* np.abs(np.sum(moment_centrifugal[:, 0]) + (np.sum(moment_aero[:,0])))/np.sum(moment_aero[:, 0]):.3f}')
print(f'r_0: {r_0:.2f}m')

print(f'')
print(f'Tether Tension')
r_0_previous = r_0
force_centrifugal_i_array, r_0, v_k_i_array, center_of_gravity = calculate_force_centrifugal_distribution_with_tether_tension(
    r_0_previous, force_aero, config.kite.mass_points, vel_app, points)
force_centrifugal = np.array(
    [np.array([0, 1, 0])*force for i, force in enumerate(force_centrifugal_i_array)])
moment_centrifugal = np.cross(points, force_centrifugal)
print(f'Fa_side: {np.sum(force_aero[:, 1]):.0f}N, Fc: {np.sum(force_centrifugal_i_array):.0f}N, diff: {np.abs(np.sum(force_centrifugal_i_array) + (np.sum(force_aero[:,1]))):.0f}N')
print(f'Ma_side: {np.sum(moment_aero[:, 0]):.0f}Nm, Mc: {np.sum(moment_centrifugal[:, 0]):.0f}Nm, diff: {np.abs(np.sum(moment_centrifugal[:, 0]) + (np.sum(moment_aero[:,0]))):.0f}Nm, %: {100* np.abs(np.sum(moment_centrifugal[:, 0]) + (np.sum(moment_aero[:,0])))/np.sum(moment_aero[:, 0]):.3f}')
print(f'r_0: {r_0:.2f}m')


def calculate_vel_app_distribution(force_aero, config.kite.mass_points, vel_app, points):
    ''' calculates the apparant wind speed at each node, given the velocity of the kite and the turning radius'''
    # calculate center of gravity
    center_of_gravity = calculate_center_of_gravity(points, config.kite.mass_points)

    # initialize some variables
    force_aero_side = np.sum(force_aero[:, 1])
    sum_of_masses = np.sum(config.kite.mass_points)
    v_k = vel_app[0]

    # calculate turning radius of the whole kite r_k
    r_k = (sum_of_masses*v_k**2) / -force_aero_side

    # calculate distance from turning-axis to x,y plane r_0
    r_0 = r_k - center_of_gravity[1]

    # Calculate kite velocity of each node v_k,i
    r_ratio_array = np.array([((r_0 + point[1])/r_k) for point in points])
    v_k_i_array = v_k * r_ratio_array

    # calculate new apparant wind speed
    vel_app_array = [(config.vel_wind - v_k_i) for v_k_i in v_k_i_array]

    return vel_app_array, r_0, v_k_i_array, center_of_gravity

# TODO: should take spanwise variation of Va into account for aero model
# calculate new apparant wind speed
# vel_app_array = [(config.vel_wind - vel_kite_i_array[i]) for i,point in enumerate(points)]


print(f' ')
print(f' --- turning variable initialiasiation --- ')
print(f'Distance to turning center: r_0: {r_0:.2f}m')
# print(f'angular velocity: {angular_velocity_cg:.2f}rad/s')
# print(f'center of gravity: {center_of_gravity}, r_i_cg: {center_of_gravity[1]:.2f}m')
print(
    f'Diff Fy: {np.sum(force_centrifugal_i_array) - (-np.sum(force_aero[:,1])) :.3f}N, Fc: {np.sum(force_centrifugal_i_array):.2f}N, Fa,y: {np.sum(force_aero[:,1]):.2f}N, Initialisation time: {time.time()-time_start_finding_vk:.2f}s')


#%% TURNING FLIGHT
# simulating turning flight

# initialisation
vel_app[0] = -vel_kite_x  # from initialisation above
functions_print.print_settings(
    config.u_p, vel_app, AEROSTRUCTURAL_SETTINGS, AERO_CONFIG, config.kite.mass_points, SOLVER_CONFIG)
it_check_turning = 1  # number of iterations between printing out results
it_delta_ls_update = 2 # every 5 iterations, increase delta_ls
# it_update_vk = 2 #numer of iterations at which Vk needs updating
tol_fx_ratio_to_fz = 0.05  # tolerance for fx, if fx>tol*fz, then update vk_x
r_0_new = r_0
# tolerance for convergence
tol_r_0 = 0.25
tol_points_turning = 0.05
tol_solver_max_error = 0.5
# stopping conditions
max_points_diff_turning = 5  # if changes are larger than 5meters, then stop
max_sideangle = 5  # if more to the side than 5degrees, then stop
max_aerostructural_iter = 300
SOLVER_CONFIG['MAXFEV'] = 200 #TODO: lowering the maxfev as a way of relaxation

is_convergence = False
for it in range(max_aerostructural_iter+1):

    ## updating delta_ls, maxing it now at .15
    if it % it_delta_ls_update == 0 and it < 60:
        delta_ls = 0.005
        bridle_rest_lengths[2] += delta_ls  # spring[1] is the left steering tape
        bridle_rest_lengths[3] -= delta_ls  # spring[2] is the right steering tape

    print(f' ')
    # AERODYNAMICS
    time_iter = time.time()  # to print out each iteration time
    # Wing Aerodynamics
    points_left_to_right = coupling_struc2aero.order_struc_nodes_right_to_left(
        points, config.kite.connectivity.plate_point_indices)
    force_aero_wing_VSM, moment_aero_wing_VSM, F_rel, ringvec, controlpoints, wingpanels, rings, coord_L, coord_refined = VSM.calculate_force_aero_wing_VSM(
        points_left_to_right, BILLOWING_ANGLES, vel_app, AERO_CONFIG, airfoil_data, config.kite.n_segments)
    force_aero_wing = coupling_aero2struc.aero2struc(points, config.kite.connectivity.wing_ci, config.kite.connectivity.wing_cj, config.kite.connectivity.plate_point_indices, force_aero_wing_VSM,
                                                     moment_aero_wing_VSM, ringvec, controlpoints)  # Put in the new positions of the points
    # Bridle Aerodynamics
    if config.is_with_aero_bridle:
        force_aero_bridle = bridle_line_system_aero.calculate_force_aero_bridle_thedens2022(
            config.vel_wind, config.kite.connectivity.bridle_ci, config.kite.connectivity.bridle_cj, points, vel_app, config.rho, AERO_CONFIG, BRIDLE_DATA, KCU_DATA)
    else:
        force_aero_bridle = [0]
    force_aero = force_aero_wing+force_aero_bridle

    #TODO: remove
    # calculating the moment of the aerodynamic forces
    moment_aero = np.cross(points, force_aero)
    # saving old force_aero for printing
    force_aero = force_aero
    # moment_aero = np.array([np.cross(point,force_aero[i]) for i,point in enumerate(points)])
    # Ma_x = np.array(force_aero[:, 2] * points[:, 1] - force_aero[:, 1] * points[:, 2])
    # Ma_y = np.array(force_aero[:, 0] * points[:, 2] - force_aero[:, 2] * points[:, 0])
    # Ma_z = np.array(force_aero[:, 1] * points[:, 0] - force_aero[:, 0] * points[:, 1])

    time_iter_aero = time.time()

    # CHECKING FOR VK UPDATE NECESSITY
    # Updating the velocity of the kite in the x-direction, to make Fx go to zero
    # run if fx > tol_fx_ratio_to_fz*Fz
    force_resultant_x = solver_utils.calculate_fx(-vel_app[0], vel_app, points, args,args_bridle)
    if np.abs(force_resultant_x) > tol_fx_ratio_to_fz*np.abs(np.sum(force_aero[:, 2])):
        print(f'i{it}, updating vel_kite_x')
        # vel_app[0] = -finding_vk_which_makes_fx_zero(initial_guess_vk_x,points,config.vel_wind,BILLOWING_ANGLES, AERO_CONFIG, airfoil_data, config.kite.n_segments, config.kite.connectivity.wing_ci, config.kite.connectivity.wing_cj, config.kite.connectivity.plate_point_indices)
        args = (config.vel_wind, BILLOWING_ANGLES, AERO_CONFIG, airfoil_data,
                config.kite.n_segments, config.kite.connectivity.wing_ci, config.kite.connectivity.wing_cj, config.kite.connectivity.plate_point_indices)
        vel_app[0] = - \
            solver_utils.optimalisation_of_vk_for_fx_0(
                -vel_app[0], vel_app, points, args, args_bridle)

    # CALCULATING CENTRIFUGAL FORCE DISTRIBUTION
    # calculate centrifugal force distribution, i.e. the force on each node Fc,i
    force_centrifugal_i_array, r_0_new, v_k_i_array, center_of_gravity = calculate_force_centrifugal_distribution_moment_based(
        moment_aero, config.kite.mass_points, vel_app, points)
    fc_moment = np.sum(force_centrifugal_i_array)
    mc_moment = np.sum(np.cross(points,np.array(
        [np.array([0, 1, 0])*force for i, force in enumerate(force_centrifugal_i_array)])),axis=0)[0]
    force_centrifugal_i_array, r_0_new, v_k_i_array, center_of_gravity = calculate_force_centrifugal_distribution(
        force_aero, config.kite.mass_points, vel_app, points)
    fc_equal = np.sum(force_centrifugal_i_array)
    mc_equal = np.sum(np.cross(points,np.array(
        [np.array([0, 1, 0])*force for i, force in enumerate(force_centrifugal_i_array)])),axis=0)[0]
    force_centrifugal_i_array, r_0_new, v_k_i_array, center_of_gravity = calculate_force_centrifugal_distribution_with_tether_tension(
        r_0_new, force_aero, config.kite.mass_points, vel_app, points)
    fc_ft= np.sum(force_centrifugal_i_array)
    mc_ft = np.sum(np.cross(points,np.array(
        [np.array([0, 1, 0])*force for i, force in enumerate(force_centrifugal_i_array)])),axis=0)[0]   
    print(f'-------Ma:{np.sum(moment_aero,axis=0)[0]:.0f}Nm----------')
    print(f'Mc_equal: {mc_equal:.0f}Nm, Mc_moment:{mc_moment:.0f}Nm, Mc_ft:{mc_ft:.0f}Nm')
    print(f'Fc_equal: {fc_equal:.0f}N, Fc_moment:{fc_moment:.0f}N, Fc_ft:{fc_ft:.0f}N')
     #TODO: this is just oposing each y force directly, essentially just setting it to zero
    # for i,fa_i in enumerate(force_aero):
    #     force_aero[i][1] = 0


    # TODO: remove, additional lines were for testing the calculation
    # calculating the moment due to the centrifugal forces
    # moment_centrifugal_x = np.array([-force_centrifugal_i_array[i]*point[2] for i,point in enumerate(points)])
    # moment_centrifugal_y = np.array([0*point[0] for i,point in enumerate(points)])
    # moment_centrifugal_z = np.array([force_centrifugal_i_array[i]*point[0] for i,point in enumerate(points)])
    force_centrifugal = np.array(
        [np.array([0, 1, 0])*force for i, force in enumerate(force_centrifugal_i_array)])
    moment_centrifugal = np.cross(points, force_centrifugal)
    # moment_centrifugal = np.array([np.cross(point,force_centrifugal[i]) for i,point in enumerate(points)])
    # Mc_x = np.array(force_centrifugal[:, 2] * points[:, 1] - force_centrifugal[:, 1] * points[:, 2])
    # Mc_y = np.array(force_centrifugal[:, 0] * points[:, 2] - force_centrifugal[:, 2] * points[:, 0])
    # Mc_z = np.array(force_centrifugal[:, 1] * points[:, 0] - force_centrifugal[:, 0] * points[:, 1])

    # TODO: unphysical correction, which also doesn't work
    # force_centrifugal = solver_utils.optimize_fc_to_get_mx_to_zero(force_centrifugal,points)

    # STRUCTURAL
    # TODO: commented below is initial attempt at bending-resistance, needed for multiple chord-wise points
    # if config.kite_name == 'V9_60C':
    #     force_strut = np.zeros(config.kite.points_ini.shape)
    #     for key in points_between_dict.keys():
    #         point = points[int(key)]
    #         line_end1 = points[points_between_dict[key][0]]
    #         line_end2 = points[points_between_dict[key][1]]
    #         distance,distance_vector,angle = structural_model.distance_point_to_line(point,[line_end1,line_end2])
    #         force_strut[int(key)] += STIFFNESS_DATA['STRUT']*distance_vector        # adding force to force_aero, as an external force
    #     # adding strut force to aero force
    #     force_aero += force_strut

    # EXTERNAL
    force_external = np.zeros(config.kite.points_ini.shape)
    moment_external = np.zeros(config.kite.points_ini.shape)
    for i in range(config.kite.points_ini.shape[0]):
        force_external[i] = force_aero[i] + force_centrifugal[i]
        moment_external[i] = moment_aero[i] + moment_centrifugal[i]

    # SOLVER (here the structural solver is called)
    points_new, solver_max_error = solver_main.calculate_points_new(
        points, solver_arguments, force_external, SOLVER_CONFIG)

    #TODO: RELAXATION APPLIED TO POINT DISPLACEMENT
    points_new = 0.75*points_new + 0.25*points
    #TODO: A type of damping applied after x number of iterations
    if it > 150 and it%10 == 0:
        STIFFNESS_DATA['BRIDLE'] *= 1.05
        STIFFNESS_DATA['TUBE'] *= 1.05
        STIFFNESS_DATA['TE'] *= 1.05

    # CALCULATE NEW APPARANT WIND SPEED
    # calculate new apparant wind speed distribution
    # TODO: include this effect
    # vel_app_array, r_0, v_k_i_array, center_of_gravity = calculate_vel_app_distribution(force_aero,config.kite.mass_points,vel_app,points_new)

    # Convergence
    # iter_difference with last iteration
    iter_diff = np.amax(np.linalg.norm(points_new-points, axis=1))
    iter_diff_r_0 = np.abs(r_0_new - r_0)
    sideangle_leading_edge_midpoint = functions_print.calculate_side_angle_leading_edge_midpoint(
        points_new[config.kite.connectivity.plate_point_indices[int(config.kite.n_segments/2)]][:2])

    if iter_diff < tol_points_turning and solver_max_error < tol_solver_max_error and iter_diff_r_0 < tol_r_0:  # if converged
        print(f'Converged at: i{it}, final update of vel_kite_x')
        # vel_app[0] = -finding_vk_which_makes_fx_zero(initial_guess_vk_x,points,config.vel_wind,BILLOWING_ANGLES, AERO_CONFIG, airfoil_data, config.kite.n_segments, config.kite.connectivity.wing_ci, config.kite.connectivity.wing_cj, config.kite.connectivity.plate_point_indices)
        args = (config.vel_wind, BILLOWING_ANGLES, AERO_CONFIG, airfoil_data,
                config.kite.n_segments, config.kite.connectivity.wing_ci, config.kite.connectivity.wing_cj, config.kite.connectivity.plate_point_indices)
        vel_app[0] = - \
            solver_utils.optimalisation_of_vk_for_fx_0(
                -vel_app[0], vel_app, points, args, args_bridle)
        is_convergence = True
        break

    # elif it == AEROSTRUCTURAL_SETTINGS['MAX_ITER']: #if reached max-iterations
    #     is_convergence = False
    #     break

    # w hen changes are bigger than 2meters
    # or when sideangle > 5deg

    elif iter_diff > max_points_diff_turning or np.abs(sideangle_leading_edge_midpoint) > max_sideangle:
        print(
            f'iter_diff: {iter_diff:.3f} > {max_points_diff_turning} and/or {np.abs(sideangle_leading_edge_midpoint)}>{max_sideangle:.3f}deg')
        print(f'--> Change is too large to be reasonable, stopping iterations')
        is_convergence = False
        break

    else:  # if not convergenced
        points = points_new
        r_0 = r_0_new
        if config.is_print_mid_results:
            if it % it_check_turning == 0:
                trim_angle = functions_print.calculate_trim_angle(
                    coord_refined, AERO_CONFIG['N_SPLITS'], points_new[config.kite.connectivity.plate_point_indices[int(config.kite.n_segments/2)]][:2])
                print(f'i{it} - diff: {1e3*iter_diff:.3f}mm - diff r_0: {1e3*iter_diff_r_0:.3f}mm- error: {solver_max_error:.3f}N? - side: {sideangle_leading_edge_midpoint:.2f}deg - trim: {trim_angle:.3f}deg - dt: {(time.time() - time_iter):.2f}s (aero: {(time_iter_aero - time_iter):.2f} - struc: {(time.time() - time_iter_aero):.2f})')

    #TODO: remove
    # , Fc: {np.sum(force_centrifugal_i_array):.2f}N, Fa,y: {np.sum(force_aero_pre_fc[:,1]):.2f}N')
    print(
        f'r_0: {r_0_new:.2f}m, Diff Fy: {np.sum(force_centrifugal_i_array) - (-np.sum(force_aero[:,1])) :.3f}N, vel_app:{vel_app}')
    print(f'Mx(roll) : {np.sum(moment_external[:,0]):.1f}, Fx: {np.sum(force_external[:,0]):.1f}, Ma_x: {np.sum(moment_aero[:,0]):.1f}, Mc_x: {np.sum(moment_centrifugal[:,0]):.1f}, Fa_x: {np.sum(force_aero[:,0]):.1f}, Fc_x: {np.sum(force_centrifugal[:,0]):.1f}')
    print(f'My(pitch): {np.sum(moment_external[:,1]):.1f}, Fy: {np.sum(force_external[:,1]):.1f}, Ma_y: {np.sum(moment_aero[:,1]):.1f}, Mc_y: {np.sum(moment_centrifugal[:,1]):.1f}, Fa_y: {np.sum(force_aero[:,1]):.1f}, Fc_y: {np.sum(force_centrifugal[:,1]):.1f}')
    print(f'Mz(yaw)  : {np.sum(moment_external[:,2]):.1f}, Fz: {np.sum(force_external[:,2]):.1f}, Ma_z: {np.sum(moment_aero[:,2]):.1f}, Mc_z: {np.sum(moment_centrifugal[:,2]):.1f}, Fa_z: {np.sum(force_aero[:,2]):.1f}, Fc_z: {np.sum(force_centrifugal[:,2]):.1f}')

np.save(f"results/V3_25/points_up_{str(int(100*config.u_p))}_turning.npy", points_new)

# post-processing
aerostructural_total_time = time.time()-aerostructural_start_time

# Computing trim-angle
if AERO_CONFIG['MODEL'] == ('VSM' or 'LLT'):
    trim_angle = trim_angle
elif AERO_CONFIG['MODEL'] == 'plate_aero':
    trim_angle = 10.  # TODO: write a function for this as well

# Printing out the solution nicely
functions_print.print_aerostructural(
    is_convergence, it, AEROSTRUCTURAL_SETTINGS['MAX_ITER'], aerostructural_total_time)
functions_print.print_structural(points_new, config.kite.points_ini, config.kite.connectivity.bridle_ci, config.kite.connectivity.bridle_cj, config.kite.connectivity.wing_ci,
                                 config.kite.connectivity.wing_cj, bridle_rest_lengths, wing_rest_lengths, config.u_p, STIFFNESS_DATA, config.kite.connectivity.tube_line_indices)

# computing the aero-again as its somehow not correctly stored ##TODO: fix-this
# Wing Aerodynamics
points_left_to_right = coupling_struc2aero.order_struc_nodes_right_to_left(
    points_new, config.kite.connectivity.plate_point_indices)
force_aero_wing_VSM, moment_aero_wing_VSM, F_rel, ringvec, controlpoints, wingpanels, rings, coord_L, coord_refined = VSM.calculate_force_aero_wing_VSM(
    points_left_to_right, BILLOWING_ANGLES, vel_app, AERO_CONFIG, airfoil_data, config.kite.n_segments)
force_aero_wing = coupling_aero2struc.aero2struc(points_new, config.kite.connectivity.wing_ci, config.kite.connectivity.wing_cj, config.kite.connectivity.plate_point_indices,
                                                 force_aero_wing_VSM, moment_aero_wing_VSM, ringvec, controlpoints)  # Put in the new positions of the points

# Bridle Aerodynamics
if config.is_with_aero_bridle:
    force_aero_bridle = bridle_line_system_aero.calculate_force_aero_bridle_thedens2022(
        config.vel_wind, config.kite.connectivity.bridle_ci, config.kite.connectivity.bridle_cj, points_new, vel_app, config.rho, AERO_CONFIG, BRIDLE_DATA, KCU_DATA)
else:
    force_aero_bridle = [0]
force_aero = force_aero_wing+force_aero_bridle
# Total force
force_total = solver_main.calculate_total_force(
    config.kite.points_ini, solver_arguments, force_aero)
print(
    f'force_total: Fx: {np.sum(force_total[:,0]):.2f}N, Fy: {np.sum(force_total[:,1]):.2f}N, Fz: {np.sum(force_total[:,2]):.2f}N')
# Spring force
force_spring = structural_model.calculate_force_spring(
    config.kite.points_ini, calculate_force_spring_arguments)
print(
    f'force_spring: Fx: {np.sum(force_spring[:,0]):.5f}N, Fy: {np.sum(force_spring[:,1]):.5f}N, Fz: {np.sum(force_spring[:,2]):.5f}N')
print(f' ---------------------------')
functions_print.print_aero(points_new, config.kite.connectivity.plate_point_indices, config.kite.n_segments, vel_app, config.vel_wind, force_aero,
                           config.kite.area_projected, AERO_CONFIG, force_aero_bridle, KCU_DATA['INDEX'], config.rho, config.kite.ref_chord, config.mu)
functions_plot.plot_aero(config.vel_wind, config.vel_kite, vel_app, wingpanels,
                         controlpoints, rings, coord_L, config.kite.n_segments, F_rel, elev=10, azim=230)
# functions_plot.plot_aero(config.vel_wind, config.vel_kite, vel_app, wingpanels, controlpoints, rings, coord_L, config.kite.n_segments, F_rel, elev=0, azim=-90)


# Plotting the solution - offline PLOTLY
plot_points_old = [2, config.kite.points_ini, 'black', 2]
plot_points_new = [2, points_new, 'green', 4]
elongation_values = post_processing_utils.calculate_elongation(
    points_new, config.kite.connectivity.wing_ci, config.kite.connectivity.wing_cj, config.kite.connectivity.bridle_ci, config.kite.connectivity.bridle_cj, wing_rest_lengths, bridle_rest_lengths, config.kite.connectivity.tube_line_indices)[2]
plot_lines_old = [config.kite.points_ini, np.concatenate(
    (config.kite.connectivity.bridle_ci, config.kite.connectivity.wing_ci)), np.concatenate((config.kite.connectivity.bridle_cj, config.kite.connectivity.wing_cj)), 'grey', 1.5, []]
plot_lines_new = [points_new, np.concatenate((config.kite.connectivity.bridle_ci, config.kite.connectivity.wing_ci)), np.concatenate(
    (config.kite.connectivity.bridle_cj, config.kite.connectivity.wing_cj)), 'purple', 4, elongation_values]
plot_surface_old = [True, config.kite.points_ini, config.kite.connectivity.plate_point_indices, 'purple', 0.0]
plot_surface_new = [True, points_new, config.kite.connectivity.plate_point_indices, 'lightgrey', 0.3]
functions_plot.plot_kite([plot_points_new, plot_points_old], [
                         plot_lines_new, plot_lines_old], f" ", [plot_surface_new, plot_surface_old])


#
#
#
#
#
#
#
#
#
#
#
#
#
#
#

# %% PLOTTING OUTPUT SEPARATELY
# plotting turning output separately

points_us_10 = np.load(
    f"results/V3_25/points_up_{str(int(100*config.u_p))}_turning_us_10.npy")
points_new = points_us_10

# Computing trim-angle
if AERO_CONFIG['MODEL'] == ('VSM' or 'LLT'):
    trim_angle = trim_angle
elif AERO_CONFIG['MODEL'] == 'plate_aero':
    trim_angle = 10.  # TODO: write a function for this as well

# Printing out the solution nicely
functions_print.print_aerostructural(
    is_convergence, it, AEROSTRUCTURAL_SETTINGS['MAX_ITER'], aerostructural_total_time)
functions_print.print_structural(points_new, config.kite.points_ini, config.kite.connectivity.bridle_ci, config.kite.connectivity.bridle_cj, config.kite.connectivity.wing_ci,
                                 config.kite.connectivity.wing_cj, bridle_rest_lengths, wing_rest_lengths, config.u_p, STIFFNESS_DATA, config.kite.connectivity.tube_line_indices)

# computing the aero-again as its somehow not correctly stored ##TODO: fix-this
# Wing Aerodynamics
points_left_to_right = coupling_struc2aero.order_struc_nodes_right_to_left(
    points_new, config.kite.connectivity.plate_point_indices)
force_aero_wing_VSM, moment_aero_wing_VSM, F_rel, ringvec, controlpoints, wingpanels, rings, coord_L, coord_refined = VSM.calculate_force_aero_wing_VSM(
    points_left_to_right, BILLOWING_ANGLES, vel_app, AERO_CONFIG, airfoil_data, config.kite.n_segments)
force_aero_wing = coupling_aero2struc.aero2struc(points_new, config.kite.connectivity.wing_ci, config.kite.connectivity.wing_cj, config.kite.connectivity.plate_point_indices,
                                                 force_aero_wing_VSM, moment_aero_wing_VSM, ringvec, controlpoints)  # Put in the new positions of the points

# Bridle Aerodynamics
if config.is_with_aero_bridle:
    force_aero_bridle = bridle_line_system_aero.calculate_force_aero_bridle_thedens2022(
        config.vel_wind, config.kite.connectivity.bridle_ci, config.kite.connectivity.bridle_cj, points_new, vel_app, config.rho, AERO_CONFIG, BRIDLE_DATA, KCU_DATA)
else:
    force_aero_bridle = [0]
force_aero = force_aero_wing+force_aero_bridle
# Total force
force_total = solver_main.calculate_total_force(
    config.kite.points_ini, solver_arguments, force_aero)
print(
    f'force_total: Fx: {np.sum(force_total[:,0]):.2f}N, Fy: {np.sum(force_total[:,1]):.2f}N, Fz: {np.sum(force_total[:,2]):.2f}N')
# Spring force
force_spring = structural_model.calculate_force_spring(
    config.kite.points_ini, calculate_force_spring_arguments)
print(
    f'force_spring: Fx: {np.sum(force_spring[:,0]):.5f}N, Fy: {np.sum(force_spring[:,1]):.5f}N, Fz: {np.sum(force_spring[:,2]):.5f}N')
print(f' ---------------------------')
functions_print.print_aero(points_new, config.kite.connectivity.plate_point_indices, config.kite.n_segments, vel_app, config.vel_wind, force_aero,
                           config.kite.area_projected, AERO_CONFIG, force_aero_bridle, KCU_DATA['INDEX'], config.rho, config.kite.ref_chord, config.mu)
functions_plot.plot_aero(config.vel_wind, config.vel_kite, vel_app, wingpanels,
                         controlpoints, rings, coord_L, config.kite.n_segments, F_rel, elev=10, azim=230)
# functions_plot.plot_aero(config.vel_wind, config.vel_kite, vel_app, wingpanels, controlpoints, rings, coord_L, config.kite.n_segments, F_rel, elev=0, azim=-90)


# Plotting the solution - offline PLOTLY
plot_points_old = [2, config.kite.points_ini, 'black', 2]
plot_points_new = [2, points_new, 'green', 4]
elongation_values = post_processing_utils.calculate_elongation(
    points_new, config.kite.connectivity.wing_ci, config.kite.connectivity.wing_cj, config.kite.connectivity.bridle_ci, config.kite.connectivity.bridle_cj, wing_rest_lengths, bridle_rest_lengths, config.kite.connectivity.tube_line_indices)[2]
plot_lines_old = [config.kite.points_ini, np.concatenate(
    (config.kite.connectivity.bridle_ci, config.kite.connectivity.wing_ci)), np.concatenate((config.kite.connectivity.bridle_cj, config.kite.connectivity.wing_cj)), 'grey', 1.5, []]
plot_lines_new = [points_new, np.concatenate((config.kite.connectivity.bridle_ci, config.kite.connectivity.wing_ci)), np.concatenate(
    (config.kite.connectivity.bridle_cj, config.kite.connectivity.wing_cj)), 'purple', 4, elongation_values]
plot_surface_old = [True, config.kite.points_ini, config.kite.connectivity.plate_point_indices, 'purple', 0.0]
plot_surface_new = [True, points_new, config.kite.connectivity.plate_point_indices, 'lightgrey', 0.3]
functions_plot.plot_kite([plot_points_new, plot_points_old], [
                         plot_lines_new, plot_lines_old], f" ", [plot_surface_new, plot_surface_old])

# %% what should be included according to python standards


def main():
    print('This is the main function, need to add the right functionalities in here')


# To show that this file is a script and not a module to be imported elsewhere
if __name__ == '__main__':  # "if we run the module directly
    main()  # run the code within this conditional"
