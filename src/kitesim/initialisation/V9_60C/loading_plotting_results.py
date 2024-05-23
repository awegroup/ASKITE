# %% IMPORTING

name_file = "V9_60C_4_53_M_lines"
path_folder = "initialisation/V9_60C/"
path_results = "results_basic_bridle_determined_shape/"

import numpy as np
import sys

sys.path.append("initialisation/V9_60C/")
from kitesim.initialisation.V9_60C.functions_plotting import (
    get_3D_plot_N_D_6times_PLUS_2lines,
)


# %% LOADING PLOTTING RESULTS

### LOADING
points_ini_model = np.load(f"{path_folder}results/points_ini_model_bridle.npy")
# rib_strut_LE                = np.load(f'{path_folder}results/rib_strut_LE.npy' )
# rib_strut_canopy            = np.load(f'{path_folder}results/rib_strut_canopy.npy',allow_pickle=True )
# rib_strut_TE                = np.load(f'{path_folder}results/rib_strut_TE.npy' )
# rib_bridle_points_list      = np.load(f'{path_folder}results/rib_bridle_points_list.npy' )
# supervector_wing            = np.load(f'{path_folder}{path_results}supervector_wing.npy' )
# supervector_bridles         = np.load(f'{path_folder}{path_results}supervector_bridles.npy' )
ci_bridle = np.load(f"{path_folder}{path_results}ci_bridle_indexed_on_bridle.npy")
cj_bridle = np.load(f"{path_folder}{path_results}cj_bridle_indexed_on_bridle.npy")
ci_wing = np.load(f"{path_folder}{path_results}conn_wing_i.npy")
cj_wing = np.load(f"{path_folder}{path_results}conn_wing_j.npy")
# plot_mid_tip_LE             = np.load(f'{path_folder}results/plot_mid_tip_LE.npy',allow_pickle=True )
# vertices_wing               = np.load(f'{path_folder}results/vertices_wing.npy' )
# pulley_idx_list             = np.load(f'{path_folder}results/pulley_idx_list.npy' )
# supervector_bridles_pulleys = np.load(f'{path_folder}results/supervector_bridles_pulleys.npy',allow_pickle=True )

### PLOTTING
## Plotting wing
# get_3D_plot_N_D(2,vertices_wing, 'blue',1)

# #%% plotting wind discretization by points
# plot_strut_LE = [3,rib_strut_LE, 'red',1]
# plot_strut_canopy = [3,rib_strut_canopy, 'green',1]
# plot_strut_TE = [3,rib_strut_TE, 'blue',1]
# plot_strut_particles = [2,supervector_wing, 'black',4]
# plot_strut_particles_bridle = [2,rib_bridle_points_list, 'orange',5]
# get_3D_plot_N_D_4times(plot_strut_particles_bridle,plot_strut_LE,plot_strut_canopy,plot_strut_particles)

# #%% plotting the wing connectivity
# supervector = np.concatenate((supervector_wing,supervector_bridles),axis=0)
# plot_line_varibles = [supervector,ci_wing,cj_wing, 'blue',2]
# get_3D_plot_N_D_4times_PLUS_lines(plot_strut_particles_bridle,plot_strut_LE,plot_mid_tip_LE,plot_strut_particles,plot_line_varibles)

# %% Plotting wing and bridle

##TODO: could add surface to the wing, might look nicer

## wing
# LE circles
plot_mid_tip_LE = [3, plot_mid_tip_LE[1], "red", 0.1]
plot_strut_LE = [3, rib_strut_LE, "red", 0.1]

# canopy
plot_strut_particles = [2, supervector_wing, "blue", 2]

# bridle line connection points
plot_strut_particles_bridle = [2, rib_bridle_points_list, "orange", 3]

# connections
plot_line_varibles_wing = [supervector, ci_wing, cj_wing, "blue", 1]

## bridles
# bridle points
plot_bridle_points = [2, supervector_bridles, "black", 2]

# pulleys
plot_pulleys = [2, supervector_bridles_pulleys, "red", 3]

# bridle lines connections
plot_line_varibles_bridle = [supervector, ci_bridle, cj_bridle, "black", 1]

get_3D_plot_N_D_6times_PLUS_2lines(
    plot_strut_particles_bridle,
    plot_strut_LE,
    plot_mid_tip_LE,
    plot_strut_particles,
    plot_bridle_points,
    plot_pulleys,
    plot_line_varibles_wing,
    plot_line_varibles_bridle,
)


# %% saving code as txt

### LOADING
points_ini_model = np.load("results/" + name_file + "_points_ini_model.npy")
supervector_wing = np.load("results/" + name_file + "_supervector_wing.npy")
supervector_bridles = np.load("results/" + name_file + "_supervector_bridles.npy")
ci_bridle = np.load("results/" + name_file + "_ci_bridle.npy")
cj_bridle = np.load("results/" + name_file + "_cj_bridle.npy")
ci_wing = np.load("results/" + name_file + "_ci_wing.npy")
cj_wing = np.load("results/" + name_file + "_cj_wing.npy")
rib_strut_LE = np.load("results/" + name_file + "_rib_strut_LE.npy")
rib_strut_canopy = np.load(
    "results/" + name_file + "_rib_strut_canopy.npy", allow_pickle=True
)
rib_strut_TE = np.load("results/" + name_file + "_rib_strut_TE.npy")
rib_bridle_points_list = np.load("results/" + name_file + "_rib_bridle_points_list.npy")
plot_mid_tip_LE = np.load(
    "results/" + name_file + "_plot_mid_tip_LE.npy", allow_pickle=True
)
vertices_wing = np.load("results/" + name_file + "_vertices_wing.npy")
pulley_idx_list = np.load("results/" + name_file + "_pulley_idx_list.npy")
supervector_bridles_pulleys = np.load(
    "results/" + name_file + "_supervector_bridles_pulleys.npy", allow_pickle=True
)
