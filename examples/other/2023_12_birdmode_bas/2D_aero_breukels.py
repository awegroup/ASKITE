import numpy as np
from functions import *
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

## Define here your variables, below are the standards for the V9_60C
wing_nodes = np.array(
    [
        [0.514901, -6.54240262, 11.11837575],
        [2.15951, -6.471868, 10.63388412],
        [0.040401, -6.21104462, 12.371314],
        [2.467556, -6.15357625, 12.22461763],
        [-0.39389788, -5.23368925, 13.60858212],
        [2.58715463, -5.18912088, 13.37017787],
        [-0.690224, -4.06784637, 14.45577938],
        [2.60166238, -4.0510165, 14.20744913],
        [-0.882376, -2.83106275, 15.0050615],
        [2.60087538, -2.82559063, 14.756022],
        [-0.99232887, -1.574333, 15.31896337],
        [2.60026538, -1.571548, 15.07032875],
        [-1.03755962, -0.315193, 15.44822762],
        [2.599966, -0.30895725, 15.19885587],
        [-1.03755962, 0.315193, 15.44822763],
        [2.599966, 0.30895725, 15.19885587],
        [-0.99232887, 1.574333, 15.31896337],
        [2.60026538, 1.571548, 15.07032875],
        [-0.882376, 2.83106275, 15.0050615],
        [2.60087538, 2.82559063, 14.756022],
        [-0.690224, 4.06784638, 14.45577938],
        [2.60166238, 4.0510165, 14.20744913],
        [-0.39389788, 5.23368925, 13.60858212],
        [2.58715463, 5.18912088, 13.37017788],
        [0.040401, 6.21104462, 12.371314],
        [2.467556, 6.15357625, 12.22461763],
        [0.514901, 6.54240262, 11.11837575],
        [2.15951, 6.471868, 10.63388412],
    ]
)

## plotting the wing_nodes
fig = plt.figure()
ax = fig.add_subplot(111, projection="3d")
ax.scatter(
    wing_nodes[:, 0],
    wing_nodes[:, 1],
    wing_nodes[:, 2],
    c="b",
    marker="o",
)
ax.set_xlabel("X")
ax.set_ylabel("Y")
ax.set_zlabel("Z")
ax.axis("equal")
ax.set_title("V9_60C Wing Nodes")
plt.show()

##TODO: I never got to correcting these! So pleas updated
# these are the angles at which the canopy leaves the strut
# compared to the line that would make for a flat-segment
billowing_angles = np.array(
    [5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0]
)
ring_geometry = "5fil"
model = "VSM"
n_segments = 13
##TODO: I never got to correcting these (I think), so please verify/update
tube_diameters = np.array(
    [
        0.2185577929737451,
        0.2504488308906025,
        0.27475599315277693,
        0.29206438873880386,
        0.3037224448203036,
        0.3106631128151186,
        0.3137430787012053,
        0.3137430787012053,
        0.3106631128151186,
        0.3037224448203036,
        0.29206438873880386,
        0.27475599315277693,
        0.2504488308906025,
        0.2185577929737451,
    ]
)
is_tube_diameter_dimensionless = False
##TODO: I never got to correcting these (I think), so please verify/update
canopy_max_heights = np.array(
    [
        0.09531306981488247,
        0.16014092745597014,
        0.18267014293931405,
        0.15352458016180195,
        0.11363127219081569,
        0.06805281069119812,
        0.01751718992699616,
        0.01751718992699616,
        0.06805281069119812,
        0.11363127219081569,
        0.15352458016180195,
        0.18267014293931405,
        0.16014092745597014,
        0.09531306981488247,
    ]
)
is_canopy_max_height_dimensionless = True

########################################################
## Tuning paramaters
vel_app = np.array([10.0, 0.0, 2.0])
n_splits = 2


## Defining some additional geometry and variables
coord = refine_LEI_mesh_ballooning(wing_nodes, billowing_angles, n_splits + 1)
n_panels_aero = int(len(coord) / 2) - 1
controlpoints, rings, wingpanels, ringvec, coord_L = create_geometry_LEI(
    coord, vel_app, n_panels_aero + 1, ring_geometry, model
)

## 2D precomputed CFD (RANS) from Breukels 2011 PhD Thesis
data_airf = calculate_polar_lookup_table(
    controlpoints,
    n_segments,
    n_panels_aero,
    n_splits,
    tube_diameters,
    is_tube_diameter_dimensionless,
    canopy_max_heights,
    is_canopy_max_height_dimensionless,
)

# Plot cl-alpha
section_number = 10  # there are 13 wingsegments and 2 splits, so 26 sections
plot_polars(data_airf, section_number)
