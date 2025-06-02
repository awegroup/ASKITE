# %% Loading necessary libraries and loading the input data

import time
import numpy as np
import pandas as pd
import logging
import os
from pathlib import Path
import json
import matplotlib.pyplot as plt

from kitesim.logging_config import *
from PSS.particleSystem import ParticleSystem
from PSS.particleSystem.SpringDamper import SpringDamperType

# Find the root directory of the repository
root_dir = os.path.abspath(os.path.dirname(__file__))
while not os.path.isfile(os.path.join(root_dir, ".gitignore")):
    root_dir = os.path.abspath(os.path.join(root_dir, ".."))
    if root_dir == "/":
        raise FileNotFoundError("Could not find the root directory of the repository.")
## defining paths
path_kite_example_input = Path(root_dir) / "processed_data" / "V3_25" / "PSS_test_input"
logging.info(f"current path: {Path.cwd()}")
logging.info(f"path_kite_example_input: {path_kite_example_input}")


def load_data_from_json(folder_path: Path) -> dict:
    """Load all JSON files from the specified folder into a dictionary.

    Args:
        folder_path (Path): Path to the folder containing JSON files.

    Returns:
        dict: Dictionary containing the loaded data.
    """
    data = {}
    for file_path in folder_path.glob("*.json"):
        with open(file_path, "r") as f:
            name = file_path.stem
            data[name] = json.load(f)
        logging.info(f"Loaded {name} from {file_path}")
    return data


sim_input = load_data_from_json(path_kite_example_input)
logging.info(f"sim_input: {sim_input.keys()}")

# %% Transforming to the right inputs and instantiating the Particle System


connectivity_matrix = sim_input["connectivity_matrix"]


# converting the string spring types to the enum
def string_to_springdampertype(link_type: str) -> SpringDamperType:
    """
    Convert a string representation of a link type to a SpringDamperType enum value.

    Args:
        link_type (str): String representation of the link type.

    Returns:
        SpringDamperType: Corresponding enum value.

    Raises:
        ValueError: If the input string doesn't match any SpringDamperType.
    """
    try:
        return SpringDamperType(link_type.lower())
    except ValueError:
        raise ValueError(f"Invalid link type: {link_type}")


# Convert the connectivity matrix
connectivity_matrix = [
    [conn[0], conn[1], conn[2], conn[3], string_to_springdampertype(conn[4])]
    for conn in connectivity_matrix
]

initial_conditions = sim_input["initial_conditions"]
params = sim_input["params"]

logging.debug(f"connectivity_matrix: {connectivity_matrix}")
logging.debug(f"initial_conditions: {initial_conditions}")
logging.debug(f"params: {params}")

## Instating the PSM
PS = ParticleSystem(
    connectivity_matrix,
    initial_conditions,
    params,
)

# Changing the params to check if the problem lies there
params["dt"] = 1e-1
params["c"] = 1

# Checking if all the right linktypes are present
for link in PS.springdampers:
    logging.debug(f"link: {link.linktype}")

initial_positions = [
    [particle.x, particle.v, particle.m, particle.fixed] for particle in PS.particles
]

# %% Running the simulation

f_ext = sim_input["f_external"]

t_vector = np.linspace(
    params["dt"], params["t_steps"] * params["dt"], params["t_steps"]
)
tol = 1e-10  # 1e-29
final_step = 0
E_kin = []
f_int = []

# And run the simulation
for step in t_vector:
    PS.kin_damp_sim(f_ext)

    final_step = step
    (
        x,
        v,
    ) = PS.x_v_current
    E_kin.append(np.linalg.norm(v * v))
    f_int.append(np.linalg.norm(PS.f_int))

    converged = False
    if step > 10:
        logging.debug(
            f"step: {step}, x: {x}, v: {v}, E_kin: {E_kin[-1]}, f_int: {f_int[-1]}"
        )

        if np.max(E_kin[-10:-1]) <= tol:
            converged = True
    if converged and step > 1:
        print(f"Kinetic damping PS converged: {step}s")
        break
    if (
        np.isnan(x).any()
        or np.isnan(v).any()
        or np.isnan(E_kin[-1])
        or np.isnan(f_int[-1])
    ):
        print(f"Kinetic damping PS diverged after: {step}s")
        break


# %% Checking the convergence plot

# Change these to zoom in on a specific region.
plotstop = len(E_kin)
plotstart = 0

plt.plot(t_vector[:plotstop], E_kin[:plotstop], label="Kinetic Energy [J]")
plt.plot(
    t_vector[plotstart:plotstop], f_int[plotstart:plotstop], label="Internal Forces [N]"
)

# Filtering out the peaks to mark where the kinetic damping algorithm kicked in.
df = pd.DataFrame(E_kin, index=t_vector[0:plotstop])
peaksonly = df[(df[0].shift(1) < df[0]) & (df[0].shift(-1) < df[0])]
plt.scatter(
    peaksonly.index,
    peaksonly[0],
    c="r",
    linewidths=1,
    marker="+",
    label="Kinetic Damping Enacted",
)

plt.legend()
plt.yscale("log")
plt.title("Convergence Plot")
plt.xlabel("Time [s]")
plt.ylabel("Quantity of interest")
plt.xlim(t_vector[0], t_vector[plotstop])
plt.show()

# %% Reviewing results

final_positions = [
    [particle.x, particle.v, particle.m, particle.fixed] for particle in PS.particles
]

# Plotting final results
fig = plt.figure()
ax = fig.add_subplot(projection="3d")

for i, node in enumerate(final_positions):
    if node[3]:
        ax.scatter(node[0][0], node[0][1], node[0][2], color="red", marker="o")
    else:
        ax.scatter(node[0][0], node[0][1], node[0][2], color="red", marker="o")

    ax.quiver(
        *node[0].tolist(), f_ext[3 * i], f_ext[3 * i + 1], f_ext[3 * i + 2], length=0.01
    )

for connection in connectivity_matrix:
    line = np.column_stack(
        [final_positions[connection[0]][0], final_positions[connection[1]][0]]
    )

    ax.plot(line[0], line[1], line[2], color="red")

ax.legend(["Fixed nodes", "Forces", "Free nodes"])

# Finding bounding box and setting aspect ratio
xyz = np.array([particle.x for particle in PS.particles])
bb = [np.ptp(i) for i in xyz.T]
ax.set_box_aspect(bb)

plt.title("Final state")
plt.show()

# %% Plotting Initial (black) and Final (red) States with Connectivity Lines

fig = plt.figure()
ax = fig.add_subplot(111, projection="3d")

# Plot initial state in black
for particle in initial_positions:
    pos = particle[0]  # particle.x is the position
    ax.scatter(pos[0], pos[1], pos[2], color="black", marker="o", s=10)

# Draw connectivity lines for initial positions (black)
for connection in connectivity_matrix:
    p1 = initial_positions[connection[0]][0]
    p2 = initial_positions[connection[1]][0]
    line = np.column_stack((p1, p2))
    ax.plot(line[0], line[1], line[2], color="black")

# Plot final state in red
for particle in final_positions:
    pos = particle[0]  # particle.x is the position
    ax.scatter(pos[0], pos[1], pos[2], color="red", marker="o", s=10)

# Draw connectivity lines for final positions (red)
for connection in connectivity_matrix:
    p1 = final_positions[connection[0]][0]
    p2 = final_positions[connection[1]][0]
    line = np.column_stack((p1, p2))
    ax.plot(line[0], line[1], line[2], color="red")

ax.set_title("Initial (black) and Final (red) States")
ax.set_xlabel("X")
ax.set_ylabel("Y")
ax.set_zlabel("Z")

# Combine positions from both states to set an appropriate aspect ratio.
all_positions = np.array(
    [p[0] for p in initial_positions] + [p[0] for p in final_positions]
)
bb = [np.ptp(axis) for axis in all_positions.T]
ax.set_box_aspect(bb)

plt.show()
