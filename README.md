<p><a target="_blank" href="https://app.eraser.io/workspace/TpUXeyQRIqink6xG9qgk" id="edit-in-eraser-github-link"><img alt="Edit in Eraser" src="https://firebasestorage.googleapis.com/v0/b/second-petal-295822.appspot.com/o/images%2Fgithub%2FOpen%20in%20Eraser.svg?alt=media&amp;token=968381c8-a7e7-472a-8ed6-4a6626da5501"></a></p>

Github should automatically render a TOC?

## Purpose
This code is a design-tool, to be used in the kite-design process for finding performance characteristics of soft-wing kites without the need for flying them.

![AI-generated-AWE_illustration](doc/images/AI_generated_AWE.png)

# Files Structure
| File/Directory | Content |
| -------------- | ------- |
| data           | Input data and output/results of the code. |
| doc            | Currently contains an attempt at making Sphinx documentation. |
| Pipfile        | Dependencies. |
| pycache        | Description of Python byte code cached for faster loading. |
| README.md      | Description of the project/repository overview. |
| pyvenv.cfg     | Virtual environment configuration. |
| requirements.txt | Description of project dependencies. |
| scripts       | Python scripts and Jupyter-notebooks. |
| settings       | config.yaml file with all settings. |
| src            | The Python Packages (see detailed list below). |
| test           | Attempt at writing test functions. |
| venv_phd       | Virtual environment files. |


## More details per folder
[data](data/doc/README.md)

## Detailed List for `src` Directory
The `src` directory contains various Python packages involved in the project:
<!-- - [coconut](doc/coconut.md) - Description or purpose of the `coconut` package. -->
- [aerodynamic](src/aerodynamic/doc/aerodynamic.md) - Description or purpose of the `aerodynamic` package.
- [coupling](src/coupling.doc/coupling.md) - Description or purpose of the `coupling` package.
- [initialisation](src/initialisation/doc/initialisation.md) - Description or purpose of the `initialisation` package.
- [particleSystem](src/particleSystem/doc/particleSystem.md) - Description or purpose of the `particleSystem` package.
- [post_processing](src/post_processing/doc/post_processing.md) - Description or purpose of the `post_processing` package.
- [solver](src/solver/doc/solver.md) - Description or purpose of the `solver` package.
- [structural](src/structural/doc/structural.md) - Description or purpose of the `structural` package.

# Performing a simulation

## Understanding the code structure
The `scripts.py` folder contains python files and jupyter-notebooks that run a simulation. Each simulation undergoes roughly the same steps. Here the steps for the `main.ipynb` are discussed:

<details>
<summary>1. Initialisation</summary>

1. **Setting up the environment**: The code starts by setting up the environment for autoreloading and defining the working directory. It then appends the current working directory to the system path.
2. **Importing necessary modules**: Various modules are imported from different packages such as `src.initialisation`, `src.particleSystem`, `src.coupling`, `src.structural`, `src.solver`, `src.post_processing`, `src.aerodynamic`, and `test`. Other necessary Python libraries are also imported.
3. **Loading configuration**: The configuration for the simulation is loaded from a YAML file using the `config` dataclass from `src.initialisation.yaml_loader`.
4. **Initializing mutable variables**: Mutable variables such as `points` and `vel_app` are initialized.
5. **Defining connectivity matrix and parameters**: The connectivity matrix and parameters for the particle system are defined using the `define_connectivity_matrix` and `define_params` functions from `input_particleSystem`.
6. **Defining initial conditions**: The initial conditions for the kite are defined using the `define_initial_conditions_kite` function from `input_particleSystem`.
7. **Setting up rotational resistance (if applicable)**: If the kite model is "V9_60C", rotational resistance is set up using the `extract_rotational_resistances_dicts` and `initialize_bending_spring` functions from `particles_with_rotational_resistance`.
8. **Creating the ParticleSystem object**: A `ParticleSystem` object is created with the connectivity matrix, initial conditions, and parameters.
9. **Printing initial dimensions**: The initial dimensions of the kite, such as the scaling factor, reference chord, wing span, wing height, wing area, and projected area, are printed.

</details>
<details>
<summary>2. Aero-Structural iterations</summary>

The simulation process involves several steps:

1. **Setting up simulation parameters**: Various parameters for the simulation are set, such as the simulation name, whether to perform VK optimization, whether it's a circular case, whether to run only one time step, whether to print intermediate results, and whether to include gravity.

2. **Running the aerostructural [solver](doc/markdown_files/src/solver.md)**: The `run_aerostructural_solver` function from `solver_main` is called with the initial points, apparent velocity, particle system, parameters dictionary, configuration, input VSM, input bridle aero, and the previously set simulation parameters. This function returns the final points, print data, plot data, and animation data.

</details>

<details>
<summary> 3. Post-processing</summary>

The process involves several steps:

1. **Loading data**: A list of data names is defined, and then each corresponding data file is loaded from the specified directory using the `dill` library. The loaded data is stored in a dictionary, and then each data item is assigned to its respective variable.

2. **Post-processing**: Several flags are set to determine whether to print results, plot data, animate the simulation, and save the results. Depending on these flags, the corresponding functions from `post_processing_main` are called.

3. **Saving data**: If the `is_with_save` flag is set to `True`, the data to be saved is defined in a list of lists, where each inner list contains the data item and its name. If the specified folder doesn't exist, it is created. Then, each data item is serialized and saved to a file in the specified folder using the `dill` library.

</details>

# Installation instructions

...
Should contain instructions about how to get this as a pip package and for the developers, how to run the pip package in developer mode.

## Authors
J.A.W. Poland

>  [!WARNING]
should show a warning 

## License is not yet clear
Apache? CC-BY-2?

>  [!NOTE]
should show a note 

## How to Cite
Poland JAW, Schmehl R. Modelling Aero-Structural Deformation of Flexible Membrane Kites. Energies. 2023; 16(14):5264. [﻿https://doi.org/10.3390/en16145264](https://doi.org/10.3390/en16145264) 

## Example usage
>  [!IMPORTANT]
should show a sort of important message 

## Installation Instructions
`﻿Here we could enter some cool looking-code` 

### Depencies and their versions
## Questions
First check out the FAQs section: -link to FAQs-
For now questions open an issue: -link to issues page-

## Contribution Guide
### Trying some stuff
Short summary 

Lorem ipsum dolor sit amet, consectetur adipiscing elit, sed do eiusmod
tempor incididunt ut labore et dolore magna aliqua. Ut enim ad minim veniam,
quis nostrud exercitation ullamco laboris nisi ut aliquip ex ea commodo
consequat. Duis aute irure dolor in reprehenderit in voluptate velit esse
cillum dolore eu fugiat nulla pariatur. Excepteur sint occaecat cupidatat non
proident, sunt in culpa qui officia deserunt mollit anim id est laborum.




<!--- Eraser file: https://app.eraser.io/workspace/TpUXeyQRIqink6xG9qgk --->
