# Solver Package

## Files Structure
| File/Directory | Content |
| -------------- | ------- |
| solver_main.py | The simulation loop is formulated and executed here |
| solver_utils.py | Contains helper functions to run the simulation-loop. |
| euler_attempts.py | Attempt at writing an Euler-Solver |
| kdr_attempt.py | Attempt at writing a KinematicDynamicRelaxation solver |

# `solver_main.py`

### First more initialization is done
- Initialise `position` data frame; rows are time-stamps, and columns represent positions like x1,y1,z1,x2,y2,z2,...
- Depending on the case, depicted by some booleans (is_with_velocity_initialization, is_with_vk_optimization,is_circular_case), more initialization occurs.   

## Simulation loop
A for loop runs over each time step, each row in the `positions` data frame

If the kite-velocity must be optimized, this is done using: `solver_utils.calculate_fx`, `tether_aero.calculate_tether_drag` and `solver_utils.optimalisation_of_vk_for_fx_0`

### External Force calculation
1. The centrifugal force is calculated if a circular case is considered.
2. The structural mesh input mesh is converted to aerodynamic mesh through `coupling_struc2aero.order_struc_nodes_right_to_left`
3. Wing aerodynamic force is calculated using `VSM.calculate_force_aero_wing_VSM`
4. The aerodynamic force is distributed back to the structural nodes, using `coupling_aero2struc.aero2struc`
5. The bridle aerodynamic force is calculated `bridle_line_system_aero.calculate_force_aero_bridle_thedens2022`
6. The total external force is calculated by summing up: force_aero_wing + force_aero_bridle + force_gravity + force_centrifugal

### Structural Calculation, finding the displacement
With the external force known on the structural nodes, the displacement is calculated by feeding the flattened external force vector to the particle-system model: `psystem.kin_damp_sim`.

### Checking-Convergence
- the residual forces are calculated
- the changes in depower- and steering tape lengths if the kite was actuated during the last loop
- A series of checks is done sequentially to ensure all conditions for convergence are satisfied.
  1. residual forces are checked to be below `params["aerostructural_tol"]`
  2. residual forces are checked to be not NaN
  3. The number of iterations is checked
- A second series of checks is done to make sure all the initialization steps have been completed
  1. the velocity initialization is finalized
  2. the depower tape is at its final length; changes are only allowed if the residual is low enough
  3. the steering tape is at its final length, changes are only allowed if the residual is low enough
  4. if all initialization steps are finalized and the residual is low enough, convergence is reached, and a `break` statement is used to stop the loop.

### Finalizing
- If necessary, the kite velocity is once more optimized
- the post-processing variables are defined
