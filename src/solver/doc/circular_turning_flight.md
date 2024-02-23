### Computing Va variation
Calculated using the `r_ratio_i = (r_0 + y_i)/r_k`. [The variables r_0 and r_k are calculated for each circular case.](https://github.com/jellepoland/kitesim/blob/1bed3a0d92e2d51a32a68daf2d4cb5fd460f9f8e/src/solver/solver_utils.py#L204C1-L230C74) 
`r_0` represents the distance from the turning-axis to the bridle point. And `r_k` represents the distance from the turning axis to the kite's CG.

Only when they are equal will the average velocity equal the nominal velocity, because otherwise the kite has rolled which makes the observed velocity deviate from that at the radial location of the bridle-point.

For y_i the 1/4chord control points for each aero-panel are used and obtained from the variable `controlpoints` from the previous iteration (therefore not computed for the first-iteration).

Using `config.is_with_varying_va` one can toggle this functionality on/off. 
 
### Feeding the variation of Va to VSM
The variation of Va is fed into the main VSM function `calculate_force_aero_wing_VSM`, which has the following sub-functions that use Va.

1. Geometry creation: `create_geometry_LEI` uses Uinf to define the orientation of the trailing vortices. Because these shouldn't cross, the **nominal value `Uinf`** is used. 

2. Solver-function, depending on whether the `r_ratio_array` is None or filled a different version of `solve_lifting_line_system_matrix_approach_semiinfinite_new` is called. Inside `Uinf` is used by 
   - Sub-function `velocity_induced_single_ring_semiinfinite` that only uses `Uinf`(=Va) to calculate the strength of the vortex core and can therefore use the **nominal value `Uinf`**.
   - [Uses Uinf to calculate the relative velocity Urel](https://github.com/jellepoland/kitesim/blob/1bed3a0d92e2d51a32a68daf2d4cb5fd460f9f8e/src/aerodynamic/VSM.py#L2588C1-L2589C59) , for which the **local value `Uinf_local`** is calculated and used.[(here as well)](https://github.com/jellepoland/kitesim/blob/1bed3a0d92e2d51a32a68daf2d4cb5fd460f9f8e/src/aerodynamic/VSM.py#L2675C1-L2675C59)