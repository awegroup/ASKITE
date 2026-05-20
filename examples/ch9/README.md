# Chapter 9 Scripts

Run commands from the ASKITE repository root. The workflows are split by topic:

- `examples/ch9/convergence/` for coupled-solver convergence diagnostics.
- `examples/ch9/effect_of_trim/` for `V_a` and trim/depower sensitivity.
- `examples/ch9/force_validation/` for ASKITE versus EKF force validation.
- `examples/ch9/shape_validation/` for final wing-shape validation plots.

All workflows use ASKITE inputs from `data/ch9/` by default. The only external
input is the EKF case/statistics folder used by force validation:

```text
/home/jellepoland/ownCloud/phd/code/EKF-AWE/data/ch9_3_2_straight_vwt/
```

Shared helpers remain in `examples/ch9/ch9_analysis_utils.py`.

## Drag Accounting

- Bridle drag is applied when `is_with_aero_bridle` is true. The cable force is
  computed in `src/kitesim/aerodynamic_bridle_line_drag.py` and added to the
  aerodynamic force totals in the coupled solvers.
- KCU drag is not currently applied as a solver force. The force-validation
  runner records `--include-kcu-drag` as requested metadata only.
- Tether drag is not directly added to ASKITE level-1 or QSM external forces.
  The QSM path creates an AWETrim tether model, but `f_tether_drag` remains a
  zero placeholder in the ASKITE force assembly.
