# Force Validation

ASKITE force-coefficient validation against the EKF straight-flight bins. The
ASKITE scripts default to:

- Config: `data/ch9/config.yaml`
- Structural geometry: `data/ch9/struc_geometry_PSM_reduced.yaml`
- Aerodynamic geometry: `data/ch9/aero_geometry.yaml`
- Results: `results/ch9/force_validation/`

The EKF input remains external:

```text
/home/jellepoland/ownCloud/phd/code/EKF-AWE/data/ch9_3_2_straight_vwt/
```

## Run ASKITE Cases

```bash
python examples/ch9/force_validation/run_ch9_3_2_vwt_cases.py \
  --cases-csv /home/jellepoland/ownCloud/phd/code/EKF-AWE/data/ch9_3_2_straight_vwt/ch9_3_2_vwt_cases_for_askite.csv \
  --mode apparent_wind_prescribed \
  --max-iter 750
```

The summary is written to:

```text
results/ch9/force_validation/ch9_3_2_askite_case_summary.csv
```

## Validation Plots

```bash
python examples/ch9/force_validation/plot_ch9_3_2_validation.py \
  --ekf-dir /home/jellepoland/ownCloud/phd/code/EKF-AWE/data/ch9_3_2_straight_vwt \
  --format pdf,png
```

## Quick Wiring Check

```bash
python examples/ch9/force_validation/run_ch9_3_2_vwt_cases.py \
  --max-cases 1 \
  --max-iter 1 \
  --output-root /tmp/askite_force_validation_check
```

## Drag Flags

- `--include-bridle-drag` enables ASKITE bridle-line drag.
- `--include-tether-drag` sets `is_with_aero_tether`, but no direct ASKITE
  tether-drag force is currently added to the coupled external force vector.
- `--include-kcu-drag` is recorded in output metadata only; it does not apply a
  KCU drag force.
