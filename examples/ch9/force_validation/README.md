# Force Validation

ASKITE force-coefficient validation against the EKF straight-flight bins. The
ASKITE scripts default to:

- Config: `data/ch9/config.yaml`
- Structural geometry: `data/ch9/struc_geometry_PSM_reduced.yaml`
- Aerodynamic geometry: `data/ch9/aero_geometry.yaml`
- Run data: `results/ch9/force_validation/processed_data/`
- Plots: `results/ch9/force_validation/`

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
python examples/ch9/force_validation/plot_ch9_3_2_validation.py \
  --ekf-dir /home/jellepoland/ownCloud/phd/code/EKF-AWE/data/ch9_3_2_straight_vwt \
  --format pdf,png
```


## Plot only

```bash
python examples/ch9/force_validation/plot_ch9_3_2_validation.py \
  --ekf-dir /home/jellepoland/ownCloud/phd/code/EKF-AWE/data/ch9_3_2_straight_vwt \
  --format pdf,png
```
