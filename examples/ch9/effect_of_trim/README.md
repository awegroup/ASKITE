# Effect Of Trim

ASKITE `V_a` and trim/depower sensitivity workflow. The scripts default to:

- Config: `data/ch9/config.yaml`
- Structural geometry: `data/ch9/struc_geometry_PSM_reduced.yaml`
- Aerodynamic geometry: `data/ch9/aero_geometry.yaml`
- Run data: `results/ch9/effect_of_trim/processed_data/`
- Plots: `results/ch9/effect_of_trim/`

```bash
python examples/ch9/effect_of_trim/run_sensitivity_analysis_va_trim.py \
  --sweep both \
  --solver-mode level_1 \
  --max-iter 750
python examples/ch9/effect_of_trim/plot_sensitivity_analysis_va_trim.py \
  --plot-format pdf,png
```

---
### Plot only

```bash
python examples/ch9/effect_of_trim/plot_sensitivity_analysis_va_trim.py \
  --plot-format pdf,png
```
