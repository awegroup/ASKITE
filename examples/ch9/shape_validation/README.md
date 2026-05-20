# Shape Validation

Final ASKITE wing-shape plots for force-validation cases. This workflow reads
case folders produced by `examples/ch9/force_validation/` and writes plots to
`results/ch9/shape_validation/plots/` by default.

## Plot All Force-Validation Cases

```bash
python examples/ch9/shape_validation/plot_ch9_3_2_wing_shape_xyz.py \
  --overlay-initial \
  --format pdf,png
```

## Plot One Case

```bash
python examples/ch9/shape_validation/plot_ch9_3_2_wing_shape_xyz.py \
  --case-dir results/ch9/force_validation/askite_vwt_001 \
  --overlay-initial \
  --format pdf,png
```
