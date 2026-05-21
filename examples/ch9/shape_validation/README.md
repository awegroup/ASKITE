# Shape Validation

Single-case ASKITE shape validation against the powered photogrammetry shape.
The run script writes simulation data to
`results/ch9/shape_validation/processed_data/`; the plot script writes figures
directly to `results/ch9/shape_validation/`.

The default case is:

- `u_dp = 0.4151`
- `V_a = 16.75 m/s`

```bash
python examples/ch9/shape_validation/run_ch9_3_2_shape_validation.py \
  --max-iter 750
python examples/ch9/shape_validation/plot_ch9_3_2_wing_shape_xyz.py \
  --format pdf,png
```

The powered photogrammetry data is read from:

```text
data/ch9/shape_validation/Torque_paper_data/powered_flight_for_depowering_plot.csv
```
