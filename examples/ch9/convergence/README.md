# Convergence

Coupled-solver convergence diagnostics for an existing ASKITE case, optionally
with a rerun.

Generated convergence summaries and rerun case data are written to
`results/ch9/convergence/processed_data/`. Figures are written directly to
`results/ch9/convergence/`.

## from 11 to 23 m/s

```bash
python examples/ch9/convergence/analyzing_coupled_convergence.py \
  --case-dir data/ch9/convergence_configs/ \
  --rerun \
  --rerun-udp 0.34 \
  --tol 0.01 \
  --rerun-max-iter 500 \
  --format pdf,png \
  --output-dir results/ch9/convergence/
```
