# Convergence

Coupled-solver convergence diagnostics for an existing ASKITE case, optionally
with a rerun. Output defaults to `results/ch9/convergence/`.

## Analyze Existing Case

```bash
python examples/ch9/convergence/analyzing_coupled_convergence.py \
  --case-dir results/ch9/effect_of_trim/va_va_0175_udp_0350 \
  --format pdf,png
```

## Rerun And Analyze

```bash
python examples/ch9/convergence/analyzing_coupled_convergence.py \
  --case-dir results/ch9/effect_of_trim/va_va_0175_udp_0350 \
  --rerun \
  --rerun-udp 0.35 \
  --tol 0.1 \
  --rerun-max-iter 750 \
  --format pdf,png
```

`--rerun-udp` rewrites the `depower_tape` `l0` directly in a copied structural
YAML for the rerun case.
