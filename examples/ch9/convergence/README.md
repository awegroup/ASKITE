# Convergence

Coupled-solver convergence diagnostics for an existing ASKITE case, optionally
with a rerun. Output defaults to `results/ch9/convergence/`.

## Run Fresh Simulation & Plot (1 Command)

<!-- ```bash
python examples/run_simulation_level_1.py && \
python examples/ch9/convergence/analyzing_coupled_convergence.py \
  --case-dir results/TUDELFT_V3_KITE/$(ls -t results/TUDELFT_V3_KITE/ | head -1) \
  --format pdf,png
```

This runs a fresh simulation from the default Chapter 9 config, then immediately analyzes and plots the convergence. -->

## Analyze Existing Case (No Rerun)

This requires that `results/<case_dir>/sim_output.h5` already exists from a prior run.

Minimal command:
```bash
python examples/ch9/convergence/analyzing_coupled_convergence.py \
  --case-dir results/ch9/effect_of_trim/va_va_0175_udp_0250
```

With output format options:
```bash
python examples/ch9/convergence/analyzing_coupled_convergence.py \
  --case-dir results/ch9/effect_of_trim/va_va_0175_udp_0250 \
  --format pdf,png
```

If the case directory does not contain `sim_output.h5`, use `--rerun` to generate it (see below).

## Rerun And Analyze

Re-runs the coupled solver from a previous result with optional parameter modifications,
then analyzes the new convergence:

```bash
python examples/ch9/convergence/analyzing_coupled_convergence.py \
  --case-dir results/ch9/effect_of_trim/va_va_0175_udp_0250 \
  --rerun \
  --rerun-udp 0.25 \
  --tol 0.1 \
  --rerun-max-iter 750 \
  --format pdf,png
```

**Note:** `--rerun-udp` rewrites the `depower_tape` `l0` directly in a copied structural YAML for the rerun case.
This approach can be sensitive to initial state mismatch; if convergence fails, try `--rerun` alone (without `--rerun-udp`)
to re-equilibrate from the recovered state at its original depower setting first.
