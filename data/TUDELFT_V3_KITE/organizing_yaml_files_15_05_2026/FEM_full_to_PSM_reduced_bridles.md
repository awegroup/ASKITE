# FEM Full To PSM Reduced Bridles

Source files:

- FEM full YAML: `ASKITE/data/TUDELFT_V3_KITE/organizing_yaml_files_15_05_2026/struc_geometry_FEM_full.yaml`
- PSM reduced YAML: `ASKITE/data/TUDELFT_V3_KITE/organizing_yaml_files_15_05_2026/struc_geometry_PSM_reduced.yaml`
- Measurement truth for the short pulley bookkeeping lines: `TUDELFT_V3_KITE/2025_experiment/bridle_line_system/bridle_line_measurements.yaml`

This compares FEM `bridle_lines.rest_length` against PSM `bridle_elements.l0`.
Name matching is exact and case-sensitive.

## Summary

- Exact line names present in both files: 16.
- Exact-name lines with equal length: 14.
- FEM-only exact names: 19.
- PSM-only exact names: 4.
- The PSM reduced YAML now uses FEM naming where the mapping is certain, and uses `_equiv` for reduced/collapsed front-tip and fan branches.
- Aside from PSM `a5_equiv` and `br_5_equiv`, exact-name PSM lengths now match the FEM full lengths.
- PSM `a5_equiv` and `br_5_equiv` now match the measurement-derived reduced pulley-branch formulas.
- Important: exact-name matches with `_equiv` are still reduced equivalents, not necessarily the same physical FEM sub-line.

## Exact Names Present In Both

| line name | FEM full length [m] | PSM reduced length [m] | equal length? | PSM - FEM [m] |
|---|---:|---:|---|---:|
| `a5_equiv` | 0.254 | 11.928 | No | +11.674 |
| `AI` | 3.670 | 3.670 | Yes | +0.000 |
| `AII` | 3.360 | 3.360 | Yes | +0.000 |
| `A_main` | 3.903 | 3.903 | Yes | +0.000 |
| `br_1` | 4.407 | 4.407 | Yes | +0.000 |
| `br_2` | 4.180 | 4.180 | Yes | +0.000 |
| `br_3` | 4.094 | 4.094 | Yes | +0.000 |
| `br_4` | 3.550 | 3.550 | Yes | +0.000 |
| `br_5_equiv` | 0.825 | 13.818 | No | +12.993 |
| `Br_I` | 2.333 | 2.333 | Yes | +0.000 |
| `Br_II` | 2.045 | 2.045 | Yes | +0.000 |
| `Br_main_1` | 0.610 | 0.610 | Yes | +0.000 |
| `Br_main_e` | 0.070 | 0.070 | Yes | +0.000 |
| `M` | 4.856 | 4.856 | Yes | +0.000 |
| `steering_tape` | 1.600 | 1.600 | Yes | +0.000 |
| `depower_tape` | 1.900 | 1.900 | Yes | +0.000 |

Notes:

- PSM `a5_equiv` is the reduced front tip branch previously named `a6`; it matches `AIII + 2*a5_e + a5 = 11.928 m`, not the short FEM `a5_equiv` line.
- PSM `br_5_equiv` is the reduced rear tip branch previously named `br6`; it matches `Br_5 + 2*br_5_e + br_5 = 13.818 m`, not the short FEM `br_5_equiv` line.

## FEM Names Missing Exactly In PSM

These missing FEM names are expected when the detailed FEM front fan is reduced into one equivalent PSM branch per index.

| FEM full line group | FEM line lengths [m] | PSM reduced representation |
|---|---:|---|
| `A1`, `ab1_cd1`, `a1_b1`, `c1`, `d1` | 1.860, 3.261, 0.330, 0.445, 0.405 | Collapsed into `A1_equiv`. |
| `A2`, `ab2_cd2`, `a2_b2`, `c2`, `d2` | 1.798, 3.180, 0.260, 0.400, 0.375 | Collapsed into `A2_equiv`. |
| `A3`, `ab3_cd3`, `a3_b3`, `c3`, `d3` | 1.791, 3.076, 0.325, 0.400, 0.370 | Collapsed into `A3_equiv`. |
| `A4`, `a4_b4` | 2.877, 0.220 | Collapsed into `A4_equiv`. |
| `AIII` | 11.420 | Folded into `a5_equiv`, together with the tip-side `a5`/`a5_e` equivalent. |
| `Br_5` | 12.168 | Folded into `br_5_equiv`, together with the tip-side `br_5`/`br_5_e` equivalent. |

## Likely Semantic Mapping

| PSM line | PSM length [m] | likely FEM full line(s) | FEM/reference length [m] | equal length? | note |
|---|---:|---|---:|---|---|
| `A1_equiv` | 3.520 | `A1` plus collapsed lower fan pieces | 1.860 for standalone `A1` | No | PSM value is a reduced branch, not standalone FEM `A1`. |
| `A2_equiv` | 3.352 | `A2` plus collapsed lower fan pieces | 1.798 for standalone `A2` | No | PSM value is a reduced branch, not standalone FEM `A2`. |
| `A3_equiv` | 3.253 | `A3` plus collapsed lower fan pieces | 1.791 for standalone `A3` | No | PSM value is a reduced branch, not standalone FEM `A3`. |
| `A4_equiv` | 2.782 | likely `A4` / outer front branch | 2.877 for standalone `A4` | No | PSM value is a reduced branch. |
| `AI` | 3.670 | `AI` | 3.670 | Yes | Same name and length. |
| `AII` | 3.360 | `AII` | 3.360 | Yes | Same name and length. |
| `A_main` | 3.903 | `A_main` | 3.903 | Yes | Same name and length. |
| `a5_equiv` | 11.928 | `AIII` + `2*a5_e` + `a5` | 11.928 using measurement `AIII + 2*a5_e + a5` | Yes | PSM `a5_equiv` is the old `a6` reduced branch; it is not the short FEM `a5_equiv`. |
| `br_1` | 4.407 | `br_1` | 4.407 | Yes | Same name and length. |
| `br_2` | 4.180 | `br_2` | 4.180 | Yes | Same name and length. |
| `br_3` | 4.094 | `br_3` | 4.094 | Yes | Same name and length. |
| `br_4` | 3.550 | `br_4` | 3.550 | Yes | Same name and length. |
| `br_5_equiv` | 13.818 | `Br_5` + `2*br_5_e` + `br_5` | 13.818 using measurement `Br_5 + 2*br_5_e + br_5` | Yes | PSM `br_5_equiv` is the old `br6` reduced branch; it is not the short FEM `br_5_equiv`. |
| `Br_I` | 2.333 | `Br_I` | 2.333 | Yes | Same name and length. |
| `Br_II` | 2.045 | `Br_II` | 2.045 | Yes | Same name and length. |
| `Br_main_1` | 0.610 | `Br_main_1` | 0.610 | Yes | Same name and length. |
| `Br_main_e` | 0.070 | `Br_main_e` | 0.070 | Yes | Added to PSM so `M` connects through `Br_main_e` instead of directly to the main junction. |
| `M` | 4.856 | `M` | 4.856 | Yes | Same name and length. |
| `steering_tape` | 1.600 | `steering_tape` | 1.600 | Yes | Same name and length. |
| `depower_tape` | 1.900 | `depower_tape` | 1.900 | Yes | Same name and length. |
