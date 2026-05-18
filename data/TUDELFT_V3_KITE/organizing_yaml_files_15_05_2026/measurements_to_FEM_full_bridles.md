# Measurements To FEM Full Bridles

Source files:

- Measurements: `/home/jellepoland/ownCloud/phd/code/TUDELFT_V3_KITE/2025_experiment/bridle_line_system/bridle_line_measurements.yaml`
- FEM full YAML: `ASKITE/data/TUDELFT_V3_KITE/organizing_yaml_files_15_05_2026/struc_geometry_FEM_full.yaml`


## Summary

All line lengths, and diameters are incorporated in the FEM full YAML and consistent with the measurement file, except for:
- br_5_e: which is excluded as a standalone FEM line because it is integrated into `br_5_equiv`
- a5_e: which is excluded as a standalone FEM line because it is integrated into `a5_equiv`
- `br_5_equiv`: is adjusted because `br_5` is a pulley branch in reality but is modelled as a single FEM line, with `br_5_e` added.
  - Adjusted from measured `br_5`: 1.420 m to (0.5*1.420) + 0.115 = 0.825 m.
- `a5_equiv`: is adjusted because `a5` is a pulley branch in reality but is modelled as a single FEM line, with `a5_e` added.
  - Adjusted from measured `a5`: 0.278 m to (0.5*0.278) + 0.115 = 0.254 m.

## Bridle Attachment Projection

The raw `struc_geometry_FEM_full.yaml` contains bridle attachment nodes at their measured/input positions. When the FEM model is built through `read_struc_geometry_yaml_level_2.py`, the intermediate nodes listed in each `strut_tubes.node_indices` entry are projected onto the straight chord line between the first and last node of that strut/tip section. The original and projected positions are stored in `simplified_bridle_points`, and bridle rest lengths are adjusted to account for this simplification before the FEM structure is created.

This means plots that read the YAML directly show the unprojected attachment nodes, while `stretching_lines.py` shows the projected FEM coordinates because it uses `read_struc_geometry_yaml_level_2.main(...)`.
