# TUDELFT V3 Kite Geometry Organization

This folder documents how the TUDELFT V3 kite structural geometry is assembled
from Surfplan/CAD wing geometry and measured bridle line data.

## Resulting geometric differences


The reference geometry values used here are: flat surface area `24.9 m2`,
projected area `19.4 m2`, maximum chord `2.6 m`, height `2.8 m`, and
width/span `8.3 m`. Projected area is computed in the `x-y` plane; flat area is
computed from 3D LE/TE quadrilateral panels.

| Geometry                    | Source                                                                                                                                       | Stations | Projected area [m2] |  Flat area [m2] |  Max chord [m] |    Height [m] |      Span [m] |
| --------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------- | -------: | ------------------: | --------------: | -------------: | ------------: | ------------: |
| Reference                   | given values                                                                                                                                 |        - |              19.400 |          24.900 |          2.600 |         2.800 |         8.300 |
| CAD/CFD aero truth          | [`aero_geometry_CFD_CAD_derived.yaml`](../2D_airfoils_polars_plots_BEST/aero_geometry_CFD_CAD_derived.yaml)                                  |       37 |      19.413 (+0.1%) |  24.231 (-2.7%) |  2.618 (+0.7%) | 2.734 (-2.3%) | 8.274 (-0.3%) |
| FEM full                    | [`final_files/struc_geometry_FEM_full.yaml`](final_files/struc_geometry_FEM_full.yaml)                                                       |       28 |      19.559 (+0.8%) |  24.422 (-1.9%) |  2.629 (+1.1%) | 2.821 (+0.7%) | 8.305 (+0.1%) |
| PSM reduced                 | [`final_files/struc_geometry_PSM_reduced.yaml`](final_files/struc_geometry_PSM_reduced.yaml)                                                 |       10 |     17.212 (-11.3%) | 21.056 (-15.4%) | 2.328 (-10.5%) | 2.728 (-2.6%) | 8.282 (-0.2%) |
| PSM photogrammetry adjusted | [`final_files/struc_geometry_PSM_reduced_photogrammetry_adjusted.yaml`](final_files/struc_geometry_PSM_reduced_photogrammetry_adjusted.yaml) |       10 |     17.212 (-11.3%) | 21.056 (-15.4%) | 2.328 (-10.5%) | 2.728 (-2.6%) | 8.282 (-0.2%) |

The PSM wing forces must therefore be scaled for the reduced planform: the PSM
has about `11%` lower projected area than the reference, corresponding to an
upward force scale of about `1.13` if matching the `19.4 m2` projected area
exactly.


## Main Pipeline

The chains below describe the intended source of each part of the geometry.
Files in `final_files/` are the current accepted structural YAMLs. Files in
`results/` are generated working outputs used to inspect, stretch, compare, or
prepare those final geometries.

### Wing Geometry

**FEM**

[data/Surfplan_files/TUDELFT_V3_KITE.txt](data/Surfplan_files/TUDELFT_V3_KITE.txt)
-> SurfplanAdapter's `scripts/process_surfplan_files.py`
-> `struc_geometry_all_in_surfplan_merged_notes.yaml`
-> manual tweaking
-> [final_files/struc_geometry_FEM_full.yaml](final_files/struc_geometry_FEM_full.yaml)

**PSM**

[final_files/struc_geometry_FEM_full.yaml](final_files/struc_geometry_FEM_full.yaml)
-> use the wing points associated with bridle line attachments
-> keep only the most forward attachment points needed for the reduced PSM wing
-> [final_files/struc_geometry_PSM_reduced.yaml](final_files/struc_geometry_PSM_reduced.yaml)

**PSM photogrammetry adjusted**
-> Reduced all TE lengths by 5%, as found through photogrammetry being a consistent difference between CAD geometry and kite shape in varying flight maneuvers. 

### Bridle Geometry

**FEM**

[data/bridle_line_measurements.yaml](data/bridle_line_measurements.yaml)
-> slightly modified from the real measurements by introducing `br_5_equiv` and
`a5_equiv`
-> added to the wing description in
[final_files/struc_geometry_FEM_full.yaml](final_files/struc_geometry_FEM_full.yaml)
-> stretched with [stretching_lines.py](stretching_lines.py)
-> final FEM bridle representation.

**PSM**

[final_files/struc_geometry_FEM_full.yaml](final_files/struc_geometry_FEM_full.yaml)
-> reduced by introducing `A1_equiv` through `A4_equiv`, plus new definitions
for `br_5_equiv` and `a5_equiv`
-> stretched with [stretching_lines.py](stretching_lines.py)
-> mapped with [map_psm_wing_nodes_to_fem_wing_nodes.py](map_psm_wing_nodes_to_fem_wing_nodes.py)
so the reduced PSM wing nodes follow the FEM attachment-derived wing geometry.
After nodal mapping, [adjust_wing_element_lenghts_to_nodal_position.py](adjust_wing_element_lenghts_to_nodal_position.py)
can be used to recompute `wing_elements.l0` from the current node positions.

**PSM photogrammetry adjusted**
-> Consistent with PSM for now.

### Mass

**Wing**

The wing mass distribution is handled by SurfplanAdapter and follows the nodal
discretisation used in each structural model. FEM and PSM therefore inherit the
mass distribution associated with their own discretisation.

**Bridles**

The bridle line density is tuned so the model matches the measured bridle line
mass. The measurement target and density bookkeeping are reported in
[data/bridle_line_measurements.yaml](data/bridle_line_measurements.yaml).
Use [check_full_kite_masses.py](check_full_kite_masses.py) to re-check masses
after changing a final YAML file.

Mass summary for the current files in `final_files/`, generated with
[check_full_kite_masses.py](check_full_kite_masses.py):

| File | Kind | Wing kg | Bridle kg | Line kg | Pulley kg | Pulleys | Kite kg | Kite+KCU kg | rho Dyneema | rho for 2.6 kg bridle |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `struc_geometry_FEM_full.yaml` | FEM full | 11.000 | 2.600 | 1.940 | 0.660 | 20 | 13.600 | 22.000 | 1499 | 1499 |
| `struc_geometry_PSM_reduced.yaml` | PSM reduced | 10.993 | 2.600 | 2.060 | 0.540 | 6 | 13.592 | 21.992 | 1717 | 1717 |
| `struc_geometry_PSM_reduced_photogrammetry_adjusted.yaml` | PSM reduced | 10.993 | 2.600 | 2.060 | 0.540 | 6 | 13.592 | 21.992 | 1717 | 1717 |

Reference values used by the checker are `11.000 kg` wing mass without
bridles, `2.600 kg` full bridle system mass, and `13.600 kg` full kite mass
excluding KCU.

---

## Important Files

| File | Purpose |
| --- | --- |
| [data/Surfplan_files/TUDELFT_V3_KITE.txt](data/Surfplan_files/TUDELFT_V3_KITE.txt) | Surfplan source file for the wing geometry. |
| [data/bridle_line_measurements.yaml](data/bridle_line_measurements.yaml) | Bridle measurements, equivalent-line definitions, density assumptions, and mass target data. |
| [final_files/struc_geometry_FEM_full.yaml](final_files/struc_geometry_FEM_full.yaml) | Current accepted full FEM structural geometry. |
| [final_files/struc_geometry_PSM_reduced.yaml](final_files/struc_geometry_PSM_reduced.yaml) | Current accepted reduced PSM structural geometry. |
| [final_files/struc_geometry_PSM_reduced_photogrammetry_adjusted.yaml](final_files/struc_geometry_PSM_reduced_photogrammetry_adjusted.yaml) | Alternative PSM file with photogrammetry-adjusted values. |
| [results/](results/) | Generated reports, plots, stretched YAML files, and inspection outputs. |

## Scripts

| Script | Purpose |
| --- | --- |
| [stretching_lines.py](stretching_lines.py) | Stretches FEM and PSM bridle systems, writes stretched YAML files, and creates the bridle comparison report. |
| [map_psm_wing_nodes_to_fem_wing_nodes.py](map_psm_wing_nodes_to_fem_wing_nodes.py) | Maps reduced PSM wing nodes to FEM attachment-derived wing positions. |
| [adjust_wing_element_lenghts_to_nodal_position.py](adjust_wing_element_lenghts_to_nodal_position.py) | Recomputes `wing_elements.l0` from current nodal distances. The filename keeps the historical `lenghts` spelling. |
| [measure_fem_bridle_attachment_chord_lengths.py](measure_fem_bridle_attachment_chord_lengths.py) | Reports forward/rearward FEM bridle attachment points and chord distances. |
| [plot_aero_geometry_and_wing_elements.py](plot_aero_geometry_and_wing_elements.py) | Plots CAD/CFD-derived aero sections and overlays structural wing nodes. |
| [compare_aero_structural_area.py](compare_aero_structural_area.py) | Compares projected and flat wing area between aero truth, FEM, and PSM geometry. |
| [check_full_kite_masses.py](check_full_kite_masses.py) | Reports wing, bridle, pulley, kite, and kite-plus-KCU masses for structural YAML files. |

## Common Commands

Regenerate stretched bridle geometries and the stretch report:

```bash
python data/TUDELFT_V3_KITE/organizing_yaml_files_15_05_2026/stretching_lines.py
```

Map the reduced PSM wing nodes to the FEM-derived wing geometry:

```bash
python data/TUDELFT_V3_KITE/organizing_yaml_files_15_05_2026/map_psm_wing_nodes_to_fem_wing_nodes.py
```

Update PSM wing-element rest lengths after changing nodal positions:

```bash
python data/TUDELFT_V3_KITE/organizing_yaml_files_15_05_2026/adjust_wing_element_lenghts_to_nodal_position.py
```

Check masses for the current final files:

```bash
python data/TUDELFT_V3_KITE/organizing_yaml_files_15_05_2026/check_full_kite_masses.py \
  data/TUDELFT_V3_KITE/organizing_yaml_files_15_05_2026/final_files/struc_geometry_FEM_full.yaml \
  data/TUDELFT_V3_KITE/organizing_yaml_files_15_05_2026/final_files/struc_geometry_PSM_reduced.yaml
```

## Notes And TODOs

- FEM bridle attachment projection is a model-building step, not a raw YAML
  edit. Direct YAML plots and FEM-built plots can therefore differ at bridle
  attachment nodes.
- Area consistency is checked by [compare_aero_structural_area.py](compare_aero_structural_area.py).
  The default projected area is the `x-y` planform area after dropping `z`;
  flat area is the sum of two 3D triangle areas per LE/TE quadrilateral.
- The PSM trailing-edge and diagonal `wing_elements.l0` values should be kept
  consistent with the current point-to-point nodal distances, or explicitly
  reported as such. This is a modeling compromise: the real canopy is curved
  and can deform over a surface, while the reduced PSM represents these members
  as straight links between nodes.
- The old CAD-derived `-1.0 deg` pitch tilt remains in the FEM file and
  therefore propagates into the mapped PSM wing geometry.
