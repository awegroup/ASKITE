# Area Consistency Report


Projected area: 19.4m2
Surface area: 24.9m2

Projection plane: `xy`. For the default `xy` plane, `z` is dropped and the result is the planform area.

Flat area is computed as the sum of two 3D triangle areas per LE/TE quadrilateral. The panel vertex order is `forward_i -> rearward_i -> rearward_j -> forward_j`.

The aero geometry is treated as the truth. The PSM area is a reduced-panel approximation because the reduced PSM wing has fewer stations than the CAD/FEM geometry.

## Summary

| Model | Source | Stations | Panels | Projected area [m2] | Projected delta vs aero | Flat area [m2] | Flat delta vs aero | Flat/projected |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Aero CAD-derived truth | `/home/jellepoland/ownCloud/phd/code/ASKITE/data/TUDELFT_V3_KITE/2D_airfoils_polars_plots_BEST/aero_geometry_CFD_CAD_derived.yaml` | 37 | 36 | 19.413150 | - | 24.231180 | - | 1.248184 |
| FEM stretched | `/home/jellepoland/ownCloud/phd/code/ASKITE/data/TUDELFT_V3_KITE/organizing_yaml_files_15_05_2026/results/struc_geometry_FEM_full_stretched.yaml` | 28 | 27 | 19.559037 | +0.145887 (+0.75%) | 24.422332 | +0.191151 (+0.79%) | 1.248647 |
| PSM reduced mapped stretched | `/home/jellepoland/ownCloud/phd/code/ASKITE/data/TUDELFT_V3_KITE/organizing_yaml_files_15_05_2026/results/struc_geometry_PSM_reduced_stretched_mapped_to_FEM_stretched.yaml` | 10 | 9 | 17.211877 | -2.201273 (-11.34%) | 21.055871 | -3.175310 (-13.10%) | 1.223334 |

## Aero CAD-derived truth Panels

| Panel | Projected area [m2] | Flat area [m2] | Flat/projected |
| --- | ---: | ---: | ---: |
| panel_01: section_01_airfoil_19 -> section_02_airfoil_18 | 0.012841 | 0.106944 | 8.328228 |
| panel_02: section_02_airfoil_18 -> section_03_airfoil_17 | 0.053482 | 0.276613 | 5.172110 |
| panel_03: section_03_airfoil_17 -> section_04_airfoil_16 | 0.070595 | 0.334993 | 4.745305 |
| panel_04: section_04_airfoil_16 -> section_05_airfoil_15 | 0.082918 | 0.370902 | 4.473136 |
| panel_05: section_05_airfoil_15 -> section_06_airfoil_14 | 0.280396 | 0.808439 | 2.883205 |
| panel_06: section_06_airfoil_14 -> section_07_airfoil_13 | 0.259552 | 0.432605 | 1.666738 |
| panel_07: section_07_airfoil_13 -> section_08_airfoil_12 | 0.337426 | 0.456609 | 1.353211 |
| panel_08: section_08_airfoil_12 -> section_09_airfoil_11 | 0.794968 | 0.961243 | 1.209160 |
| panel_09: section_09_airfoil_11 -> section_10_airfoil_10 | 0.766728 | 1.020368 | 1.330809 |
| panel_10: section_10_airfoil_10 -> section_11_airfoil_9 | 0.446872 | 0.526269 | 1.177671 |
| panel_11: section_11_airfoil_9 -> section_12_airfoil_8 | 0.504878 | 0.538106 | 1.065814 |
| panel_12: section_12_airfoil_8 -> section_13_airfoil_7 | 1.063374 | 1.094735 | 1.029492 |
| panel_13: section_13_airfoil_7 -> section_14_airfoil_6 | 1.018672 | 1.123728 | 1.103130 |
| panel_14: section_14_airfoil_6 -> section_15_airfoil_5 | 0.547220 | 0.572707 | 1.046575 |
| panel_15: section_15_airfoil_5 -> section_16_airfoil_4 | 0.576581 | 0.578504 | 1.003335 |
| panel_16: section_16_airfoil_4 -> section_17_airfoil_3 | 1.158803 | 1.159529 | 1.000626 |
| panel_17: section_17_airfoil_3 -> section_18_airfoil_2 | 1.146683 | 1.167176 | 1.017872 |
| panel_18: section_18_airfoil_2 -> section_19_airfoil_1 | 0.584587 | 0.586023 | 1.002457 |
| panel_19: section_19_airfoil_1 -> section_20_airfoil_2 | 0.584587 | 0.586023 | 1.002457 |
| panel_20: section_20_airfoil_2 -> section_21_airfoil_3 | 1.146683 | 1.167177 | 1.017872 |
| panel_21: section_21_airfoil_3 -> section_22_airfoil_4 | 1.158803 | 1.159528 | 1.000626 |
| panel_22: section_22_airfoil_4 -> section_23_airfoil_5 | 0.576581 | 0.578504 | 1.003335 |
| panel_23: section_23_airfoil_5 -> section_24_airfoil_6 | 0.547220 | 0.572707 | 1.046575 |
| panel_24: section_24_airfoil_6 -> section_25_airfoil_7 | 1.018672 | 1.123728 | 1.103130 |
| panel_25: section_25_airfoil_7 -> section_26_airfoil_8 | 1.063374 | 1.094734 | 1.029491 |
| panel_26: section_26_airfoil_8 -> section_27_airfoil_9 | 0.504878 | 0.538106 | 1.065814 |
| panel_27: section_27_airfoil_9 -> section_28_airfoil_10 | 0.446872 | 0.526269 | 1.177672 |
| panel_28: section_28_airfoil_10 -> section_29_airfoil_11 | 0.766728 | 1.020368 | 1.330809 |
| panel_29: section_29_airfoil_11 -> section_30_airfoil_12 | 0.794968 | 0.961240 | 1.209157 |
| panel_30: section_30_airfoil_12 -> section_31_airfoil_13 | 0.337426 | 0.456609 | 1.353211 |
| panel_31: section_31_airfoil_13 -> section_32_airfoil_14 | 0.259552 | 0.432605 | 1.666737 |
| panel_32: section_32_airfoil_14 -> section_33_airfoil_15 | 0.280396 | 0.808442 | 2.883214 |
| panel_33: section_33_airfoil_15 -> section_34_airfoil_16 | 0.082918 | 0.370905 | 4.473171 |
| panel_34: section_34_airfoil_16 -> section_35_airfoil_17 | 0.070595 | 0.335020 | 4.745679 |
| panel_35: section_35_airfoil_17 -> section_36_airfoil_18 | 0.053482 | 0.276738 | 5.174441 |
| panel_36: section_36_airfoil_18 -> section_37_airfoil_19 | 0.012841 | 0.106985 | 8.331368 |

## FEM stretched Panels

| Panel | Projected area [m2] | Flat area [m2] | Flat/projected |
| --- | ---: | ---: | ---: |
| panel_01: nodes_1_2 -> nodes_3_4 | 0.012717 | 0.167536 | 13.173762 |
| panel_02: nodes_3_4 -> nodes_5_6 | 0.069485 | 0.362122 | 5.211518 |
| panel_03: nodes_5_6 -> nodes_7_8 | 0.122922 | 0.582971 | 4.742594 |
| panel_04: nodes_7_8 -> nodes_9_10 | 0.303654 | 0.826310 | 2.721221 |
| panel_05: nodes_9_10 -> nodes_11_12 | 0.623306 | 0.915004 | 1.467985 |
| panel_06: nodes_11_12 -> nodes_13_14 | 0.802161 | 0.974933 | 1.215383 |
| panel_07: nodes_13_14 -> nodes_15_16 | 0.770319 | 1.023345 | 1.328469 |
| panel_08: nodes_15_16 -> nodes_17_18 | 0.958696 | 1.062546 | 1.108324 |
| panel_09: nodes_17_18 -> nodes_19_20 | 1.067166 | 1.100084 | 1.030846 |
| panel_10: nodes_19_20 -> nodes_21_22 | 1.023091 | 1.126924 | 1.101489 |
| panel_11: nodes_21_22 -> nodes_23_24 | 1.124544 | 1.144416 | 1.017671 |
| panel_12: nodes_23_24 -> nodes_25_26 | 1.163088 | 1.163792 | 1.000605 |
| panel_13: nodes_25_26 -> nodes_27_28 | 1.152543 | 1.172927 | 1.017686 |
| panel_14: nodes_27_28 -> nodes_29_30 | 1.171651 | 1.171828 | 1.000151 |
| panel_15: nodes_29_30 -> nodes_31_32 | 1.152543 | 1.172927 | 1.017686 |
| panel_16: nodes_31_32 -> nodes_33_34 | 1.163088 | 1.163792 | 1.000605 |
| panel_17: nodes_33_34 -> nodes_35_36 | 1.124544 | 1.144416 | 1.017671 |
| panel_18: nodes_35_36 -> nodes_37_38 | 1.023091 | 1.126924 | 1.101489 |
| panel_19: nodes_37_38 -> nodes_39_40 | 1.067166 | 1.100084 | 1.030846 |
| panel_20: nodes_39_40 -> nodes_41_42 | 0.958696 | 1.062546 | 1.108324 |
| panel_21: nodes_41_42 -> nodes_43_44 | 0.770319 | 1.023345 | 1.328469 |
| panel_22: nodes_43_44 -> nodes_45_46 | 0.802161 | 0.974933 | 1.215384 |
| panel_23: nodes_45_46 -> nodes_47_48 | 0.623306 | 0.915004 | 1.467985 |
| panel_24: nodes_47_48 -> nodes_49_50 | 0.303654 | 0.826309 | 2.721219 |
| panel_25: nodes_49_50 -> nodes_51_52 | 0.122922 | 0.582979 | 4.742662 |
| panel_26: nodes_51_52 -> nodes_53_54 | 0.069485 | 0.362221 | 5.212944 |
| panel_27: nodes_53_54 -> nodes_55_56 | 0.012717 | 0.172117 | 13.533964 |

## PSM reduced mapped stretched Panels

| Panel | Projected area [m2] | Flat area [m2] | Flat/projected |
| --- | ---: | ---: | ---: |
| panel_01: nodes_1_2 -> nodes_3_4 | 0.184127 | 0.928202 | 5.041110 |
| panel_02: nodes_3_4 -> nodes_5_6 | 1.473407 | 2.293182 | 1.556381 |
| panel_03: nodes_5_6 -> nodes_7_8 | 2.470664 | 2.770824 | 1.121489 |
| panel_04: nodes_7_8 -> nodes_9_10 | 2.935462 | 2.993155 | 1.019654 |
| panel_05: nodes_9_10 -> nodes_11_12 | 3.084558 | 3.085016 | 1.000148 |
| panel_06: nodes_11_12 -> nodes_13_14 | 2.935462 | 2.993155 | 1.019654 |
| panel_07: nodes_13_14 -> nodes_15_16 | 2.470664 | 2.770824 | 1.121490 |
| panel_08: nodes_15_16 -> nodes_17_18 | 1.473407 | 2.293182 | 1.556381 |
| panel_09: nodes_17_18 -> nodes_19_20 | 0.184127 | 0.928331 | 5.041809 |
