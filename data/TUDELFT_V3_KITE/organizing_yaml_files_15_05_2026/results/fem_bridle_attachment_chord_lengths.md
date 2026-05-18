# FEM Bridle Attachment Chord Lengths

Input geometry: `struc_geometry_FEM_full_stretched.yaml`

Attachment nodes are nodes that appear in `bridle_connections` and also belong to a `strut_tubes.node_indices` station. Direct `wing_particles` attachments are assigned to the nearest station by y-coordinate.

| Station | Attachment nodes sorted by x | Most forward node [x,y,z] | Most rearward node [x,y,z] | Distance [m] |
| --- | --- | --- | --- | ---: |
| le_tip_start | 3, 120 | 3 [-0.183243, 4.085799, 8.602814] | 120 [0.789240, 4.140976, 8.378346] | 0.999577 |
| strut_1 | 103, 109, 121 | 103 [-0.692027, 3.967853, 9.239621] | 121 [0.839743, 3.956753, 9.236237] | 1.531813 |
| strut_2 | 96, 99, 133, 112, 122 | 96 [-0.906897, 3.134159, 10.244571] | 122 [1.086881, 3.118511, 10.224542] | 1.993940 |
| strut_3 | 92, 97, 132, 111, 123 | 92 [-1.050409, 1.959598, 10.846579] | 123 [1.169142, 1.947901, 10.813494] | 2.219829 |
| strut_4 | 91, 94, 131, 110, 124 | 91 [-1.122026, 0.664217, 11.106687] | 124 [1.205432, 0.661074, 11.066604] | 2.327805 |
| strut_5 | 57, 60, 128, 76, 90 | 57 [-1.122026, -0.664217, 11.106687] | 90 [1.205432, -0.661074, 11.066604] | 2.327805 |
| strut_6 | 58, 63, 129, 77, 89 | 58 [-1.050409, -1.959598, 10.846579] | 89 [1.169142, -1.947901, 10.813494] | 2.219829 |
| strut_7 | 62, 65, 130, 78, 88 | 62 [-0.906897, -3.134159, 10.244571] | 88 [1.086881, -3.118511, 10.224542] | 1.993940 |
| strut_8 | 69, 75, 87 | 69 [-0.692027, -3.967853, 9.239621] | 87 [0.839743, -3.956753, 9.236237] | 1.531813 |
| le_tip_end | 53, 86 | 53 [-0.183243, -4.085799, 8.602814] | 86 [0.789240, -4.140976, 8.378346] | 0.999577 |

## Attachment Details

| Station | Node | Source table | Position [x,y,z] | Bridle connections |
| --- | ---: | --- | --- | --- |
| le_tip_start | 3 | wing_particles | [-0.183243, 4.085799, 8.602814] | a5_equiv |
| le_tip_start | 120 | bridle_particles | [0.789240, 4.140976, 8.378346] | br_5_equiv |
| strut_1 | 103 | bridle_particles | [-0.692027, 3.967853, 9.239621] | a4_b4 |
| strut_1 | 109 | bridle_particles | [-0.588761, 3.967105, 9.239393] | a4_b4 |
| strut_1 | 121 | bridle_particles | [0.839743, 3.956753, 9.236237] | br_4 |
| strut_2 | 96 | bridle_particles | [-0.906897, 3.134159, 10.244571] | a3_b3 |
| strut_2 | 99 | bridle_particles | [-0.772484, 3.133104, 10.243221] | a3_b3 |
| strut_2 | 133 | bridle_particles | [-0.557434, 3.131416, 10.241060] | c3 |
| strut_2 | 112 | bridle_particles | [-0.257239, 3.129060, 10.238045] | d3 |
| strut_2 | 122 | bridle_particles | [1.086881, 3.118511, 10.224542] | br_3 |
| strut_3 | 92 | bridle_particles | [-1.050409, 1.959598, 10.846579] | a2_b2 |
| strut_3 | 97 | bridle_particles | [-0.900776, 1.958809, 10.844349] | a2_b2 |
| strut_3 | 132 | bridle_particles | [-0.701544, 1.957759, 10.841379] | c2 |
| strut_3 | 111 | bridle_particles | [-0.327184, 1.955786, 10.835799] | d2 |
| strut_3 | 123 | bridle_particles | [1.169142, 1.947901, 10.813494] | br_2 |
| strut_4 | 91 | bridle_particles | [-1.122026, 0.664217, 11.106687] | a1_b1 |
| strut_4 | 94 | bridle_particles | [-0.965119, 0.664005, 11.103985] | a1_b1 |
| strut_4 | 131 | bridle_particles | [-0.743317, 0.663705, 11.100165] | c1 |
| strut_4 | 110 | bridle_particles | [-0.363641, 0.663193, 11.093626] | d1 |
| strut_4 | 124 | bridle_particles | [1.205432, 0.661074, 11.066604] | br_1 |
| strut_5 | 57 | bridle_particles | [-1.122026, -0.664217, 11.106687] | a1_b1 |
| strut_5 | 60 | bridle_particles | [-0.965119, -0.664005, 11.103985] | a1_b1 |
| strut_5 | 128 | bridle_particles | [-0.743317, -0.663705, 11.100165] | c1 |
| strut_5 | 76 | bridle_particles | [-0.363641, -0.663193, 11.093626] | d1 |
| strut_5 | 90 | bridle_particles | [1.205432, -0.661074, 11.066604] | br_1 |
| strut_6 | 58 | bridle_particles | [-1.050409, -1.959598, 10.846579] | a2_b2 |
| strut_6 | 63 | bridle_particles | [-0.900776, -1.958809, 10.844349] | a2_b2 |
| strut_6 | 129 | bridle_particles | [-0.701544, -1.957759, 10.841379] | c2 |
| strut_6 | 77 | bridle_particles | [-0.327184, -1.955786, 10.835799] | d2 |
| strut_6 | 89 | bridle_particles | [1.169142, -1.947901, 10.813494] | br_2 |
| strut_7 | 62 | bridle_particles | [-0.906897, -3.134159, 10.244571] | a3_b3 |
| strut_7 | 65 | bridle_particles | [-0.772484, -3.133104, 10.243221] | a3_b3 |
| strut_7 | 130 | bridle_particles | [-0.557434, -3.131416, 10.241060] | c3 |
| strut_7 | 78 | bridle_particles | [-0.257239, -3.129060, 10.238045] | d3 |
| strut_7 | 88 | bridle_particles | [1.086881, -3.118511, 10.224542] | br_3 |
| strut_8 | 69 | bridle_particles | [-0.692027, -3.967853, 9.239621] | a4_b4 |
| strut_8 | 75 | bridle_particles | [-0.588761, -3.967105, 9.239393] | a4_b4 |
| strut_8 | 87 | bridle_particles | [0.839743, -3.956753, 9.236237] | br_4 |
| le_tip_end | 53 | wing_particles | [-0.183243, -4.085799, 8.602814] | a5_equiv |
| le_tip_end | 86 | bridle_particles | [0.789240, -4.140976, 8.378346] | br_5_equiv |
