## Data content
The data consists of input and output.
In input `paths.yaml` describes the paths to the necessary **input and output** files.

## V3 data structure ##
For the V3 kite the `points.npy` file is read which contains all the structural nodes that are used for the calculation procedure.
These nodes are stored as np.array in x,y,z and are structured as follows:

| indices | python | Description | 
| --- | --- | --- |  
| 0 | points[0] | Bridle Point |
| (1,2..,20) | points[1:21] | Wing Nodes |
| (21,22,..37) | points[21:-1] | Bridle Nodes |

</details>
<details>
<summary>2. Wing Nodes </summary>

The wing nodes consists of the bridle Line Attachment Points(LAPs), of which one is located at the Leading-edge and one at the Trailing-edge of the kite. The indices are shown here.

<img src="wing_nodes_indices.png" alt="wing_nodes_indicese" width="500">

</details>

</details>
<details>
<summary>2. Bridle Nodes </summary>

The indices of the bridle line nodes are shown here.

<img src="bridle_nodes_leading_edge_indices.png" alt="bridle_nodes_leading_edge_indices" width="500">
<img src="bridle_nodes_trailing_edge_indices.png" alt="bridle_nodes_trailing_edge_indices" width="500">

</details>

