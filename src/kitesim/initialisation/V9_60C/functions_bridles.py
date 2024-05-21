# libraries
import numpy as np
import plotly.graph_objects as go
import plotly.offline as offline
import re
import pandas as pd
from sklearn.cluster import DBSCAN  # for clustering


# %% Define the Functions


def read_obj_file_bridles_faces(file_path):

    # file_path = 'V9_10_3d.obj'
    vertices = []
    normals = []
    faces = []

    #  v 2669.676921 661.911471 25.857700
    vertex_pattern = re.compile(r"v\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)")
    normal_pattern = re.compile(r"vn\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)")
    # f 18003/11004/21081 18009/11011/21088 18250/11012/21089 18244/11005/21082
    face_pattern = re.compile(
        r"f\s + (-?\d+)/(-?\d+)/(-?\d+)\s + (-?\d+)/(-?\d+)/(-?\d+)\s + (-?\d+)/(-?\d+)/(-?\d+)\s "
    )

    bridle_pattern = re.compile(r"g\s+B4_bridle")
    bridle_pattern = re.compile(
        r"usemtl\s+mtl_18_bridle"
    )  # TODO: can probably write this as r'.bridle.'
    vt_pattern = re.compile(r"vt\s+0\s0")

    f_pattern = re.compile(r"f\s+.")

    flag = False
    with open(file_path, "r") as f:
        # for line in range(len(f)):
        for index, line in enumerate(f):

            vertex_match = vertex_pattern.match(line)
            if vertex_match and flag == True:
                # This line describes a vertex
                x, y, z = map(float, vertex_match.groups())
                vertices.append((x, y, z))

            elif flag == True:
                normal_match = normal_pattern.match(line)
                if normal_match:
                    # This line describes a vertex normal
                    x, y, z = map(float, normal_match.groups())
                    normals.append((x, y, z))

                elif f_pattern.match(line):
                    # elif not vt_pattern.match(line): #and not 'B4_bridle' and not 'mtl_18_bridle':
                    # face_match = face_pattern.match(line)
                    faces.append([num.split("/")[0] for num in line.split()[1:]])

            # Setting the flag = True once the bridle line pattern is found
            # This is done at last, to skip the first line
            # done by checking if the line contains 'bridle'
            if "bridle" in re.split(r"-|,|_|\n", line):
                flag = True
                first_vertex_index = index + 2

    # the faces list doesn't start at the correct vertex index
    # so some transformation is necessary
    # the first_vertex_index is known, which should correspond with the first faces item
    # i.e. the int(faces[0][0]) = first_vertex_index
    alteration_of_index = -int(faces[0][0])

    # getting each string item to become an integer and transforming the index
    new_faces = []
    for face in faces:
        new_face = []
        for vertex in face:
            new_face.append(int(vertex) + alteration_of_index)
        new_faces.append(new_face)

    # this the vertices connectivity matrix
    # looks like {0 : (1,3,2)}
    #       { vertex_index: indices of vertices that it is connected to}
    bridle_vertices_conn = {}
    for i in range(len(new_faces)):  # looping through each face
        bridle_vertices_conn[i] = new_faces[i][1:0]  # making a connectivity dictionairy

    return vertices, normals, new_faces


def split_for_proximity(points, eps, min_samples):
    """
    eps = (float) distance tolerance
    min_samples = (ithing, non-float) number of points that you want to form the cluster
    points = (list) list of points
    """

    db = DBSCAN(eps=eps, min_samples=min_samples, metric="euclidean").fit(points)
    labels = db.labels_
    unique_labels = set(labels)
    clusters = []
    for label in unique_labels:
        if label != -1:
            indices = np.where(labels == label)
            cluster = np.take(points, indices, axis=0)
            clusters.append(cluster[0])

    return clusters


def get_indices_from_proximity_list(proximity_list, vertices):
    proximity_list_index = []

    for pointcloud in proximity_list:
        pointcloud_index = []

        for point in pointcloud:
            for vertex_index, vertex in enumerate(vertices):
                if (
                    np.allclose(point, vertex, rtol=1e-6)
                    and vertex_index not in pointcloud_index
                ):
                    pointcloud_index.append(vertex_index)

        proximity_list_index.append(pointcloud_index)
    return proximity_list_index


def get_bridle_connectivity(faces_raw, proximity_list_index):
    bridle_connectivity = []
    bridle_connectivity_i = []
    bridle_connectivity_j = []
    pointcloud_a, pointcloud_b = 0, 0  # initialising the pointcloud indices

    for face in faces_raw:  # looping through all the faces
        if len(face) == 4:  # if the face has 4 vertices
            flag_a, flag_b = False, False  # initialising flags as false

            # finding the pointcloud index that contains the first vertex of the face
            for i, pointcloud_index_list in enumerate(proximity_list_index):
                if face[0] in pointcloud_index_list:
                    pointcloud_a = i
                    flag_a = True
                    break

            # finding another pointcloud index that contains the other vertices of the face
            for i, pointcloud_index_list in enumerate(proximity_list_index):
                if (
                    any(item in pointcloud_index_list for item in face[1:])
                    and i != pointcloud_a
                    and flag_a == True
                ):
                    pointcloud_b = i
                    flag_b = True
                    break

            if flag_b and [pointcloud_a, pointcloud_b] not in bridle_connectivity:
                bridle_connectivity_i.append(pointcloud_a)
                bridle_connectivity_j.append(pointcloud_b)
                bridle_connectivity.append([pointcloud_a, pointcloud_b])

    return bridle_connectivity, bridle_connectivity_i, bridle_connectivity_j


def get_midpoints_of_proximity_list(proximity_list):

    midpoints = []
    for pointcloud in proximity_list:
        midpoints.append(np.mean(pointcloud, axis=0))
    return midpoints


def transform_coordinate_system_Surfplan_to_Model(point_list):
    point_list_new = []
    for point in point_list:
        x_new = -point[2]
        y_new = -point[0]
        z_new = point[1]
        point_list_new.append([x_new, y_new, z_new])
    return point_list_new


def get_dataframe_bridle(
    midpoints_surfplan, midpoints_model, bridle_connectivity_i, bridle_connectivity_j
):

    bridle_midpoints_Surfplan = midpoints_surfplan
    bridle_midpoints_Model = midpoints_model
    df_surfplan = pd.DataFrame(
        bridle_midpoints_Surfplan, columns=["x_Surfplan", "y_Surfplan", "z_Surfplan"]
    )
    df_model = pd.DataFrame(
        bridle_midpoints_Model, columns=["x_Model", "y_Model", "z_Model"]
    )
    df_bridle = pd.concat([df_surfplan, df_model], axis=1)

    df_conn_i = pd.DataFrame(bridle_connectivity_i, columns=["conn_i"])
    df_conn_j = pd.DataFrame(bridle_connectivity_j, columns=["conn_j"])
    df_bridle_conn = pd.concat([df_conn_i, df_conn_j], axis=1)

    return df_bridle, df_bridle_conn


# %% Maybe usefull for later
def main_bridles(filename, proximity_tol=10, proximity_min_points=5):

    # (1) Reading out the bridle part of the obj file
    obj_raw = read_obj_file_bridles(filename)
    vertices = np.array(obj_raw[0])
    faces_raw = obj_raw[2]  # extracting the vertices that form a face
    # get_3D_plot(vertices, 'red',True,1)

    # (2) Splitting the vertices on proximity
    proximity_list = split_for_proximity(vertices, proximity_tol, proximity_min_points)

    # (3) Getting the indices of the proximity list
    proximity_list_index = get_indices_from_proximity_list(proximity_list, vertices)

    # (4) Getting the connectivity of the bridle
    bridle_connectivity, bridle_connectivity_i, bridle_connectivity_j = (
        get_bridle_connectivity(faces_raw, proximity_list_index)
    )

    # (5) Getting the midpoitns of eahc bridle
    midpoints_surfplan = get_midpoints_of_proximity_list(proximity_list)
    # get_3D_plot_with_lines(midpoints_surfplan,bridle_connectivity_i,bridle_connectivity_j,show_legend=False)

    # (6) Transforming the coordinate system from Surfplan to Model
    midpoints_model = transform_coordinate_system_Surfplan_to_Model(midpoints_surfplan)
    get_3D_plot_with_lines(
        midpoints_model, bridle_connectivity_i, bridle_connectivity_j, show_legend=False
    )

    # (7) Storing the data in a dataframe
    df_bridle, df_bridle_conn = get_dataframe_bridle(
        midpoints_surfplan,
        midpoints_model,
        bridle_connectivity_i,
        bridle_connectivity_j,
    )

    return df_bridle, df_bridle_conn
