# libraries
import re
import numpy as np

def read_obj_file_part(file_path,part_1,part_2):
    ''' Reads the .obj file and returns the vertices, normals, and faces'''
    vertices = []
    normals = []
    faces = []
                                
    vertex_pattern = re.compile(r'v\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)')  #  v 2669.676921 661.911471 25.857700
    normal_pattern = re.compile(r'vn\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)')
    #face_pattern = re.compile(r'f\s + (-?\d+)/(-?\d+)/(-?\d+)\s + (-?\d+)/(-?\d+)/(-?\d+)\s + (-?\d+)/(-?\d+)/(-?\d+)\s ') # f 18003/11004/21081 18009/11011/21088 18250/11012/21089 18244/11005/21082
    face_pattern = re.compile(r'f\s+.')

    flag = False
    with open(file_path, 'r') as f:
        #for line in range(len(f)):
        for index, line in enumerate(f):

            ## reading out each vertex, normal and face when flag is true.    
            if flag == True:
                ## reading out each vertex, normal and face.
                vertex_match = vertex_pattern.match(line)
                normal_match = normal_pattern.match(line)
                face_match = face_pattern.match(line)
                if vertex_match:
                    x, y, z = map(float, vertex_match.groups())
                    vertices.append((x, y, z))
                elif normal_match:
                    x, y, z = map(float, normal_match.groups())
                    normals.append((x, y, z))
                elif face_match:
                    faces.append([num.split('/')[0] for num in line.split()[1:]])

            ## figuring out when to start reading out the vertices, normals and faces
            if part_1 in re.split(r"-|,|_|\n", line):    
                flag = True 

            ## figuring out when to stop reading out the vertices, normals and faces
            if part_2 in re.split(r"-|,|_|\n", line):    
                flag = False
                break

    return vertices, normals, faces

def read_obj_file_rib(file_path,part_1,part_2):
    ''' Reads the .obj file and returns the vertices, normals, and faces'''

    vertices = []
    normals = []
    faces = []
    ribs = []
    idx = 0
                                
    vertex_pattern = re.compile(r'v\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)')  #  v 2669.676921 661.911471 25.857700
    normal_pattern = re.compile(r'vn\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)')
    #face_pattern = re.compile(r'f\s + (-?\d+)/(-?\d+)/(-?\d+)\s + (-?\d+)/(-?\d+)/(-?\d+)\s + (-?\d+)/(-?\d+)/(-?\d+)\s ') # f 18003/11004/21081 18009/11011/21088 18250/11012/21089 18244/11005/21082
    face_pattern = re.compile(r'f\s+.')

    flag = False
    with open(file_path, 'r') as f:
        #for line in range(len(f)):
        for index, line in enumerate(f):

            ## reading out each vertex, normal and face when flag is true.    
            if flag == True:
                ## reading out each vertex, normal and face.
                vertex_match = vertex_pattern.match(line)
                normal_match = normal_pattern.match(line)
                face_match = face_pattern.match(line)
                if vertex_match:
                    x, y, z = map(float, vertex_match.groups())
                    vertices.append((x, y, z))
                elif normal_match:
                    x, y, z = map(float, normal_match.groups())
                    normals.append((x, y, z))
                elif face_match:
                    faces.append([num.split('/')[0] for num in line.split()[1:]])

            ## figuring out when to start reading out the vertices, normals and faces
            if part_1 in re.split(r"-|,|_|\n", line):    
                flag = True
                
                if idx >0: #if not the first occurence, such that vertices is filled already 
                    ribs.append(np.array(vertices)) #append the vertices of the rib to the ribs list
                    vertices = []

                idx = idx + 1 #increase the index of the rib

            ## figuring out when to stop reading out the vertices, normals and faces
            if part_2 in re.split(r"-|,|_|\n", line):
                ## also getting the last rib
                ribs.append(vertices) #append the vertices of the rib to the ribs list
                vertices = []
                flag = False
                break

    return vertices, normals, faces, ribs


def analyze_obj_file(file_path):
    ''' Reads the .obj file and returns the vertices, normals, and faces'''
    
    

    flag = False
    with open(file_path, 'r') as f:
        #for line in range(len(f)):
        for index, line in enumerate(f):
            
            if re.compile(r'g\s+.').match(line) or re.compile(r'usemtl\s+.').match(line):
                print('index: ',index,' | text: ',line)
    return