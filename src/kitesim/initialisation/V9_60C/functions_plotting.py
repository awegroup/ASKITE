#%% libraries
import plotly.graph_objects as go
import plotly.offline as offline
import numpy as np


def get_xyz_list_MinMax_range(N_D,points_list): #use sinternally for making nice plots
    x_lst, y_lst, z_lst = [], [], []

    if N_D == 4: #if it has 3 dimensions, i.e. list of list of lists of points

        flat_point_list = []
        for list_3 in points_list:
            for list in list_3:
                for point in list:
                    # Extract the x, y, and z coordinates from the array
                    x_lst.append(point[0])
                    y_lst.append(point[1])
                    z_lst.append(point[2])
                    flat_point_list.append(point[0])
                    flat_point_list.append(point[1])
                    flat_point_list.append(point[2])

        # Get the min and max range of the points, for aspect ratio
        min_range = min(flat_point_list)
        max_range = max(flat_point_list)

    elif N_D == 3: #if it has 3 dimensions, i.e. list of lists of points

        flat_point_list = []
        for list in points_list:
            for point in list:
                # Extract the x, y, and z coordinates from the array
                x_lst.append(point[0])
                y_lst.append(point[1])
                z_lst.append(point[2])
                flat_point_list.append([point[0],point[1],point[2]])

        # Get the min and max range of the points, for aspect ratio
        min_range = min(np.concatenate(flat_point_list))
        max_range = max(np.concatenate(flat_point_list))

    elif N_D == 2: #it has 2 dimensions, i.e. a list of point

        flat_point_list = []
        for point in points_list:
            # Extract the x, y, and z coordinates from the array
            x_lst.append(point[0])
            y_lst.append(point[1])
            z_lst.append(point[2])
            flat_point_list.append([point[0],point[1],point[2]])
        
        # Get the min and max range of the points, for aspect ratio
        min_range = min(np.concatenate(flat_point_list))
        max_range = max(np.concatenate(flat_point_list))

    return x_lst,y_lst,z_lst,min_range,max_range

def get_3D_plot_N_D(N_D,points_list, color, size=3, show_legend=False,plot_title='3D Plot'):

    x_lst,y_lst,z_lst,min_range,max_range = get_xyz_list_MinMax_range(N_D,points_list)

    ## Getting the plot
    min_range = min_range - 0.1*max_range
    max_range = max_range + 0.1*max_range

    # Create a figure, set aspect ratio and plot data
    fig = go.Figure(
        layout=go.Layout(   scene=dict(aspectmode='cube',
                            xaxis=dict(range=[min_range,max_range]), 
                            yaxis=dict(range=[min_range,max_range]), 
                            zaxis=dict(range=[min_range,max_range])),
                            title=dict( text = plot_title,
                                        x = 0.5,
                                        y = 0.95,
                                        font =dict( family= 'Arial, sans-serif',
                                                    size=30,
                                                    color='black'))),
                data=[go.Scatter3d(x=x_lst, y=y_lst, z=z_lst, 
                                    mode='markers',
                                    marker=dict(
                                        color=color,
                                            size=size))])
    
    # Show the plot
    offline.plot(fig, filename='filename.html', auto_open=True)
    #fig.savefig("figure.svg")
    #fig.show()

    return()

def get_3D_plot_N_D_double(variables1, variables2,plot_title='3D Plot'):

    N_D1,points_list1, color1, size1 = variables1
    N_D2,points_list2, color2, size2 = variables2
    x_lst1,y_lst1,z_lst1,min_range1,max_range1 = get_xyz_list_MinMax_range(N_D1,points_list1)
    x_lst2,y_lst2,z_lst2,min_range2,max_range2 = get_xyz_list_MinMax_range(N_D2,points_list2)

    ## Getting the plot
    min_range = min(min_range1,min_range2) - 0.1*max(max_range1,max_range2)
    max_range = max(max_range1,max_range2) + 0.1*max(max_range1,max_range2)

    # Create a figure, set aspect ratio and plot data
    fig = go.Figure(
        layout=go.Layout(   scene=dict(aspectmode='cube',
                            xaxis=dict(range=[min_range,max_range]), 
                            yaxis=dict(range=[min_range,max_range]), 
                            zaxis=dict(range=[min_range,max_range])),
                            title=dict( text = plot_title,
                                        x = 0.5,
                                        y = 0.95,
                                        font =dict( family= 'Arial, sans-serif',
                                                    size=30,
                                                    color='black'))),
                data=[go.Scatter3d(x=x_lst1, y=y_lst1, z=z_lst1, 
                                    mode='markers',
                                    marker=dict(color=color1,size=size1)),
                        go.Scatter3d(x=x_lst2, y=y_lst2, z=z_lst2, 
                                    mode='markers',
                                    marker=dict(color=color2,size=size2)),
                        ])
    
    # Show the plot
    offline.plot(fig, filename='filename.html', auto_open=True)
    #fig.savefig("figure.svg")
    #fig.show()

    return()

def get_3D_plot_N_D_triple(variables1, variables2,variables3,plot_title='3D Plot'):

    N_D1,points_list1, color1, size1 = variables1
    N_D2,points_list2, color2, size2 = variables2
    N_D3,points_list3, color3, size3 = variables3
    x_lst1,y_lst1,z_lst1,min_range1,max_range1 = get_xyz_list_MinMax_range(N_D1,points_list1)
    x_lst2,y_lst2,z_lst2,min_range2,max_range2 = get_xyz_list_MinMax_range(N_D2,points_list2)
    x_lst3,y_lst3,z_lst3,min_range3,max_range3 = get_xyz_list_MinMax_range(N_D3,points_list3)

    ## Getting the plot
    min_range = min(min_range1,min_range2,min_range3) - 0.1*max(max_range1,max_range2,max_range3)
    max_range = max(max_range1,max_range2,max_range3) + 0.1*max(max_range1,max_range2,max_range3)

    # Create a figure, set aspect ratio and plot data
    fig = go.Figure(
        layout=go.Layout(   scene=dict(aspectmode='cube',
                            xaxis=dict(range=[min_range,max_range]), 
                            yaxis=dict(range=[min_range,max_range]), 
                            zaxis=dict(range=[min_range,max_range])),
                            title=dict( text = plot_title,
                                        x = 0.5,
                                        y = 0.95,
                                        font =dict( family= 'Arial, sans-serif',
                                                    size=30,
                                                    color='black'))),
                data=[go.Scatter3d(x=x_lst1, y=y_lst1, z=z_lst1, 
                                    mode='markers',
                                    marker=dict(color=color1,size=size1)),
                        go.Scatter3d(x=x_lst2, y=y_lst2, z=z_lst2, 
                                    mode='markers',
                                    marker=dict(color=color2,size=size2)),
                        go.Scatter3d(x=x_lst3, y=y_lst3, z=z_lst3,
                                    mode='markers',
                                    marker=dict(color=color3,size=size3)),
                        ])
    
    # Show the plot
    offline.plot(fig, filename='3D_plot.html', auto_open=True)
    #fig.savefig("figure.svg")
    #fig.show()

    return()

def get_3D_plot_N_D_4times(variables1, variables2,variables3,variables4,plot_title='3D Plot'):

    N_D1,points_list1, color1, size1 = variables1
    N_D2,points_list2, color2, size2 = variables2
    N_D3,points_list3, color3, size3 = variables3
    N_D4,points_list4, color4, size4 = variables4
    x_lst1,y_lst1,z_lst1,min_range1,max_range1 = get_xyz_list_MinMax_range(N_D1,points_list1)
    x_lst2,y_lst2,z_lst2,min_range2,max_range2 = get_xyz_list_MinMax_range(N_D2,points_list2)
    x_lst3,y_lst3,z_lst3,min_range3,max_range3 = get_xyz_list_MinMax_range(N_D3,points_list3)
    x_lst4,y_lst4,z_lst4,min_range4,max_range4 = get_xyz_list_MinMax_range(N_D4,points_list4)

    ## Getting the plot
    min_range = min(min_range1,min_range2,min_range3,min_range4) - 0.1*max(max_range1,max_range2,max_range3,max_range4)
    max_range = max(max_range1,max_range2,max_range3,max_range4) + 0.1*max(max_range1,max_range2,max_range3,max_range4)

    # Create a figure, set aspect ratio and plot data
    fig = go.Figure(
        layout=go.Layout(   scene=dict(aspectmode='cube',
                            xaxis=dict(range=[min_range,max_range]), 
                            yaxis=dict(range=[min_range,max_range]), 
                            zaxis=dict(range=[min_range,max_range])),
                            title=dict( text = plot_title,
                                        x = 0.5,
                                        y = 0.95,
                                        font =dict( family= 'Arial, sans-serif',
                                                    size=30,
                                                    color='black'))),
                data=[go.Scatter3d(x=x_lst1, y=y_lst1, z=z_lst1, 
                                    mode='markers',
                                    marker=dict(color=color1,size=size1)),
                        go.Scatter3d(x=x_lst2, y=y_lst2, z=z_lst2, 
                                    mode='markers',
                                    marker=dict(color=color2,size=size2)),
                        go.Scatter3d(x=x_lst3, y=y_lst3, z=z_lst3,
                                    mode='markers',
                                    marker=dict(color=color3,size=size3)),
                        go.Scatter3d(x=x_lst4, y=y_lst4, z=z_lst4,
                                    mode='markers',
                                    marker=dict(color=color4,size=size4)),
                        ])
    
    # Show the plot
    offline.plot(fig, filename='filename.html', auto_open=True)
    #fig.savefig("figure.svg")
    #fig.show()

    return()

def get_3D_plot_N_D_5times(variables1, variables2,variables3,variables4,variables5,plot_title='3D Plot'):


    N_D1,points_list1, color1, size1 = variables1
    N_D2,points_list2, color2, size2 = variables2
    N_D3,points_list3, color3, size3 = variables3
    N_D4,points_list4, color4, size4 = variables4
    N_D5,points_list5, color5, size5 = variables5
    x_lst1,y_lst1,z_lst1,min_range1,max_range1 = get_xyz_list_MinMax_range(N_D1,points_list1)
    x_lst2,y_lst2,z_lst2,min_range2,max_range2 = get_xyz_list_MinMax_range(N_D2,points_list2)
    x_lst3,y_lst3,z_lst3,min_range3,max_range3 = get_xyz_list_MinMax_range(N_D3,points_list3)
    x_lst4,y_lst4,z_lst4,min_range4,max_range4 = get_xyz_list_MinMax_range(N_D4,points_list4)
    x_lst5,y_lst5,z_lst5,min_range5,max_range5 = get_xyz_list_MinMax_range(N_D5,points_list5)

    ## Getting the plot
    min_range = min(min_range1,min_range2,min_range3,min_range4,min_range5) - 0.1*max(max_range1,max_range2,max_range3,max_range4,max_range5)
    max_range = max(max_range1,max_range2,max_range3,max_range4,max_range5) + 0.1*max(max_range1,max_range2,max_range3,max_range4,max_range5)

    # Create a figure, set aspect ratio and plot data
    fig = go.Figure(
        layout=go.Layout(   scene=dict(aspectmode='cube',
                            xaxis=dict(range=[min_range,max_range]), 
                            yaxis=dict(range=[min_range,max_range]), 
                            zaxis=dict(range=[min_range,max_range])),
                            title=dict( text = plot_title,
                                        x = 0.5,
                                        y = 0.95,
                                        font =dict( family= 'Arial, sans-serif',
                                                    size=30,
                                                    color='black'))),
                data=[go.Scatter3d(x=x_lst1, y=y_lst1, z=z_lst1, 
                                    mode='markers',
                                    marker=dict(color=color1,size=size1)),
                        go.Scatter3d(x=x_lst2, y=y_lst2, z=z_lst2, 
                                    mode='markers',
                                    marker=dict(color=color2,size=size2)),
                        go.Scatter3d(x=x_lst3, y=y_lst3, z=z_lst3,
                                    mode='markers',
                                    marker=dict(color=color3,size=size3)),
                        go.Scatter3d(x=x_lst4, y=y_lst4, z=z_lst4,
                                    mode='markers',
                                    marker=dict(color=color4,size=size4)),
                        go.Scatter3d(x=x_lst5, y=y_lst5, z=z_lst5,
                                    mode='markers',
                                    marker=dict(color=color5,size=size5)),            
                        ])
    
    # Show the plot
    offline.plot(fig, filename='filename.html', auto_open=True)
    #fig.savefig("figure.svg")
    #fig.show()

    return()


def get_3D_plot_N_D_4times_PLUS_lines(variables1, variables2,variables3,variables4,line_variables,plot_title='3D Plot'):

    N_D1,points_list1, color1, size1 = variables1
    N_D2,points_list2, color2, size2 = variables2
    N_D3,points_list3, color3, size3 = variables3
    N_D4,points_list4, color4, size4 = variables4
    x_lst1,y_lst1,z_lst1,min_range1,max_range1 = get_xyz_list_MinMax_range(N_D1,points_list1)
    x_lst2,y_lst2,z_lst2,min_range2,max_range2 = get_xyz_list_MinMax_range(N_D2,points_list2)
    x_lst3,y_lst3,z_lst3,min_range3,max_range3 = get_xyz_list_MinMax_range(N_D3,points_list3)
    x_lst4,y_lst4,z_lst4,min_range4,max_range4 = get_xyz_list_MinMax_range(N_D4,points_list4)

    points,conn_i,conn_j,line_color,line_width = line_variables

    ## Getting the plot
    min_range = min(min_range1,min_range2,min_range3,min_range4) - 0.1*max(max_range1,max_range2,max_range3,max_range4)
    max_range = max(max_range1,max_range2,max_range3,max_range4) + 0.1*max(max_range1,max_range2,max_range3,max_range4)

    # Create a figure, set aspect ratio and plot data
    fig = go.Figure(
        layout=go.Layout(   scene=dict(aspectmode='cube',
                            xaxis=dict(range=[min_range,max_range]), 
                            yaxis=dict(range=[min_range,max_range]), 
                            zaxis=dict(range=[min_range,max_range])),
                            title=dict( text = plot_title,
                                        x = 0.5,
                                        y = 0.95,
                                        font =dict( family= 'Arial, sans-serif',
                                                    size=30,
                                                    color='black'))),
                data=[go.Scatter3d(x=x_lst1, y=y_lst1, z=z_lst1, 
                                    mode='markers',
                                    marker=dict(color=color1,size=size1)),
                        go.Scatter3d(x=x_lst2, y=y_lst2, z=z_lst2, 
                                    mode='markers',
                                    marker=dict(color=color2,size=size2)),
                        go.Scatter3d(x=x_lst3, y=y_lst3, z=z_lst3,
                                    mode='markers',
                                    marker=dict(color=color3,size=size3)),
                        go.Scatter3d(x=x_lst4, y=y_lst4, z=z_lst4,
                                    mode='markers',
                                    marker=dict(color=color4,size=size4)),
                        ])
    
    # Create the edges using `Scatter3d` and set `mode` to `lines`
    for i in range(len(conn_i)):
        fig.add_trace(go.Scatter3d(
            x=[points[conn_i[i]][0], points[conn_j[i]][0]],
            y=[points[conn_i[i]][1], points[conn_j[i]][1]],
            z=[points[conn_i[i]][2], points[conn_j[i]][2]],
            mode='lines', 
            line=dict(width=line_width, color=line_color),
            showlegend= False
        ))

    # Show the plot
    offline.plot(fig, filename='filename.html', auto_open=True)
    #fig.savefig("figure.svg")
    #fig.show()

    return()

def get_3D_plot_N_D_5times_PLUS_2lines(variables1, variables2,variables3,variables4,variables5,line_variables1,line_variables2,plot_title='3D Plot'):

    N_D1,points_list1, color1, size1 = variables1
    N_D2,points_list2, color2, size2 = variables2
    N_D3,points_list3, color3, size3 = variables3
    N_D4,points_list4, color4, size4 = variables4
    N_D5,points_list5, color5, size5 = variables5
    x_lst1,y_lst1,z_lst1,min_range1,max_range1 = get_xyz_list_MinMax_range(N_D1,points_list1)
    x_lst2,y_lst2,z_lst2,min_range2,max_range2 = get_xyz_list_MinMax_range(N_D2,points_list2)
    x_lst3,y_lst3,z_lst3,min_range3,max_range3 = get_xyz_list_MinMax_range(N_D3,points_list3)
    x_lst4,y_lst4,z_lst4,min_range4,max_range4 = get_xyz_list_MinMax_range(N_D4,points_list4)
    x_lst5,y_lst5,z_lst5,min_range5,max_range5 = get_xyz_list_MinMax_range(N_D5,points_list5)

    points1,conn_i1,conn_j1,line_color1,line_width1 = line_variables1
    points2,conn_i2,conn_j2,line_color2,line_width2 = line_variables2

    ## Getting the plot
    min_range = min(min_range1,min_range2,min_range3,min_range4,min_range5) - 0.1*max(max_range1,max_range2,max_range3,max_range4,max_range5)
    max_range = max(max_range1,max_range2,max_range3,max_range4,max_range5) + 0.1*max(max_range1,max_range2,max_range3,max_range4,max_range5)

    # Create a figure, set aspect ratio and plot data
    fig = go.Figure(
        layout=go.Layout(   scene=dict(aspectmode='cube',
                            xaxis=dict(range=[min_range,max_range]), 
                            yaxis=dict(range=[min_range,max_range]), 
                            zaxis=dict(range=[min_range,max_range])),
                            title=dict( text = plot_title,
                                        x = 0.5,
                                        y = 0.95,
                                        font =dict( family= 'Arial, sans-serif',
                                                    size=30,
                                                    color='black'))),
                data=[go.Scatter3d(x=x_lst1, y=y_lst1, z=z_lst1, 
                                    mode='markers',
                                    marker=dict(color=color1,size=size1)),
                        go.Scatter3d(x=x_lst2, y=y_lst2, z=z_lst2, 
                                    mode='markers',
                                    marker=dict(color=color2,size=size2)),
                        go.Scatter3d(x=x_lst3, y=y_lst3, z=z_lst3,
                                    mode='markers',
                                    marker=dict(color=color3,size=size3)),
                        go.Scatter3d(x=x_lst4, y=y_lst4, z=z_lst4,
                                    mode='markers',
                                    marker=dict(color=color4,size=size4)),
                        go.Scatter3d(x=x_lst5, y=y_lst5, z=z_lst5,
                                    mode='markers',
                                    marker=dict(color=color5,size=size5)),
                        ])
    
    # Create the edges using `Scatter3d` and set `mode` to `lines`
    for i in range(len(conn_i1)):
        fig.add_trace(go.Scatter3d(
            x=[points1[conn_i1[i]][0], points1[conn_j1[i]][0]],
            y=[points1[conn_i1[i]][1], points1[conn_j1[i]][1]],
            z=[points1[conn_i1[i]][2], points1[conn_j1[i]][2]],
            mode='lines', 
            line=dict(width=line_width1, color=line_color1),
            showlegend= False
        ))
    
    for i in range(len(conn_i2)):
        fig.add_trace(go.Scatter3d(
            x=[points2[conn_i2[i]][0], points2[conn_j2[i]][0]],
            y=[points2[conn_i2[i]][1], points2[conn_j2[i]][1]],
            z=[points2[conn_i2[i]][2], points2[conn_j2[i]][2]],
            mode='lines', 
            line=dict(width=line_width2, color=line_color2),
            showlegend= False
        ))

    # Show the plot
    offline.plot(fig, filename='filename.html', auto_open=True)
    #fig.savefig("figure.svg")
    #fig.show()

    return()


def get_3D_plot_N_D_6times_PLUS_2lines(variables1, variables2,variables3,variables4,variables5,variables6,line_variables1,line_variables2,plot_title='3D Plot'):

    N_D1,points_list1, color1, size1 = variables1
    N_D2,points_list2, color2, size2 = variables2
    N_D3,points_list3, color3, size3 = variables3
    N_D4,points_list4, color4, size4 = variables4
    N_D5,points_list5, color5, size5 = variables5
    N_D6,points_list6, color6, size6 = variables6
    x_lst1,y_lst1,z_lst1,min_range1,max_range1 = get_xyz_list_MinMax_range(N_D1,points_list1)
    x_lst2,y_lst2,z_lst2,min_range2,max_range2 = get_xyz_list_MinMax_range(N_D2,points_list2)
    x_lst3,y_lst3,z_lst3,min_range3,max_range3 = get_xyz_list_MinMax_range(N_D3,points_list3)
    x_lst4,y_lst4,z_lst4,min_range4,max_range4 = get_xyz_list_MinMax_range(N_D4,points_list4)
    x_lst5,y_lst5,z_lst5,min_range5,max_range5 = get_xyz_list_MinMax_range(N_D5,points_list5)
    x_lst6,y_lst6,z_lst6,min_range6,max_range6 = get_xyz_list_MinMax_range(N_D6,points_list6)

    points1,conn_i1,conn_j1,line_color1,line_width1 = line_variables1
    points2,conn_i2,conn_j2,line_color2,line_width2 = line_variables2

    ## Getting the plot
    min_range = min(min_range1,min_range2,min_range3,min_range4,min_range5,min_range6) - 0.1*max(max_range1,max_range2,max_range3,max_range4,max_range5,max_range6)
    max_range = max(max_range1,max_range2,max_range3,max_range4,max_range5,max_range6) + 0.1*max(max_range1,max_range2,max_range3,max_range4,max_range5,max_range6)

    # Create a figure, set aspect ratio and plot data
    fig = go.Figure(
        layout=go.Layout(   scene=dict(aspectmode='cube',
                            xaxis=dict(range=[min_range,max_range]), 
                            yaxis=dict(range=[min_range,max_range]), 
                            zaxis=dict(range=[min_range,max_range])),
                            title=dict( text = plot_title,
                                        x = 0.5,
                                        y = 0.95,
                                        font =dict( family= 'Arial, sans-serif',
                                                    size=30,
                                                    color='black'))),
                data=[go.Scatter3d(x=x_lst1, y=y_lst1, z=z_lst1, 
                                    mode='markers',
                                    marker=dict(color=color1,size=size1)),
                        go.Scatter3d(x=x_lst2, y=y_lst2, z=z_lst2, 
                                    mode='markers',
                                    marker=dict(color=color2,size=size2)),
                        go.Scatter3d(x=x_lst3, y=y_lst3, z=z_lst3,
                                    mode='markers',
                                    marker=dict(color=color3,size=size3)),
                        go.Scatter3d(x=x_lst4, y=y_lst4, z=z_lst4,
                                    mode='markers',
                                    marker=dict(color=color4,size=size4)),
                        go.Scatter3d(x=x_lst5, y=y_lst5, z=z_lst5,
                                    mode='markers',
                                    marker=dict(color=color5,size=size5)),
                        go.Scatter3d(x=x_lst6, y=y_lst6, z=z_lst6,
                                    mode='markers',
                                    marker=dict(color=color6,size=size6)),
                        ])
    
    # Create the edges using `Scatter3d` and set `mode` to `lines`
    for i in range(len(conn_i1)):
        fig.add_trace(go.Scatter3d(
            x=[points1[conn_i1[i]][0], points1[conn_j1[i]][0]],
            y=[points1[conn_i1[i]][1], points1[conn_j1[i]][1]],
            z=[points1[conn_i1[i]][2], points1[conn_j1[i]][2]],
            mode='lines', 
            line=dict(width=line_width1, color=line_color1),
            showlegend= False
        ))
    
    for i in range(len(conn_i2)):
        fig.add_trace(go.Scatter3d(
            x=[points2[conn_i2[i]][0], points2[conn_j2[i]][0]],
            y=[points2[conn_i2[i]][1], points2[conn_j2[i]][1]],
            z=[points2[conn_i2[i]][2], points2[conn_j2[i]][2]],
            mode='lines', 
            line=dict(width=line_width2, color=line_color2),
            showlegend= False
        ))

    # Show the plot
    offline.plot(fig, filename='filename.html', auto_open=True)
    #fig.savefig("figure.svg")
    #fig.show()

    return()




def get_3D_plot(points_list, color, size=3, show_legend=False):


    x_lst, y_lst, z_lst = [], [], []

    if np.ndim(np.concatenate(points_list)) == 2: #if it has 3 dimensions, i.e. list of lists of points
        
        for list in points_list:
            for point in list:
                # Extract the x, y, and z coordinates from the array
                x_lst.append(point[0])
                y_lst.append(point[1])
                z_lst.append(point[2])

        # Get the min and max range of the points, for aspect ratio
        min_range = min(np.concatenate(np.concatenate(points_list)))
        max_range = max(np.concatenate(np.concatenate(points_list)))
        min_range = min_range - 0.1*max_range
        max_range = max_range + 0.1*max_range

        # Create a figure, set aspect ratio and plot data
        fig = go.Figure(layout=go.Layout(scene=dict(aspectmode='cube',
                    xaxis=dict(range=[min_range,max_range]), 
                    yaxis=dict(range=[min_range,max_range]), 
                    zaxis=dict(range=[min_range,max_range]))),
                    data=[go.Scatter3d(x=x_lst, y=y_lst, z=z_lst, 
                                        mode='markers',
                                        marker=dict(
                                            color=color,
                                                size=size))])

    else: #it has 2 dimensions, i.e. a list of point

        # Get the min and max range of the points, for aspect ratio
        min_range = min(np.concatenate(points_list))
        max_range = max(np.concatenate(points_list))
        min_range = min_range - 0.1*max_range
        max_range = max_range + 0.1*max_range

        for point in points_list:
            # Extract the x, y, and z coordinates from the array
            x_lst.append(point[0])
            y_lst.append(point[1])
            z_lst.append(point[2])

        # Create a figure, set aspect ratio and plot data
        fig = go.Figure(layout=go.Layout(scene=dict(aspectmode='cube',
                    xaxis=dict(range=[min_range,max_range]), 
                    yaxis=dict(range=[min_range,max_range]), 
                    zaxis=dict(range=[min_range,max_range]))),
                    data=[go.Scatter3d(x=x_lst, y=y_lst, z=z_lst, 
                                        mode='markers',
                                        marker=dict(
                                            color=color,
                                                size=size))])


    # Show the plot
    offline.plot(fig, filename='filename.html', auto_open=True)
    #fig.savefig("figure.svg")
    #fig.show()

    return()


def get_3D_plot_with_lines(points,conn_i,conn_j,show_legend=False):

    # Get the min and max range of the points, for aspect ratio
    min_range = min(np.concatenate(points))
    max_range = max(np.concatenate(points))
    min_range = min_range - 0.1*max_range
    max_range = max_range + 0.1*max_range

    # Create an empty figure and set aspect ratio
    fig = go.Figure(layout=go.Layout(scene=dict(aspectmode='cube',
                xaxis=dict(range=[min_range,max_range]), 
                yaxis=dict(range=[min_range,max_range]), 
                zaxis=dict(range=[min_range,max_range]))) )

    # Create the nodes using `Scatter3d
    for point in points:
        fig.add_trace(go.Scatter3d(
            x = [point[0]],
            y =[point[1]],
            z =[point[2]],
            mode='markers',
            marker=dict(
                size=2,
                color='rgb(0, 0, 255)',
                symbol='circle'),
            showlegend=show_legend
        ))

    # Create the edges using `Scatter3d` and set `mode` to `lines`
    for i in range(len(conn_i)):
        fig.add_trace(go.Scatter3d(
            x=[points[conn_i[i]][0], points[conn_j[i]][0]],
            y=[points[conn_i[i]][1], points[conn_j[i]][1]],
            z=[points[conn_i[i]][2], points[conn_j[i]][2]],
            mode='lines', 
            line=dict(width=1, color='blue'),
            showlegend=show_legend
        ))

    #plotting the figure offline
    offline.plot(fig, filename='filename.html', auto_open=True)

    return


def get_3D_plot_wing_and_bridle(wing_data,bridle_data,show_legend=False):
    
    points_wing,conn_i_wing,conn_j_wing = wing_data
    points_bridle,conn_i_bridle,conn_j_bridle = bridle_data

    # Get the min and max range of the points, for aspect ratio
    min_range_wing = min(np.concatenate(points_wing))
    max_range_wing = max(np.concatenate(points_wing))
    min_range_bridle = min(np.concatenate(points_bridle))
    max_range_bridle = min(np.concatenate(points_bridle))

    min_range = min(min_range_wing,min_range_bridle)
    max_range = max(max_range_wing,max_range_bridle)
    min_range = min_range - 0.1*max_range
    max_range = max_range + 0.1*max_range

    # Create an empty figure and set aspect ratio
    fig = go.Figure(layout=go.Layout(scene=dict(aspectmode='cube',
                xaxis=dict(range=[min_range,max_range]), 
                yaxis=dict(range=[min_range,max_range]), 
                zaxis=dict(range=[min_range,max_range]))) )

    # Create the nodes using `Scatter3d
    for point in points_wing:
        fig.add_trace(go.Scatter3d(
            x = [point[0]],
            y =[point[1]],
            z =[point[2]],
            mode='markers',
            marker=dict(
                size=2,
                color='rgb(0, 0, 255)',
                symbol='circle'),
            showlegend=show_legend
        ))

    for point in points_bridle:
        fig.add_trace(go.Scatter3d(
            x = [point[0]],
            y =[point[1]],
            z =[point[2]],
            mode='markers',
            marker=dict(
                size=2,
                color='rgb(0, 0, 255)',
                symbol='circle'),
            showlegend=show_legend
        ))

    # Create the edges using `Scatter3d` and set `mode` to `lines`
    for i in range(len(conn_i_wing)):
        fig.add_trace(go.Scatter3d(
            x=[points_wing[conn_i_wing[i]][0], points_wing[conn_j_wing[i]][0]],
            y=[points_wing[conn_i_wing[i]][1], points_wing[conn_j_wing[i]][1]],
            z=[points_wing[conn_i_wing[i]][2], points_wing[conn_j_wing[i]][2]],
            mode='lines', 
            line=dict(width=1, color='blue'),
            showlegend=show_legend
        ))
    for i in range(len(conn_i_bridle)):
        fig.add_trace(go.Scatter3d(
            x=[points_bridle[conn_i_bridle[i]][0], points_bridle[conn_j_bridle[i]][0]],
            y=[points_bridle[conn_i_bridle[i]][1], points_bridle[conn_j_bridle[i]][1]],
            z=[points_bridle[conn_i_bridle[i]][2], points_bridle[conn_j_bridle[i]][2]],
            mode='lines', 
            line=dict(width=1, color='blue'),
            showlegend=show_legend
        ))

    #plotting the figure offline
    offline.plot(fig, filename='filename.html', auto_open=True)

    return




#%% Old functions maybe useful later

def get_3D_plot_double(list_1,list_2,size_1=1,size_2=1,color_1='blue',color_2='red',aspect=True):
    

    x_lst_1, y_lst_1, z_lst_1 = [], [], []
    if len(list_1[0])==3: #points_list is a list of points [[]]
        for point in list_1:
            # Extract the x, y, and z coordinates from the array
            x_lst_1.append(point[0])
            y_lst_1.append(point[1])
            z_lst_1.append(point[2])

    else :   #points_list is a list with lists of points [[[]]]
        for list in list_1:
            for point in list:
            # Extract the x, y, and z coordinates from the array
                x_lst_1.append(point[0])
                y_lst_1.append(point[1])
                z_lst_1.append(point[2])


    x_lst_2, y_lst_2, z_lst_2 = [], [], []
    if len(list_2[0])==3: #points_list is a list of points [[]]
        for point in list_2:
            # Extract the x, y, and z coordinates from the array
            x_lst_2.append(point[0])
            y_lst_2.append(point[1])
            z_lst_2.append(point[2])

    else:   #points_list is a list with lists of points [[[]]]
        for list in list_2:
            for point in list:
            # Extract the x, y, and z coordinates from the array
                x_lst_2.append(point[0])
                y_lst_2.append(point[1])
                z_lst_2.append(point[2])

    
    # Create a scatter plot
    fig = go.Figure(data=[go.Scatter3d(x=x_lst_1, y=y_lst_1, z=z_lst_1, 
                                    mode='markers',
                                    marker=dict(color=color_1,size=size_1)),
                        go.Scatter3d(x=x_lst_2, y=y_lst_2, z=z_lst_2, 
                                    mode='markers',
                                    marker=dict(color=color_2,size=size_2)),
                        ])
    fig.update_layout(scene=dict(aspectmode="manual", aspectratio=dict(x=1, y=1, z=1)))


    if aspect:

        fig.update_layout(scene=dict(aspectmode="manual", aspectratio=dict(x=1, y=1, z=1)))   
        fig.update_xaxes(range=[min(y_lst_1,), max(y_lst_1)])
        #fig.update_yaxes(range=[min(y_lst), max(y_lst)])
        #fig.update_zaxes(range=[min(y_lst), max(y_lst)])
        fig.update_layout(scene=dict(zaxis=dict(range=[min(y_lst_1,), max(y_lst_1)])))

    offline.plot(fig, filename='filename.html', auto_open=True)
    return

# %%
