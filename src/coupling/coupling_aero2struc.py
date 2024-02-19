#%%
import numpy as np

def force2nodes(F,Fpoint,nodes,tangential): 
    
    P1 = line_intersect(nodes[0,:],nodes[1,:],Fpoint, Fpoint+tangential)
    d1 = Fpoint-P1
    
    M1 = np.cross(d1, F)
    
    P2 = line_intersect(nodes[3,:],nodes[2,:],Fpoint, Fpoint+tangential)
    
    d2 = P2-P1
    Fp2 = np.cross(M1,d2)/np.linalg.norm(d2)**2
    
    Fp1 = F-Fp2
    
    M3 = np.cross(P1-nodes[0,:], Fp1)
    d3 = nodes[1,:]-nodes[0,:]
    F3 = np.cross(M3,d3)/np.linalg.norm(d3)**2
    
    node1 = Fp1-F3
    node2 = F3
    
    M4 = np.cross(P2-nodes[2,:], Fp2)
    d4 = nodes[3,:]-nodes[2,:]
    F4 = np.cross(M4,d4)/np.linalg.norm(d4)**2
    
    node4 = F4
    node3 = Fp2-F4
    
    Fnode = np.array([node1,
             node2,
             node3,
             node4])
    
    return Fnode

def moment2nodes(M,Mpoint,nodes,tangential):  
    d = tangential*0.05
    
    dF = np.cross(M,d)
    dF = dF/np.linalg.norm(dF)
    Fmag = np.linalg.norm(M)/np.linalg.norm(np.cross(dF,d))
    F = dF*Fmag
    
    P1 = Mpoint+d
    
    Fnode1 = force2nodes(F, P1, nodes,tangential)
    Fnode2 = force2nodes(-F, Mpoint, nodes,tangential)

    Fnode = np.array(Fnode1)+np.array(Fnode2)
    
    return Fnode

def line_intersect(p1,p2,p3,p4):
    
    p13 = np.empty(3)
    p43 = np.empty(3)
    p21 = np.empty(3)
    pa = np.empty(3)
    pb = np.empty(3)
    
    p13[0] = p1[0] - p3[0];
    p13[1] = p1[1] - p3[1];
    p13[2] = p1[2] - p3[2];
    
    p43[0] = p4[0] - p3[0];
    p43[1] = p4[1] - p3[1];
    p43[2] = p4[2] - p3[2];
    
    
    p21[0] = p2[0] - p1[0];
    p21[1] = p2[1] - p1[1];
    p21[2] = p2[2] - p1[2];
    
    
    d1343 = p13[0] * p43[0] + p13[1] * p43[1] + p13[2] * p43[2];
    d4321 = p43[0] * p21[0] + p43[1] * p21[1] + p43[2] * p21[2];
    d1321 = p13[0] * p21[0] + p13[1] * p21[1] + p13[2] * p21[2];
    d4343 = p43[0] * p43[0] + p43[1] * p43[1] + p43[2] * p43[2];
    d2121 = p21[0] * p21[0] + p21[1] * p21[1] + p21[2] * p21[2];
    
    denom = d2121 * d4343 - d4321 * d4321;
    
    
    
    numer = d1343 * d4321 - d1321 * d4343;
    
    mua = numer / denom;
    mub = (d1343 + d4321 * mua) / d4343;
    
    pa[0] = p1[0] + mua * p21[0];
    pa[1] = p1[1] + mua * p21[1];
    pa[2] = p1[2] + mua * p21[2];
    pb[0] = p3[0] + mub * p43[0];
    pb[1] = p3[1] + mub * p43[1];
    pb[2] = p3[2] + mub * p43[2];

    return pa

def aero2struc(pos,ci,cj,plates,F, M,ringvec,controlpoints):
    
    lift_force = np.zeros(pos.shape)
    N_struct = len(plates)
    N_split = int(len(controlpoints)/N_struct)
    for i in np.arange(0,len(controlpoints)): #looping through each panel
        sec = (N_struct-1)-int((i+1)/N_split-0.01)
        # Fi = (F[i][0] +F[i][1])*np.linalg.norm(ringvec[i]['r0'])
        # Fi = (F[i])*np.linalg.norm(ringvec[i]['r0'])
        # Mi = M[i]*np.linalg.norm(ringvec[i]['r0'])
        Fi = F[i]
        Mi = M[i]
        Mi = Mi*controlpoints[i]['airf_coord'][:,2]  
        
        if sec>4:
            Pnodes = np.array([pos[plates[sec][0],:],
                           pos[plates[sec][1],:],
                           pos[plates[sec][2],:],
                           pos[plates[sec][3],:]])
        else:
            Pnodes = np.array([pos[plates[sec][1],:],
                           pos[plates[sec][0],:],
                           pos[plates[sec][3],:],
                           pos[plates[sec][2],:]])
        Fnode = force2nodes(Fi, controlpoints[i]['coordinates_aoa'], Pnodes, controlpoints[i]['tangential'])
        # print(sum(Fnode)-Fi)
        # M1 = np.cross(Pnodes[0]-controlpoints[i]['coordinates_aoa'], Fnode[0,:])
        # M2 = np.cross(Pnodes[1]-controlpoints[i]['coordinates_aoa'], Fnode[1,:])
        # M3 = np.cross(Pnodes[2]-controlpoints[i]['coordinates_aoa'], Fnode[2,:])
        # M4 = np.cross(Pnodes[3]-controlpoints[i]['coordinates_aoa'], Fnode[3,:])
        # MT = M1+M2+M3+M4
        # print(MT)
        Fnode += moment2nodes(Mi, controlpoints[i]['coordinates_aoa'], Pnodes, controlpoints[i]['tangential'])
        
        
    
        if sec>4:
            
            lift_force[plates[sec][0],:] += Fnode[0]
            lift_force[plates[sec][1],:] += Fnode[1]
            lift_force[plates[sec][2],:] += Fnode[2]
            lift_force[plates[sec][3],:] += Fnode[3]
        else: 
            lift_force[plates[sec][1],:] += Fnode[0]
            lift_force[plates[sec][0],:] += Fnode[1]
            lift_force[plates[sec][3],:] += Fnode[2]
            lift_force[plates[sec][2],:] += Fnode[3]
            
    return lift_force


#%% V9_60 trials?

def aero2struc_V9_60(lift_force,plate_point_indices):
    """
    Takes lift-force, applied only at the plate corner points, 
    and distributes it to the structural wing nodes
    """
    ##  print(lift_force.shape)
    ##  print(points_CAD_discretized.shape)
    ##  print(plate_point_indices.shape)


    ## 1) Instead of taking the points_wing from ci_wing,cj_wing -which is not ordered
    #     we can take the points_wing from rib_db_whole_model_struts, which is ordered

    ## 2) Lift is applied at each plate_point_indices index, and we know its order
    #     thus, we can take the struts front and rear force and distribute it
    #     by taking the average of the two strut forces

    wing_node_indices_list = []

    for i,strut in enumerate(rib_db_whole_model_struts):
        wing_node_indices = strut[5][7:]
        wing_node_indices_list.append(wing_node_indices)
        ##  print("i",i, "| wing_node_indices",wing_node_indices)

        ## build in a check, to make sure we are only taking the bridle-attachment nodes
        for j in wing_node_indices: # loop through each index
            #print('j',j)
            if j not in np.hstack((ci_wing,cj_wing)): # verify that it is in the ci_wing or cj_wing
                print('error') # if not print an error

        ## if it is not the last strut, take the left-points of the plate
        if i < len(rib_db_whole_model_struts)-1:
            ## left points of the plate
            idx_plate = i 
            idx_plate_front = 0
            idx_plate_rear  = 3

        ## if it is the last strut, take the right_points of the plate
        if i == len(rib_db_whole_model_struts):
            ## right points of the plate
            idx_plate = i-1
            idx_plate_front = 1
            idx_plate_rear  = 2

        ## get the lift at the strut
        lift_at_strut_i_front = lift_force[plate_point_indices[idx_plate][idx_plate_front]] # get the lift at the strut
        lift_at_strut_i_rear  = lift_force[plate_point_indices[idx_plate][idx_plate_rear]] # get the lift at the strut
        lift_at_strut_i       = (lift_at_strut_i_front + lift_at_strut_i_rear)/2 # get the average lift at the strut

        ## distribute the lift EQUALLY over the present wing node ##TODO: equally is not physically correct
        lift_per_wing_node_i = lift_at_strut_i/len(wing_node_indices) # get the lift per wing node
        for idx in wing_node_indices: # loop through each index
            lift_force[idx] = lift_per_wing_node_i # distribute the lift over the wing node

    return lift_force,wing_node_indices_list
