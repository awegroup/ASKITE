import numpy as np

def airfoil_coeffs(alpha, coeffs):
    cl = np.interp(np.rad2deg(alpha),coeffs[:,0],coeffs[:,1])
    cd = np.interp(np.rad2deg(alpha),coeffs[:,0],coeffs[:,2])
    cm = np.interp(np.rad2deg(alpha),coeffs[:,0],coeffs[:,3])
           
    return cl,cd,cm

def LEI_airf_coeff(t,k,alpha):
    """
    ----------
    t : TYPE
        DESCRIPTION.
    k : TYPE
        DESCRIPTION.
    alpha : TYPE
        DESCRIPTION.

    Returns
    -------
    Cl : TYPE
        DESCRIPTION.
    Cd : TYPE
        DESCRIPTION.
    Cm : TYPE
        DESCRIPTION.

    """
    C20 = -0.008011
    C21 = -0.000336
    C22 = 0.000992
    C23 = 0.013936
    C24 = -0.003838
    C25 = -0.000161
    C26 = 0.001243
    C27 = -0.009288
    C28 = -0.002124
    C29 = 0.012267
    C30 = -0.002398
    C31 = -0.000274
    C32 = 0
    C33 = 0
    C34 = 0
    C35 = -3.371000
    C36 = 0.858039
    C37 = 0.141600
    C38 = 7.201140
    C39 = -0.676007
    C40 = 0.806629
    C41 = 0.170454
    C42 = -0.390563
    C43 = 0.101966
    C44 = 0.546094
    C45 = 0.022247
    C46 = -0.071462
    C47 = -0.006527
    C48 = 0.002733
    C49 = 0.000686
    C50 = 0.123685
    C51 = 0.143755
    C52 = 0.495159
    C53 = -0.105362
    C54 = 0.033468
    C55 = -0.284793
    C56 = -0.026199
    C57 = -0.024060
    C58 = 0.000559
    C59 = -1.787703
    C60 = 0.352443
    C61 = -0.839323
    C62 = 0.137932
    
    S9 = C20*t**2+C21*t+C22
    S10 = C23*t**2+C24*t+C25
    S11 = C26*t**2+C27*t+C28
    S12 = C29*t**2+C30*t+C31
    S13 = C32*t**2+C33*t+C34
    S14 = C35*t**2+C36*t+C37
    S15 = C38*t**2+C39*t+C40
    S16 = C41*t**2+C42*t+C43
    
    lambda5 = S9*k+S10
    lambda6 = S11*k+S12
    lambda7 = S13*k+S14
    lambda8 = S15*k+S16
    
    Cl = lambda5*alpha**3+lambda6*alpha**2+lambda7*alpha+lambda8
    Cd = ((C44*t+C45)*k**2+(C46*t+C47)*k+(C48*t+C49))*alpha**2+\
        (C50*t+C51)*k+(C52*t**2+C53*t+C54)
    Cm = ((C55*t+C56)*k+(C57*t+C58))*alpha**2+(C59*t+C60)*k+(C61*t+C62)
    
    if alpha > 20 or alpha<-20:
        Cl = 2*np.cos(np.deg2rad(alpha))*np.sin(np.deg2rad(alpha))**2
        Cd = 2*np.sin(np.deg2rad(alpha))**3
        
    return Cl,Cd,Cm

def calculate_2D_coefficients(controlpoints,N_segments_refined,N_segments,N_split,tubediams,canopyheights):

    ''' Returns the 2D coefficients of the kite
        using the interpolated correlations from [Breukels2011]
        
        input:
            controlpoints       : list of control points of each refined segment 
            N_segments_refined  : number of refined segments (aero sections)
            N_segments          : number of segments (canopy pieces between struts
            N_split             : number of times each segment is splitted for the aero-refinement

        output:
            data_airf           : array [[alpha,Cl,Cd,Cm],..] of 2D coefficients for each refined segment              
    '''

    ## finding the tube diameters of each refined element
    tubediams_refined = np.array([])
    for i in range(N_segments): #loop through each segments

        # define a linspace for the chordwise aero lines in between each chordwise struc line
        # linspace goes from thickness 'i' to 'i+1',and has 'N_split+1' points
        tubediam_at_lines_i = np.linspace(tubediams[i],tubediams[i+1],N_split+1)
        tubediam_i = np.array([])
        # tubediam_i = []
        for j in range(len(tubediam_at_lines_i)-1):
            tubediam_i = np.append(tubediam_i,(tubediam_at_lines_i[j] +tubediam_at_lines_i[j+1])/2)
            # tubediam_i.append((tubediam_at_lines_i[j] +tubediam_at_lines_i[j+1])/2)
        tubediams_refined = np.append(tubediams_refined,tubediam_i)

    ## finding the canopy height of each refined element    
    canopyheights_refined = np.array([])
    for i in range(N_segments):
        canopyheight_at_lines_i = np.linspace(canopyheights[i],canopyheights[i+1],N_split+1)
        canopyheight_i = np.array([])
        # canopyheight_i = []
        for j in range(len(canopyheight_at_lines_i)-1):
            canopyheight_i = np.append(canopyheight_i,(canopyheight_at_lines_i[j] +canopyheight_at_lines_i[j+1])/2)
            # canopyheight_i.append((canopyheight_at_lines_i[j] +canopyheight_at_lines_i[j+1])/2)
        canopyheights_refined = np.append(canopyheights_refined,canopyheight_i)


    # defining the angle of attack range
    aoas = np.arange(-20,21,1)
    
    # defining the t (= tubediam/chord) of each aero section
    t_list_refined = np.empty(N_segments_refined)
    # defining the k (= MaxCanopyHeight/chord) of each aero section
    k_list_refined = np.empty(N_segments_refined)
    #defining the data_airf array struc [len(aoas),[alpha,Cl,Cd,Cm],N_chordwise_aero_lines-1]
    data_airf = np.empty((len(aoas),4,N_segments_refined))

    ## loop through each aero section
    for i in range(N_segments_refined):
        for j,alpha in enumerate(aoas): # looping through each angle
            ## Find the chord
            ### --- DEFINITION OF CHORD ---
            ### Is defined as the longest distance possible measured from the TE to the LE
            ### which due to the "vertical alignment of TE and LE tube centre" [Surfplan], 
            ### is initially always the horizontal

            chord_i = controlpoints[i]['chord']
            t_list_refined[i] = tubediams_refined[i]/chord_i
            k_list_refined[i] = canopyheights_refined[i]/chord_i

            alpha = aoas[j] 

            # using [Breukels2011] correlation model to find the 2D airfoil data
            Cl,Cd,Cm = LEI_airf_coeff(t_list_refined[i], k_list_refined[i], alpha)
            data_airf[j,0,i] = alpha
            data_airf[j,1,i] = Cl
            data_airf[j,2,i] = Cd
            data_airf[j,3,i] = Cm
        
    return data_airf





#%% new 2D polar things

def LEI_airf_coeff(t,k,alpha):
    """
    ----------
    t : TYPE
        DESCRIPTION.
    k : TYPE
        DESCRIPTION.
    alpha : TYPE
        DESCRIPTION.

    Returns
    -------
    Cl : TYPE
        DESCRIPTION.
    Cd : TYPE
        DESCRIPTION.
    Cm : TYPE
        DESCRIPTION.

    """
    C20 = -0.008011
    C21 = -0.000336
    C22 = 0.000992
    C23 = 0.013936
    C24 = -0.003838
    C25 = -0.000161
    C26 = 0.001243
    C27 = -0.009288
    C28 = -0.002124
    C29 = 0.012267
    C30 = -0.002398
    C31 = -0.000274
    C32 = 0
    C33 = 0
    C34 = 0
    C35 = -3.371000
    C36 = 0.858039
    C37 = 0.141600
    C38 = 7.201140
    C39 = -0.676007
    C40 = 0.806629
    C41 = 0.170454
    C42 = -0.390563
    C43 = 0.101966
    C44 = 0.546094
    C45 = 0.022247
    C46 = -0.071462
    C47 = -0.006527
    C48 = 0.002733
    C49 = 0.000686
    C50 = 0.123685
    C51 = 0.143755
    C52 = 0.495159
    C53 = -0.105362
    C54 = 0.033468
    C55 = -0.284793
    C56 = -0.026199
    C57 = -0.024060
    C58 = 0.000559
    C59 = -1.787703
    C60 = 0.352443
    C61 = -0.839323
    C62 = 0.137932
    
    S9 = C20*t**2+C21*t+C22
    S10 = C23*t**2+C24*t+C25
    S11 = C26*t**2+C27*t+C28
    S12 = C29*t**2+C30*t+C31
    S13 = C32*t**2+C33*t+C34
    S14 = C35*t**2+C36*t+C37
    S15 = C38*t**2+C39*t+C40
    S16 = C41*t**2+C42*t+C43
    
    lambda5 = S9*k+S10
    lambda6 = S11*k+S12
    lambda7 = S13*k+S14
    lambda8 = S15*k+S16
    
    Cl = lambda5*alpha**3+lambda6*alpha**2+lambda7*alpha+lambda8
    Cd = ((C44*t+C45)*k**2+(C46*t+C47)*k+(C48*t+C49))*alpha**2+\
        (C50*t+C51)*k+(C52*t**2+C53*t+C54)
    Cm = ((C55*t+C56)*k+(C57*t+C58))*alpha**2+(C59*t+C60)*k+(C61*t+C62)
    
    if alpha > 20 or alpha<-20:
        Cl = 2*np.cos(np.deg2rad(alpha))*np.sin(np.deg2rad(alpha))**2
        Cd = 2*np.sin(np.deg2rad(alpha))**3
        
    return Cl,Cd,Cm

def calculate_polar_lookup_table(controlpoints,n_segments,n_panels_aero,N_split,airfoil_data):

    ''' Returns the 2D coefficients of the kite
        using the interpolated correlations from [Breukels2011]
        
        input:
            controlpoints       : list of control points of each refined segment 
            N_segments_refined  : number of refined segments (aero sections)
            N_segments          : number of segments (canopy pieces between struts
            N_split             : number of times each segment is splitted for the aero-refinement

        output:
            data_airf           : array [[alpha,Cl,Cd,Cm],..] of 2D coefficients for each refined segment              
    '''
    TUBE_DIAMETERS = airfoil_data['TUBE_DIAMETERS']
    IS_TUBE_DIAMETER_DIMENSIONLESS = airfoil_data['IS_TUBE_DIAMETER_DIMENSIONLESS']
    CANOPY_MAX_HEIGHTS = airfoil_data['CANOPY_MAX_HEIGHTS']
    IS_MAX_CANOPY_HEIGHT_DIMENSIONLESS = airfoil_data['IS_MAX_CANOPY_HEIGHT_DIMENSIONLESS']

    thicc = np.array([])
    for i in range(n_segments):
        temp = np.linspace(TUBE_DIAMETERS[i],TUBE_DIAMETERS[i+1],N_split+1)
        temp1 = []
        for a in range(len(temp)-1):
            temp1 = np.append(temp1,(temp[a] +temp[a+1])/2)
            # temp1.append((temp[a] +temp[a+1])/2)
        thicc = np.append(thicc,temp1)
        
    camb = np.array([])
    for i in range(n_segments):
        temp = np.linspace(CANOPY_MAX_HEIGHTS[i],CANOPY_MAX_HEIGHTS[i+1],N_split+1)
        temp1 = []
        for a in range(len(temp)-1):
            temp1 = np.append(temp1,(temp[a] +temp[a+1])/2)
            # temp1.append((temp[a] +temp[a+1])/2)
        camb = np.append(camb,temp1)

    aoas = np.arange(-20,21,1)
    data_airf = np.empty((len(aoas),4,n_panels_aero))
    thicc_c = np.empty(n_panels_aero)
    camb_c = np.empty(n_panels_aero)
    for i in range(n_panels_aero):
        for j in range(len(aoas)):

            if IS_TUBE_DIAMETER_DIMENSIONLESS:
                thicc_c[i] = thicc[i]
            else:
                thicc_c[i] = thicc[i]/controlpoints[i]['chord']
            
            if IS_MAX_CANOPY_HEIGHT_DIMENSIONLESS:
                camb_c[i] = camb[i]
            else:
                camb_c[i] = camb[i]/controlpoints[i]['chord']
            
            alpha = aoas[j]
            Cl,Cd,Cm = LEI_airf_coeff(thicc_c[i], camb_c[i], alpha)
            data_airf[j,0,i] = alpha
            data_airf[j,1,i] = Cl
            data_airf[j,2,i] = Cd
            data_airf[j,3,i] = Cm

    return data_airf
