#%% LIFT-FORCE
import numpy as np

##TODO: used
def distance(A, B):
	return np.linalg.norm(B - A)
	#return np.sqrt( ((B[0] - A[0]) ** 2) + ((B[1] - A[1]) ** 2) + ((B[2] - A[2]) ** 2) )


def get_lift_force_orientation(plates,pos,lift_force,alpha_0,Vw,equal_boolean):
	# Specifying tubular frame + plate diagonals.
	# Reasons for incorporating plate diagonals is that otherwise the quadrilaterals can deform into other shapes
	# This phenomena is best explained by looking at a rectangle that turns into a diamond/kite shape

	for i in np.arange(0,len(plates)): #looping through each panel

		### Calculating the area of the plate
		# calculating the sides
		side_1 = np.linalg.norm(pos[plates[i][0]]-pos[plates[i][1]])
		side_2 = np.linalg.norm(pos[plates[i][1]]-pos[plates[i][2]])
		side_3 = np.linalg.norm(pos[plates[i][2]]-pos[plates[i][3]])
		side_4 = np.linalg.norm(pos[plates[i][3]]-pos[plates[i][0]])

		# Calculating the semi-perimeter (s) of the given quadilateral
		s = (side_1 + side_2 + side_3 + side_4) / 2

		# Applying Brahmagupta's formula to #https://en.wikipedia.org/wiki/Brahmagupta%27s_formula
		# get maximum area of quadrilateral
		import math
		area_p =  math.sqrt((s - side_1) *(s - side_2) * (s - side_3) * (s - side_4))
		area_w = 27708749.72261707 ##TODO: should adapt to this actual area, calculated per time-step

		### Calculating the angle of attack
		from numpy.linalg import norm

		middle_LE_point = (np.subtract(pos[plates[i][1]], pos[plates[i][0]]) / 2) + pos[plates[i][0]]
		middle_TE_point = (np.subtract(pos[plates[i][3]], pos[plates[i][2]]) / 2) + pos[plates[i][2]]
		middle_vec = np.subtract(middle_TE_point,middle_LE_point)
		middle_vec_unit = middle_vec / norm(middle_vec)
		vw_vec_unit = np.array([1,0,0]) / norm(np.array([1,0,0]))
		aoa_p = np.arccos(np.dot(middle_vec_unit,vw_vec_unit)/(norm(middle_vec_unit)*norm(vw_vec_unit)))

		###  Define the perpendicular unit vector to each segment
		# Defining 2 vectors tangential to the plane, by using the diagonals
		diagonal_1 = np.subtract(pos[plates[i][0]], pos[plates[i][2]])
		diagonal_2 = np.subtract(pos[plates[i][1]], pos[plates[i][3]])

		# Finding the cross product of these two vectors, should give a perpendicular vector / two solutions actually.
		perpendicular_vector_1 = np.cross(diagonal_1,diagonal_2)
		unit_perp_vec = perpendicular_vector_1 / np.linalg.norm(perpendicular_vector_1)

		### Find the direction of lift and drag
		#L_vec = [0,unit_perp_vec[1],unit_perp_vec[2]] #remove x-component from perpendicular vector
		#L_vec_unit = L_vec / np.linalg.norm(L_vec) #normalize the vector
		L_vec_unit = unit_perp_vec

		### Calculating Lift
		rho = 1.225
		
		# correcting the angle of attack for the orientation of the plate
		if unit_perp_vec[0] >0: #the plate is tilted backwards
			aoa_p = aoa_p
		elif unit_perp_vec[0] <0: #the plate is tilted forwards
			aoa_p = -aoa_p

		#print('middle_vec_unit',middle_vec_unit)
		#print('aoa_p',np.rad2deg(aoa_p), 'alpha_0',alpha_0)

		# AoA is equal to aoa_p + alpha_0
		Cl_p = 2*np.pi*np.sin(aoa_p+np.deg2rad(alpha_0))
		L_p = 0.5*rho*(Vw**2)*Cl_p*(area_p*1e-6)

		#print('aoa_p',np.rad2deg(aoa_p), 'alpha_0',alpha_0)
		#print('aoa_p+np.deg2rad(alpha_0)',aoa_p+np.deg2rad(alpha_0))
		#print('Cl_p',Cl_p,'L_p',L_p)
		#print(' ')

		### Splitting lift into Fz_p and Fy_p
		# increase in x-direction due to the initial alpha_0
		Fx_p = L_p*L_vec_unit[0] ##TODO:+L_p*L_vec_unit[0]*np.sin(np.deg2rad(alpha_0))
		Fy_p = L_p*L_vec_unit[1] # side-force
		# reduction in z-direction due to the alpha_0
		Fz_p = L_p*L_vec_unit[2] ##TODO: *np.cos(np.deg2rad(alpha_0)) # up-force ('lift')

		# ### Calculating Drag
		#Cd_w = 0.07 + np.rad2deg((aoa_p+np.deg2rad(alpha_0)))*(0.23/20)
		#Cd_p = Cd_w*(area_p/area_w)
		#D_p = 0.5*rho*(Vw**2)*Cd_p*(area_p*1e-6)
		D_p = 0

		#print('Cl_p',Cl_p,'Cd_p',Cd_p,'L_p',L_p,'D_p',D_p)

		### Apply 1/4 of perpendicular unit vector to each node of the respective panel
		for k,j in enumerate(plates[i]): #loop Through all the indexes for the defined plate
			if equal_boolean == True:
				lift_force[j, 0] += 0.25*D_p+0.25*Fx_p  # Drag only works parallel to the wind
				lift_force[j, 1] += 0.25*Fy_p
				lift_force[j, 2] += 0.25*Fz_p				
			else:
				if k == 0 or k == 1: #LEADING EDGE
					lift_force[j, 0] += 0.25*D_p+0.25*Fx_p  # Drag only works parallel to the wind
					lift_force[j, 1] += 0.5*0.75*Fy_p
					lift_force[j, 2] += 0.5*0.75*Fz_p
				elif k == 2 or k == 3: #TRAILING EDGE
					lift_force[j, 0] += 0.25*D_p+0.25*Fx_p  # Drag only works parallel to the wind
					lift_force[j, 1] += 0.5*0.25*Fy_p
					lift_force[j, 2] += 0.5*0.25*Fz_p				

	# # giving the bridle points also drag
	# for i in range(0,len(pos)):
	# 	if i not in np.concatenate(plates):
	# 		lift_force[i, 0] = 0 #2
			
	return lift_force

'''

###for testing
# points_CAD = np.loadtxt('geometry_and_connectivity/Geometry_modified_kcu.csv',delimiter = ',')
# pos = points_CAD
u_p = 1
date = '17jun'
filename = '../AerostructuralModel_v'+date+'/run_results/pos_up_'+str(int(100*u_p))+'.csv'
pos = np.loadtxt(filename,delimiter = ',')
alpha_0 = 0

import sys
sys.path.insert(0, '../AerostructuralModel_v'+date+'/geometry_and_connectivity/') 
import functions_connectivity as conn
ci_kite,cj_kite, plates = conn.get_kite_plate_connectivity()
Vw = 20
equal_boolean = False
lift_force = np.zeros((len(pos), 3))
lift_force = get_lift_force_orientation(plates,pos,lift_force,alpha_0,Vw,equal_boolean)
print('lift_force',lift_force)
print(np.max(lift_force[:,0]))
print(np.max(lift_force[:,1]))
print(np.max(lift_force[:,2]))
print("SUM:",np.sum(lift_force[:,0]))

'''
