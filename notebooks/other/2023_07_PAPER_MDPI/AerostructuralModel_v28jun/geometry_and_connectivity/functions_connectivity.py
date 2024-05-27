import numpy as np


def get_bridle_line_system_connectivity():

    ci_bridle_TE = [ 0, 0, 0,  21,21,  22,22,  26,26,  23,23,  27,27,  24,24,  28,28,  25,25,  29,29]
    cj_bridle_TE = [21,22,26,  23,27,  23, 1,  27,10,  24,25,  28,29,  18,17,  11,12,  16,15,  13,14]
    ci_bridle_LE = [  0, 0,  30,30,  33,33,  31,31,  34,34,  32,32,  35, 35,  23,  27,  30,33]
    cj_bridle_LE = [ 30,33,  31,32,  34,35,   2, 3,   9, 8,   4, 5,   7, 6,   19,  20,  19,20]

    ci_bridle = np.append(ci_bridle_TE,ci_bridle_LE)
    cj_bridle = np.append(cj_bridle_TE,cj_bridle_LE)

    return ci_bridle,cj_bridle

def get_bridle_line_system_connectivity_KCU():
    
    ci_bridle_TE = [ 0, 21, 21, 21,  22,22,  23,23,  27,27,  24,24,  28,28,  25,25,  29,29,  26,26,  30,30]
    cj_bridle_TE = [21,22,23,27,  24,28,  24, 1,  28,10,  25,26,  29,30,  18,17,  11,12,  16,15,  13,14]
    ci_bridle_LE = [  0, 0,  31,31,  34,34,  32,32,  35,35,  33,33,  36, 36,  24,  28,  31,34]
    cj_bridle_LE = [ 31,34,  32,33,  35,36,   2, 3,   9, 8,   4, 5,   7, 6,   19,  20,  19,20]
    
    ci_bridle = np.append(ci_bridle_TE,ci_bridle_LE)
    cj_bridle = np.append(cj_bridle_TE,cj_bridle_LE)
    
    return ci_bridle,cj_bridle

def inflatable_tubes_idx():
    return [0,1,3,6,7,9,12,13,15,18,19,21,24,25,27,30,31,33,36,37,39,42,43,45,48,49,51]

def get_kite_plate_connectivity():
	### Plate connectivity
	plate_1 = [19, 2, 18, 1]
	plate_2 = [2, 3, 17, 18]
	plate_3 = [3, 4, 16, 17]
	plate_4 = [4, 5, 15, 16]
	plate_5 = [5, 6, 14, 15]
	plate_6 = [6, 7, 13, 14]
	plate_7 = [7, 8, 12, 13]
	plate_8 = [8, 9, 11, 12]
	plate_9 = [9, 20, 10, 11]
	plates = [plate_1, plate_2, plate_3, plate_4, plate_5, plate_6, plate_7, plate_8, plate_9]

	ci_kite, cj_kite = [], []
	for i in np.arange(0, len(plates)):
		# The 4 lines describing the tubular frame
		ci_kite.append(plates[i][0]) #LE
		cj_kite.append(plates[i][1]) #LE

		ci_kite.append(plates[i][1]) #Strut right
		cj_kite.append(plates[i][2]) #Strut right

		ci_kite.append(plates[i][2]) #TE
		cj_kite.append(plates[i][3]) #TE

		ci_kite.append(plates[i][3]) #Strut left
		cj_kite.append(plates[i][0]) #Strut left

		# Them diagonals
		ci_kite.append(plates[i][0])
		cj_kite.append(plates[i][2])

		ci_kite.append(plates[i][1])
		cj_kite.append(plates[i][3])
	
        
	ci_kite = np.reshape(ci_kite,len(ci_kite))
	cj_kite = np.reshape(cj_kite, len(cj_kite))

	return ci_kite, cj_kite, plates