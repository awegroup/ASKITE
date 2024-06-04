#%%
import numpy as np
filename = './run_results/pos_up_100_plate_aero_dynamic_blieb.csv'
pos_initial_guess = np.loadtxt(filename,delimiter = ',')
width_lst = [np.round(((pos_initial_guess[1][1]-pos_initial_guess[10][1])),2)]

for up_it in np.arange(0,1.1,0.1):
    u_p = round(1-up_it,1)
    filename = './run_results/pos_up_'+str(int(100*(u_p)))+'_plate_aero_dynamic.csv'
    pos_width = np.loadtxt(filename,delimiter = ',')
    width_lst.append(np.round(((pos_width[1][1]-pos_width[10][1])),2))

for width in width_lst[1:]:
    print(width)
# %%
