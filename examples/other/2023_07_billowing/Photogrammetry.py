#%%
import openpyxl as opn
# from openpyxl import load_workbook
#from Variables_Literature import *

#%% Retrieving data from excel

# Opening the workbook
wb = opn.load_workbook(filename='billowing/Photogrammetry_data.xlsx', read_only=True)

# Getting the widths
kcu_pow_w,kcu_dep_w,launch_pow_w,launch_dep_w = [],[],[],[]
for row in [2,3,4,5,6]:
     kcu_pow_w.append(wb['kcu_pow']['B'+str(row)].value)

for row in [21,22,23,24,25]:
     kcu_dep_w.append(wb['kcu_dep']['C'+str(row)].value)

for row in [2,3,4,5,6]:
     launch_pow_w.append(wb['launch_pow']['B'+str(row)].value)

for row in [21,22,23,24,25]:
     launch_dep_w.append(wb['launch_dep']['C'+str(row)].value)

total_averaged_width = wb['kcu_dep']['C51'].value

#Getting the Ballooning
# ballooning_plate_1 = wb['kcu_dep']['D46'].value
# ballooning_plate_2 = wb['kcu_dep']['F46'].value
# ballooning_plate_3 = wb['kcu_dep']['H46'].value
# ballooning_plate_4 = 0 #output is an outlier, thus 0
ballooning_plate_1 = -0.056
ballooning_plate_2 = -0.0216  
ballooning_plate_3 = +0.01
ballooning_plate_4 = 0#-0.1258 #output is an outlier, thus 0

# Getting pulley & knots
pulley_pow_percentage = wb['kcu_dep']['O51'].value
pulley_dep_percentage = wb['kcu_dep']['O52'].value
knot_pow_percentage = wb['kcu_dep']['R51'].value
knot_dep_percentage = wb['kcu_dep']['R52'].value

#Closing the workbook after reading
wb.close()

#%% Usefull functions based on measurements, 'correlation based'

#don't forget need to add a factor of 2 to all models
#u_p = u_p  # Since the robe splits up, effect will be halved for each pulley

def ballooning_up(u_p):
     ballooning_plate_1 = -0.056
     ballooning_plate_2 = -0.0216  
     ballooning_plate_3 = +0.01
     ballooning_plate_4 = 0#-0.1258 #output is an outlier, thus 0
     
     return [1+ ballooning_plate_1*(1-u_p), 1+ ballooning_plate_2*(1-u_p), 1+ ballooning_plate_3*(1-u_p),1 + ballooning_plate_4*(1-u_p)]

def ballooning_up_perc(u_p):
    return [ballooning_plate_1*(1-u_p)*100,ballooning_plate_2*(1-u_p)*100, ballooning_plate_3*(1-u_p)*100]

def plotting_photogrammetry_straight(line_font_size,dot_font_size,color):
    import matplotlib.pyplot as plt
    y_kcu = np.append(kcu_pow_w,kcu_dep_w)
    y_launch = np.append(launch_pow_w,launch_dep_w)
    #y = y_kcu
    y = np.append(y_kcu,y_launch)

    x_half = np.append(np.ones(5),np.zeros(5))
    #x = x_half
    x = np.append(x_half,x_half)

    y = y * (w_powered/ (total_averaged_width * 1000) )

    # fit a linear curve an estimate its y-values and their error.
    a, b = np.polyfit(x, y, deg=1)
    y_est = a * x + b
    #y_err = y_err / np.sqrt(len(x))
    y_err = y.std() / np.sqrt(len(x)/2)

    plt.plot(x, y_est,linestyle='-',linewidth=line_font_size,  color=color, label = 'Photogrammetry linearized mean')
    plt.fill_between(x, y_est - y_err, y_est + y_err, alpha=0.2, color=color, label='Photogrammetry confidence band')
    #plt.plot(x, y, 'o', color=color,linewidth = dot_font_size,label = 'Photogrammetry measurements')
    #print('std: ',y.std())
    #print('SE: ', y_err)
    print(y_est)
    return


#  #%% Cornering for Future Reference
#  
#  def us_to_ls(u_s):
#      #u_s = -0.5 -> min_steer and u_s = 0.5 -> max_steer
#      # negative u_s means steering to the right!
#      min_steer, max_steer = -40,40  # Defining the min/max steering numbers observed from the measurements (averaged)
#      ratio = (max_steer-min_steer)*u_s #This will be number between 0 and (-)80
#      ls = delta_ls_max*(ratio/100)
#      ls_left,ls_right = ls_0 - ls, ls_0 + ls
#      return ls_left,ls_right
#  
#  
#  # Validation data - turning
#  
#  #Should case about a 4.9% difference, which is does now.
#  w_turn_us_0_5 = w_powered/2  *0.975
#  w_non_turn_us_0_5 = w_powered/2 *1.025
#  
#  #print(1-w_turn_us_0_5/w_non_turn_us_0_5)
#  
#  def ballooning_us(u_s):
#      maxplate_1 = 1.043 - 1
#  
#      maxplate_2_non = 1.023 -1
#      maxplate_3_non = 0.997 -1
#      maxplate_4_non = 1.009 -1
#  
#      maxplate_2_turn = 1.071 -1
#      maxplate_3_turn = 1.050 -1
#      maxplate_4_turn = 1.054 -1
#  
#      return [1 + maxplate_4_non*u_s*2 , 1 + maxplate_3_non*u_s*2,  1 + maxplate_2_non*u_s*2,
#      1 + maxplate_1*u_s*2,
#      1 + maxplate_4_turn*u_s*2 , 1 + maxplate_3_turn*u_s*2,  1 + maxplate_2_turn*u_s*2]
#      
#  
