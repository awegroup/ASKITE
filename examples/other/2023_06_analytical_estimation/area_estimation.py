#%%
import numpy as np

L_to_D = 8.5 #[-] [Bryan]
Va_design = 44.7 #m/s [SP80]
F_fwd = 6800 #N [SP80]
L_design = 25000 #N [SP80]

Sw          = 20  #m^2 [Jelle]
CD_bridle   = 1.2  #[Literature]
R_Sw_Sb     = 1.37/30 #30m2 [Bryan]
D_bridle    = 0.5*1.225*(Va_design**2)*CD_bridle*Sw*R_Sw_Sb #m [SP80]
print("D_bridle :", D_bridle)

D_wing = L_design/L_to_D #N
print("D_wing (from L/D) :", D_wing)

D_kite = D_bridle + D_wing #N
print("D_kite   :", D_kite)

Fa = np.sqrt(L_design**2+D_kite**2) #N
print("Fa       :", Fa)

angle_Fa_to_Fwd = 15.5 #deg
print("angle_Fa_to_Fwd :", angle_Fa_to_Fwd)

F_fwd_from_Fa = Fa*np.sin(np.deg2rad(angle_Fa_to_Fwd)) #N
print("F_fwd    :", F_fwd_from_Fa)

print(' ')

CL = 1     # [V3 data]
AR = 6.7   # [Bryan]
e = 0.8    # [Vlugt, Flysurfer]
CD0 = 0.06 # [Vlugt, Flysurfer]

CD_i = CL**2/(np.pi*AR*e)
print("CD_i     :", CD_i)

CD = CD0 + CD_i
print("CD       :", CD)
D_wing_from_CL = 0.5*1.225*(Va_design**2)*CD*Sw
print("D_wing (calc)  :", D_wing_from_CL)
L_D_calc = CL/CD
print("L/D (calc)  :" , CL/CD)
L_from_drag = D_wing_from_CL*L_to_D
print("L (from drag)  :" , L_from_drag)

A = L_from_drag / (0.5*1.225*(Va_design)**2*CL)

print("A        :", A)
# %%
