# %%
import numpy as np
import matplotlib.pyplot as plt

date = "28jun"

w_lst_8, w_lst_13, u_p_lst = [], [], []
w_lst_8_billow, w_lst_13_billow = [], []
for i in np.arange(0.0, 1.1, 0.1):
    u_p = np.round(i, 1)
    u_p_lst.append(u_p)

    delta_ld_used = 8
    billowing_boolean = False
    filename = (
        "../AerostructuralModel_v"
        + date
        + "/run_results/pos_up_"
        + str(int(100 * (u_p)))
        + "_"
        + str(delta_ld_used)
        + "_"
        + str(billowing_boolean)
        + ".csv"
    )
    # filename = '../AerostructuralModel_v'+date+'/run_results_final_v2/pos_up_'+str(int(100*u_p))+'_8.csv'
    pos = np.loadtxt(filename, delimiter=",")
    w_lst_8.append(pos[1][1] - pos[10][1])

    delta_ld_used = 13
    billowing_boolean = False
    filename = (
        "../AerostructuralModel_v"
        + date
        + "/run_results/pos_up_"
        + str(int(100 * (u_p)))
        + "_"
        + str(delta_ld_used)
        + "_"
        + str(billowing_boolean)
        + ".csv"
    )
    # filename = '../AerostructuralModel_v'+date+'/run_results_final_v2/pos_up_'+str(int(100*u_p))+'_13.csv'
    # filename = '/run_results_final/pos_up_'+str(int(100*u_p))+'_13.csv'
    pos = np.loadtxt(filename, delimiter=",")
    w_lst_13.append(pos[1][1] - pos[10][1])

    delta_ld_used = 8
    billowing_boolean = True
    filename = (
        "../AerostructuralModel_v"
        + date
        + "/run_results/pos_up_"
        + str(int(100 * (u_p)))
        + "_"
        + str(delta_ld_used)
        + "_"
        + str(billowing_boolean)
        + ".csv"
    )
    # filename = '../AerostructuralModel_v'+date+'/run_results_final_v2/pos_up_'+str(int(100*u_p))+'_8_billowing.csv'
    # filename = '/run_results_final/pos_up_'+str(int(100*u_p))+'_13.csv'
    pos = np.loadtxt(filename, delimiter=",")
    w_lst_8_billow.append(pos[1][1] - pos[10][1])

    delta_ld_used = 13
    billowing_boolean = True
    filename = (
        "../AerostructuralModel_v"
        + date
        + "/run_results/pos_up_"
        + str(int(100 * (u_p)))
        + "_"
        + str(delta_ld_used)
        + "_"
        + str(billowing_boolean)
        + ".csv"
    )
    # filename = '../AerostructuralModel_v'+date+'/run_results_final_v2/pos_up_'+str(int(100*u_p))+'_13_billowing.csv'
    # filename = '/run_results_final/pos_up_'+str(int(100*u_p))+'_13.csv'
    pos = np.loadtxt(filename, delimiter=",")
    w_lst_13_billow.append(pos[1][1] - pos[10][1])


### PLOTTING THE RESULTS

# defining settings
label_font_size = 12
legend_font_size = 11
dot_size = 2
line_font_size = 2.5
factor_white_dashes = 0.8

u_p_lst = np.array(u_p_lst)
w_lst_8 = 1e-3 * np.array(w_lst_8).T
w_lst_13 = 1e-3 * np.array(w_lst_13).T
w_lst_8_billow = 1e-3 * np.array(w_lst_8_billow).T
w_lst_13_billow = 1e-3 * np.array(w_lst_13_billow).T

## Re-adjusted the initial value as it somehow went wrong there
# w_lst_8[-1] = 8.203240342440058
# w_lst_13[-1] = 8.203240342439997
# w_lst_8_billow[-1] = 8.203240342440058
# w_lst_13_billow[-1] = 8.203240342439997

w_lst_tetrahedon_8 = [
    8.005754349090248,
    8.034376497759158,
    8.06223228758389,
    8.089345199872435,
    8.115737556829757,
    8.141430587137126,
    8.166444486711464,
    8.190798475050514,
    8.214510847526352,
    8.237599023955054,
    8.260079593739977,
]
w_lst_tetrahedon_13 = [
    7.756243235479863,
    7.809213443728077,
    7.860090637834112,
    7.9089835408530105,
    7.955992578911484,
    8.001210649154606,
    8.044723799023792,
    8.086611828758926,
    8.126948827167995,
    8.165803649181093,
    8.203240342439997,
]
# plt.plot(u_p_lst, w_lst_tetrahedon_8 , label=r'Tetrahedon ($\delta_{\rm{d}}$ = 8%), from presimulated' , color='green' , linewidth=line_font_size*.8, linestyle ='-')
# plt.plot(u_p_lst, w_lst_tetrahedon_13 , label=r'Tetrahedon ($\delta l_d$ = 13%), from presimulated' , color='green' , linewidth=line_font_size*.5, linestyle ='dashed')

# plt.plot(u_p_lst, w_lst_8 , label=r'PSM ($\delta_{\rm{d}}$ = 8%), excl. billowing' , color='blue' , linewidth=line_font_size, linestyle ='-')
# plt.plot(u_p_lst, w_lst_8_billow , label=r'_PSM ($\delta_{\rm{d}}$ = 8%), incl. billowing' , color='blue' , linewidth=line_font_size*factor_white_dashes, linestyle ='dashed')
plt.plot(
    u_p_lst,
    w_lst_8_billow,
    label=r"Particle System Model",
    color="blue",
    linewidth=line_font_size * factor_white_dashes,
    linestyle="dashed",
)

# plt.plot(u_p_lst, w_lst_13, label=r'PSM ($\delta_{\rm{d}}$ = 13%), excl. billowing', color='red'  , linewidth=line_font_size, linestyle ='-')
# plt.plot(u_p_lst, w_lst_13_billow, label=r'_PSM ($\delta_{\rm{d}}$ = 13%), incl. billowing', color='red'  , linewidth=line_font_size*factor_white_dashes, linestyle ='dashed')

### Photogrammetry measurements
data_photogrammetry_width = np.array(
    [
        8.2377444,
        8.2114677,
        8.25745192,
        8.25088275,
        8.2377444,
        7.86407394,
        7.84413124,
        8.0435583,
        7.96378748,
        8.00367289,
        8.34285119,
        8.39540458,
        8.2377444,
        8.17205266,
        8.25745192,
        7.7922407,
        7.91707424,
        7.68054754,
        7.68054754,
        7.67397736,
    ]
)
data_photogrammetry_up = [
    1.0,
    1.0,
    1.0,
    1.0,
    1.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    1.0,
    1.0,
    1.0,
    1.0,
    1.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
]

plt.plot(
    data_photogrammetry_up,
    data_photogrammetry_width,
    "o",
    color="white",
    linewidth=dot_size,
    mec="black",
    label="Photogrammetry measurements",
)

# configuring the plot and saving it
plt.grid(True)
plt.xlabel(r"Power Setting $(u_{\rm{p}})\,[-]$", fontsize=label_font_size)
plt.ylabel(r"Kite Width $\,[m]$", fontsize=label_font_size)
plt.ylim(7.4, 8.6)
plt.legend(fontsize=legend_font_size, loc="lower right", bbox_to_anchor=(1.0, 0.0))
plt.rcParams["pdf.fonttype"] = 42  # Output Type 3 (Type3) or Type 42 (TrueType)
plt.savefig("results_PSM_graph.svg")
plt.savefig("results_PSM_graph.pdf")
plt.show()

# %%
