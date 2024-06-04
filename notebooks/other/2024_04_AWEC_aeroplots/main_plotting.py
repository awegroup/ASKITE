import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Reading out .csv files
# df = pd.read_csv("scripts/other/2024_04_AWEC_aeroplots/EKF_results_uri.csv")
df = pd.read_csv("EKF_results_uri.csv")

# filtering the data
df = df[df["turn_straight"] == "straight"]
df = df[df["powered"] == "powered"]

# Assuming df is your DataFrame with 'CL' and 'aoa' columns
# Generate a scatter plot with KDE overlay
sns.scatterplot(
    data=df, x="CL", y="aoa", alpha=0.1, color="darkblue", label="EKF"
)  # alpha controls transparency
sns.kdeplot(
    data=df,
    x="aoa",
    y="CL",
    cmap="viridis",
    levels=25,
    alpha=1.0,
    linewidths=0,
    fill=True,
    # thresh=0.01,
)  # cmap specifies colormap

# 10.7 and 8.4
aoa_sim = [
    7.38,
    7.53,
    7.69,
    7.90,
    8.07,
    8.11,
    9.10,
    9.40,
    9.63,
    10.78,
    11.79,
]
CL_simm = [0.88, 0.89, 0.90, 0.91, 0.92, 0.93, 0.98, 0.99, 1.02, 1.088, 1.154]
plt.plot(aoa_sim, CL_simm, "ro", label="Simulation")

plt.legend()
plt.grid()
plt.xlim([0, 15])
plt.ylim([0.0, 1.8])
plt.xlabel(f"$\\alpha$ [deg]")
plt.ylabel(f"$C_L$ [-]")
plt.title("Density-based scatter plot")
plt.show()
