# Python program to analyze MJD results with and without drift
# We use the statistical data contained in the repository for different time steps and plot the results.

from matplotlib import pyplot as plt, rcParams
import numpy as np
import pandas as pd
from math import sqrt, exp, pi, log,factorial

# Conversion factor for inches to centimeters (used for figure sizing)
cm = 1 / 2.54
tw = 4.777      # text width in inches
fs = 10         # font size

plt.rcParams.update({
    "text.usetex": True,  # Enable LaTeX rendering for text
    "font.family": "serif",  # Set font family
    "font.serif": ["Times New Roman"],  # Specific serif font
})

# Function to calculate the normal distribution (Gaussian)
def Norm(x, mu, var):
    if var == 0:
        return 1 if x == 0 else 0
    return exp(-(x - mu) ** 2 / (2 * var)) / sqrt(2 * pi * var)

# Define Merton Jump-Diffusion (MJD) PDF
def MJD(S, t, lamda, std1, std2, mu1, mu2, S0, order):
    """Merton Jump-Diffusion model"""
    tot = 0
    for n in range(order + 1):
        term = Norm(np.log(S), np.log(S0) + (mu1 - std1 ** 2 / 2) * t + n * mu2, std1 ** 2 * t + n * std2 ** 2)
        tot += term * (lamda * t) ** n / (S * factorial(n))
    return tot * exp(-lamda * t)

# Parameters
S0 = 10
mu1_values = [-0.4, 0.0]  #
mu2_values = [-0.4/3.0, 0.0]  # Drift values
std = 0.1
lamda=3
order=10
times = [0.1, 0.5, 1.0]  # Time steps
nBins = 24               # Number of bins for histograms

# File paths for data (drift and no drift cases)
path1 = "data_vec/"
drift_file = path1 + "MJD_n120_M10000_v0.5_mu1-0.4_mu2-0.133333333333333_lambda3_S010_std10.1_std20.1_tmax1_vectors.dat"
no_drift_file = path1 + "MJD_n120_M10000_v0.5_mu10_mu20_lambda3_S010_std10.1_std20.1_tmax1_vectors.dat"

# Read data for drift case
df_drift = pd.read_csv(drift_file, header=1, sep=" ")
x01_d, x05_d, x10_d = df_drift["x1"], df_drift["x2"], df_drift["x3"]

# Read data for no drift case
df_no_drift = pd.read_csv(no_drift_file, header=1, sep=" ")
x01_nd, x05_nd, x10_nd = df_no_drift["x1"], df_no_drift["x2"], df_no_drift["x3"]

# Generate x values for plotting
M = 300
xmin = min(min(x01_nd), min(x01_d), min(x05_nd), min(x05_d), min(x10_nd), min(x10_d))
xmax = max(max(x01_nd), max(x01_d), max(x05_nd), max(x05_d), max(x10_nd), max(x10_d))
x = np.linspace(xmin, xmax, M)

# Generate data for macroscopic GBM solution
y_values_d = []
y_values_nd = []

for i in [0,1]:
    y_data = []
    for t in times:
        mu1,mu2 = mu1_values[i],mu2_values[i]
        y = np.array([MJD(x_, t, lamda, std, std, mu1, mu2, S0, order) for x_ in x])
        y_data.append(y)
    if i == 0:
        y_values_d = y_data
    else:
        y_values_nd = y_data
        
# Prepare plot dimensions
largo = 9
alto = 3
dimensioni = np.array([largo, alto]) / largo * tw
fig, axes = plt.subplots(nrows=1, ncols=3, figsize=dimensioni)
axes = axes.flatten()  # Flatten axes for easier access

lww = 1  # Line width for plots
lwww = 0.2  # Line width for histogram borders

# Plotting for each time step

for i, (y_d, y_nd, t) in enumerate(zip(y_values_d, y_values_nd, times)):
    ax = axes[i]
    
    # Plot histograms for drift and no drift
    ax.hist(x01_nd if i == 0 else (x05_nd if i == 1 else x10_nd), bins=nBins, density=True, alpha=0.6, color='tab:blue', edgecolor='black', lw=lwww)
    ax.hist(x01_d if i == 0 else (x05_d if i == 1 else x10_d), bins=nBins, density=True, alpha=0.6, color='tab:red', edgecolor='black', lw=lwww)
    
    # Plot macroscopic GBM curves
    ax.plot(x, y_nd, ls="--", color="darkblue", lw=lww)
    ax.plot(x, y_d, ls="--", color="darkred", lw=lww)
    
    # Set axis labels and limits
    ax.set_xlabel('Price', fontsize=fs)
    ax.set_title(f'$t= {t}$', fontsize=fs)
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(0, 1.45)

# Remove y-axis labels for second and third plots
axes[1].set_yticklabels([])
axes[2].set_yticklabels([])

# Set custom ticks for each plot
axes[0].set_yticks([0, 1])
axes[0].set_yticklabels(["0", "1"], fontsize=10)
axes[0].set_xticks([9.2, 10])
axes[0].set_xticklabels(["9.2", "10"], fontsize=10)
axes[1].set_xticks([6.7, 10])
axes[1].set_xticklabels(["6.7", "10"], fontsize=10)
axes[2].set_xticks([4.5, 10])
axes[2].set_xticklabels(["4.5", "10"], fontsize=10)

# Add legend
#fig.legend(loc='upper center', bbox_to_anchor=(0.5, -0.12), ncol=4, fontsize=7)

# Tight layout and save the figure
plt.tight_layout()
plt.savefig("figures/MJD_with_without_drift.pdf")
plt.show()
