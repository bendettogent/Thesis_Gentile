# Python program to make the analysis of Donsken results.
# We use the statistical data contained in repository data_vec and the cumulants in repository data.

from matplotlib import pyplot as plt, rcParams
import numpy as np
import pandas as pd
from math import factorial, sqrt, exp, pi, log

# Conversion factor for inches to centimeters (used for figure sizing)
cm = 1 / 2.54
tw = 4.777  # text width in inches
fs = 10  # font size

plt.rcParams.update({
    "text.usetex": True,  # Enable LaTeX rendering for text
    "font.family": "serif",  # Set font family
    "font.serif": ["Times New Roman"],  # Specific serif font
})

# Parameters
M = 10000  # Number of Monte Carlo runs
times = [0.1, 0.5, 1.0]  # Time steps
S0 = 50  # Initial stock price
nBins = 24  # Number of bins for histograms
means = [0, 5]  # Mean values

# DataFrames to store results
dF = pd.DataFrame()
dFdata = pd.DataFrame()

# Variables for tracking min and max values across all data
MAXMAX = S0
MINMIN = S0

path1 = "data_vec/"
# Loop over mean values to read data and compute statistics
for i, mean in enumerate(means):
    # Read data for each mean value
    path_file = path1 + "/Donsken_stat_n40_M10000_mean" + str(mean) + "_S050_vectors.dat"
    df = pd.read_csv(path_file, header=1, sep=" ")
    
    # Extract columns and stack them
    x1, x2, x3 = np.array(df["x1"]), np.array(df["x2"]), np.array(df["x3"])
    stacked_matrices = np.stack((x1, x2, x3))
    
    # Update global min and max values
    MAX = max(max(x1), max(x3), max(x2))
    MIN = min(min(x1), min(x3), min(x2))
    MAXMAX = max(MAX, MAXMAX)
    MINMIN = min(MIN, MINMIN)

    # Loop to calculate histograms for each dataset (x1, x2, x3)
    for j, data in enumerate([x1, x2, x3], start=1):
        counts, bins = np.histogram(data, bins=nBins, density=True)
        bin_centers = 0.5 * (bins[:-1] + bins[1:])
        
        dFdata[f"data{j}_mean{mean}"] = data
        dF[f"counts{j}_mean{mean}"] = counts
        dF[f"bins{j}_mean{mean}"] = bin_centers

# Define a function to calculate normal distribution (Gaussian)
def Norm(x_, a, var):
    return ((2 * pi * var) ** (-0.5) * np.exp(-(x_ - a) ** 2 / (var * 2)))


# Read data with no drift
df = pd.read_csv(path1 + "Donsken_stat_n40_M10000_mean0_S050_vectors.dat", header=1, sep=" ")
x01_nd, x05_nd, x10_nd = np.array(df["x1"]), np.array(df["x2"]), np.array(df["x3"])

# Read data with drift
df = pd.read_csv(path1 + "Donsken_stat_n40_M10000_mean5_S050_vectors.dat", header=1, sep=" ")
x01_d, x05_d, x10_d = np.array(df["x1"]), np.array(df["x2"]), np.array(df["x3"])

# Generate linearly spaced x-values for the normal distribution curve
M = 300  # Number of points
xmin = min(min(x01_nd), min(x01_d), min(x05_nd), min(x05_d), min(x10_nd), min(x10_d))
xmax = max(max(x01_nd), max(x01_d), max(x05_nd), max(x05_d), max(x10_nd), max(x10_d))
x = np.linspace(xmin, xmax, M)

# Generate normal distributions for different time points and drifts using a loop
y_d = []
y_nd = []
for mu, drift in zip([5, 0], ['d', 'nd']):  # mu=5 with drift, mu=0 without drift
    for t in times:  # Loop over times
        var = 1 * t
        y = np.array([Norm(x_, S0 + mu * t, var) for x_ in x])
        if drift == 'd':
            y_d.append(y)
        else:
            y_nd.append(y)

# Prepare plot dimensions
largo = 9
alto = 4
dimensioni = np.array([largo, alto]) / largo * tw
fig, axes = plt.subplots(nrows=1, ncols=3, figsize=dimensioni)
axes = axes.flatten()  # Flatten axes for easier access

lww = 1  # Line width for plots
lwww = 0.2  # Line width for histogram borders

# Plotting in a loop to handle the three time points
for i, t in enumerate(times):
    ax = axes[i]
    
    # Plot histograms for no drift and drift
    ax.hist(x01_nd if i == 0 else (x05_nd if i == 1 else x10_nd), bins=25, density=True, alpha=0.6, color='tab:blue', edgecolor='black', lw=lwww)
    ax.hist(x01_d if i == 0 else (x05_d if i == 1 else x10_d), bins=nBins, density=True, alpha=0.6, color='tab:red', edgecolor='black', lw=lwww)
    
    # Plot normal distribution curves
    ax.plot(x, y_nd[i], ls="--", color="blue", lw=lww)
    ax.plot(x, y_d[i], ls="--", color="darkred", lw=lww)
    
    # Set axis labels, title, and limits
    ax.set_xlabel('Price', fontsize=fs)
    ax.set_title(f'$t= {t}$', fontsize=fs)
    ax.set_xlim(min(x), max(x))
    ax.set_ylim(0, 1.45)

# Remove y-axis tick labels for the second and third plots
axes[1].set_yticklabels([])
axes[2].set_yticklabels([])

# Tight layout for clean figure formatting
plt.tight_layout()

# Save figure to file
plt.savefig("figures/Donsker_with_without_drift.pdf")
plt.show()