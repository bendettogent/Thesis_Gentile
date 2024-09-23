import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.integrate import quad
from math import sqrt, exp, pi

# Impostazioni grafiche
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.serif": ["Times New Roman"],
})

def Norm(x, mu=0, var=1):
    """Gaussian PDF"""
    return 1 / sqrt(2 * pi * var) * np.exp(-(x - mu)**2 / (2 * var))

def cdf_empirica(x, data):
    """Compute empirical CDF at x"""
    return np.mean(data <= x)

def cdf_from_pdf(pdf, x, mean=0, factor_=1, lower_limit=-np.inf):
    """Compute theoretical CDF by integrating the given PDF"""
    result, _ = quad(pdf, lower_limit, x)
    return result

# File paths
path1 = "data_vec/"
files = {
    "n40_nd": path1 + "Donsken_stat_n40_M10000_mean0_S050_vectors.dat",
    "n40_d": path1 + "Donsken_stat_n40_M10000_mean5_S050_vectors.dat",
    "n120_nd": path1 + "Donsken_stat_n120_M10000_mean0_S050_vectors.dat",
    "n120_d": path1 + "Donsken_stat_n120_M10000_mean5_S050_vectors.dat"
}

# Read data for both n=40 and n=120 cases
data_n40_nd = pd.read_csv(files["n40_nd"], header=1, sep=" ")
data_n40_d = pd.read_csv(files["n40_d"], header=1, sep=" ")
data_n120_nd = pd.read_csv(files["n120_nd"], header=1, sep=" ")
data_n120_d = pd.read_csv(files["n120_d"], header=1, sep=" ")

# Extract x1, x2, x3 for n=40 and n=120
x_n40_nd = [data_n40_nd["x1"], data_n40_nd["x2"], data_n40_nd["x3"]]
x_n40_d = [data_n40_d["x1"], data_n40_d["x2"], data_n40_d["x3"]]
x_n120_nd = [data_n120_nd["x1"], data_n120_nd["x2"], data_n120_nd["x3"]]
x_n120_d = [data_n120_d["x1"], data_n120_d["x2"], data_n120_d["x3"]]

# Set time values and corresponding variances
t_values = [0.1, 0.5, 1.0]
var_values = [t for t in t_values]

# Function to plot empirical vs theoretical CDF
def plot_cdf_comparison(vec, t, mean, title, ax):
    varia = t * 1  # Variance as a function of time
    def Norm_(x):
        return Norm(x, mu=mean, var=varia)

    data = np.sort(vec)
    cdf_emp = np.array([cdf_empirica(x, data) for x in data])
    cdf_teo = np.array([cdf_from_pdf(Norm_, x) for x in data])
    D_ = np.max(np.abs(cdf_emp - cdf_teo))

    ax.plot(data, cdf_emp, c="red", lw=1, label="Empirical")
    ax.plot(data, cdf_teo, c="blue", label="Theoretical")
    ax.set_title(f"{title} (D = {D_:.5f})")
    ax.legend()

# Create subplots
fig1, axs1 = plt.subplots(3, 2, figsize=(12, 10))  # First 6 plots
fig2, axs2 = plt.subplots(3, 2, figsize=(12, 10))  # Last 6 plots

# Titles for each plot
titles = [
    r"$n=40, \mu=0, t=0.1$", r"$n=40, \mu=0, t=0.5$", r"$n=40, \mu=0, t=1.0$",
    r"$n=40, \mu=5, t=0.1$", r"$n=40, \mu=5, t=0.5$", r"$n=40, \mu=5, t=1.0$",
    r"$n=120, \mu=0, t=0.1$", r"$n=120, \mu=0, t=0.5$", r"$n=120, \mu=0, t=1.0$",
    r"$n=120, \mu=5, t=0.1$", r"$n=120, \mu=5, t=0.5$", r"$n=120, \mu=5, t=1.0$"
]

# Plot the first 6 cases in fig1
for i in range(3):  # n=40, mu=0
    plot_cdf_comparison(x_n40_nd[i], t_values[i], 50, titles[i], axs1[i, 0])
for i in range(3):  # n=40, mu=5
    plot_cdf_comparison(x_n40_d[i], t_values[i], 50 + 5 * t_values[i], titles[i + 3], axs1[i, 1])

# Plot the next 6 cases in fig2
for i in range(3):  # n=120, mu=0
    plot_cdf_comparison(x_n120_nd[i], t_values[i], 50, titles[i + 6], axs2[i, 0])
for i in range(3):  # n=120, mu=5
    plot_cdf_comparison(x_n120_d[i], t_values[i], 50 + 5 * t_values[i], titles[i + 9], axs2[i, 1])

# Adjust layout
fig1.tight_layout()
fig2.tight_layout()

fig1.savefig("figures/Donsken_KStest_n40.pdf")
fig2.savefig("figures/Donsken_KStest_n120.pdf")

# Show both figures
plt.show()

