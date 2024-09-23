import numpy as np
from scipy.integrate import quad
import matplotlib.pyplot as plt
import pandas as pd
from math import factorial, sqrt, exp, pi

# Set plot configuration for LaTeX style
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.serif": ["Times New Roman"],
})

# Define the Gaussian PDF
def Norm(x, mu, var):
    """Normal (Gaussian) distribution PDF"""
    return exp(-(x - mu) ** 2 / (2 * var)) / sqrt(2 * pi * var)

# Compute empirical CDF for a given x and data vector
def cdf_empirica(x, data):
    """Compute empirical CDF for a given value x"""
    return np.mean(data <= x)

# Compute theoretical CDF by integrating a PDF
def cdf_from_pdf(pdf, x, discontinuous=False, mean=0, factor_=1, lower_limit=0.001):
    """Calculate the CDF from the given PDF by numerical integration"""
    if discontinuous:
        if x < mean:
            result, _ = quad(pdf, lower_limit, x)
        else:
            result, _ = quad(pdf, lower_limit, x)
            result += np.exp(-factor_)
    else:
        result, _ = quad(pdf, lower_limit, x)
    return result

# Define Merton Jump-Diffusion (MJD) PDF
def MJD(S, t, lamda, std1, std2, mu1, mu2, S0, order):
    """Merton Jump-Diffusion model"""
    tot = 0
    for n in range(order + 1):
        term = Norm(np.log(S), np.log(S0) + (mu1 - std1 ** 2 / 2) * t + n * mu2, std1 ** 2 * t + n * std2 ** 2)
        tot += term * (lamda * t) ** n / (S * factorial(n))
    return tot * exp(-lamda * t)


# File with drift
path1 = "data_vec/"
path_file_drift = path1 + "MJD_n120_M10000_v0.5_mu1-0.4_mu2-0.133333333333333_lambda3_S010_std10.1_std20.1_tmax1_vectors.dat"
df_drift = pd.read_csv(path_file_drift, header=1, sep=" ")

x01_d = df_drift["x1"]
x05_d = df_drift["x2"]
x10_d = df_drift["x3"]

# File no drift
path_file_no_drift = path1 + "MJD_n120_M10000_v0.5_mu10_mu20_lambda3_S010_std10.1_std20.1_tmax1_vectors.dat"
df_no_drift = pd.read_csv(path_file_no_drift, header=1, sep=" ")

x01_nd = df_no_drift["x1"]
x05_nd = df_no_drift["x2"]
x10_nd = df_no_drift["x3"]


# Function to plot empirical vs theoretical CDF
def plot_cdf_comparison(vec, t, lamda, std1, std2, mu1, mu2, S0, order, title, ax):
    # Definisco la pdf teorica
    def MJD_(x):
        return MJD(x, t, lamda, std1, std2, mu1, mu2, S0, order)

    # Dati ordinati
    data = np.sort(vec)

    # Calcolo CDF empirica
    cdf_emp = np.array([cdf_empirica(x, data) for x in data])

    # Calcolo CDF teorica
    cdf_teo = np.array([cdf_from_pdf(MJD_, x) for x in data])

    # Calcolo della distanza di Kolmogorov-Smirnov (D)
    D_ = np.max(np.abs(cdf_emp - cdf_teo))

    # Plot
    ax.plot(data, cdf_emp, c="red", lw=1, label="Empirical")
    ax.plot(data, cdf_teo, c="blue", label="Theoretical")
    ax.set_title(f"{title} (D = {D_:.5f})")
    ax.legend()
    
# Parametters
mu1_d = -0.4
mu2_d = -0.4/3.0
mu1_nd = 0
mu2_nd = 0
lamda=3
order=10

S0 = 10
std1 = 0.1
std2 = 0.1
t_values = [0.1, 0.5, 1.0]

# Creazione dei sottoplot
fig, axs = plt.subplots(3, 2, figsize=(12, 10))

# Titoli dei plot
titles = [
    r"$b,\mu =0, t=0.1$", r"$b,\mu=0, t=0.5$", r"$b,\mu=0, t=1.0$",
    r"$b,\mu \neq 0, t=0.1$", r"$b,\mu \neq 0, t=0.5$", r"$b,\mu \neq 0, t=1.0$"
]

# Confronto per il caso senza drift
vecs_no_drift = [x01_nd, x05_nd, x10_nd]
for i in range(3):
    plot_cdf_comparison(vecs_no_drift[i], t_values[i], lamda, std1, std2, mu1_nd, mu2_nd, S0, order, titles[i], axs[i, 0])
# Confronto per il caso con drift
vecs_drift = [x01_d, x05_d, x10_d]
for i in range(3):
    plot_cdf_comparison(vecs_drift[i], t_values[i], lamda, std1, std2, mu1_d, mu2_d, S0, order, titles[i + 3], axs[i, 1])

# Regolazione layout e salvataggio
fig.tight_layout()
fig.savefig("figures/MJD_KStest_n120.pdf")

# Mostra il plot
plt.show()
