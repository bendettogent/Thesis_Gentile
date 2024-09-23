import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.integrate import quad
from math import sqrt, exp, pi, log

# Impostazioni grafiche
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.serif": ["Times New Roman"],
})

# Funzione Gaussiana
def Norm(x, mu, var):
    """Funzione normale"""
    if var == 0:
        return 1 if x == 0 else 0
    return exp(-(x - mu)**2 / (2 * var)) / sqrt(2 * pi * var)

# Funzione macro per il moto browniano geometrico (GBM)
def macro_GBM(x, t, mu, std, S0):
    var = std ** 2
    return Norm(log(x), t * (mu - var * 0.5) + log(S0), var * t) / x

# Funzione per calcolare la cdf empirica
def cdf_empirica(x, data):
    return np.mean(data <= x)

# Funzione per calcolare la cdf teorica a partire da una pdf
def cdf_from_pdf(pdf, x, lower_limit=0.00001):
    result, _ = quad(pdf, lower_limit, x)
    return result

# Dati forniti
S0 = 10

# File con drift
path1 = "data_vec/"
path_file_drift = path1 + "BSdrift_n120_M10000_v0.5_mu-0.8_S010_std0.1_tmax1_vectors.dat"
df_drift = pd.read_csv(path_file_drift, header=1, sep=" ")

x01_d = df_drift["x1"]
x05_d = df_drift["x2"]
x10_d = df_drift["x3"]

# File senza drift
path_file_no_drift = path1 + "BSdrift_n120_M10000_v0.5_mu0_S010_std0.1_tmax1_vectors.dat"
df_no_drift = pd.read_csv(path_file_no_drift, header=1, sep=" ")

x01_nd = df_no_drift["x1"]
x05_nd = df_no_drift["x2"]
x10_nd = df_no_drift["x3"]

# Funzione per il confronto e plot tra CDF empirica e teorica
def plot_cdf_comparison(vec, t, mu, std, S0, title, ax):
    # Definisco la pdf teorica
    def macro_GBM_(x):
        return macro_GBM(x, t, mu, std, S0)

    # Dati ordinati
    data = np.sort(vec)

    # Calcolo CDF empirica
    cdf_emp = np.array([cdf_empirica(x, data) for x in data])

    # Calcolo CDF teorica
    cdf_teo = np.array([cdf_from_pdf(macro_GBM_, x) for x in data])

    # Calcolo della distanza di Kolmogorov-Smirnov (D)
    D_ = np.max(np.abs(cdf_emp - cdf_teo))

    # Plot
    ax.plot(data, cdf_emp, c="red", lw=1, label="Empirical")
    ax.plot(data, cdf_teo, c="blue", label="Theoretical")
    ax.set_title(f"{title} (D = {D_:.5f})")
    ax.legend()

# Parametri per i casi con e senza drift
mu_drift = -0.8
mu_no_drift = 0
S0 = 10
std = 0.1
t_values = [0.1, 0.5, 1.0]

# Creazione dei sottoplot
fig, axs = plt.subplots(3, 2, figsize=(12, 10))

# Titoli dei plot
titles = [
    r"$\mu=0, t=0.1$", r"$\mu=0, t=0.5$", r"$\mu=0, t=1.0$",
    r"$\mu=-0.8, t=0.1$", r"$\mu=-0.8, t=0.5$", r"$\mu=-0.8, t=1.0$"
]

# Confronto per il caso senza drift
vecs_no_drift = [x01_nd, x05_nd, x10_nd]
for i in range(3):
    plot_cdf_comparison(vecs_no_drift[i], t_values[i], mu_no_drift, std, S0, titles[i], axs[i, 0])

# Confronto per il caso con drift
vecs_drift = [x01_d, x05_d, x10_d]
for i in range(3):
    plot_cdf_comparison(vecs_drift[i], t_values[i], mu_drift, std, S0, titles[i + 3], axs[i, 1])

# Regolazione layout e salvataggio
fig.tight_layout()
fig.savefig("figures/BS_KStest_n120.pdf")

# Mostra il plot
plt.show()
