# Price Models in Regulated Financial Markets: A Multiscale Approach

**Author:** Benedetto Gentile  
**Supervisors:** E. Scalas, A. Vulpiani  
**University:** Università La Sapienza di Roma  
**Year:** 2024

---

## Description

This repository contains the code, data, and documentation related to my thesis, titled **Price Models in Regulated Financial Markets: A Multiscale Approach**. The project explores various price modeling techniques in financial markets, with a focus on multiscale approaches.

## Repository Structure

The repository is organized as follows:
├── src/                # Source code \\
│   ├── BS              # Black-Scholes models \\
│   ├── Donsken         # Donsker models \\
│   └── MJD             # Merton Jump Diffusion models \\
└── README              # This file \\ 

In each of the three model directories (`BS`, `Donsken`, `MJD`), the following files and folders are present:

├── [model]/
│   ├── [model]_base.R       # R script to run a single simulation and generate plots
│   ├── [model]_stat.R       # R script for statistical simulations
│   ├── data/                # Folder containing cumulants used for statistical simulations 
│   │                        # (these results were used in the thesis)
│   ├── data_vec/            # Folder containing individual position data at specific time points 
│   │                        # for each statistical simulation (these results were used in the thesis)
│   ├── [model]_analisi.py   # Python script to plot data with and without the mean
│   ├── [model]_kstest.py    # Python script to perform the KS test and plot cumulative graphs
│   ├── figures/             # Folder containing Python-generated figures


The files in the `data` and `data_vec` directories contain the results used in my thesis.

**CAVEAT:** If you wish to run a new statistical simulation, you must change the file names in the `.py` scripts accordingly.
