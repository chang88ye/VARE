# VARE
Vector Autoregressive Evolution for Dynamic Multi-Objective Optimisation

---------------------------------------
A MATLAB implementation of four evolutionary algorithms for Dynamic Multi-objective Optimization

- VARE/VARE_MOEAD: Vector AutoRegressive Evolution

- SGEA: Steay and Generational Evolutionary Algorithm

- PPS: Population Prediction Strategy

- MOEASVR: Support Vector Regression with MOEA/D

- TrRMMEDA: Transfer Learning with Regularity Model based Estimation of Distribution

It includes the following files:

- four algorithms mentioned above

- A public folder that contains public functions used by four algorithms

- Indicators: performance assessment related functions

- main.m: the main entry of running experiments

- run_VARE_10_10.sh: batch file to run experiments in a distributed computing system such HPC

Running the code will create new folders (e.g., Results_nt10_taut10) to save experimental results.

Please Note: 

The current implementation does not work for time-varying number of objectives/variables, i.e. SDP12, SDP13. 

The code was primarily tested in Linux/Mac with MATLAB 2020a/2022a. However, bugs may occur in a different operating system or other MATLAB versions. Any report on bugs is welcome and suggestions are much appreciated.

# Cite Us
Shouyong Jiang, Yong Wang, Yaru Hu, Qingyang Zhang, andShengxiang Yang (2023). Vector Autoregressive Evolution for Dynamic
Multi-Objective Optimisation.  https://doi.org/10.48550/arXiv.2305.12752
