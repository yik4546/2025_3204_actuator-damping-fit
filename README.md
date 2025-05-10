# Actuator Fitting Project

This project contains MATLAB code for fitting actuator dynamics based on experimental data. The main goal is to optimize the damping coefficient `c` of the actuator model by comparing simulated outputs with experimental data.

## Description

This repository contains a MATLAB script for fitting an actuator model to experimental data. The script performs the following tasks:

- Reads CSV files containing experimental time-series data.
- Uses an optimization algorithm (`fminsearch`) to estimate the best damping coefficient for the actuator model.
- Simulates the actuator's dynamics and compares the simulated output to the experimental data.
- Calculates a matching percentage between the experimental and simulated data, and visualizes the fitting results.

### Key Features

- **Actuator Model Simulation**: Simulates the dynamics of an actuator based on physical parameters.
- **Optimization**: Uses the `fminsearch` function to find the optimal damping coefficient (`c`).
- **Data Processing**: Automatically processes experimental data to remove duplicates and normalize the time series.
- **Visualization**: Plots the experimental and fitted data for visual comparison.
- **Matching Percentage**: Calculates a matching percentage to assess the quality of the fit.

## Getting Started

To run the script, follow these steps:

1. Clone the repository:
   ```bash
   git clone https://github.com/your-username/actuator-fitting.git
   cd actuator-fitting