#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  9 13:20:49 2025

@author: shar_sp
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 8 09:17:36 2025
@author: shar_sp
"""

import numpy as np
import pickle
import os
import matplotlib.pyplot as plt
from matplotlib.mlab import psd
from multiprocessing import Pool, cpu_count

# Configure Matplotlib to use LaTeX for text rendering
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.serif": ["sans"],
    "axes.labelsize": 22,
    "axes.titlesize": 22,
    "xtick.labelsize": 22,
    "ytick.labelsize": 22,
    "legend.fontsize": 22,
    "figure.titlesize": 16
})

# Create 'postAcoustics' folder if it doesn't exist
output_directory = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'postAcoustics')
os.makedirs(output_directory, exist_ok=True)

# Create 'figures' folder inside 'postAcoustics' if it doesn't exist
figures_directory = os.path.join(output_directory, 'figures')
os.makedirs(figures_directory, exist_ok=True)

# ----------------
# Read deltaT from controlDict and calculate fs
with open('system/controlDict', 'r') as controlDict_file:
    for line in controlDict_file:
        if 'deltaT' in line:
            deltaT_str = line.split()[1].rstrip(';')
            deltaT = float(deltaT_str)
            fs = int(1 / deltaT)
            break

c = 343  # speed of sound
rho0 = 1.2  # density of air
dt = deltaT  # time step
# ---------------------

# Load geometry data
with open('data_geo.pckl', 'rb') as f:
    obj = pickle.load(f)

A = obj[0]  # area of each cell on the surface (face)
n = obj[1]  # normal vector of face
c_xyz = obj[2]  # centroid of A
p_cxyz = obj[3]

# Additional Parameters
num_angles = 101  # Number of angles to consider (0 to 360 degrees)
theta_values = np.linspace(0, 360, num_angles)  # Angles in degrees

# Specify the frequencies of interest
freq_of_interest = [0.001, 0.002, 0.003, 0.004, 0.005, 0.006, 0.008, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.08, 0.1]
freq_of_interest_scaled = [100.0, 200.0, 300.0, 400.0, 500.0, 600.0, 800.0, 1000.0, 2000.0, 3000.0, 4000.0, 5000.0, 6000.0, 8000.0, 10000.0]

# Observer positions (radii)
x0_values = [0.2, 0.4, 0.6, 0.8, 1.]  # Add the x0 values you want to consider

# Progress bar
def print_progress_bar(iteration, total, prefix='', suffix='', decimals=1, length=50, fill='█'):
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filled_length = int(length * iteration // total)
    bar = fill * filled_length + '-' * (length - filled_length)
    print(f'\r{prefix} |{bar}| {percent}% {suffix}', end='\r')
    if iteration == total:
        print()

# Function to process SPL at each angle and radius
def process_angle(params):
    x0, theta = params
    x0_theta = np.array([x0 * np.cos(np.radians(theta)), x0 * np.sin(np.radians(theta)), 0.07])

    F_theta = (-rho0 * A * p_cxyz.T).T
    r_theta = np.sqrt(((c_xyz - x0_theta)[:, 0])**2 + ((c_xyz - x0_theta)[:, 1])**2 + ((c_xyz - x0_theta)[:, 2])**2)

    e_r_theta = np.empty(c_xyz.shape, dtype=float)
    for comp in range(len(x0_theta)):
        e_r_theta[:, comp] = np.divide((x0_theta - c_xyz)[:, comp], r_theta)

    tau_theta = np.divide(r_theta - np.min(r_theta), c)
    p_L_theta = np.ones(F_theta.shape[1] - 6, dtype=float) - 1

    for i, Frow in enumerate(F_theta):
        p_tau_theta = (((45 * Frow[4:-2] - 45 * Frow[2:-4] - 9 * Frow[5:-1] + 9 * Frow[1:-5] - Frow[0:-6] + Frow[6:]).T) *
                       np.dot(e_r_theta[i, :], n[i, :])) * (4 * np.pi * c * r_theta[i] * 60 * dt)**-1
        pL_tau_theta = np.zeros(len(p_tau_theta))
        pL_tau_theta[:len(p_tau_theta)] = p_tau_theta[:len(pL_tau_theta)]
        p_L_theta = p_L_theta + pL_tau_theta

    H_NC_theta = ((np.abs(p_L_theta)) / (2.e-5))**2
    val_NC_theta, freq_NC_theta = psd(H_NC_theta, fs, detrend='mean')

    indices_of_interest = [np.argmin(np.abs(freq_NC_theta - f)) for f in freq_of_interest]
    SPL_values = 5 * np.log10(val_NC_theta[indices_of_interest] + 1e-10)

    return (x0, theta, SPL_values)

# Generate SPL Heatmap (2D Polar Visualization) for all frequencies
def generate_heatmaps(theta_values, x0_values, SPL_results, freq_of_interest_scaled):
    theta_radians = np.radians(theta_values)

    # Calculate global SPL min and max across all data
    SPL_min = np.inf
    SPL_max = -np.inf
    for result in SPL_results:
        # SPL_min = min(SPL_min, np.min(result[2]))
        # SPL_max = max(SPL_max, np.max(result[2]))
        SPL_min = 20
        SPL_max = 70

    print(f"Global SPL Range: {SPL_min:.2f} dB to {SPL_max:.2f} dB")

    for freq_idx, freq_scaled in enumerate(freq_of_interest_scaled):
        SPL_heatmap = np.zeros((len(x0_values), len(theta_values)))

        for x0_idx, x0 in enumerate(x0_values):
            for result in SPL_results:
                if result[0] == x0:
                    theta, SPL = result[1], result[2]
                    idx = np.where(theta_values == theta)[0][0]
                    SPL_heatmap[x0_idx, idx] = SPL[freq_idx]

        # Create polar heatmap
        theta_grid, radius_grid = np.meshgrid(theta_radians, x0_values)

        fig, ax = plt.subplots(figsize=(8, 8), subplot_kw={'projection': 'polar'})
        heatmap = ax.contourf(
            theta_grid, radius_grid, SPL_heatmap, levels=np.linspace(SPL_min, SPL_max, 100), cmap='viridis', vmin=SPL_min, vmax=SPL_max
        )

        # Add colorbar
        cbar = plt.colorbar(heatmap, ax=ax, orientation='horizontal', pad=0.1, shrink=0.8)
        cbar.set_label('SPL (dB)', fontsize=14)
        cbar.ax.tick_params(labelsize=12)

        # Set colorbar ticks and labels to match SPL range
        cbar.set_ticks(np.linspace(SPL_min, SPL_max, 5))  # 5 evenly spaced ticks

        # Customize plot
        ax.set_title(f'SPL Heatmap (Frequency: {freq_scaled:.1f} Hz)', va='bottom', fontsize=16)
        ax.set_theta_zero_location('E')
        ax.set_theta_direction(-1)
        ax.set_rticks(x0_values)

        # Save plot
        filename = f'SPL_Heatmap_Polar_{freq_scaled:.1f}Hz.pdf'
        plt.savefig(os.path.join(figures_directory, filename), bbox_inches='tight', dpi=300)
        plt.show()


if __name__ == '__main__':
    tasks = [(x0, theta) for x0 in x0_values for theta in theta_values]
    total_tasks = len(tasks)

    with Pool(cpu_count()) as pool:
        results = []
        for i, result in enumerate(pool.imap_unordered(process_angle, tasks)):
            results.append(result)
            print_progress_bar(i + 1, total_tasks, prefix='Progress:', suffix='Complete', length=50)

    # Generate heatmaps
    generate_heatmaps(theta_values, x0_values, results, freq_of_interest_scaled)
