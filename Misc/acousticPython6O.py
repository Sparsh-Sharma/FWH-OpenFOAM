#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 13 16:33:40 2019

@author: Sparsh Sharma
"""
import numpy as np
import pickle
import os
from matplotlib.mlab import psd
import matplotlib.pyplot as plt
from matplotlib.pyplot import grid
from matplotlib.ticker import ScalarFormatter


# ----------------
fs = 100000  # sampling frequency
c = 343  # speed of sound
x0 = np.array([0., 1., -0.4])  # real observer position
rho0 = 1.2  # density of air
dt = fs**-1  # timesep
# ---------------------

f = open('data_geo.pckl', 'rb')  # wurde so gespeichert: pickle.dump([A, n, c_xyz, p_cxyz], f)
obj = pickle.load(f)
f.close()
p_cxyz = obj[3]

A = obj[0]  # area of each cell on the surface (face)
n = obj[1]  # normal vetor of face
c_xyz = obj[2]  # centroid of A

# calculate the force for every face
F = (-rho0 * A * p_cxyz.T).T

# calculate distance r from centroid to receiver location x0
r = np.sqrt(((c_xyz - x0)[:, 0])**2 + ((c_xyz - x0)[:, 1])**2 + ((c_xyz - x0)[:, 2])**2)

# unity vector in the direction r
e_r = np.empty(c_xyz.shape, dtype=float)
for comp in range(len(x0)):
    e_r[:, comp] = np.divide((x0 - c_xyz)[:, comp], r)

# retarded time tau0, find the difference to the minimum distance of all points an divided by propagation speed
tau0 = np.divide(r - np.min(r), c)

# for every time step add all points in the observer point x0
p_L = np.ones(F.shape[1] - 6, dtype=float) - 1  # -6 in shape because sixth order

for i, Frow in enumerate(F):
    p_tau = (((45 * Frow[4:-2] - 45 * Frow[2:-4] - 9 * Frow[5:-1] + 9 * Frow[1:-5] - Frow[0:-6] + Frow[6:]).T) *
             np.dot(e_r[i, :], n[i, :])) * (4 * np.pi * c * r[i] * 60 * dt)**-1  # 6order

    pL_tau = np.zeros(len(p_tau))
    pL_tau[:len(p_tau)] = p_tau[:len(pL_tau)]  # shift the pressure signal by tau0
    p_L = p_L + pL_tau  # add the radiated sound from each face at the observer point x0

# save the acoustic pressure time signal p_L at observerpoint x0
f = open('data_pL.pckl', 'wb')
pickle.dump([p_L, x0, dt], f)
f.close()
print('Finished')

H_NC = ((np.abs(p_L)) / (2.e-5))**2
val_NC, freq_NC = psd(H_NC, fs, detrend='mean')
SPL = 10 * np.log10(val_NC + 1e-10)  # Add a small value to avoid zero

plt.figure(1)
plt.semilogx(freq_NC * 100000, 5 * np.log10(val_NC), label='NACA0012', color='blue')
plt.title('Sound Pressure Level (dB)')
plt.xlabel('Frequency')
plt.ylabel('SPL, dB', labelpad=1.5)
grid(color='0.5', linestyle=':', linewidth=0.2)
plt.legend()

plt.savefig('SPL_plot.png')  # Save the plot as an image file
plt.show()

# Plotting
# Set font to sans serif
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'Verdana', 'DejaVu Sans']

plt.figure(figsize=(6, 6))
plt.semilogx(freq_NC * 100000, 5 * np.log10(val_NC), label='SPL', color='blue', linewidth=2)

plt.title('Sound Pressure Level (SPL)', fontsize=16)
plt.xlabel('Frequency (Hz)', fontsize=14)
plt.ylabel('SPL (dB)', fontsize=14)
plt.grid(True, linestyle='--', alpha=0.7)
plt.legend(fontsize=12)

# Customize tick label size
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)

# Use a more readable format for y-axis labels
plt.gca().yaxis.set_major_formatter(ScalarFormatter())

# Adding a horizontal line at 0 dB for reference
plt.axhline(0, color='black', linestyle='--', linewidth=1)

# Set x-axis and y-axis limits
plt.xlim(10, 10000)  # Set your desired x-axis limits
plt.ylim(-20, 100)     # Set your desired y-axis limits

# Save the plot as an image file
plt.savefig('SPL_plot_stylish.png', bbox_inches='tight', dpi=300)

# Display the plot
plt.show()
