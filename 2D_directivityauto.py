import numpy as np
import pickle
import os
from matplotlib.mlab import psd
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter

# Create 'postAcoustics' folder if it doesn't exist
output_directory = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'postAcoustics')
os.makedirs(output_directory, exist_ok=True)

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
dt = deltaT  # timesep
# ---------------------

# Get the directory name where the script is stored
script_directory = os.path.basename(os.path.dirname(os.path.realpath(__file__)))

f = open('data_geo.pckl', 'rb')
obj = pickle.load(f)
f.close()
p_cxyz = obj[3]

A = obj[0]  # area of each cell on the surface (face)
n = obj[1]  # normal vector of face
c_xyz = obj[2]  # centroid of A

# Additional Parameters
num_angles = 5  # Number of angles to consider (0 to 360 degrees)
theta_values = np.linspace(0, 360, num_angles)  # Angles in degrees

# Specify the frequencies of interest
freq_of_interest = [0.001, 0.002, 0.003, 0.004, 0.005, 0.006, 0.008, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.08, 0.1]  # Add the frequencies you are interested in

# Specify the observer positions (x0 values)
x0_values = [0.2, 0.4]  # Add the x0 values you want to consider

# Loop through different observer positions (x0 values)
for x0 in x0_values:
    # Initialize an array to store SPL values for each angle and frequency
    SPL_values = np.zeros((num_angles, len(freq_of_interest)))

    # Loop through different angles
    for idx, theta in enumerate(theta_values):
        # Calculate observer position for the current angle
        x0_theta = np.array([x0 * np.cos(np.radians(theta)), x0 * np.sin(np.radians(theta)), 0.07])

        # Calculate the force for every face at the new observer position x0_theta
        F_theta = (-rho0 * A * p_cxyz.T).T

        # Calculate distance r_theta from centroid to the new observer position x0_theta
        r_theta = np.sqrt(((c_xyz - x0_theta)[:, 0])**2 + ((c_xyz - x0_theta)[:, 1])**2 + ((c_xyz - x0_theta)[:, 2])**2)

        # Calculate unity vector in the direction r_theta
        e_r_theta = np.empty(c_xyz.shape, dtype=float)
        for comp in range(len(x0_theta)):
            e_r_theta[:, comp] = np.divide((x0_theta - c_xyz)[:, comp], r_theta)

        # Calculate retarded time tau_theta
        tau_theta = np.divide(r_theta - np.min(r_theta), c)

        # Calculate pressure at each face for the new observer position x0_theta
        p_L_theta = np.ones(F_theta.shape[1] - 6, dtype=float) - 1  # -6 in shape because sixth order

        for i, Frow in enumerate(F_theta):
            p_tau_theta = (((45 * Frow[4:-2] - 45 * Frow[2:-4] - 9 * Frow[5:-1] + 9 * Frow[1:-5] - Frow[0:-6] + Frow[6:]).T) *
                           np.dot(e_r_theta[i, :], n[i, :])) * (4 * np.pi * c * r_theta[i] * 60 * dt)**-1  # 6order

            pL_tau_theta = np.zeros(len(p_tau_theta))
            pL_tau_theta[:len(p_tau_theta)] = p_tau_theta[:len(pL_tau_theta)]  # shift the pressure signal by tau_theta
            p_L_theta = p_L_theta + pL_tau_theta

        # Calculate Power Spectral Density (PSD) for the new observer position x0_theta
        H_NC_theta = ((np.abs(p_L_theta)) / (2.e-5))**2
        val_NC_theta, freq_NC_theta = psd(H_NC_theta, fs, detrend='mean')

        # Find indices corresponding to frequencies of interest
        indices_of_interest = [np.argmin(np.abs(freq_NC_theta - f)) for f in freq_of_interest]

        # Extract SPL values at frequencies of interest
        SPL_values[idx, :] = 5 * np.log10(val_NC_theta[indices_of_interest] + 1e-10)  # Add a small value to avoid zero

        # Print progress
        print(f'x0: {x0}, Iteration {idx + 1}/{num_angles} - Angle: {theta}, SPL at {freq_of_interest} Hz: {SPL_values[idx, :]} dB')

    # Save data to CSV file
    csv_data = np.column_stack((np.radians(theta_values), SPL_values))
    csv_filename = f'Directivity_Data_x0_{x0}_{script_directory}.csv'
    csv_path = os.path.join(output_directory, csv_filename)
    np.savetxt(csv_path, csv_data, delimiter=',', header='Theta (radians), ' + ', '.join(f'SPL_{freq}Hz' for freq in freq_of_interest), comments='')
    print(f'Directivity data saved to {csv_path}')

    # Set font to sans serif
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'Verdana', 'DejaVu Sans']

    # Loop through different frequencies
    for i, freq in enumerate(freq_of_interest):
        plt.figure(figsize=(8, 8), num=f"x0_{x0}_freq_{freq}Hz")
        ax = plt.subplot(111, projection='polar')

        ax.plot(np.radians(theta_values), SPL_values[:, i], color='blue', linewidth=2)

        ax.set_rlabel_position(0)
        ax.set_theta_zero_location('E')
        ax.set_rmax(np.max(SPL_values) + 5)  # Set maximum radius for better visibility

        ax.set_title(f'Sound Pressure Level Directivity at {freq} Hz (x0={x0})', va='bottom', fontsize=14)
        ax.grid(True)

        ax.legend(fontsize=12)

        # Customize tick label size
        ax.tick_params(axis='both', labelsize=12)

        # Save the polar plot as an image file
        figure_filename = f'Directivity_Polar_plot_x0_{x0}_freq_{freq}Hz_{script_directory}.pdf'
        figure_path = os.path.join(output_directory, 'figures', figure_filename)
        plt.savefig(figure_path, bbox_inches='tight', dpi=300)
        print(f'Directivity polar plot saved to {figure_path}')

# Display the polar plot
plt.show()
