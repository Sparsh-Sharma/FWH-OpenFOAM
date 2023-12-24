import numpy as np
import pickle
import os
from matplotlib.mlab import psd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# ----------------
fs = 100000  # sampling frequency
c = 343  # speed of sound
rho0 = 1.2  # density of air
dt = fs**-1  # timesep
# ---------------------

f = open('data_geo.pckl', 'rb')
obj = pickle.load(f)
f.close()
p_cxyz = obj[3]

A = obj[0]  # area of each cell on the surface (face)
n = obj[1]  # normal vector of face
c_xyz = obj[2]  # centroid of A

# Additional Parameters
num_angles = 21  # Number of angles to consider (0 to 360 degrees)
num_phi = 11  # Number of phi angles to consider (0 to 180 degrees)
theta_values = np.linspace(0, 360, num_angles)  # Angles in degrees
phi_values = np.linspace(0, 180, num_phi)  # Phi angles in degrees

# Specify the frequencies of interest
freq_of_interest = [0.08]

# Specify the observer positions (x0 values)
x0_values = [0.4]

# Initialize an array to store SPL values for each angle, phi, and frequency
SPL_values_3D = np.zeros((num_angles, num_phi, len(freq_of_interest)))

# Loop through different observer positions (x0 values)
for x0 in x0_values:
    # Loop through different angles
    for idx_theta, theta in enumerate(theta_values):
        # Loop through different phi angles
        for idx_phi, phi in enumerate(phi_values):
            # Calculate observer position for the current angle and phi
            x0_theta_phi = np.array([x0 * np.sin(np.radians(phi)) * np.cos(np.radians(theta)),
                                     x0 * np.sin(np.radians(phi)) * np.sin(np.radians(theta)),
                                     x0 * np.cos(np.radians(phi))])

            # Calculate the force for every face at the new observer position x0_theta_phi
            F_theta_phi = (-rho0 * A * p_cxyz.T).T

            # Calculate distance r from centroid to the new observer position x0_theta_phi
            r_theta_phi = np.sqrt(((c_xyz - x0_theta_phi)[:, 0])**2 + ((c_xyz - x0_theta_phi)[:, 1])**2 +
                                 ((c_xyz - x0_theta_phi)[:, 2])**2)

            # Calculate unity vector in the direction r_theta_phi
            e_r_theta_phi = (x0_theta_phi - c_xyz) / r_theta_phi[:, np.newaxis]

            # Calculate retarded time tau_theta_phi
            tau_theta_phi = (r_theta_phi - np.min(r_theta_phi)) / c

            # Calculate pressure at each face for the new observer position x0_theta_phi
            p_L_theta_phi = np.zeros(F_theta_phi.shape[1] - 6, dtype=float)  # -6 in shape because sixth order

            for i, Frow in enumerate(F_theta_phi):
                p_tau_theta_phi = (((45 * Frow[4:-2] - 45 * Frow[2:-4] - 9 * Frow[5:-1] + 9 * Frow[1:-5] -
                                     Frow[0:-6] + Frow[6:]).T) *
                                   np.dot(e_r_theta_phi[i, :], n[i, :])) * \
                                  (4 * np.pi * c * r_theta_phi[i] * 60 * dt)**-1  # 6order

                pL_tau_theta_phi = np.zeros(len(p_tau_theta_phi))
                pL_tau_theta_phi[:len(p_tau_theta_phi)] = p_tau_theta_phi[:len(pL_tau_theta_phi)]
                p_L_theta_phi = p_L_theta_phi + pL_tau_theta_phi

            # Calculate Power Spectral Density (PSD) for the new observer position x0_theta_phi
            H_NC_theta_phi = ((np.abs(p_L_theta_phi)) / (2.e-5))**2
            val_NC_theta_phi, freq_NC_theta_phi = psd(H_NC_theta_phi, fs, detrend='mean')

            # Find indices corresponding to frequencies of interest
            indices_of_interest = [np.argmin(np.abs(freq_NC_theta_phi - f)) for f in freq_of_interest]

            # Extract SPL values at frequencies of interest
            SPL_values_3D[idx_theta, idx_phi, :] = 5 * np.log10(val_NC_theta_phi[indices_of_interest] + 1e-10)

            # Print progress
            print(f'x0: {x0}, Iteration {idx_theta + 1}/{num_angles}, {idx_phi + 1}/{num_phi} - '
                  f'Angle: {theta}, Phi: {phi}, SPL at {freq_of_interest} Hz: {SPL_values_3D[idx_theta, idx_phi, :]} dB')
#%%
import numpy as np
import pickle
import os
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

def interp_array(N1):  # add interpolated rows and columns to array
    N2 = np.empty([int(N1.shape[0]), int(2*N1.shape[1] - 1)])  # insert interpolated columns
    N2[:, 0] = N1[:, 0]  # original column
    for k in range(N1.shape[1] - 1):  # loop through columns
        N2[:, 2*k+1] = np.mean(N1[:, [k, k + 1]], axis=1)  # interpolated column
        N2[:, 2*k+2] = N1[:, k+1]  # original column
    N3 = np.empty([int(2*N2.shape[0]-1), int(N2.shape[1])])  # insert interpolated columns
    N3[0] = N2[0]  # original row
    for k in range(N2.shape[0] - 1):  # loop through rows
        N3[2*k+1] = np.mean(N2[[k, k + 1]], axis=0)  # interpolated row
        N3[2*k+2] = N2[k+1]  # original row
    return N3

# Set up spherical coordinates
theta_mesh, phi_mesh = np.meshgrid(np.linspace(0, 360, num_phi), np.linspace(0, 180, num_angles))
# R = SPL_values_3D.mean(axis=1)  # Use mean SPL for radius
R = np.mean(SPL_values_3D, axis=1)  # Use mean SPL for radius


# Convert spherical coordinates to Cartesian coordinates
X = R * np.sin(np.radians(phi_mesh)) * np.cos(np.radians(theta_mesh))
Y = R * np.sin(np.radians(phi_mesh)) * np.sin(np.radians(theta_mesh))
Z = R * np.cos(np.radians(phi_mesh))

# Interpolate between points to increase the number of faces
for _ in range(3):  # Interpolate three times (adjust as needed)
    X = interp_array(X)
    Y = interp_array(Y)
    Z = interp_array(Z)

fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(1, 1, 1, projection='3d')
# ax = fig.add_subplot(111, projection='3d')
ax.grid(True)
ax.axis('off')
ax.set_xticks([])
ax.set_yticks([])
ax.set_zticks([])

N = np.sqrt(X**2 + Y**2 + Z**2)
Rmax = np.max(N)
N = N / Rmax

# axes_length = 1.5
# ax.plot([0, axes_length * Rmax], [0, 0], [0, 0], linewidth=2, color='red', label='X-axis')
# ax.plot([0, 0], [0, axes_length * Rmax], [0, 0], linewidth=2, color='green', label='Y-axis')
# ax.plot([0, 0], [0, 0], [0, axes_length * Rmax], linewidth=2, color='blue', label='Z-axis')

# Find middle points between values for face colors
N = interp_array(N)[1::2, 1::2]

mycol = cm.viridis(N)

surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, facecolors=mycol, linewidth=0.5, antialiased=True, shade=False)
# surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, facecolors=mycol, cmap='viridis', edgecolor='k')
# ax.set_xlim([-axes_length * Rmax, axes_length * Rmax])
# ax.set_ylim([-axes_length * Rmax, axes_length * Rmax])
# ax.set_zlim([-axes_length * Rmax, axes_length * Rmax])

m = cm.ScalarMappable(cmap=cm.viridis)
m.set_array(R)
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')

# Add color bar
cbar = fig.colorbar(m, ax=ax, shrink=0.8)
cbar.set_label('SPL (dB)', rotation=270, labelpad=15)

# Set legend
ax.legend(loc='upper right', fontsize='medium')

# Set title
ax.set_title('SPL Directivity in 3D')

# Set view angle
ax.view_init(azim=300, elev=30)

plt.show()
