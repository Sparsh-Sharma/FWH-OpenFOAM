import numpy as np
import pickle
from matplotlib.mlab import psd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# ----------------
fs = 100000  # sampling frequency
c = 343  # speed of sound
rho0 = 1.2  # density of air
dt = fs**-1  # timesep
# ---------------------

# Load data_geo.pckl
f = open('data_geo.pckl', 'rb')
obj = pickle.load(f)
f.close()
p_cxyz = obj[3]

A = obj[0]  # area of each cell on the surface (face)
n = obj[1]  # normal vector of face
c_xyz = obj[2]  # centroid of A

# Specify the specific frequency you're interested in (in Hz)
specific_frequency = 0.01  # Replace this with your desired frequency

# Define the observer positions on the surface of the sphere
theta = np.linspace(0, 2 * np.pi, 10)
phi = np.linspace(0, np.pi, 10)
observer_positions = []

for t in theta:
    for p in phi:
        x = 0.4 * np.sin(p) * np.cos(t)
        y = 0.4 * np.sin(p) * np.sin(t)
        z = 0.4 * np.cos(p)
        observer_positions.append(np.array([x, y, z]))

# Initialize an array to store SPL values
SPL_field = np.zeros(len(observer_positions))

# Calculate SPL at each observer position
print("Calculating SPL at observer positions...")
total_positions = len(observer_positions)

for idx, x0 in enumerate(observer_positions):
    # Calculate the force for every face
    F = (-rho0 * A * p_cxyz.T).T

    # Calculate distance r from centroid to receiver location x0
    r = np.sqrt(((c_xyz - x0)[:, 0])**2 + ((c_xyz - x0)[:, 1])**2 + ((c_xyz - x0)[:, 2])**2)

    # Unity vector in the direction r
    e_r = np.empty(c_xyz.shape, dtype=float)
    for comp in range(len(x0)):
        e_r[:, comp] = np.divide((x0 - c_xyz)[:, comp], r)

    # Retarded time tau0, find the difference to the minimum distance of all points and divide by propagation speed
    tau0 = np.divide(r - np.min(r), c)

    # For every time step, add all points in the observer point x0
    p_L = np.ones(F.shape[1] - 6, dtype=float) - 1  # -6 in shape because sixth order

    total_faces = len(F)
    for i, Frow in enumerate(F):
        p_tau = (((45 * Frow[4:-2] - 45 * Frow[2:-4] - 9 * Frow[5:-1] + 9 * Frow[1:-5] - Frow[0:-6] + Frow[6:]).T) *
                 np.dot(e_r[i, :], n[i, :])) * (4 * np.pi * c * r[i] * 60 * dt)**-1  # 6order

        pL_tau = np.zeros(len(p_tau))
        pL_tau[:len(p_tau)] = p_tau[:len(pL_tau)]  # Shift the pressure signal by tau0
        p_L = p_L + pL_tau  # Add the radiated sound from each face at the observer point x0

        # Print progress for each position
        progress = (i + 1) / total_faces * 100
        remaining_iterations = total_faces - (i + 1)
        print(f"\rPosition Progress: [{int(progress) * '=' + (100 - int(progress)) * ' '}] {int(progress)}% | Remaining Iterations: {remaining_iterations}", end='', flush=True)
        

    # Calculate SPL
    H_NC = ((np.abs(p_L)) / (2.e-5))**2
    val_NC, freq_NC = psd(H_NC, fs, detrend='mean')
    
    # Find the index corresponding to the frequency closest to the specified value
    freq_index = np.argmin(np.abs(freq_NC - specific_frequency))
    SPL = 10 * np.log10(val_NC[freq_index] + 1e-10)  # Add a small value to avoid zero

    # Save SPL at this position
    SPL_field[idx] = SPL

    # Print overall progress
    overall_progress = (idx + 1) / total_positions * 100
    remaining_iterations = total_positions - (idx + 1)
    # print(f"\rOverall Progress: [{'=' * int(overall_progress) + ' ' * (100 - int(overall_progress))}] {int(overall_progress)}% | Remaining Iterations: {remaining_iterations}", end='', flush=True)
    print(f"\rOverall Progress: [{int(overall_progress) * '=' + (100 - int(overall_progress)) * ' '}] {int(overall_progress)}% | Remaining Iterations: {remaining_iterations}", end='', flush=True)

print("\nCalculations complete. Plotting SPL field...")
#%%
# Plotting the SPL field as a surface without grids
fig = plt.figure(figsize=(12, 10))
ax = fig.add_subplot(111, projection='3d')

# Extract x, y, z coordinates from observer_positions
x = [pos[0] for pos in observer_positions]
y = [pos[1] for pos in observer_positions]
z = [pos[2] for pos in observer_positions]

# Create a 3D scatter plot of SPL_field
sc = ax.scatter(x, y, z, c=SPL_field, cmap='viridis', s=80, edgecolor='k', linewidth=0.5)

# Create a 3D surface plot from the scatter plot
surf = ax.plot_trisurf(x, y, z, cmap='viridis', linewidth=0, antialiased=False, alpha=0.8)

# Set labels and title
ax.set_xlabel('X', fontsize=14)
ax.set_ylabel('Y', fontsize=14)
ax.set_zlabel('Z', fontsize=14)
ax.set_title(f'SPL Field at {specific_frequency} Hz', fontsize=16)

# Add a color bar
norm = plt.Normalize(SPL_field.min(), SPL_field.max())
sm = plt.cm.ScalarMappable(cmap='viridis', norm=norm)
sm.set_array([])  # This line is necessary for ScalarMappable to work properly

# Display the color bar
cbar = plt.colorbar(sm, ax=ax, pad=0.05)
cbar.set_label('SPL (dB)', fontsize=14)

# Remove grid lines
ax.grid(False)

# Customize tick label size
ax.tick_params(axis='x', labelsize=12)
ax.tick_params(axis='y', labelsize=12)
ax.tick_params(axis='z', labelsize=12)

# Set background color to pure white
fig.patch.set_facecolor('white')

# Save the plot as an image file
plt.savefig(f'SPL_Field_{specific_frequency}Hz_surface_professional_white_background.png', bbox_inches='tight', dpi=300, facecolor='white')

# Display the plot
plt.show()
