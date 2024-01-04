import os
import glob
import re
import pandas as pd
import matplotlib.pyplot as plt

# Get the directory where the script is located
script_directory = os.path.dirname(os.path.realpath(__file__))

# Get the parent directory name
parent_directory_name = os.path.basename(os.path.dirname(script_directory))

# Extract information from the parent directory name using regular expressions
match = re.match(r'(\d+)_(\d+)_(\d+)', parent_directory_name)
if match:
    # Extract the values from the regular expression match
    diameter, thickness, nozzle = match.groups()

    # Create the label for the plot dynamically
    label = f'{diameter}mm|{thickness}mm|{"Nozzle" if nozzle == "1000" else "HW"}'
else:
    print(f"Warning: Unable to extract information from the parent directory name: {parent_directory_name}. Using default label.")
    label = 'Default Label'

# Get a list of all CSV files in the script's directory
csv_files = glob.glob(os.path.join(script_directory, '*.csv'))

# Check if there are matching files
if not csv_files:
    print("No CSV files found in the script's directory")
    exit()

# Initialize variables to store time-averaged data
total_tu = None
total_time = 0

# Loop through each CSV file
for csv_file in csv_files:
    try:
        # Read CSV file into a Pandas DataFrame
        df = pd.read_csv(csv_file)

        # Check if 'Tu' column exists
        if 'Tu' not in df.columns:
            print(f"Warning: 'Tu' column not found in {csv_file}. Skipping...")
            continue

        # Extract time and Tu columns
        time = df['Time'].values
        tu = df['Tu'].values

        # Accumulate time-averaged Tu
        total_tu = total_tu + tu if total_tu is not None else tu
        total_time += 1
    except Exception as e:
        print(f"Error reading {csv_file}: {e}")

# Check if any valid files were processed
if total_time == 0:
    print("No valid CSV files found.")
    exit()

# Calculate time-averaged Tu
average_tu = total_tu / total_time

# Create a DataFrame for the averaged data
averaged_data = pd.DataFrame({'X-coordinate': df['Points:0'].values,
                              'Average_Tu': average_tu})

# Construct filenames with the parent directory name as a prefix
pdf_filename = f"{parent_directory_name}_averaged_turbulence_intensity.pdf"
csv_filename = f"{parent_directory_name}_averaged_turbulence_intensity.csv"

# Plot and save the figure
fig, ax = plt.subplots()
line, = ax.plot(averaged_data['X-coordinate'] / 0.2, averaged_data['Average_Tu'] * 100, color='blue', linewidth=2,
                label=label)  # Use the dynamically created label
ax.set_xlabel('x/c')
ax.set_ylabel('Tu (%)')
# ax.set_title('Time-Averaged Turbulence Intensity')
ax.set_xlim([0, 10])  # Replace with actual x-axis limits
ax.set_ylim([0, 10])  # Replace with actual y-axis limits
ax.legend(loc='upper right', fontsize=14)
ax.grid(True)

# Set sans-serif font
font = {'family': 'sans-serif', 'size': 18}
plt.rc('font', **font)

# Adjust layout to prevent label cutoff
fig.tight_layout()

# Save the figure with the parent directory name as a prefix
fig.savefig(os.path.join(script_directory, pdf_filename))

# Save the averaged data to a new CSV file with the parent directory name as a prefix
averaged_data.to_csv(os.path.join(script_directory, csv_filename), index=False)

# Show the plot (optional)
plt.show()
