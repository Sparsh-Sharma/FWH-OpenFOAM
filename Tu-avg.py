import os
import glob
import pandas as pd
import matplotlib.pyplot as plt

# Specify the directory where your CSV fils are located
csv_directory = '/scratch/ws/m0/shar_sp-LES/LEN/80_2_1000Mean'

# Get a list of all CSV files in the directory
csv_files = glob.glob(os.path.join(csv_directory, '*.csv'))

# Create an empty DataFrame to store the aggregated data
aggregate_data = pd.DataFrame()

# Iterate through each CSV file and append the data to the DataFrame
for csv_file in csv_files:
    df = pd.read_csv(csv_file)
    aggregate_data = aggregate_data.append(df, ignore_index=True)

# Calculate the time-averaged turbulence intensity
time_averaged_data = aggregate_data.groupby('Time').mean().reset_index()

# Plot the time-averaged turbulence intensity
plt.plot(time_averaged_data['Time'], time_averaged_data['Tu'])
plt.xlabel('Time')
plt.ylabel('Tu (Turbulence Intensity)')
plt.title('Time-Averaged Turbulence Intensity along 2D Line')
plt.grid(True)
plt.savefig('/scratch/ws/m0/shar_sp-LES/LEN/80_2_1000Mean/Tufigure.png')
plt.show()

# Save the time-averaged data to a new CSV file
time_averaged_data.to_csv('/scratch/ws/m0/shar_sp-LES/LEN/80_2_1000Mean/averaged_data.csv', index=False)
