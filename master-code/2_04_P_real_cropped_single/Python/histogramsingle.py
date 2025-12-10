import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Load the CSV files (assuming the structure is [x, y, energy])
file1 = "../../matRad-master/x_y_energy_values_at_80_percent.csv"
file2 = "../../matRad-master/x_y_energy_values_at_80_percent_MATLAB.csv"
#file2 = "../../matRad-master/x_y_energy_values_at_80_percent_MATLAB.csv"

# Read the CSV files (only first two columns: x and y)
data1 = pd.read_csv(file1, usecols=[0, 1], header=None, names=['x', 'y'])
data2 = pd.read_csv(file2, usecols=[0, 1], header=None, names=['x', 'y'])

# Calculate Euclidean distance between corresponding (x, y) pairs
distances = np.sqrt((data1['x'] - data2['x'])**2 + (data1['y'] - data2['y'])**2)

# Calculate the median and mean (expected value) of the distances
median_distance = np.median(distances)
mean_distance = np.mean(distances)

# Set the figure size to (15, 9)
plt.figure(figsize=(15, 9))

# Plot histogram of the distances in royal blue
plt.hist(distances, bins=20, edgecolor='black', color='royalblue')

# Mark the median and mean on the plot
plt.axvline(median_distance, color='red', linestyle='--', linewidth=2, label=f'Median = {median_distance:.2f} mm')
plt.axvline(mean_distance, color='lightgreen', linestyle='--', linewidth=2, label=f'Mean = {mean_distance:.2f} mm')

# Add titles and labels
plt.title('Distances between unoptimized ($B_z$ = 0 T) and matRad proton stopping positions', fontsize=24)
plt.xlabel('Euclidean distance [mm]', fontsize=24)
plt.ylabel('Frequency [N]', fontsize=24)

# Remove the grid from the plot
plt.grid(False)

# Add a legend
plt.legend(fontsize=20)

# Set font size for tick labels on both axes
plt.tick_params(axis='x', labelsize=24)  # Adjust tick size for x-axis
plt.tick_params(axis='y', labelsize=24)  # Adjust tick size for y-axis

# Show the plot
plt.show()
