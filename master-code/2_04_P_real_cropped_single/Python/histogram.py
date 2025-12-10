import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Load the CSV files (assuming the structure is [x, y, energy])
file1 = "../../matRad-master/x_y_energy_values_at_80_percent.csv"
file2 = "../../matRad-master/x_y_energy_values_at_80_percent_gradient.csv"
file3 = "../../matRad-master/x_y_energy_values_at_80_percent_MATLAB.csv"

# Read the CSV files (only first two columns: x and y)
data1 = pd.read_csv(file1, usecols=[0, 1], header=None, names=['x', 'y'])
data2 = pd.read_csv(file2, usecols=[0, 1], header=None, names=['x', 'y'])
data3 = pd.read_csv(file3, usecols=[0, 1], header=None, names=['x', 'y'])

# Calculate Euclidean distance between corresponding (x, y) pairs for the first and second datasets
distances1 = np.sqrt((data1['x'] - data2['x'])**2 + (data1['y'] - data2['y'])**2)
distances2 = np.sqrt((data1['x'] - data3['x'])**2 + (data1['y'] - data3['y'])**2)

# Calculate the median and mean (expected value) of the distances for both datasets
median_distance1 = np.median(distances1)
mean_distance1 = np.mean(distances1)
median_distance2 = np.median(distances2)
mean_distance2 = np.mean(distances2)

# Set the figure size to (15, 9)
plt.figure(figsize=(15, 9))

# Determine the range of distances to cover both datasets
min_distance = min(distances1.min(), distances2.min())
max_distance = max(distances1.max(), distances2.max())

# Create common bin edges using np.linspace
num_bins = 23  # Specify the number of bins
common_bins = np.linspace(min_distance, max_distance, num_bins + 1)

# Plot histogram of the second dataset (MATLAB) in light gray, placed behind the main one
plt.hist(distances2, bins=common_bins, edgecolor='black', color='gray', alpha=0.6, label='Developed algorithm ($B_z$ = 0 T)')

# Plot histogram of the first dataset in royal blue, same opacity as the gray one
plt.hist(distances1, bins=common_bins, edgecolor='black', color='royalblue', alpha=0.6, label='Developed algorithm ($B_z$ = 1.5 T) optimized')

# Add titles and labels
plt.title('Distances between optimized ($B_z = 1.5 \, \, \mathrm{T}$) and matRad proton stopping positions', fontsize=24)
plt.xlabel('Euclidean distance [mm]', fontsize=24)
plt.ylabel('Frequency [N]', fontsize=24)

# Remove the grid from the plot
plt.grid(False)

# Add a legend to label the datasets and statistical lines
plt.legend(fontsize=20)

# Set font size for tick labels on both axes
plt.tick_params(axis='x', labelsize=24)  # Adjust tick size for x-axis
plt.tick_params(axis='y', labelsize=24)  # Adjust tick size for y-axis

# Add a table with statistics
table_data = [
    ["Dataset", "Mean [mm]", "Median [mm]"],
    ["$B_z$ = 0 T", f"{mean_distance2:.2f}", f"{median_distance2:.2f}"],
    ["Opt. $B_z$ = 1.5 T", f"{mean_distance1:.2f}", f"{median_distance1:.2f}"]
]

# Create the table
table = plt.table(
    cellText=table_data, 
    loc='upper right', 
    colWidths=[0.3, 0.3, 0.3],  # Increased column width
    cellLoc='center', 
    colLabels=None, 
    bbox=[0.45, 0.4, 0.5, 0.4]  # Adjusted table size (increased height)
)

# Adjust the font size for the table
for key, cell in table.get_celld().items():
    cell.set_fontsize(18)  # Set font size for each cell

# Adjust cell size (ensure they have more space by adjusting the bbox)
for key, cell in table.get_celld().items():
    cell.set_height(0.12)  # Increase height of the cells
    cell.set_width(0.3)    # Ensure width is enough for content

# Show the plot
plt.show()