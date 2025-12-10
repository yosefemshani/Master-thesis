import pandas as pd
import numpy as np
from matplotlib import pyplot as plt



# Read in the csv-file
outputfile_topas = '../TOPAS/TOPAS_Beam.csv'



# Skip header lines and load the data
df = pd.read_csv(outputfile_topas, comment='#', header=None)


# Convert the dataframe to a numpy array
topas_datamatrix = np.array(df)


# Extract depth (x) and dose values
depth = topas_datamatrix[:, 0]  # depth (x)
y_positions = topas_datamatrix[:, 1]  # y positions (height)
dose = topas_datamatrix[:, 3]  # dose values



# Reshape the dose data into a 2D grid based on 341 x bins, 225 y bins, and 7 z bins
n_x_bins = 341
n_y_bins = 225
n_z_bins = 62


dose_grid = dose.reshape(n_x_bins, n_y_bins, n_z_bins)


# Define x axis for the plot (convert from cm to mm)
x_bin_size = 0.109375 * 10  # mm per bin in x
x_axis = np.arange(n_x_bins) * x_bin_size  # depth in mm



# Plotting only the last z-bin (z_bin 7, index 6)
z_bin = 20  # The last z-bin (index 6)



# Integrate the dose along the y-axis for the last z-bin (collapse into 1D dose profile along x)
integrated_dose = np.sum(dose_grid[:, :, z_bin], axis=1)



# Normalize the integrated dose by its maximum value
integrated_dose_norm = (integrated_dose / np.max(integrated_dose)) * 100


# Find the peak dose index (where max dose occurs)
peak_idx = np.argmax(integrated_dose_norm)



# Search for the position after the peak where dose crosses below 80% of the peak dose
r80_value = 80  # 80% of the maximum dose
post_peak_dose = integrated_dose_norm[peak_idx:]
crossing_idx = np.where(post_peak_dose <= r80_value)[0][0] + peak_idx


# Get the values at the two points surrounding the R80 crossing
x1, x2 = x_axis[crossing_idx - 1], x_axis[crossing_idx]
y1, y2 = integrated_dose_norm[crossing_idx - 1], integrated_dose_norm[crossing_idx]



# Perform linear interpolation to find the exact x position for R80
r80_x = x1 + (r80_value - y1) * (x2 - x1) / (y2 - y1)


# Analyze where the dose remains around 0 after the peak
threshold = 0.0005  # Define a threshold close to zero
consecutive_length = 5  # Number of consecutive x values to check
zero_indices = np.where(integrated_dose_norm[peak_idx:] < threshold)[0]  # Find indices where dose < threshold


# Find the index where we have 50 consecutive values
cutoff_index = None
for i in range(len(zero_indices) - consecutive_length + 1):
    if zero_indices[i + consecutive_length - 1] - zero_indices[i] == consecutive_length - 1:
        cutoff_index = zero_indices[i] + peak_idx
        break


# If cutoff_index is found, limit x_axis and integrated_dose_norm
if cutoff_index is not None:
    x_axis = x_axis[:cutoff_index + 1]
    integrated_dose_norm = integrated_dose_norm[:cutoff_index + 1]


# Create a figure with specified size
plt.figure(figsize=(12, 8))


# Plot the normalized Bragg peak for the last z-bin
plt.plot(x_axis, integrated_dose_norm)


# Highlight the interpolated R80 point with a red dot and add it to the legend
plt.plot(r80_x, r80_value, 'ro', label=f'R$_{{80}}$ = {r80_x:.2f} mm')


# Add title, labels, and legend
plt.title('Depth dose curve for protons (N = $10^5$, E = 200 MeV) in water')
plt.xlabel('Depth [mm]')
plt.ylabel('Relative Dose [%]')
plt.grid(True)
plt.legend()

# Show the plot
plt.show()

