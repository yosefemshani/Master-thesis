import numpy as np
import pandas as pd
import pydicom
import os
import matplotlib.pyplot as plt
import matplotlib.cm as cm

# Read the DICOM file
ct = pydicom.dcmread("../DICOM/IMG0036_cropped.dcm")

# Extract pixel data
pixel_data = ct.pixel_array

# Save HU values to a CSV file for this slice
np.savetxt("IMG0036_cropped.csv", pixel_data, delimiter=",")

# Load HU values from the CSV file
hu_values_csv_file = "IMG0036_cropped.csv"
hu_values_df = pd.read_csv(hu_values_csv_file, header=None)

# Define the conversion functions
def DEMono79opt_conversion(HU):
    DEMono79opt_5k = (5000 - 3071) * (2.86 - 1.07) / (3071 - 130) + 2.86
    return np.interp(HU, [-1024, -990, -900, -150, -120, -70, -20, 70, 130, 3071, 5000],
                     [0.01, 0.01, 0.12, 0.83, 0.92, 0.95, 0.98, 1.06, 1.07, 2.86, DEMono79opt_5k])

def proCT_conversion(HU):
    proCT_5k = (5000 - 3071) * (2.32 - 1.15) / (3071 - 184) + 2.32
    return np.interp(HU, [-1024, -1000, -71, 184, 3071, 5000],
                     [0.01, 0.01, 0.96, 1.15, 2.32, proCT_5k])

# Choose the default conversion method
default_conversion = "proCT"

# Apply the chosen conversion method
if default_conversion == "DEMono79opt":
    SPR = DEMono79opt_conversion(hu_values_df)
elif default_conversion == "proCT":
    SPR = proCT_conversion(hu_values_df)
else:
    print("Invalid default conversion method specified.")

# Save SPR values to a CSV file
spr_values_csv_file = "SPR_values.csv"
np.savetxt(spr_values_csv_file, SPR, delimiter=",")

print("SPR values saved to", spr_values_csv_file)

# Plot the colormap
fig, ax = plt.subplots(figsize=(15, 9))
im = ax.imshow(SPR, cmap=cm.jet, aspect='auto')
plt.colorbar(im, ax=ax, label='SPR Value')
plt.title('SPR Value Colormap')

# Load the trajectory data
trajectory_data = np.loadtxt('trajectory_TOPAS_filtered.csv', delimiter=',')

# Calculate the center y-coordinate of the CT image
ct_center_y = pixel_data.shape[0] // 2

# Calculate the offset needed to align the trajectory start point to the center of the CT image
trajectory_start_y = trajectory_data[0, 1]
y_offset = ct_center_y - trajectory_start_y

# Apply the offset to the y-coordinates to align the trajectory start point to the center of the CT image
trajectory_data[:, 1] += y_offset

# Invert the y-coordinates to correct the downward shift
trajectory_data[:, 1] = pixel_data.shape[0] - trajectory_data[:, 1]

# Plot the adjusted trajectory
plt.plot(trajectory_data[:, 0], trajectory_data[:, 1], 'w.-', linewidth=1)
plt.show()
plt.close()
