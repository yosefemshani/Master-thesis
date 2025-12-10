import numpy as np
import pandas as pd
import pydicom
import os
import matplotlib.pyplot as plt
import matplotlib.cm as cm

# Create new DICOM file with rescaled pixel size

ct = pydicom.dcmread("../DICOM/IMG0036_cropped.dcm")

# Pixel spacing
original_pixel_spacing = ct.PixelSpacing
print("The original pixel spacing is ", original_pixel_spacing)

# Calculate the factor to change the resolution
factor = original_pixel_spacing[0] / 0.001

# Change the resolution of pixel spacing to 0.001mm
new_pixel_spacing = [x / factor for x in original_pixel_spacing]
print("The new pixel spacing is ", new_pixel_spacing)

# Update the pixel spacing in the DICOM metadata
ct.PixelSpacing = new_pixel_spacing

# Save the modified DICOM file
ct.save_as("IMG0036_cropped_updated.dcm")

ct_updated = pydicom.dcmread("IMG0036_cropped_updated.dcm")

# Work with new DICOM file!

# Extract pixel data
pixel_data = ct_updated.pixel_array
# Save HU values to a CSV file for this slice
np.savetxt("IMG0036_cropped_updated.csv", pixel_data, delimiter=",")
        
print(f"HU values for slice IMG0036_cropped_updated exported successfully.")



# Load HU values from the CSV file
hu_values_csv_file = "IMG0036_cropped_updated.csv"
hu_values_df = pd.read_csv(hu_values_csv_file, header=None)

# Flatten the DataFrame to get a 1D array of HU values
#hu_values = hu_values_df.values.flatten()

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
    SPR = DEMono79opt_conversion(hu_values_df) #change to hu_values for flattened! next line as well
elif default_conversion == "proCT":
    SPR = proCT_conversion(hu_values_df)
else:
    print("Invalid default conversion method specified.")

# Save SPR values to a CSV file
spr_values_csv_file = "SPR_values_updated.csv"
np.savetxt(spr_values_csv_file, SPR, delimiter=",")

print("SPR values saved to", spr_values_csv_file)


# Pixel spacing
pixel_spacing = ct_updated.PixelSpacing
print("The pixel spacing is ", pixel_spacing)



# Plot the colormap
fig, ax = plt.subplots()
im = ax.imshow(SPR, cmap=cm.jet, aspect='auto') #hu_values_df for HU, SPR for SPR
plt.colorbar(im, ax=ax, label='SPR Value')
plt.title('SPR Value Colormap')

# Load the trajectory data
trajectory_data = np.loadtxt('trajectory_mm.csv', delimiter=',')  # Load the trajectory data

# Plot the trajectory
plt.plot(trajectory_data[:, 0], trajectory_data[:, 1], 'w.-', linewidth=1)  # Plot the trajectory

plt.show()
plt.close()

# Load and average SPR matrices from all specified CSV files
spr_files = [f"IMG{str(i).zfill(4)}_cropped_SPR.csv" for i in range(10, 72)]
spr_matrices = [np.loadtxt(file, delimiter=',') for file in spr_files]

# Calculate the average matrix
average_spr_matrix = np.mean(spr_matrices, axis=0)

# Save the averaged matrix to a new CSV file
average_spr_csv_file = "DICOM_SPR.csv"
np.savetxt(average_spr_csv_file, average_spr_matrix, delimiter=",")

print("Averaged SPR values saved to", average_spr_csv_file)