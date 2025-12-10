import numpy as np
import pandas as pd
import pydicom
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import csv

# Function to read CSV values
def read_csv_values(filepath):
    with open(filepath, 'r') as infile:
        reader = csv.reader(infile)
        return [float(row[0]) for row in reader]

def read_isoCenter(filepath):
    with open(filepath, 'r') as infile:
        reader = csv.reader(infile)
        isoCenter = [float(row[0]) for row in reader]
    return isoCenter[0], isoCenter[1]  # Return x, y

def read_target(filepath):
    with open(filepath, 'r') as infile:
        reader = csv.reader(infile)
        target = [list(map(float, row)) for row in csv.reader(infile)]
    return target

def transform_target(target, isoCenter):
    transformed_target = [[x + isoCenter[0], y + isoCenter[1]] for x, _, y in target]
    return transformed_target

# Load the DICOM file
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

# Calculate extent for imshow (extent = [xmin, xmax, ymin, ymax])
resolution_mm = 1.09375  # Pixel spacing in mm
num_rows, num_cols = SPR.shape
x_extent = num_cols * resolution_mm
y_extent = num_rows * resolution_mm
extent = [0, x_extent, 0, y_extent]

# Plot the colormap with the correct resolution for the CT image
fig, ax = plt.subplots(figsize=(15, 9))
im = ax.imshow(SPR, cmap=cm.jet, extent=extent, aspect='auto')
plt.colorbar(im, ax=ax, label='SPR Value')
plt.title('SPR Value Colormap with Correct Resolution')
plt.xlabel('X (mm)')
plt.ylabel('Y (mm)')

# Plot coordsX and coordsY on the colormap without transforming them
coordsX_filepath = './coordsX_transformed.csv'
coordsY_filepath = './coordsY_transformed.csv'

coordsX = read_csv_values(coordsX_filepath)
coordsY = read_csv_values(coordsY_filepath)

# Plot coordsX and coordsY on the SPR colormap without transformation
ax.plot(coordsX, coordsY, 'ro-', label='coordsX vs coordsY')
ax.legend()

# Load isoCenter and target values
isoCenter_filepath = '../../matRad-master/isoCenter.csv'
target_filepath = '../../matRad-master/target.csv'

isoCenter = read_isoCenter(isoCenter_filepath)
target = read_target(target_filepath)

# Transform the target coordinates using isoCenter
transformed_target = transform_target(target, isoCenter)

# Extract transformed x and y coordinates for plotting
transformed_x = [point[0] for point in transformed_target]
transformed_y = [point[1] for point in transformed_target]

# Scatter plot the transformed target values onto the SPR colormap without scaling them
ax.scatter(transformed_x, transformed_y, color='blue', label='Transformed Target Points', s=10)

plt.legend()
plt.show()
plt.close()
