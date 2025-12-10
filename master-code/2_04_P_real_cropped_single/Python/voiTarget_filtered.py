import numpy as np
import pandas as pd
import pydicom
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import csv

# Function to read matrix from CSV
def read_matrix(filepath):
    return pd.read_csv(filepath, header=None).values

# Function to read isoCenter from CSV
def read_isoCenter(filepath):
    with open(filepath, 'r') as infile:
        reader = csv.reader(infile)
        isoCenter = [float(row[0]) for row in reader]
    return isoCenter[0], isoCenter[1]  # Return x, y (ignore z)

# Function to read rayPos from CSV
def read_rayPos(filepath):
    rayPos_df = pd.read_csv(filepath, header=None)
    x_positions = rayPos_df.iloc[0].to_numpy()
    y_positions = rayPos_df.iloc[2].to_numpy()
    return x_positions, y_positions

# Function to transform rayPos coordinates using isoCenter
def transform_rayPos(x_positions, y_positions, isoCenter):
    transformed_x = x_positions + isoCenter[0]
    transformed_y = y_positions + isoCenter[1]
    return transformed_x, transformed_y

# Function to read and filter maxRows3D
def read_and_filter_maxRows3D(filepath, threshold):
    max_rows_df = pd.read_csv(filepath, header=None)
    x_positions = max_rows_df.iloc[:, 0].to_numpy(dtype=float)
    y_positions = max_rows_df.iloc[:, 1].to_numpy(dtype=float)
    z_positions = max_rows_df.iloc[:, 2].to_numpy(dtype=float)

    std_dev = np.std(x_positions)
    threshold_value = std_dev * threshold
    mask = np.abs(x_positions - np.mean(x_positions)) <= threshold_value
    filtered_x = x_positions[mask]
    filtered_y = y_positions[mask]
    filtered_z = z_positions[mask]
    
    return filtered_x, filtered_y, filtered_z

# Load the DICOM file
ct = pydicom.dcmread("../DICOM/IMG0036_cropped.dcm")

# Extract pixel data
pixel_data = ct.pixel_array

# Save HU values to a CSV file
np.savetxt("IMG0036_cropped.csv", pixel_data, delimiter=",")

# Load HU values
hu_values_csv_file = "IMG0036_cropped.csv"
hu_values_df = pd.read_csv(hu_values_csv_file, header=None)

# Define conversion functions
def DEMono79opt_conversion(HU):
    DEMono79opt_5k = (5000 - 3071) * (2.86 - 1.07) / (3071 - 130) + 2.86
    return np.interp(HU, [-1024, -990, -900, -150, -120, -70, -20, 70, 130, 3071, 5000],
                     [0.01, 0.01, 0.12, 0.83, 0.92, 0.95, 0.98, 1.06, 1.07, 2.86, DEMono79opt_5k])

def proCT_conversion(HU):
    proCT_5k = (5000 - 3071) * (2.32 - 1.15) / (3071 - 184) + 2.32
    return np.interp(HU, [-1024, -1000, -71, 184, 3071, 5000],
                     [0.01, 0.01, 0.96, 1.15, 2.32, proCT_5k])

# Conversion method
default_conversion = "proCT"

# Apply conversion
if default_conversion == "DEMono79opt":
    SPR = DEMono79opt_conversion(hu_values_df)
elif default_conversion == "proCT":
    SPR = proCT_conversion(hu_values_df)

# Load the voiTarget matrix
voi_target_filepath = '../../matRad-master/voiTarget_page_2.csv'
voi_target_matrix = read_matrix(voi_target_filepath)
target_positions = np.argwhere(voi_target_matrix == 1)

# Calculate extent for imshow
resolution_mm = 1.09375
num_rows, num_cols = SPR.shape
x_extent = num_cols * resolution_mm
y_extent = num_rows * resolution_mm
extent = [0, x_extent, 0, y_extent]

# Plot the SPR colormap
fig, ax = plt.subplots(figsize=(15, 9))
im = ax.imshow(SPR, cmap=cm.jet, extent=extent, aspect='auto')
plt.colorbar(im, ax=ax, label='SPR Value')
plt.title('SPR Value Colormap with Correct Resolution')
plt.xlabel('X (mm)')
plt.ylabel('Y (mm)')

# Adjust and plot target volume positions
x_positions_target = target_positions[:, 1] * resolution_mm
y_positions_target = target_positions[:, 0] * resolution_mm
ax.scatter(x_positions_target, y_positions_target, color='red', label='Target Volume', s=200)

# Load isoCenter and rayPos
isoCenter_filepath = '../../matRad-master/isoCenter.csv'
rayPos_filepath = '../../matRad-master/rayPos.csv'

isoCenter = read_isoCenter(isoCenter_filepath)
x_positions, y_positions = read_rayPos(rayPos_filepath)

# Transform and plot rayPos
transformed_x, transformed_y = transform_rayPos(x_positions, y_positions, isoCenter)

# Filter and plot maxRows3D from original matRad file
maxRows3D_filepath = '../../matRad-master/x_y_energy_values_at_80_percent.csv'
x_positions_3d, y_positions_3d, energy_values_3d = read_and_filter_maxRows3D(maxRows3D_filepath, threshold=0.7)

x_positions_3d *= 1  # No resolution scaling needed
y_positions_3d *= 1  # No resolution scaling needed
ax.scatter(x_positions_3d, y_positions_3d, color='green', label='matRad (green)', s=20, marker='x')  # Crosses for matRad

# Load filtered MATLAB data (no filtering, direct loading)
maxRows_MATLAB_filepath = '../../matRad-master/x_y_energy_values_at_80_percent_MATLAB.csv' 
# '../../matRad-master/maxRows_MATLAB_with_energy_filtered_0.csv' '../../matRad-master/maxRows_MATLAB_with_energy_filtered_1_5.csv' '../../matRad-master/maxRows_MATLAB_with_energy_filtered_3.csv'
matlab_data_df = pd.read_csv(maxRows_MATLAB_filepath, header=None)
x_positions_MATLAB = matlab_data_df.iloc[:, 0].to_numpy(dtype=float)
y_positions_MATLAB = matlab_data_df.iloc[:, 1].to_numpy(dtype=float)

# Scatter plot for MATLAB data using dots
ax.scatter(x_positions_MATLAB, y_positions_MATLAB, color='blue', label='MATLAB (blue)', s=20, marker='.')  # Dots for MATLAB

# Add legend and display plot
plt.legend()
plt.show()
plt.close()
