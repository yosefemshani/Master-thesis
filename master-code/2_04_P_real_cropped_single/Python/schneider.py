import numpy as np
import pandas as pd
import pydicom
import os
import matplotlib.pyplot as plt
import matplotlib.cm as cm

# Define the conversion functions
def DEMono79opt_conversion(HU):
    DEMono79opt_5k = (5000 - 3071) * (2.86 - 1.07) / (3071 - 130) + 2.86
    return np.interp(HU, [-1024, -990, -900, -150, -120, -70, -20, 70, 130, 3071, 5000],
                     [0.01, 0.01, 0.12, 0.83, 0.92, 0.95, 0.98, 1.06, 1.07, 2.86, DEMono79opt_5k])

def proCT_conversion(HU):
    proCT_5k = (5000 - 3071) * (2.32 - 1.15) / (3071 - 184) + 2.32
    return np.interp(HU, [-1024, -1000, -71, 184, 3071, 5000],
                     [0.01, 0.01, 0.96, 1.15, 2.32, proCT_5k])

# Define the range of DICOM files
dicom_files = [f"IMG{i:04d}_cropped.dcm" for i in range(10, 72)]

for dicom_file in dicom_files:
    # Read DICOM file
    ct = pydicom.dcmread(os.path.join("../DICOM", dicom_file))
    
    # Extract pixel data
    pixel_data = ct.pixel_array
    
    # Save HU values to a CSV file for this slice
    hu_values_csv_file = f"{dicom_file.split('.')[0]}.csv"
    np.savetxt(hu_values_csv_file, pixel_data, delimiter=",")
    print(f"HU values for slice {dicom_file} exported successfully.")
    
    # Load HU values from the CSV file
    hu_values_df = pd.read_csv(hu_values_csv_file, header=None)
    
    # Apply the chosen conversion method
    default_conversion = "proCT"  # Change this to "DEMono79opt" if needed
    if default_conversion == "DEMono79opt":
        SPR = DEMono79opt_conversion(hu_values_df)
    elif default_conversion == "proCT":
        SPR = proCT_conversion(hu_values_df)
    else:
        print("Invalid default conversion method specified.")
    
    # Save SPR values to a CSV file
    spr_values_csv_file = f"{dicom_file.split('.')[0]}_SPR.csv"
    np.savetxt(spr_values_csv_file, SPR, delimiter=",")
    print(f"SPR values for slice {dicom_file} saved to {spr_values_csv_file}.")
    
    # Pixel spacing
    original_pixel_spacing = ct.PixelSpacing
    print(f"The original pixel spacing for {dicom_file} is {original_pixel_spacing}.")
    
   # Plot the colormap only for IMG0035_cropped.dcm
    if dicom_file == "IMG0035_cropped.dcm":
        # Plot the colormap
        fig, ax = plt.subplots()
        im = ax.imshow(SPR, cmap=cm.jet, aspect='auto')
        plt.colorbar(im, ax=ax, label='SPR Value')
        plt.title(f'SPR Value Colormap for {dicom_file}')
        
        # Load the trajectory data
        trajectory_data = np.loadtxt('trajectory_mm.csv', delimiter=',')
        
        # Plot the trajectory
        plt.plot(trajectory_data[:, 0], trajectory_data[:, 1], 'w.-', linewidth=1)
        
        plt.show()
        plt.close()