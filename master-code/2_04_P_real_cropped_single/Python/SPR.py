import numpy as np
import pandas as pd
import pydicom
import os
import matplotlib.pyplot as plt
import matplotlib.cm as cm

def dicom_to_spr(input_dicom_file, default_conversion="matRad"):
    """
    Processes a DICOM file and generates HU and SPR values, saving them to CSV files.
    
    Args:
        input_dicom_file (str): The path to the DICOM file (e.g., "../DICOM_DoubleSlice/IMG0036_cropped.dcm").
        default_conversion (str): The conversion method to use (options: "DEMono79opt", "proCT", "matRad"). Default is "matRad".
    
    Output:
        Two CSV files:
            - IMG(...)_HU.csv: The HU values for the input DICOM slice.
            - SPR_values.csv: The SPR values based on the chosen conversion method.
    """
    # Load the DICOM file
    ct = pydicom.dcmread(input_dicom_file)

    # Extract pixel spacing and pixel data
    original_pixel_spacing = ct.PixelSpacing
    print("The original pixel spacing is ", original_pixel_spacing)
    
    pixel_data = ct.pixel_array

    # Generate the HU values CSV file name based on the input DICOM file name
    dicom_filename = os.path.basename(input_dicom_file)
    hu_csv_filename = dicom_filename.replace(".dcm", "_HU.csv")

    # Save HU values to the CSV file
    np.savetxt(hu_csv_filename, pixel_data, delimiter=",")
    print(f"HU values for slice {dicom_filename} exported successfully to {hu_csv_filename}.")

    # Load HU values from the CSV file to process SPR
    hu_values_df = pd.read_csv(hu_csv_filename, header=None)

    # Define the conversion functions
    def DEMono79opt_conversion(HU):
        DEMono79opt_5k = (5000 - 3071) * (2.86 - 1.07) / (3071 - 130) + 2.86
        return np.interp(HU, [-1024, -990, -900, -150, -120, -70, -20, 70, 130, 3071, 5000],
                         [0.01, 0.01, 0.12, 0.83, 0.92, 0.95, 0.98, 1.06, 1.07, 2.86, DEMono79opt_5k])

    def proCT_conversion(HU):
        proCT_5k = (5000 - 3071) * (2.32 - 1.15) / (3071 - 184) + 2.32
        return np.interp(HU, [-1024, -1000, -71, 184, 3071, 5000],
                         [0.01, 0.01, 0.96, 1.15, 2.32, proCT_5k])

    def matRad_conversion(HU):
        return np.interp(HU, [-1024, 200, 449, 2000, 2048, 3071],
                        [0.00324, 1.2, 1.20001, 2.49066, 2.5306, 2.53061])

    # Apply the chosen conversion method
    if default_conversion == "DEMono79opt":
        SPR = DEMono79opt_conversion(hu_values_df)
    elif default_conversion == "proCT":
        SPR = proCT_conversion(hu_values_df)
    elif default_conversion == "matRad":
        SPR = matRad_conversion(hu_values_df)
    else:
        raise ValueError("Invalid default conversion method specified.")

    # Save SPR values to a CSV file
    spr_values_csv_file = "SPR_values.csv"
    np.savetxt(spr_values_csv_file, SPR, delimiter=",")
    print(f"SPR values saved to {spr_values_csv_file}.")

# Example usage of the function
dicom_to_spr("../DICOM_bone/IMG0036_cropped.dcm", default_conversion="proCT")
