import numpy as np
import pandas as pd
import pydicom
import os
import re
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import csv
import numpy as np
import pandas as pd
import pydicom
import os

def dicom_to_spr(input_dicom_file, default_conversion="matRad"):
    """
    Processes a DICOM file and generates HU and SPR values, saving them to CSV files.
    
    Args:
        input_dicom_file (str): The path to the DICOM file (e.g., "../DICOM_DoubleSlice/IMG0036_cropped.dcm").
        default_conversion (str): The conversion method to use (options: "DEMono79opt", "proCT", "matRad"). Default is "matRad".
    
    Output:
        Two CSV files:
            - IMG(...)_HU.csv: The HU values for the input DICOM slice.
            - SPR_values_XX_conversion.csv: The SPR values based on the chosen conversion method, where XX corresponds to the slice number and conversion refers to the method used.
    """
    # Load the DICOM file
    ct = pydicom.dcmread(input_dicom_file)

    # Extract pixel spacing and pixel data
    original_pixel_spacing = ct.PixelSpacing
    print("The original pixel spacing is ", original_pixel_spacing)
    
    pixel_data = ct.pixel_array

    # Generate the HU values CSV file name based on the input DICOM file name
    dicom_filename = os.path.basename(input_dicom_file)
    slice_number = re.search(r'IMG00(\d+)_cropped.dcm', dicom_filename).group(1)
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

    # Save SPR values to a CSV file with the dynamic filename that includes the slice number and conversion method
    spr_values_csv_file = f"SPR_values_{slice_number}_{default_conversion}.csv"
    np.savetxt(spr_values_csv_file, SPR, delimiter=",")
    print(f"SPR values saved to {spr_values_csv_file}.")


def fake_pixel_data(slice_index, directory="../DICOM/", output_directory="../DICOM_autofake/"):
    """
    Kopiert die PixelData von der ausgewählten DICOM-Datei auf ihre benachbarten Dateien.
    
    Args:
        slice_index (int): Index der Datei, die gefälscht wird, z.B. 46 für IMG0046_cropped.dcm.
        directory (str): Verzeichnis der Original-DICOM-Dateien.
        output_directory (str): Verzeichnis, in das die modifizierten DICOM-Dateien gespeichert werden.
        
    Returns:
        str: The file path to the fake DICOM file.
    """
    # Erstellen des Dateipfads der ausgewählten DICOM-Datei basierend auf dem Index
    fake_file = f"{directory}IMG00{slice_index}_cropped.dcm"
    
    # Prüfen, ob die ausgewählte Fake-Datei existiert
    if not os.path.exists(fake_file):
        raise FileNotFoundError(f"{fake_file} existiert nicht.")
    
    # Laden der Fake-DICOM-Datei
    fake_dataset = pydicom.dcmread(fake_file)

    # Schleife über den Bereich der Originaldateien um die ausgewählte Datei und ihre Nachbarn
    for i in range(slice_index - 1, slice_index + 2):
        original_file = f"{directory}IMG00{i}_cropped.dcm"
        
        # Prüfen, ob die Datei existiert
        if not os.path.exists(original_file):
            print(f"{original_file} existiert nicht. Überspringen...")
            continue
        
        # Laden der Original-DICOM-Datei
        original_dataset = pydicom.dcmread(original_file)
        
        # Ersetzen der PixelData im originalen Dataset mit der aus dem Fake-Dataset
        original_dataset.PixelData = fake_dataset.PixelData
        
        # Festlegen des Ausgabedateinamens: 
        # Wenn es die ausgewählte Fake-Datei ist, den ursprünglichen Namen behalten
        # Sonst den Namen mit "_fake" speichern
        if original_file == fake_file:
            output_file = f"{output_directory}IMG00{i}_cropped.dcm"
        else:
            output_file = f"{output_directory}IMG00{i}_cropped_fake.dcm"
        
        # Speichern der modifizierten Datei
        original_dataset.save_as(output_file)
        
        print(f"PixelData erfolgreich von {fake_file} nach {original_file} kopiert und gespeichert als {output_file}.")
    
    print("Die ausgewählte Datei und ihre nächsten Nachbarn wurden erfolgreich verarbeitet.")
    
    # Return the path to the fake file (original fake slice)
    return f"{output_directory}IMG00{slice_index}_cropped.dcm"

# Function to plot SPR data with all required file paths as inputs
def plot_CT_matRad_MATLAB_r80(dicom_file, voi_target_file, max_rows_3d_file, max_rows_matlab_file, resolution_mm=1.09375, default_conversion="proCT"):
    
    # Function to read matrix from CSV
    def read_matrix(filepath):
        return pd.read_csv(filepath, header=None).values

    # Function to read all x and y positions from maxRows3D file
    def read_maxRows3D(filepath):
        max_rows_df = pd.read_csv(filepath, header=None)
        # Extract x, y values
        x_positions = max_rows_df.iloc[:, 0].to_numpy(dtype=float)  # Convert to float
        y_positions = max_rows_df.iloc[:, 1].to_numpy(dtype=float)  # Convert to float
        return x_positions, y_positions

    # Load the DICOM file and extract pixel data
    ct = pydicom.dcmread(dicom_file)
    pixel_data = ct.pixel_array

    # Save HU values to a CSV file for this slice
    hu_csv_file = dicom_file.replace(".dcm", "_HU.csv")
    np.savetxt(hu_csv_file, pixel_data, delimiter=",")
    
    # Load HU values from the CSV file
    hu_values_df = pd.read_csv(hu_csv_file, header=None).values

    # Define the conversion functions
    def DEMono79opt_conversion(HU):
        DEMono79opt_5k = (5000 - 3071) * (2.86 - 1.07) / (3071 - 130) + 2.86
        return np.interp(HU, [-1024, -990, -900, -150, -120, -70, -20, 70, 130, 3071, 5000],
                         [0.01, 0.01, 0.12, 0.83, 0.92, 0.95, 0.98, 1.06, 1.07, 2.86, DEMono79opt_5k])

    def proCT_conversion(HU):
        proCT_5k = (5000 - 3071) * (2.32 - 1.15) / (3071 - 184) + 2.32
        return np.interp(HU, [-1024, -1000, -71, 184, 3071, 5000],
                         [0.01, 0.01, 0.96, 1.15, 2.32, proCT_5k])

    # Apply the chosen conversion method
    if default_conversion == "DEMono79opt":
        SPR = DEMono79opt_conversion(hu_values_df)
    elif default_conversion == "proCT":
        SPR = proCT_conversion(hu_values_df)
    else:
        raise ValueError("Invalid default conversion method specified.")
    
    # Load the voiTarget matrix
    voi_target_matrix = read_matrix(voi_target_file)
    # Find positions of 1's in the voiTarget matrix
    target_positions = np.argwhere(voi_target_matrix == 1)
    
    # Calculate extent for imshow (extent = [xmin, xmax, ymin, ymax])
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

    # Adjust target positions to the correct scale using resolution
    x_positions_target = target_positions[:, 1] * resolution_mm
    y_positions_target = target_positions[:, 0] * resolution_mm

    # Scatter plot of target volume positions on the SPR colormap
    ax.scatter(x_positions_target, y_positions_target, color='white', label='Target Volume', s=200)

    # Load all maxRows3D values
    x_positions_3d, y_positions_3d = read_maxRows3D(max_rows_3d_file)

    # Scatter plot for all maxRows3D positions with color green
    ax.scatter(x_positions_3d, y_positions_3d, color='green', label='matRad (green)', s=10)

    # Load MATLAB data
    matlab_data_df = pd.read_csv(max_rows_matlab_file, header=None)
    x_positions_MATLAB = matlab_data_df.iloc[:, 0].to_numpy(dtype=float)
    y_positions_MATLAB = matlab_data_df.iloc[:, 1].to_numpy(dtype=float)

    # Scatter plot for MATLAB data using dots
    ax.scatter(x_positions_MATLAB, y_positions_MATLAB, color='red', label='MATLAB (blue)', s=20, marker='.')

    # Show the plot with legends
    plt.legend()
    plt.show()
    plt.close()


# Example usage: first run fake_pixel_data, then use its output in dicom_to_spr
# Comment out functions that are not relevant to use right now!

#fake_file_path = fake_pixel_data(44) # Fakes HU values from input slice onto neighbor slices and saves in DICOM_autofake
dicom_to_spr("../DICOM_custom/IMG0036_cropped.dcm", default_conversion="proCT") 

#plot_CT_matRad_MATLAB_r80(
#    dicom_file="../DICOM/IMG0036_cropped.dcm",
#    voi_target_file='../../matRad-master/voiTarget_page_2.csv',
#    max_rows_3d_file='../../matRad-master/x_y_energy_values_at_80_percent.csv', # cihangir.m 
#    max_rows_matlab_file='../../matRad-master/x_y_energy_values_at_80_percent_MATLAB.csv' # run_Loop.m
#)