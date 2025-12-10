import numpy as np
import pandas as pd
import pydicom
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.patches import FancyArrowPatch

# Function to plot SPR data with optional simplification
def plot_CT_matRad_MATLAB_r80(dicom_file, voi_target_file, max_rows_3d_file, max_rows_matlab_file, matlab_shifted_file=None, 
                               resolution_mm=1.09375, default_conversion="proCT", simplify_spr=False, show_matlab_shifted=False):
    
    # Function to read matrix from CSV
    def read_matrix(filepath):
        return pd.read_csv(filepath, header=None).values

    # Function to read all x and y positions from maxRows3D file
    def read_maxRows3D(filepath):
        max_rows_df = pd.read_csv(filepath, header=None)
        x_positions = max_rows_df.iloc[:, 0].to_numpy(dtype=float)
        y_positions = max_rows_df.iloc[:, 1].to_numpy(dtype=float)
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
    
    # Capture the original min and max SPR values (before simplification)
    original_min = np.min(SPR)
    original_max = np.max(SPR)
    
    # Simplify the SPR values if simplify_spr is True
    if simplify_spr:
        outer_layer_value = 0.2
        body_contour_value = 0.9
        
        # Create a mask for the outer air layer
        outer_layer_mask = (SPR <= outer_layer_value + 0.05)

        # Everything inside the outer layer will be treated as part of the body contour
        body_contour_mask = ~outer_layer_mask

        # Apply the simplification: set all internal anatomy to the body contour value
        simplified_spr = np.copy(SPR)
        simplified_spr[body_contour_mask] = body_contour_value

        SPR = simplified_spr

    # Load the voiTarget matrix
    voi_target_matrix = read_matrix(voi_target_file)
    target_positions = np.argwhere(voi_target_matrix == 1)
    
    # Calculate extent for imshow
    num_rows, num_cols = SPR.shape
    x_extent = num_cols * resolution_mm
    y_extent = num_rows * resolution_mm
    extent = [0, x_extent, 0, y_extent]

    # Plot the colormap with the fixed vmin and vmax (from original data)
    fig, ax = plt.subplots(figsize=(16, 9))
    im = ax.imshow(SPR, cmap=cm.jet, extent=extent, aspect='auto', vmin=original_min, vmax=original_max)
    plt.colorbar(im, ax=ax, label='SPR Value')
    plt.title('Proton stopping positions extracted from matRad', fontsize=24)
    #plt.xlim(155,220)
    #plt.ylim(90,170)
    plt.xlabel('x [mm]', fontsize=24)
    plt.ylabel('y [mm]', fontsize=24)

    # Set font size for x and y tick labels
    ax.tick_params(axis='x', labelsize=24)
    ax.tick_params(axis='y', labelsize=24)

    # Adjust target positions to the correct scale using resolution
    x_positions_target = target_positions[:, 1] * resolution_mm
    y_positions_target = target_positions[:, 0] * resolution_mm

    # Scatter plot of target volume positions on the SPR colormap
    ax.scatter(x_positions_target, y_positions_target, color='white', label='Target', s=250)

    # Load all maxRows3D values
    x_positions_3d, y_positions_3d = read_maxRows3D(max_rows_3d_file)

    # Scatter plot for all maxRows3D positions with color green
    ax.scatter(x_positions_3d, y_positions_3d, color='royalblue', label='matRad stopping positions', s=40)

    # Load MATLAB data
    matlab_data_df = pd.read_csv(max_rows_matlab_file, header=None)
    x_positions_MATLAB = matlab_data_df.iloc[:, 0].to_numpy(dtype=float)
    y_positions_MATLAB = matlab_data_df.iloc[:, 1].to_numpy(dtype=float)

    # Scatter plot for MATLAB data using dots
    #ax.scatter(x_positions_MATLAB, y_positions_MATLAB, color='red', label='Developed algorithm stopping positions ($B_z$ = 1.5 T) optimized', s=40)

    # Load MATLAB shifted positions if requested
    if show_matlab_shifted and matlab_shifted_file is not None:
        matlab_shifted_data_df = pd.read_csv(matlab_shifted_file, header=None)
        x_positions_shifted = matlab_shifted_data_df.iloc[:, 0].to_numpy(dtype=float)
        y_positions_shifted = matlab_shifted_data_df.iloc[:, 1].to_numpy(dtype=float)

        # Scatter plot for shifted MATLAB data using light red
        #ax.scatter(x_positions_shifted, y_positions_shifted, color='orange', label='Developed algorithm stopping positions ($B_z$ = 0 T)', s=40, alpha=1.0) #0.4 for all

        # List of indices where arrows should be
        #  drawn
        #arrow_indices = [-66, -79, -52, -38, -25, -13, -4, -1] # all last bixels of each row
        arrow_indices = [-1] # all last bixels of each row


        ## Loop through the indices to draw arrows
        #for index in arrow_indices:
        #    if len(x_positions_shifted) > index and len(x_positions_MATLAB) > index:
        #        # Access the corresponding spots for both shifted and MATLAB data
        #        x_shifted_pos = x_positions_shifted[index]
        #        y_shifted_pos = y_positions_shifted[index]
        #        x_matlab_pos = x_positions_MATLAB[index]
        #        y_matlab_pos = y_positions_MATLAB[index]
        #
        #        # Create a curved arrow
        #        arrow = FancyArrowPatch(
        #            posA=(x_shifted_pos, y_shifted_pos), 
        #            posB=(x_matlab_pos, y_matlab_pos),
        #            arrowstyle='->',
        #            mutation_scale=15,
        #            color='orange',
        #            lw=2.5,
        #            linestyle='-',
        #            alpha=1.0,
        #            connectionstyle="arc3,rad=-0.3"  # 'rad' controls the curvature
        #        )
        #        # Add the arrow to the plot
        #        ax.add_patch(arrow)

    # Show the plot with legends
    plt.legend(loc='best', fontsize=20)
    plt.show()
    plt.close()

# Example usage
plot_CT_matRad_MATLAB_r80(
    dicom_file="../DICOM/IMG0036_cropped.dcm",
    voi_target_file='../../matRad-master/voiTarget_page_2.csv',
    max_rows_3d_file='../../matRad-master/x_y_energy_values_at_80_percent.csv', 
    max_rows_matlab_file='../../matRad-master/x_y_energy_values_at_80_percent_gradient.csv',
    matlab_shifted_file='../../matRad-master/x_y_energy_values_at_80_percent_MATLAB.csv', 
    simplify_spr=True,
    show_matlab_shifted=True
)
