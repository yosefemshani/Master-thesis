import numpy as np
import pandas as pd
import pydicom
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.lines as mlines

# Define the conversion functions
def DEMono79opt_conversion(HU):
    DEMono79opt_5k = (5000 - 3071) * (2.86 - 1.07) / (3071 - 130) + 2.86
    return np.interp(HU, [-1024, -990, -900, -150, -120, -70, -20, 70, 130, 3071, 5000],
                     [0.01, 0.01, 0.12, 0.83, 0.92, 0.95, 0.98, 1.06, 1.07, 2.86, DEMono79opt_5k])

def proCT_conversion(HU):
    proCT_5k = (5000 - 3071) * (2.32 - 1.15) / (3071 - 184) + 2.32
    return np.interp(HU, [-1024, -1000, -71, 184, 3071, 5000],
                     [0.01, 0.01, 0.96, 1.15, 2.32, proCT_5k])

# Function to read the CT image and convert to SPR values
def get_spr_values(ct_image_file, conversion_method="proCT"):
    ct = pydicom.dcmread(ct_image_file)
    pixel_data = ct.pixel_array

    if conversion_method == "DEMono79opt":
        return DEMono79opt_conversion(pixel_data)
    elif conversion_method == "proCT":
        return proCT_conversion(pixel_data)
    else:
        raise ValueError("Invalid default conversion method specified.")

# Function to plot both trajectories
def plot_combined_trajectories(matlab_trajectory_file, topas_trajectory_file, ct_image_file, conversion_method="proCT"):
    # Get SPR values
    SPR = get_spr_values(ct_image_file, conversion_method)

    # Read MATLAB trajectory data
    matlab_trajectory_data = np.loadtxt(matlab_trajectory_file, delimiter=',')

    # Extract the last row of MATLAB trajectory data
    depth_in_matlab = matlab_trajectory_data[-1, 0]  # Last x entry (depth)
    height_in_matlab = matlab_trajectory_data[-1, 1] * 1.09375  # Last y entry (height)

    # Read TOPAS trajectory data
    topas_df = pd.read_csv(topas_trajectory_file, comment='#', header=None)
    topas_data_matrix = np.array(topas_df)
    topas_x = topas_data_matrix[:, 0] * 1.09375
    topas_y = topas_data_matrix[:, 1] * 1.09375
    topas_dose = topas_data_matrix[:, 3] 

    # Mirror the TOPAS y values
    topas_y = np.max(topas_y) - topas_y

    # Ensure that topas_dose is normalized between 0 and 100%
    topas_norm = (topas_dose / np.max(topas_dose)) * 100
    topas_norm = np.clip(topas_norm, 0, 100)  # Clipping to make sure values are within the 0–100 range

    # Step 1: Get the count of points in each bin (for averaging purposes)
    count, _, _ = np.histogram2d(topas_x, topas_y, bins=(341, 225))

    # Step 2: Create the weighted sum as before
    weighted_sum, xedges, yedges = np.histogram2d(topas_x, topas_y, bins=(341, 225), weights=topas_norm)

    # Step 3: Compute the average dose by dividing the weighted sum by the count
    heatmap = np.divide(weighted_sum, count, out=np.zeros_like(weighted_sum), where=(count != 0))

    # Normalize the heatmap to 0-100%
    heatmap = (heatmap / np.max(heatmap)) * 100

    # Plotting
    fig, ax = plt.subplots(figsize=(15, 9))
    X, Y = np.meshgrid(xedges, yedges)
    c = ax.pcolormesh(X, Y, heatmap.T, shading='auto', cmap='jet', alpha=0.6)
    cb = fig.colorbar(c, ax=ax, label='Relative Dose [%]')
    cb.ax.tick_params(labelsize=24)  # Set font size for colorbar ticks
    cb.set_label('Relative Dose [%]', fontsize=24)  # Set font size for colorbar label

    # Plot MATLAB trajectory
    ax.plot(matlab_trajectory_data[:, 0], matlab_trajectory_data[:, 1] * 1.09375, 'w-', linewidth=15, alpha=1, label='MATLAB')

    # Create custom legend entries
    matlab_legend = mlines.Line2D([], [], color='white', linewidth=13, label='MATLAB')

    # Add the custom legend with both entries
    plt.legend(handles=[matlab_legend],loc='upper left',fontsize=20)

    plt.title('Proton beam trajectories resulting from \ndeveloped algorithm and Monte Carlo', fontsize=24)
    plt.xlabel('x [mm]', fontsize=24)
    plt.ylabel('y [mm]', fontsize=24)
    # Set font size for tick labels on both axes
    plt.tick_params(axis='x', labelsize=24)  # Adjust tick size for x-axis
    plt.tick_params(axis='y', labelsize=24)  # Adjust tick size for y-axis

    plt.show()
    plt.close()

def save_quantitative_trajectory_analysis(matlab_trajectory_file, topas_trajectory_file, ct_image_file, output_file, conversion_method="proCT"):
    # Get SPR values for plotting purposes
    SPR = get_spr_values(ct_image_file, conversion_method)

    # Read MATLAB trajectory data (x, y positions)
    matlab_trajectory_data = np.loadtxt(matlab_trajectory_file, delimiter=',')

    # Read TOPAS trajectory data (x, y, z, dose)
    topas_df = pd.read_csv(topas_trajectory_file, comment='#', header=None)
    topas_data_matrix = np.array(topas_df)

    # Extract x, y, and dose columns from TOPAS data
    topas_x = topas_data_matrix[:, 0] * 1.09375
    topas_y = np.max(topas_data_matrix[:, 1]) - topas_data_matrix[:, 1]  # Mirror the y-values
    topas_dose = topas_data_matrix[:, 3] * 100

    # Filter out positions with low dose (e.g., dose > 1% of max dose)
    dose_threshold = 0.8 * np.max(topas_dose)
    valid_positions = topas_dose > dose_threshold
    topas_x_filtered = topas_x[valid_positions]
    topas_y_filtered = topas_y[valid_positions]

    # Determine the max size of both trajectories for padding
    max_size = max(len(matlab_trajectory_data), len(topas_x_filtered))

    # Pad the MATLAB data and TOPAS data with NaN to match sizes
    matlab_x_padded = np.pad(matlab_trajectory_data[:, 0], (0, max_size - len(matlab_trajectory_data)), constant_values=np.nan)
    matlab_y_padded = np.pad(matlab_trajectory_data[:, 1], (0, max_size - len(matlab_trajectory_data)), constant_values=np.nan)
    topas_x_padded = np.pad(topas_x_filtered, (0, max_size - len(topas_x_filtered)), constant_values=np.nan)
    topas_y_padded = np.pad(topas_y_filtered, (0, max_size - len(topas_y_filtered)), constant_values=np.nan)

    # Combine MATLAB and TOPAS trajectories into a single array
    combined_trajectory_data = np.column_stack((
        matlab_x_padded, matlab_y_padded,  # MATLAB x, y
        topas_x_padded, topas_y_padded     # TOPAS x, y
    ))

    # Save to CSV file for quantitative analysis
    np.savetxt(output_file, combined_trajectory_data, delimiter=',', 
               header='MATLAB_X,MATLAB_Y,TOPAS_X,TOPAS_Y', comments='')
    
    print(f"Quantitative trajectory data saved to {output_file}")


def plot_padded_trajectories(comparison_file):
    # Load the saved trajectory data (MATLAB and TOPAS x, y)
    trajectory_data = np.loadtxt(comparison_file, delimiter=',', skiprows=1)
    
    # Extract MATLAB and TOPAS data from the file
    matlab_x = trajectory_data[:, 0]
    matlab_y = trajectory_data[:, 1]
    topas_x = trajectory_data[:, 2]
    topas_y = trajectory_data[:, 3]

    # Create a new plot
    plt.figure(figsize=(10, 6))
    
    # Plot MATLAB trajectory
    plt.plot(matlab_x, matlab_y, 'k-', label='MATLAB Trajectory', linewidth=6, alpha=0.8)
    
    # Plot TOPAS trajectory
    plt.plot(topas_x, topas_y, 'g-', label='TOPAS Trajectory', linewidth=1.5, alpha=0.8)

    # Set plot title and labels
    plt.title('Padded MATLAB and TOPAS Trajectories')
    plt.xlabel('X Position [mm]')
    plt.ylabel('Y Position [mm]')
    
    # Add a legend
    plt.legend()

    # Show the plot
    plt.grid(True)
    plt.show()


# Call the function with the appropriate file paths
# save_quantitative_trajectory_analysis('trajectory_mm.csv', '../TOPAS/TOPAS_Beam.csv', '../DICOM/IMG0036_cropped.dcm', 'quantitative_trajectory_comparison.csv')

# Call the function with the appropriate comparison file path
# plot_padded_trajectories('quantitative_trajectory_comparison.csv')

# Call the function with the appropriate file paths
#plot_combined_trajectories('../Data/trajectory_mm_water_B0_E200MeV.csv', '../TOPAS/TOPAS_Beam.csv', '../DICOM/IMG0036_cropped.dcm')
plot_combined_trajectories('../Data/trajectory_mm_water_B0_E200MeV.csv', '../Data/TOPAS_water_B0_200MeV.csv', '../DICOM/IMG0036_cropped.dcm')