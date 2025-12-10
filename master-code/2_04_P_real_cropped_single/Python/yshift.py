import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d

# Define Gaussian function for fitting
def gaussian(y, a, mean, sigma):
    return a * np.exp(-0.5 * ((y - mean) / sigma) ** 2)

# Read the csv-file
#outputfile_topas = '../Data/TOPAS_water_B0_200MeV.csv'
outputfile_topas = '../TOPAS/TOPAS_Beam.csv'
df = pd.read_csv(outputfile_topas, comment='#', header=None)

# Convert the dataframe to a numpy array
topas_datamatrix = np.array(df)

# Extract depth (x), height (y), and dose values
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

# Sum dose values across all z-bins (integrating across slices)
integrated_dose_z = np.sum(dose_grid, axis=2)

# Integrate the dose along the y-axis (collapse into 1D dose profile along x)
integrated_dose = np.sum(integrated_dose_z, axis=1)

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

print(f"R80 value at depth: {r80_x:.2f} mm")

# Define y axis (height) in mm for the plot
y_bin_size = 0.109375 * 10  # mm per bin in y
y_axis = np.arange(n_y_bins) * y_bin_size  # height in mm

# Initialize a list to store the (x, y, dose) values
yshift_array = []

# Boolean flag to track when the first Gaussian fit fails
fitting_failed = False

# Iterate over all x_bins
for x_bin in range(n_x_bins):
    if fitting_failed:
        x_value = x_bin * x_bin_size  # Convert x_bin to depth in mm
        yshift_array.append([x_value, 0.0, 0.0])
        continue

    # Extract dose values along the y-axis for the current x_bin (summed over z)
    y_dose_distribution = np.sum(dose_grid[x_bin, :, :], axis=1)  # Sum over all z-bins

    # Check if the dose is all zeros or close to zero
    if np.max(y_dose_distribution) <= 0:
        x_value = x_bin * x_bin_size  # Convert x_bin to depth in mm
        yshift_array.append([x_value, 0.0, 0.0])
        continue

    # Normalize the dose distribution
    y_dose_distribution_norm = (y_dose_distribution / np.max(y_dose_distribution)) * 100

    if np.isnan(y_dose_distribution_norm).any() or np.isinf(y_dose_distribution_norm).any():
        x_value = x_bin * x_bin_size  # Convert x_bin to depth in mm
        yshift_array.append([x_value, 0.0, 0.0])
        continue

    try:
        # Fit a Gaussian function to the normalized y-dose distribution
        popt, _ = curve_fit(gaussian, y_axis, y_dose_distribution_norm, 
                            p0=[np.max(y_dose_distribution_norm), np.mean(y_axis), np.std(y_axis)], 
                            maxfev=4000)

        # Extract the mean (which corresponds to the y position of the Gaussian peak)
        mean_fit = popt[1]

        # Get the absolute dose value at the position where the Gaussian peak occurs
        closest_y_idx = np.argmin(np.abs(y_axis - mean_fit))
        absolute_dose = y_dose_distribution[closest_y_idx]

        # Append the current x_bin (converted to depth), y value (mean_fit), and absolute dose to the array
        x_value = x_bin * x_bin_size  # Convert x_bin to depth in mm
        yshift_array.append([x_value, mean_fit, absolute_dose])

    except RuntimeError:
        # If curve fitting fails, store 0.0 values for this x_bin
        fitting_failed = True
        x_value = x_bin * x_bin_size  # Convert x_bin to depth in mm
        yshift_array.append([x_value, 0.0, 0.0])

# Convert the list to a DataFrame
yshift_df = pd.DataFrame(yshift_array, columns=['x (mm)', 'y (mm)', 'dose (absolute)'])

# Filter out rows where both y and dose are zero
yshift_df_filtered = yshift_df[(yshift_df['y (mm)'] > 0) & (yshift_df['dose (absolute)'] > 0)].copy()

# Calculate cumulative distance 's'
s_values = [0.0]  # Start with s = 0 for the first point
for i in range(1, len(yshift_df_filtered)):
    delta_x = yshift_df_filtered['x (mm)'].iloc[i] - yshift_df_filtered['x (mm)'].iloc[i - 1]
    delta_y = yshift_df_filtered['y (mm)'].iloc[i] - yshift_df_filtered['y (mm)'].iloc[i - 1]
    delta_s = np.sqrt(delta_x**2 + delta_y**2)
    s_values.append(s_values[-1] + delta_s)

yshift_df_filtered.loc[:, 's (mm)'] = s_values

# Normalize the dose
yshift_df_filtered.loc[:, 'dose_norm'] = (yshift_df_filtered['dose (absolute)'] / np.max(yshift_df_filtered['dose (absolute)'])) * 100

# Apply a threshold to remove very small dose values (e.g., less than 1%)
dose_threshold = 1.0
threshold_idx = np.where(yshift_df_filtered['dose_norm'] < dose_threshold)[0]

# Keep the next 4 points after the dose drops below the threshold
if len(threshold_idx) > 0:
    first_threshold_idx = threshold_idx[0]
    end_idx = first_threshold_idx + 4  # Keep 4 more points
    yshift_df_filtered = yshift_df_filtered.iloc[:min(end_idx, len(yshift_df_filtered))]


# Toggle to show Gaussian fit for a specific bin
show_gaussian_fit = True  # Set to False to disable

# Specify the bin for Gaussian example (this is where the Gaussian will be plotted)
example_bin = 40  # You can change this bin as needed

# Plot the Gaussian fit for a specific bin if requested
if show_gaussian_fit:
    # Extract dose values along the y-axis for the chosen x_bin (summed over z)
    y_dose_distribution_example = np.sum(dose_grid[example_bin, :, :], axis=1)  # Sum over z-bins
    
    # Normalize the dose distribution
    y_dose_distribution_example_norm = (y_dose_distribution_example / np.max(y_dose_distribution_example)) * 100
    
    # Fit a Gaussian function
    popt, _ = curve_fit(gaussian, y_axis, y_dose_distribution_example_norm, 
                        p0=[np.max(y_dose_distribution_example_norm), np.mean(y_axis), np.std(y_axis)], 
                        maxfev=4000)
    
    # Extract fitted Gaussian
    gaussian_fit_example = gaussian(y_axis, *popt)

    # Plot Gaussian fit
    plt.figure(figsize=(15, 9))
    plt.plot(y_axis, y_dose_distribution_example_norm, label='Normalized Dose Distribution', color='blue')
    plt.plot(y_axis, gaussian_fit_example, label='Gaussian Fit', color='orange', linestyle='--')
    plt.xlabel('y [mm]', fontsize=30)
    plt.ylabel('Relative Dose [%]', fontsize=30)
    plt.title(f'Gaussian fit for x = 40 mm',fontsize=30)
    #plt.title(f'Gaussian fit for x bin {example_bin}')
    plt.legend(loc="best",fontsize=21)
    plt.grid(False)
    # Set font size for tick labels on both axes
    plt.tick_params(axis='x', labelsize=30)  # Adjust tick size for x-axis
    plt.tick_params(axis='y', labelsize=30)  # Adjust tick size for y-axis
    plt.show()


# Plot the depth-dose distribution for all z-bins (summarized result)
plt.figure(figsize=(15, 9))

# Plot the depth-dose distribution averaged across all z-bins
plt.plot(x_axis, integrated_dose_norm, color='blue')

# Add a horizontal line for 80% of the max dose
plt.axhline(y=r80_value, color='gray', linestyle='--', label='80% of Max. Dose')

# Add a red dot at the R80 position
plt.plot(r80_x, r80_value, 'ro', label=f'$R_{{80}}$ = {r80_x:.2f} mm')  # Add label to legend

# Add labels and title
plt.xlabel('x [mm]', fontsize=24)
plt.ylabel('Relative Dose [%]', fontsize=24)
plt.title('Percentage depth dose curve for protons in water', fontsize=24)
plt.legend(loc='best',fontsize=20)
plt.grid(True)

plt.tick_params(axis='x', labelsize=24)  # Adjust tick size for x-axis
plt.tick_params(axis='y', labelsize=24)  # Adjust tick size for y-axis
# Show the depth-dose distribution plot
plt.show()




# Interpolation to find the dose at the calculated R80 value
r80_s = r80_x  # Use the calculated R80 value
s_values = yshift_df_filtered['s (mm)']
dose_values = yshift_df_filtered['dose_norm']

# Create interpolation function
interp_function = interp1d(s_values, dose_values, bounds_error=False, fill_value="extrapolate")

# Get interpolated dose at r80_s
interpolated_dose = interp_function(r80_s)

# Plot the dose vs distance 's' (mm)
plt.figure(figsize=(15, 9))
plt.plot(yshift_df_filtered['s (mm)'], yshift_df_filtered['dose_norm'],color='blue')
plt.xlabel(r'${s_{\mathrm{L2}}}$ [mm]', fontsize=24)
plt.ylabel('Relative Dose [%]',fontsize=24)
plt.grid(True)

# Plot the interpolated R80 point
plt.plot(r80_s, interpolated_dose, 'ro', label=f'Dose = {interpolated_dose:.2f}% ($R_{{80}}$ = {r80_s:.2f} mm)')
#plt.plot(r80_s, interpolated_dose, 'ro', label=f'R$_{{80}}$ = {r80_s:.2f} mm (Dose = {interpolated_dose:.2f}%)')
plt.title(r'Analyzing rel. dose for $s_{\mathrm{L1}} = s_{\mathrm{L2}}$ in water without magnetic field', fontsize=24)
plt.legend(loc='best',fontsize=20)
# Set font size for tick labels on both axes
plt.tick_params(axis='x', labelsize=24)  # Adjust tick size for x-axis
plt.tick_params(axis='y', labelsize=24)  # Adjust tick size for y-axis
# Show the plot
plt.show()

# Search for `s` value based on input relative dose
search_rel_dose = True  # Toggle on or off

# Allow user to input the desired relative dose % (default: 83.69%)
target_rel_dose = 83.69 #interpolated_dose # ENTER HERE HOW MUCH r80 we had for B = 0 T !!!

# Find the maximum dose and its corresponding index
max_dose_idx = np.argmax(yshift_df_filtered['dose_norm'])
max_dose_value = yshift_df_filtered['dose_norm'].iloc[max_dose_idx]

# Now, slice the data after the maximum dose
falling_side_df = yshift_df_filtered.iloc[max_dose_idx:]  # Data after reaching 100%

# Check if the target relative dose exists on the falling side
if target_rel_dose <= max_dose_value:
    # Interpolation function to find `s` for the given relative dose, only after the peak
    interp_s_for_dose_falling = interp1d(falling_side_df['dose_norm'], falling_side_df['s (mm)'], 
                                         kind='linear', bounds_error=False, fill_value="extrapolate")
    
    # Get the corresponding `s` value for the target relative dose (83.69%)
    target_s_value_falling = interp_s_for_dose_falling(target_rel_dose)
    
    # Plot the dose vs distance 's' without the R80 point, but with the searched target dose point
    plt.figure(figsize=(15, 9))
    plt.plot(yshift_df_filtered['s (mm)'], yshift_df_filtered['dose_norm'], color='blue')
    
    # Plot the interpolated point at the target dose
    plt.plot(target_s_value_falling, target_rel_dose, 'ro', label=f's = {target_s_value_falling:.2f} mm for rel. dose = {target_rel_dose:.2f} %')
    
    plt.xlabel(r'${s_{\mathrm{L1}}}$ [mm]', fontsize=24)
    plt.ylabel('Relative Dose [%]',fontsize=24)
    plt.title(r'Analyzing ${s_{\mathrm{L1}}}$ for ${s_{\mathrm{L1}}} = {s_{\mathrm{L2}}}$ with given rel. dose in water with magnetic field', fontsize=24)
    plt.legend(loc='best',fontsize=20)
    plt.grid(True)
    # Set font size for tick labels on both axes
    plt.tick_params(axis='x', labelsize=24)  # Adjust tick size for x-axis
    plt.tick_params(axis='y', labelsize=24)  # Adjust tick size for y-axis
    plt.show()
else:
    print(f"Target relative dose {target_rel_dose}% exceeds the maximum dose.")