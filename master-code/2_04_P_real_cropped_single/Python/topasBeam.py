import pandas as pd
import numpy as np
import pydicom
import matplotlib.pyplot as plt
import matplotlib.cm as cm

def plot_dose_profile_with_ct(dose_profile_file, ct_image_file, conversion_method="proCT"):

    # Read the CT image
    ct = pydicom.dcmread(ct_image_file)
    pixel_data = ct.pixel_array

    # Convert HU values to SPR values
    def DEMono79opt_conversion(HU):
        DEMono79opt_5k = (5000 - 3071) * (2.86 - 1.07) / (3071 - 130) + 2.86
        return np.interp(HU, [-1024, -990, -900, -150, -120, -70, -20, 70, 130, 3071, 5000],
                         [0.01, 0.01, 0.12, 0.83, 0.92, 0.95, 0.98, 1.06, 1.07, 2.86, DEMono79opt_5k])

    def proCT_conversion(HU):
        proCT_5k = (5000 - 3071) * (2.32 - 1.15) / (3071 - 184) + 2.32
        return np.interp(HU, [-1024, -1000, -71, 184, 3071, 5000],
                         [0.01, 0.01, 0.96, 1.15, 2.32, proCT_5k])

    if conversion_method == "DEMono79opt":
        SPR = DEMono79opt_conversion(pixel_data)
    elif conversion_method == "proCT":
        SPR = proCT_conversion(pixel_data)
    else:
        raise ValueError("Invalid default conversion method specified.")

    # Read the dose profile data
    df = pd.read_csv(dose_profile_file, comment='#', header=None)
    dp_datamatrix = np.array(df)
    x = dp_datamatrix[:, 0] * 1.09375
    y = dp_datamatrix[:, 1]
    dose = dp_datamatrix[:, 3] 

    # Mirror the y values
    y = np.max(y) - y

    # Normalize the dose values to get relative dose (0% to 100%)
    dose = (dose / np.max(dose)) * 100

    # Create a grid for the heatmap
    heatmap, xedges, yedges = np.histogram2d(x, y, bins=(341, 225), weights=dose)

    # Normalize the heatmap to 0-100%
    heatmap = (heatmap / np.max(heatmap)) * 100

    # Plotting
    fig, ax = plt.subplots(figsize=(15, 9))
    #im = ax.imshow(SPR, cmap=cm.gray, aspect='auto', extent=[xedges.min(), xedges.max(), yedges.min(), yedges.max()])
    X, Y = np.meshgrid(xedges, yedges)
    c = ax.pcolormesh(X, Y, heatmap.T, shading='auto', cmap='jet', alpha=0.6)
    cb = fig.colorbar(c, ax=ax, label='Relative Dose [%]')
    cb.ax.tick_params(labelsize=30)  # Set font size for colorbar ticks
    cb.set_label('Relative Dose [%]', fontsize=30)  # Set font size for colorbar label

    # Add a vertical dashed line at x = 40
    ax.axvline(x=40, color='black', linestyle='--', linewidth=1.5)

    plt.title('Proton trajectory with E = 100 MeV, B = 3 T', fontsize=30)
    plt.xlabel('x [mm]', fontsize=30)
    plt.ylabel('y [mm]', fontsize=30)
    # Set font size for tick labels on both axes
    plt.tick_params(axis='x', labelsize=30)  # Adjust tick size for x-axis
    plt.tick_params(axis='y', labelsize=30)  # Adjust tick size for y-axis
    plt.show()
    plt.close()

# Example call
#plot_dose_profile_with_ct('trajectory_TOPAS.csv', '../DICOM/IMG0036_cropped.dcm')
plot_dose_profile_with_ct('../TOPAS/TOPAS_Beam.csv', '../DICOM/IMG0036_cropped.dcm')
