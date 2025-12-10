import pandas as pd
import numpy as np
import pydicom
import matplotlib.pyplot as plt
import matplotlib.cm as cm

def plot_dose_profiles_with_ct(ct_image_file, dose_profile_files, conversion_method="proCT"):
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
        raise ValueError("Invalid conversion method specified.")

    # Prepare to combine dose profiles
    heatmap_combined = None
    xedges_combined = None
    yedges_combined = None

    for dose_profile_file in dose_profile_files:
        # Read the dose profile data
        df = pd.read_csv(dose_profile_file, comment='#', header=None)
        dp_datamatrix = np.array(df)
        x = dp_datamatrix[:, 0] * 1.09375
        y = dp_datamatrix[:, 1] * 1.09375
        dose = dp_datamatrix[:, 3]

        # Mirror the y values
        y = np.max(y) - y

        # Normalize the dose values to get relative dose (0% to 100%)
        dose = (dose / np.max(dose)) * 100

        # Create a grid for the heatmap
        heatmap, xedges, yedges = np.histogram2d(x, y, bins=(341, 225), weights=dose)

        # Normalize the heatmap to 0-100%
        heatmap = (heatmap / np.max(heatmap)) * 100

        if heatmap_combined is None:
            heatmap_combined = heatmap
            xedges_combined = xedges
            yedges_combined = yedges
        else:
            heatmap_combined += heatmap

    # Normalize combined heatmap to 0-100%
    heatmap_combined = (heatmap_combined / np.max(heatmap_combined)) * 100

    # Plotting
    fig, ax = plt.subplots(figsize=(15, 9))
    im = ax.imshow(SPR, cmap='gray', aspect='auto',
                   extent=[xedges_combined.min(), xedges_combined.max(),
                           yedges_combined.min(), yedges_combined.max()],
                   vmin=np.min(SPR), vmax=np.max(SPR))
    X, Y = np.meshgrid(xedges_combined, yedges_combined)
    c = ax.pcolormesh(X, Y, heatmap_combined.T, shading='auto', cmap='magma', alpha=0.6)
    cb = fig.colorbar(c, ax=ax, label='Relative Dose [%]')
    cb.ax.tick_params(labelsize=24)
    cb.set_label('Relative Dose [%]', fontsize=24)

    ax.scatter(152.73, 134.26, color='greenyellow', s=400, alpha=1, marker='.', label='Stopping position for $E$ = 150 MeV and $B_z$ = 1.5 T')
    ax.scatter(254.33, 150.72, color='cyan', s=400, alpha=1, marker='.', label='Stopping position for $E$ = 200 MeV and $B_z$ = 1.5 T')


    ax.legend(fontsize=20)

    plt.title('Proton trajectories for prostate patient with magnetic field', fontsize=24)
    plt.xlabel('x [mm]', fontsize=24)
    plt.ylabel('y [mm]', fontsize=24)
    plt.tick_params(axis='x', labelsize=24)
    plt.tick_params(axis='y', labelsize=24)
    plt.show()
    plt.close()

# Example call with two dose profiles
plot_dose_profiles_with_ct(
    '../DICOM/IMG0036_cropped.dcm',
    ['../Data/TOPAS_prostate_B15_200MeV.csv', '../Data/TOPAS_prostate_B15_150MeV.csv']
#    ['../Data/TOPAS_prostate_B15_150MeV.csv']
)
