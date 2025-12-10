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
    x = dp_datamatrix[:, 0] * 0.2
    y = dp_datamatrix[:, 1] * 0.2
    dose = dp_datamatrix[:, 3]

    # Normalize the dose values to get relative dose (0% to 100%)
    dose = (dose / np.max(dose)) * 100

    # Mirror the y values
    y = np.max(y) - y

    # Threshold for dose
    dose_threshold = 30.0  # 70% of the max dose

    # Filter x, y where dose >= dose_threshold
    mask_high_dose = dose >= dose_threshold
    x_high_dose = x[mask_high_dose]
    y_high_dose = y[mask_high_dose]

    # Find y values where x crosses zero with high dose
    epsilon = 1e-6  # tolerance for considering x as zero
    y_crossings = y_high_dose[np.abs(x_high_dose) < epsilon]

    print(f"The beam crosses x=0 at the following y values (with dose >= {dose_threshold}%): {y_crossings}")
    radius = y_crossings / 2
    print(f"Resulting radius is {radius} cm")

    # Create a grid for the heatmap
    heatmap, xedges, yedges = np.histogram2d(x, y, bins=(700, 700), weights=dose)
    #heatmap, xedges, yedges = np.histogram2d(x, y, weights=dose)

    # Plotting
    fig, ax = plt.subplots(figsize=(15, 9))
    # im = ax.imshow(SPR, cmap=cm.gray, aspect='auto', extent=[xedges.min(), xedges.max(), yedges.min(), yedges.max()])
    X, Y = np.meshgrid(xedges, yedges)
    c = ax.pcolormesh(X, Y, heatmap.T, shading='auto', cmap='jet', alpha=0.6)
    fig.colorbar(c, ax=ax, label='Relative Dose (%)')

    plt.title('Proton Beam Trajectory inside a vacuum')
    plt.xlabel('Depth [mm]')
    plt.ylabel('Height [mm]')
    plt.show()
    plt.close()

# plot_dose_profile_with_ct('trajectory_TOPAS.csv', '../DICOM/IMG0036_cropped.dcm')
plot_dose_profile_with_ct('../TOPAS/TOPAS_Beam.csv', '../DICOM/IMG0036_cropped.dcm')
