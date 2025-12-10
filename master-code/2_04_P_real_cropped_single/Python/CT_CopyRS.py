import pydicom

# Load the RS structure file
rt_struct = pydicom.dcmread('../DICOM/RS1.2.752.243.1.1.20170503093409104.7100.87740_transformed.dcm', force=True)

# Initialize variables to store x and y values from IMG0036_cropped.dcm
target_z = -26.54271
img0036_contour_data = []

# Find and store the x and y values for IMG0036
for roi in rt_struct.ROIContourSequence:
    for contour in roi.ContourSequence:
        contour_data = contour.ContourData
        # Check if the z value matches target_z
        if contour_data[2] == target_z:
            img0036_contour_data = contour_data
            break
    if img0036_contour_data:
        break

# Separate the x and y values
img0036_x_y_values = [(img0036_contour_data[i], img0036_contour_data[i+1]) for i in range(0, len(img0036_contour_data), 3)]

# Replace contour data in all slices
for roi in rt_struct.ROIContourSequence:
    for contour in roi.ContourSequence:
        contour_data = contour.ContourData
        # Create a new contour list with the x and y values from IMG0036 and the original z values
        new_contour_data = []
        for i in range(0, len(contour_data), 3):
            z_value = contour_data[i+2]
            x_y_index = (i // 3) % len(img0036_x_y_values)  # Loop over the img0036_x_y_values if needed
            new_contour_data.extend([img0036_x_y_values[x_y_index][0], img0036_x_y_values[x_y_index][1], z_value])
        # Update the contour data
        contour.ContourData = new_contour_data

# Save the modified RS structure file
rt_struct.save_as('../DICOM_DoubleSlice/RSfake.dcm')
