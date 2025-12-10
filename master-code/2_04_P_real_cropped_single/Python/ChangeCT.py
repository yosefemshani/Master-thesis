import os
import glob
import pydicom

# Path to the directory containing DICOM files
path_dicom = '../DICOM_fake/'

# Find all DICOM files in the directory
dicom_files = glob.glob(os.path.join(path_dicom, '*.dcm'))

# Sort the DICOM files based on their names
dicom_files.sort()

# Determine the base index for IMG0036_cropped_fake.dcm
base_index = None
for i, file in enumerate(dicom_files):
    if 'IMG0036_cropped.dcm' in file:
        base_index = i
        break

if base_index is None:
    raise FileNotFoundError("Base file IMG0036_cropped.dcm not found")

# Loop through all files and modify the ImagePositionPatient attribute
for i, file in enumerate(dicom_files):
    ds = pydicom.dcmread(file)
    x, y = 0, 0
    z = (i - base_index) * 2
    ds.ImagePositionPatient = [x, y, z]
    ds.save_as(file)

print("DICOM files updated successfully.")
