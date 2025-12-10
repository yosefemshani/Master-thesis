import nibabel as nib
import matplotlib.pylab as plt
import matplotlib as mpl
import SimpleITK as sitk
import random
import shutil
import os, glob
from PIL import Image
import numpy as np
import pandas as pd
from scipy import stats
from scipy import interpolate
import gzip
import shutil
from scipy import ndimage
from scipy.ndimage import zoom
import h5py
import datetime
from pathlib import Path
import tempfile
import pydicom as dicom
from pydicom.dataset import Dataset, FileMetaDataset, FileDataset
from pydicom.uid import UID, ExplicitVRLittleEndian
import datetime

# Path to the directory containing original DICOM files
path_dicom = '../DICOM_fake/'

# Path to the directory where modified DICOM files will be saved
path_dicom_water = '../DICOM_water/'

# Create the new directory if it doesn't exist
os.makedirs(path_dicom_water, exist_ok=True)

# Find all DICOM files in the original directory
dicom_files = glob.glob(path_dicom + '*.dcm')

# Ensure there are exactly 62 slices
assert len(dicom_files) == 62, f"Expected 62 slices, but found {len(dicom_files)}"

# Loop through each DICOM file, modify the pixel array, and save to the new directory
for i, file_this in enumerate(dicom_files):
    # Read the DICOM file
    exp_dcm = dicom.dcmread(file_this)
    
    # Extract the pixel array (image data)
    image = exp_dcm.pixel_array
    
    # Set all values in the image array to 0
    image[:] = 0
    
    # Convert the image to int16 (which is common for CT scans)
    CT_i = image.astype(np.int16)
    
    # Replace the original pixel data with the modified pixel data
    exp_dcm.PixelData = CT_i.tobytes()
    
    # Save the modified DICOM file to the new directory
    new_file_path = os.path.join(path_dicom_water, os.path.basename(file_this))
    exp_dcm.save_as(new_file_path)
    
    # Optionally, print the file name of the saved DICOM
    if i == 0:
        print(f"First modified DICOM saved to: {new_file_path}")
