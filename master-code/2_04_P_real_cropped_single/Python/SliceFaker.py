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

# Path to the directory containing DICOM files
path_dicom = '../DICOM/'

# Find all DICOM files in the directory
dicom_files = glob.glob(path_dicom + '*.dcm')

file_this = dicom_files[26]

# Read the DICOM file
exp_dcm = dicom.dcmread(file_this)

# Extract the pixel array (image data)
image = exp_dcm.pixel_array

# Plot the image
#plt.imshow(image)
#plt.title(f"Image from {file_this}")
#plt.show()

print('Original file image position:', exp_dcm.ImagePositionPatient)

path_fake = '../DICOM_fake'
os.makedirs(path_fake, exist_ok=True)

Thickness_ = exp_dcm.SliceThickness
PixelSpace_ = exp_dcm.PixelSpacing
CT_i = image.astype(np.int16)

for i in range(10, 72):
    # Copy the original DICOM file
    dup_dicom = dicom.dcmread(file_this)
    
    # Modify the DICOM file properties
    dup_dicom.PatientName = 'TestName'
    dup_dicom.Rows, dup_dicom.Columns = CT_i.shape
    dup_dicom.SliceThickness = Thickness_
    dup_dicom.PixelSpacing = PixelSpace_
    dup_dicom.InstanceNumber = i
    
    # Calculate the new ImagePositionPatient
    position_z = (i - 36) * Thickness_
    dup_dicom.ImagePositionPatient = [0, 0, position_z]
    dup_dicom.PixelData = CT_i.tobytes()
    
    # Create the file name
    file_name = f'../DICOM_fake/Fake{str(i).zfill(4)}.dcm'
    
    # Save the new DICOM file
    dup_dicom.save_as(file_name)
    print(f'Saved {file_name} with image position: {dup_dicom.ImagePositionPatient}')

    # Plot the image to verify
    check_example = dicom.dcmread(file_name)
    #plt.imshow(check_example.pixel_array)
    #plt.title(f"Fake Image from {file_this}")
    #plt.show()
