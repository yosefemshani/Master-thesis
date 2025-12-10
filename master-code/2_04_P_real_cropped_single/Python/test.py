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

ct = dicom.dcmread("../DICOM_fake/Fake0034.dcm")

print(ct.ImagePositionPatient)

