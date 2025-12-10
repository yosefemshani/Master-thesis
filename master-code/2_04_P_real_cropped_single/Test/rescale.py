import numpy as np
import pydicom
from scipy.interpolate import interpn

# Function to get grid ticks for 2D systems
def getGridTicks(fVP:list, voxelSize:np.array, imageDataShape:np.array):
    '''
    Return 2D ticks (i.e. coordinate pattern along x and y axes) of grid with size "imageDataShape", 
    taking the voxelSize and firstVoxelPosition of certain structure into account

    Parameters:
    -----------
    fVP (list): first voxel position of structure for which grid ticks should be calculated
    voxelSize (np.array): voxel size of structure for which grid ticks should be calculated
    imageDataShape (np.array): shape of required coordinate grid
    '''

    # calculate ticks of existing grid with phantom voxel size 
    x = fVP[0] + np.linspace(0, imageDataShape[0]-1, imageDataShape[0]) * voxelSize[0]
    y = fVP[1] + np.linspace(0, imageDataShape[1]-1, imageDataShape[1]) * voxelSize[1]
    return (x, y)

# Function to get coordinate meshgrid for 2D systems
def getCoordinateMeshgrid(minCoords:np.array, maxCoords:np.array, voxelSize:np.array):
    '''
    Return full coordinate list of 2D regular grid between minCoords and maxCoords with voxelSize as spacing

    Parameters:
    -----------
    minCoords (np.array): minimal 2D coordinates of new grid
    maxCoords (np.array): maximal 2D coordinates of new grid
    voxelSize (np.array): voxelSize in both dimensions

    Returns:
    ---------
    coordinateList (np.array): full list of all coordinate points of new grid
    gridshape (np.array): size of new grid
    '''

    # calculate list of coordinates in new grid with patient voxel size
    coord_x = np.arange(minCoords[0], maxCoords[0], voxelSize[0])
    coord_y = np.arange(minCoords[1], maxCoords[1], voxelSize[1])

    gridshape = np.array([coord_x.shape[0], coord_y.shape[0]])
    xarr, yarr = np.meshgrid(coord_x, coord_y, indexing='ij')
    xarr = np.ravel(xarr)
    yarr = np.ravel(yarr)
    coordinateList = np.array([xarr, yarr]).T

    return coordinateList, gridshape


# Read DICOM file
ct = pydicom.dcmread("IMG0035_cropped.dcm")

# Define original pixel spacing
original_pixel_spacing = ct.PixelSpacing[0]

# Define desired pixel spacing
new_pixel_spacing = 0.001

# Calculate scaling factor
scaling = original_pixel_spacing / new_pixel_spacing

# Upscale pixel spacing
ct.PixelSpacing = [new_pixel_spacing, new_pixel_spacing]

# Extract rows and columns
rows = ct.Rows
columns = ct.Columns

# Calculate upscaled rows and columns
new_rows = round(rows * scaling)
new_columns = round(columns * scaling)

# Generate new grid for upscaling
grid_rows = np.linspace(0, rows - 1, new_rows)
grid_columns = np.linspace(0, columns - 1, new_columns)

# Extract pixel data
pixel_array = ct.pixel_array


# Define the grid for interpolation
existCoordinateGrid1 = getGridTicks([0, 0], [1, 1], pixel_array.shape)

# Define the coordinates for interpolation
coordList1_split, gridshape1_split = getCoordinateMeshgrid(
    [0, 0],
    [new_rows - 1, new_columns - 1],
    [1, 1]
)

# Perform interpolation
interpSection = interpn(
    existCoordinateGrid1[:2],  # Use only the first two dimensions of the grid
    pixel_array,
    coordList1_split[:, :2],  # Use only the first two dimensions of the coordinate list
    method='linear',
    bounds_error=False,
    fill_value=-1024
)

# Reshape the interpolated section to match the grid shape
interpSection = np.reshape(interpSection, (new_rows, new_columns))

# Update pixel data
ct.Rows = new_rows
ct.Columns = new_columns
ct.PixelData = interpSection.astype(np.uint16).tobytes()

# Save the DICOM file with upscaled pixel data
ct.save_as("upscaled_dicom_file.dcm")

