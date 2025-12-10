import numpy as np
import pydicom
import matplotlib.pyplot as plt
from scipy.interpolate import RegularGridInterpolator, interpn

import yaml, os, glob
import time

# read yaml file
with open('params.yaml', 'r') as paramfile:
    params = yaml.safe_load(paramfile)

# read DICOM-related metadata from YAML
dicomDir = params['DICOMdir']
dicomPrefix = params['DICOMprefix']
patientName = params['PatientName']
patientID = params['PatientID']
studyInstanceUID = params['StudyInstanceUID']
studyID = params['StudyID']
studyDate = params['StudyDate']
studyDescription = params['StudyDescription']
seriesDescription = params['SeriesDescription']
referringPhysician = params['ReferringPhysician']

# amount of splits to be performed on lower and upper phantom part for interpolation
splitCountLower = params['splitsLower']
splitCountUpper = params['splitsUpper']

# threshold for air voxels
HUAirThresh = params['HUAirThresh']


def dicom_fusion(Phantom:object, Patient:object, finalTrans:np.array, clearDICOMdir:bool):
    '''
    Write joint DICOM dataset from patient and phantom image data, 
    taking finalTrans into consideration as neccessary deviation in voxel positions.
    Bring all arrays to same x/y-dimension and save them z-slicewise as DICOM file with valid metadata.

    Parameters:
    -----------
    Phantom (obj): Phantom object, storing all phantom-related data
    Patient (obj): Patient object, storing all patient-related data
    finalTrans (np.array): deviation in voxel positions between patient and phantom
    '''

    # consider finalTrans directly to Phantom object
    Phantom.initial_firstVoxelPosition = Phantom.firstVoxelPosition
    Phantom.initial_skeletonPositions = Phantom.skeletonPositions
    Phantom.firstVoxelPosition -= finalTrans
    Phantom.skeletonPositions = (Phantom.skeletonPositions.T - finalTrans).T
    
    # add last voxel position coordinates to structure objects
    for struct in [Phantom, Patient]:
        struct.lastVoxelPosition = struct.firstVoxelPosition + (np.array(struct.data.shape)-1)*struct.voxelSize

    # cut out phantom voxel points in range of patient CT by z coordinate
    # define min/max z coordinates from Patient
    zMinPatient = Patient.firstVoxelPosition[2]
    zMaxPatient = Patient.lastVoxelPosition[2]

    # cut out z slices in phantom data array
    # split phantom in part below CT and part above CT
    zMinPatientIndexInPhantomGrid = np.rint((zMinPatient - Phantom.firstVoxelPosition[2])/Phantom.voxelSize[2]).astype(int)
    zMaxPatientIndexInPhantomGrid = np.rint((zMaxPatient - Phantom.firstVoxelPosition[2])/Phantom.voxelSize[2]).astype(int)
    phantomPart1 = Phantom.HU[:,:,:zMinPatientIndexInPhantomGrid+1] # keep extra slice as puffer to prevent interpolated phantom
    phantomPart2 = Phantom.HU[:,:,zMaxPatientIndexInPhantomGrid-1:] # from having an empty slice between phantom and CT
    print(f'gridshapes of phantom parts: {phantomPart1.shape}, {phantomPart2.shape}')

    # interpolate lower part of phantom to patient grid
    # calculate ticks of existing grid with phantom voxel size 
    existCoordinateGrid1 = getGridTicks(Phantom.firstVoxelPosition, Phantom.voxelSize, phantomPart1.shape)
    # modify z coordinate of firstVoxelPosition for upper part of phantom
    modZ_fVP = Phantom.firstVoxelPosition[2] + (zMaxPatientIndexInPhantomGrid-1) * Phantom.voxelSize[2]
    existCoordinateGrid2 = getGridTicks([*Phantom.firstVoxelPosition[0:2],modZ_fVP], Phantom.voxelSize, phantomPart2.shape)
    
    # find range of voxel coordinates
    # in x/y: use one unit of voxelSize as extension to get enough values back
    minVoxelCoords = np.min([Patient.firstVoxelPosition, Phantom.firstVoxelPosition], axis=0)
    maxVoxelCoords = np.max([Patient.lastVoxelPosition,  Phantom.lastVoxelPosition], axis=0) + Patient.voxelSize/2
    print(f'min/max voxel coords: {minVoxelCoords}, {maxVoxelCoords}')

    ### perform interpolation from phantom data to patient voxel grid
    # prepare min/max coordinate lists with patient voxel size and modified z range: 
    # split lower/upper part more times to make program more stable
    # define approximate borders for sub-interpolations as zSplitList
    start_interpolate = time.time()
    # LOWER part
    if phantomPart1.shape[2] > 0:
        phantom1Exist = True
        # for lower part of phantom use global minCoords and zMinPatient as maximum z value, invert order
        zSplitList = np.linspace(zMinPatient-Patient.voxelSize[2],minVoxelCoords[2],splitCountLower+1)
        for i in range(splitCountLower): # splitCountLower from YAML
            # set min/max voxel coordinates for interpolations, use elements of zSplitList as z
            minVoxelCoords1_split = np.array([*minVoxelCoords[0:2], zSplitList[i+1]])
            maxVoxelCoords1_split = np.array([*maxVoxelCoords[0:2], zSplitList[i]])
            coordList1_split, gridshape1_split = getCoordinateMeshgrid(minVoxelCoords1_split, maxVoxelCoords1_split, 
                                                                        Patient.voxelSize, zDirection='decrease')
            # update next maximal z coordinate to minimum from current section, decreased by one voxelSize unit
            zSplitList[i+1] = np.min(coordList1_split, axis=0)[2] - Patient.voxelSize[2]
            #print(coordList1_split.shape, np.min(coordList1_split, axis=0), np.max(coordList1_split, axis=0))
            interpSection = interpn(existCoordinateGrid1, phantomPart1, coordList1_split, 
                                    method='linear', bounds_error=False, fill_value=-1024)
            interpSection = np.reshape(interpSection, gridshape1_split)
            # initialize resp. concatenate resulting array
            if i == 0:
                interpPhantom1 = interpSection
            else:
                interpPhantom1 = np.concatenate((interpPhantom1, interpSection), axis=2)
            print(interpPhantom1.shape)
            print(zSplitList)
        # flip lower part because of downwards interpolation, resulting in saving lower slices at the top of interpPhantom1
        interpPhantom1 = np.flip(interpPhantom1, axis=2)
    else:
        phantom1Exist = False
    
    # UPPER part
    if phantomPart2.shape[2] > 0:
        phantom2Exist = True
        # for upper part of phantom use global maxCoords and zMaxPatient as minimum z value
        zSplitList = np.linspace(zMaxPatient+Patient.voxelSize[2],maxVoxelCoords[2],splitCountUpper+1)
        print(zSplitList)
        for i in range(splitCountUpper): # splitCountUpper from YAML
            # set min/max voxel coordinates for interpolations, use elements of zSplitList as z
            minVoxelCoords2_split = np.array([*minVoxelCoords[0:2], zSplitList[i]])
            maxVoxelCoords2_split = np.array([*maxVoxelCoords[0:2], zSplitList[i+1]])
            coordList2_split, gridshape2_split = getCoordinateMeshgrid(minVoxelCoords2_split, maxVoxelCoords2_split, 
                                                                        Patient.voxelSize, zDirection='increase')
            # update next minimal z coordinate to maximum from current section, increased by one voxelSize unit
            zSplitList[i+1] = np.max(coordList2_split, axis=0)[2] + Patient.voxelSize[2]
            interpSection = interpn(existCoordinateGrid2, phantomPart2, coordList2_split, 
                                    method='linear', bounds_error=False, fill_value=-1024)
            interpSection = np.reshape(interpSection, gridshape2_split)
            # initialize resp. concatenate resulting array
            if i == 0:
                interpPhantom2 = interpSection
            else:
                interpPhantom2 = np.concatenate((interpPhantom2, interpSection), axis=2)
            print(interpPhantom2.shape)
            print(zSplitList)
    else:
        phantom2Exist = False
    print('interpolating time: ', time.time() - start_interpolate)
    print(f'shapes interpolated phantom parts:')
    print(interpPhantom1.shape if phantom1Exist else None, interpPhantom2.shape if phantom2Exist else None)

    # calculate new z slice amounts
    CTSlices = Patient.data.shape[2]
    intPhantom1Slices = interpPhantom1.shape[2]
    intPhantom2Offset = intPhantom1Slices + CTSlices

    # pad CT to phantom gridsize
    if phantom1Exist:
        padCT = np.asarray(interpPhantom1.shape)[0:2] - Patient.data.shape[0:2]
    elif phantom2Exist:
        padCT = np.asarray(interpPhantom2.shape)[0:2] - Patient.data.shape[0:2]
    else:
        raise ValueError('phantom interpolation to patient grid has failed, both parts are empty')
    padCT[padCT < 0] = 0
    print(f'required padding x/y: {padCT}')
    paddedPatientData = np.pad(Patient.data, 
        ((int(np.ceil(padCT[0]/2)),int(np.floor(padCT[0]/2))), 
         (int(np.ceil(padCT[1]/2)),int(np.floor(padCT[1]/2))),(0,0)), 'constant', constant_values=-1024)
    print(f'shape padded CT: {paddedPatientData.shape}')
    # number of voxels being added at lower coordinate borders
    addedLowerVoxels = np.array([np.ceil(padCT[0]/2), np.ceil(padCT[1]/2), 0]).astype(int)
    print('prelim addedLowerVoxels', addedLowerVoxels)

    if ((not phantom1Exist or paddedPatientData.shape[0:2] == interpPhantom1.shape[0:2])
        and (not phantom2Exist or paddedPatientData.shape[0:2] == interpPhantom2.shape[0:2])):
        print('Reshaping successful!')
    else:
        print('WARNING! Problem in padding to same size.')
    
    # find minimal and maximal non-air voxel indices
    start_minmax = time.time()
    if phantom1Exist and phantom2Exist: arrayList = [interpPhantom1, interpPhantom2, paddedPatientData]
    elif phantom1Exist: arrayList = [interpPhantom1, paddedPatientData]
    elif phantom2Exist: arrayList = [interpPhantom2, paddedPatientData]
    minIdx, maxIdx = [], []
    for dataArray in arrayList:
        minIdx.append(np.min(np.where(dataArray > HUAirThresh), axis=1))
        maxIdx.append(np.max(np.where(dataArray > HUAirThresh), axis=1))
    minGlobalIdx = np.min(np.asarray(minIdx), axis=0)
    maxGlobalIdx = np.max(np.asarray(maxIdx), axis=0)
    #minGlobalIdx = np.array([0,3,0])
    #maxGlobalIdx = np.array([511,702,457])
    print(f'min/max global non-air indices: {minGlobalIdx}, {maxGlobalIdx}')
    print(f'Time min/max global non-air indices: {time.time() - start_minmax} s')
    # subtract voxels that are cut off at minimum side to get number of voxels that are finally added
    addedLowerVoxels -= minGlobalIdx
    print('final addedLowerVoxels', addedLowerVoxels)

    # cut down all three parts of phantom and CT to HU-relevant data by previously calculated non-air indices
    for i in range(len(arrayList)):
        arrayList[i] = arrayList[i][minGlobalIdx[0]:maxGlobalIdx[0]+1, minGlobalIdx[1]:maxGlobalIdx[1]+1,:]
    if phantom1Exist and phantom2Exist: 
        [interpPhantom1, interpPhantom2, paddedPatientData] = arrayList
        totalSlices = CTSlices + interpPhantom1.shape[2] + interpPhantom2.shape[2] 
        print('cut shapes phantom1/phantom2/CT: ', interpPhantom1.shape, interpPhantom2.shape, paddedPatientData.shape)
    elif phantom1Exist: 
        [interpPhantom1, paddedPatientData] = arrayList
        totalSlices = CTSlices + interpPhantom1.shape[2]
        print('cut shapes phantom1/CT: ', interpPhantom1.shape, paddedPatientData.shape)
    elif phantom2Exist: 
        [interpPhantom2, paddedPatientData] = arrayList
        totalSlices = CTSlices + interpPhantom2.shape[2]
        print('cut shapes phantom2/CT: ', interpPhantom2.shape, paddedPatientData.shape)

    # generate new SeriesInstanceUID based on existing to ensure that series is loaded independent from original CT
    splitSeriesInstanceUID = Patient.dicomInfo[0][0x0020,0x000e].value.split('.')
    splitSeriesInstanceUID[-1] = str(int(splitSeriesInstanceUID[-1]) + 1)
    newSeriesInstanceUID = '.'.join(splitSeriesInstanceUID)
    
    # save original slices with modified metadata, especially UIDs
    # check if aimed directory exists
    if not os.path.isdir(dicomDir):
        os.makedirs(dicomDir, exist_ok=True)
    # delete existing slice files to avoid mixed datasets with incompatible metadata
    elif clearDICOMdir:
        DICOMlist = glob.glob(os.path.join(dicomDir, '*.dcm'))
        i = 0
        for file in DICOMlist[totalSlices:]:
            os.remove(file)
            i += 1
        print('removed slices: ', i)


    start_saving = time.time()
    for i in range(CTSlices):
        sliceIndex = i + intPhantom1Slices
        Patient.dicomInfo[i].InstanceNumber = sliceIndex
        #Patient.dicomInfo[i][0x0010,0x0010].value = patientName # patient name
        #Patient.dicomInfo[i][0x0010,0x0020].value = patientID # patient ID
        #splitSOPInstanceUID = Patient.dicomInfo[i][0x0008,0x0018].value.split('.')
        #splitSOPInstanceUID[-1] = str(int(splitSOPInstanceUID[-1]) + 10**5) # high number to ensure slice-individual new UID
        #Patient.dicomInfo[i][0x0008,0x0018].value = '.'.join(splitSOPInstanceUID) # SOP instance UID
        #Patient.dicomInfo[i][0x0020,0x000D].value = studyInstanceUID # study instance UID
        #Patient.dicomInfo[i].StudyID              = studyID # study ID
        Patient.dicomInfo[i].StudyDate            = studyDate # study date
        Patient.dicomInfo[i].StudyDescription     = studyDescription # study description
        #Patient.dicomInfo[i][0x0020,0x000e].value = newSeriesInstanceUID
        Patient.dicomInfo[i].SeriesDescription    = seriesDescription # series description
        Patient.dicomInfo[i].ReferringPhysicianName = referringPhysician # referring physician name
        Patient.dicomInfo[i][0x0028,0x0010].value = paddedPatientData.shape[0]
        Patient.dicomInfo[i][0x0028,0x0011].value = paddedPatientData.shape[1]
        # consider added or cut-off voxels in ImagePositionPatient to keep origin at constant location
        # swap order of x/y because medical x corresponds to technical y and vice versa
        Patient.dicomInfo[i].ImagePositionPatient[0:2] -= np.flip(addedLowerVoxels[0:2] * Patient.voxelSize[0:2])
        pixel_array = np.array(
            (paddedPatientData[:,:,i] - Patient.dicomInfo[i].RescaleIntercept)/Patient.dicomInfo[i].RescaleSlope, dtype=np.uint16)
        Patient.dicomInfo[i].PixelData = pixel_array.tostring()
        Patient.dicomInfo[i].save_as(f'{dicomDir}{dicomPrefix}{str(sliceIndex).zfill(4)}.dcm')
    refSOPInstanceUID = Patient.dicomInfo[0].SOPInstanceUID
    refSliceLocation = Patient.dicomInfo[0].SliceLocation

    # add slices below CT data
    if phantom1Exist:
        saveSlices(interpPhantom1, 0, -intPhantom1Slices, 
            refSliceLocation, refSOPInstanceUID, Patient.dicomInfo[0], Patient.voxelSize)

    # add slices above CT data
    if phantom2Exist:
        saveSlices(interpPhantom2, intPhantom2Offset, CTSlices, 
            refSliceLocation, refSOPInstanceUID, Patient.dicomInfo[0], Patient.voxelSize)
    print(f'Time DICOM saving: {time.time() - start_saving} s')


def saveSlices(phantomData:np.array, indexOffset:int, coordinateOffset:int, 
    refSliceLocation, refSOPInstanceUID, sampleMetadata, voxelSize):
    '''
    Short function to save interpolated phantom slices as DICOM files

    Parameters:
    -----------
    phantomData (np.array): array-like image data to be saved
    indexOffset (int): offset for all indexing parameters (index non-zero only for upper part)
    coordinateOffset (int): offset for all coordinate-related parameters that should take first CT slice as origin data
    refSliceLocation: reference z position of first CT slice
    refSOPInstanceUID: reference SOPInstanceUID for generation of new unique SOPInstanceUIDs
    sampleMetadata: sample DICOM metadata from original CT dataset as basis for artificial phantom metadata
    voxelSize: new applied voxel size within interpolation (here: patient voxel size)
    '''
    
    sampleMetadata[0x0028,0x0010].value = phantomData.shape[0]
    sampleMetadata[0x0028,0x0011].value = phantomData.shape[1]
    for i in range(phantomData.shape[2]):
        splitSOPInstanceUID = refSOPInstanceUID.split('.')
        splitSOPInstanceUID[-1] = str(int(splitSOPInstanceUID[-1]) + 10**7 + (indexOffset+i)*2)
        sampleMetadata[0x0020,0x1041].value = refSliceLocation + (coordinateOffset+i)*voxelSize[2] # slice location
        sampleMetadata[0x0008,0x0018].value = '.'.join(splitSOPInstanceUID) # SOP instance UID
        sampleMetadata[0x0020,0x0013].value = indexOffset + i # instance number
        sampleMetadata.ImagePositionPatient[2] = sampleMetadata[0x0020,0x1041].value
        pixel_array = np.array(
            (phantomData[:,:,i] - sampleMetadata.RescaleIntercept)/sampleMetadata.RescaleSlope, dtype=np.uint16)
        sampleMetadata[0x7fe0,0x0010].value = pixel_array.tostring()
        sampleMetadata.save_as(f'{dicomDir}{dicomPrefix}{str(indexOffset+i).zfill(4)}.dcm')


def getGridTicks(fVP:list, voxelSize:np.array, imageDataShape:np.array):
    '''
    Return 3D ticks (i.e. coordinate pattern along every axis) of grid with size "imageDataShape", 
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
    z = fVP[2] + np.linspace(0, imageDataShape[2]-1, imageDataShape[2]) * voxelSize[2]
    return (x,y,z)


def getCoordinateMeshgrid(minCoords:np.array, maxCoords:np.array, voxelSize:np.array, zDirection:str='increase'):
    '''
    Return full coordinate list of 3D regular grid between minCoords and maxCoords with voxelSize as spacing

    Parameters:
    -----------
    minCoords (np.array): minimal 3D coordinates of new grid
    maxCoords (np.array): maximal 3D coordinates of new grid
    voxelSize (np.array): voxelSize in all 3 dimensions
    zDirection (str): starting point for z coordinate in linspace, either minCoords or maxCoords

    Returns:
    ---------
    coordinateList (np.array): full list of all coordinate points of new grid
    gridshape (np.array): size of new grid
    '''

    # calculate list of coordinates in new grid with patient voxel size
    coord_x = np.arange(minCoords[0], maxCoords[0], voxelSize[0])
    coord_y = np.arange(minCoords[1], maxCoords[1], voxelSize[1])
    if zDirection == 'increase': # take minCoords in z as starting coordinate for meshgrid
        coord_z = np.arange(minCoords[2], maxCoords[2], voxelSize[2])
    elif zDirection == 'decrease': # take maxCoords in z as starting coordinate for meshgrid
        coord_z = -np.arange(-maxCoords[2], -minCoords[2], voxelSize[2])
    else:
        raise ValueError("zDirection type could not be resolved, choose 'increase' or 'decrease'")

    gridshape = np.array([coord_x.shape[0], coord_y.shape[0], coord_z.shape[0]])
    xarr, yarr, zarr = np.meshgrid(coord_x, coord_y, coord_z, indexing='ij')
    xarr = np.ravel(xarr)
    yarr = np.ravel(yarr)
    zarr = np.ravel(zarr)
    coordinateList = np.array([xarr, yarr, zarr]).T

    return coordinateList, gridshape