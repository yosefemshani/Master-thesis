import pydicom
from pydicom.dataset import Dataset, FileDataset
from pydicom.uid import generate_uid

# File paths
rs_file = "../DICOM/RS1.2.752.243.1.1.20170503093409104.7100.87740_transformed.dcm"
original_ct_file = "../DICOM/IMG0035_cropped.dcm"
fake_ct_file = "../DICOM/IMG0036_cropped.dcm"
output_rs_file = "../DICOM_DoubleSlice/RSfake.dcm"

# Load the original and fake CT slices
original_ct = pydicom.dcmread(original_ct_file)
fake_ct = pydicom.dcmread(fake_ct_file)

# Load the RS structure file
rs_dataset = pydicom.dcmread(rs_file, force=True)

# Find the contour sequence for the fake CT slice
for roi_contour_seq in rs_dataset.ROIContourSequence:
    for contour_seq in roi_contour_seq.ContourSequence:
        ref_sop_instance_uid = contour_seq.ContourImageSequence[0].ReferencedSOPInstanceUID
        if ref_sop_instance_uid == fake_ct.SOPInstanceUID:
            # Duplicate the contour sequence
            new_contour_seq = Dataset()
            new_contour_seq.ContourGeometricType = contour_seq.ContourGeometricType
            new_contour_seq.NumberOfContourPoints = contour_seq.NumberOfContourPoints
            new_contour_seq.ContourNumber = contour_seq.ContourNumber
            new_contour_seq.ContourData = contour_seq.ContourData
            
            # Create a new ContourImageSequence
            new_contour_image_seq = Dataset()
            new_contour_image_seq.ReferencedSOPClassUID = original_ct.SOPClassUID
            new_contour_image_seq.ReferencedSOPInstanceUID = original_ct.SOPInstanceUID
            new_contour_seq.ContourImageSequence = [new_contour_image_seq]
            
            # Append the new contour sequence to the ROIContourSequence
            roi_contour_seq.ContourSequence.append(new_contour_seq)
            break

# Save the updated RS structure file
rs_dataset.save_as(output_rs_file)

print(f"Contours successfully copied from {fake_ct_file} to {original_ct_file} in the RS file, saved as {output_rs_file}.")
