# Coded version of DICOM file '/gpfs_projects/prabhat.kc/lowdosect/clin_LDCT_lib/flull_3mm_sharp_sorted/000153.IMA'
# Produced by pydicom codify utility script
# pydicom /gpfs_projects/prabhat.kc/lowdosect/clin_LDCT_lib/flull_3mm_sharp_sorted/000153.IMA /home/brandon.nelson/Dev/XCIST_demo/pydicom_code.py
import pydicom
from pydicom.dataset import Dataset, FileMetaDataset
from pydicom.sequence import Sequence

# File meta info data elements
file_meta = FileMetaDataset()
file_meta.FileMetaInformationGroupLength = 192
file_meta.FileMetaInformationVersion = b'\x00\x01'
file_meta.MediaStorageSOPClassUID = '1.2.840.10008.5.1.4.1.1.2'
file_meta.MediaStorageSOPInstanceUID = '1.3.12.2.1107.5.1.4.64291.30000016012200111078900002895'
file_meta.TransferSyntaxUID = '1.2.840.10008.1.2.1'
file_meta.ImplementationClassUID = '1.3.12.2.1107.5.1.4'
file_meta.ImplementationVersionName = 'SIEMENS_S7VA44A'

# Main data elements
ds = Dataset()
ds.SpecificCharacterSet = 'ISO_IR 100'
ds.ImageType = ['ORIGINAL', 'PRIMARY', 'AXIAL', 'CT_SOM5 SPI']
ds.SOPClassUID = '1.2.840.10008.5.1.4.1.1.2'
ds.SOPInstanceUID = '1.3.12.2.1107.5.1.4.64291.30000016012200111078900002895'
ds.StudyDate = '20130111'
ds.SeriesDate = '20160120'
ds.AcquisitionDate = '20130111'
ds.ContentDate = '20130111'
ds.AcquisitionDateTime = '20130111101259.867000'
ds.StudyTime = '100953.500000'
ds.SeriesTime = '181856.851000'
ds.AcquisitionTime = '101259.867000'
ds.ContentTime = '101259.867000'
ds.AccessionNumber = ''
ds.Modality = 'CT'
ds.Manufacturer = 'SIEMENS'
ds.ReferringPhysicianName = ''
ds.ManufacturerModelName = 'SOMATOM Definition AS+'
ds.DerivationDescription = 'Force Anonymity'

# Source Image Sequence
source_image_sequence = Sequence()
ds.SourceImageSequence = source_image_sequence

# Source Image Sequence: Source Image 1
source_image1 = Dataset()
source_image1.ReferencedSOPClassUID = '1.3.12.2.1107.5.9.1'
source_image1.ReferencedSOPInstanceUID = '1.3.12.2.1107.5.1.4.64291.30000013011114102204600000047'
source_image_sequence.append(source_image1)

ds.IrradiationEventUID = '1.3.12.2.1107.5.1.4.64291.30000013011114102204600000047'
ds.PatientName = 'L096_FD_3_SHARP'
ds.PatientID = 'Anonymous'
ds.PatientBirthDate = ''
ds.PatientSex = 'O'
ds.PregnancyStatus = 4
ds.BodyPartExamined = 'CHEST'
ds.SliceThickness = '3.0'
ds.KVP = '120.0'
ds.DataCollectionDiameter = '500.0'
ds.SoftwareVersions = 'syngo CT 2011A'
ds.ReconstructionDiameter = '380.0'
ds.DistanceSourceToDetector = '1085.6'
ds.DistanceSourceToPatient = '595.0'
ds.GantryDetectorTilt = '0.0'
ds.TableHeight = '153.5'
ds.RotationDirection = 'CW'
ds.ExposureTime = '500'
ds.XRayTubeCurrent = '325'
ds.Exposure = '232'
ds.FilterType = 'FLAT'
ds.GeneratorPower = '52'
ds.FocalSpots = '1.2'
ds.DateOfLastCalibration = '20130111'
ds.TimeOfLastCalibration = ''
ds.ConvolutionKernel = 'D45f'
ds.PatientPosition = 'FFS'
ds.SingleCollimationWidth = 0.6
ds.TotalCollimationWidth = 38.4
ds.TableSpeed = 53.6
ds.TableFeedPerRotation = 26.8
ds.SpiralPitchFactor = 0.7
ds.DataCollectionCenterPatient = [0.37109375, -153.12890625, -375.0]
ds.ReconstructionTargetCenterPatient = [4.37109375, -153.12890625, -375.0]
ds.ExposureModulationType = 'XYZ_EC'
ds.EstimatedDoseSaving = 45.625
ds.CTDIvol = 15.68198149253731

# CTDI Phantom Type Code Sequence
ctdi_phantom_type_code_sequence = Sequence()
ds.CTDIPhantomTypeCodeSequence = ctdi_phantom_type_code_sequence

# CTDI Phantom Type Code Sequence: CTDI Phantom Type Code 1
ctdi_phantom_type_code1 = Dataset()
ctdi_phantom_type_code1.CodeValue = '113691'
ctdi_phantom_type_code1.CodingSchemeDesignator = 'DCM'
ctdi_phantom_type_code1.CodeMeaning = 'IEC Body Dosimetry Phantom'
ctdi_phantom_type_code_sequence.append(ctdi_phantom_type_code1)

ds.CalciumScoringMassFactorDevice = [0.7900000214576721, 0.8330000042915344, 0.871999979019165]
ds.StudyInstanceUID = '1.3.12.2.1107.5.1.4.64291.30000016012200111078900002740'
ds.SeriesInstanceUID = '1.3.12.2.1107.5.1.4.64291.30000016012200111078900002741'
ds.StudyID = ''
ds.SeriesNumber = '5'
ds.AcquisitionNumber = '2'
ds.InstanceNumber = '154'
ds.ImagePositionPatient = [-185.62890625, -343.12890625, -375]
ds.ImageOrientationPatient = [1, 0, 0, 0, 1, 0]
ds.FrameOfReferenceUID = '1.3.12.2.1107.5.1.4.64291.30000016012200111078900002739'
ds.PositionReferenceIndicator = ''
ds.SliceLocation = '375.0'
ds.SamplesPerPixel = 1
ds.PhotometricInterpretation = 'MONOCHROME2'
ds.Rows = 512
ds.Columns = 512
ds.PixelSpacing = [0.7421875, 0.7421875]
ds.BitsAllocated = 16
ds.BitsStored = 12
ds.HighBit = 11
ds.PixelRepresentation = 0
ds.SmallestImagePixelValue = 0
ds.LargestImagePixelValue = 2310
ds.WindowCenter = [40, -600]
ds.WindowWidth = [400, 1500]
ds.RescaleIntercept = '-1024.0'
ds.RescaleSlope = '1.0'
ds.RescaleType = 'HU'
ds.WindowCenterWidthExplanation = ['WINDOW1', 'WINDOW2']
ds.RequestingPhysician = 'BLOCK^MATTHEW^S'
ds.RequestedProcedureDescription = 'CT CHEST w + 3D Depend WS'
ds.PixelData = # XXX Array of 524288 bytes excluded

ds.file_meta = file_meta
ds.is_implicit_VR = False
ds.is_little_endian = True
ds.save_as(r'000153_from_codify.dcm', write_like_original=False)