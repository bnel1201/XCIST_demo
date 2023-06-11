# %%
# https://pydicom.github.io/pydicom/stable/auto_examples/input_output/plot_write_dicom.html#sphx-glr-auto-examples-input-output-plot-write-dicom-py/
# see for info on how to make new DICOMS from raw
import re
import numpy as np

raw_file = '/gpfs_projects/brandon.nelson/DLIR_Ped_Generalizability/geometric_phantom_studies/main/anthropomorphic/phantoms/adaptive_fov/male_infant_ref_atn_1.bin'
log_file = '/gpfs_projects/brandon.nelson/DLIR_Ped_Generalizability/geometric_phantom_studies/main/anthropomorphic/phantoms/adaptive_fov/male_infant_ref_log'

def get_pixel_width(log_filename, units='mm'):
    """
    reads XCAT log file and returns pixel width in specified `units`
    """
    with open(log_filename, 'r') as f:
        text = f.read()
    m = re.search('(pixel width =  )(\d.\d+)', text)
    pixel_width = float(m.group(2))
    if units == 'mm':
        pixel_width *= 10
    return pixel_width

def get_slice_width(log_filename, units='mm'):
    """
    reads XCAT log file and returns pixel width in specified `units`
    """
    with open(log_filename, 'r') as f:
        text = f.read()
    m = re.search('(slice width =  )(\d.\d+)', text)
    pixel_width = float(m.group(2))
    if units == 'mm':
        pixel_width *= 10
    return pixel_width

def get_array_size(log_filename):
    """
    reads XCAT log file and returns array size in pixels
    """
    with open(log_filename, 'r') as f:
        text = f.read()
    m = re.search('(array_size =\s+)(\d.\d+)', text)
    return int(m.group(2))

def get_water_lac(log_filename):
    """
    reads XCAT log file and returns array size in pixels
    """
    with open(log_filename, 'r') as f:
        text = f.read()
    m = re.search('(Body \(water\)\s+=\s+)(\d.\d+)', text)
    return float(m.group(2))



def mu2hu(mu, mu_water): return 1000*(mu - mu_water)/mu_water
# %% [markdown]
# need to grab water LAC from log file to convert to HU to save into DICOM for the DICOM_to_voxelized function
# then keep that water LAC to give to the DICOM_to_voxelized to convert *back* to mu
import catsim as xc
import matplotlib.pyplot as plt

data_mu = xc.rawread(raw_file, 2*[get_array_size(log_file)], 'float')
water_mu = get_water_lac(log_file)/10
pixel_width = get_pixel_width(log_file)
slice_width = get_slice_width(log_file)
data_hu = mu2hu(data_mu/pixel_width, water_mu).astype(np.int16) + 1024
# %%
# authors : Guillaume Lemaitre <g.lemaitre58@gmail.com>
# license : MIT

import datetime

import pydicom
from pydicom.dataset import FileDataset, FileMetaDataset
from pydicom.uid import UID
import numpy as np



# Create some temporary filenames
suffix = '.dcm'
filename = f'test_0{suffix}'

print("Setting file meta information...")
# Populate required values for file meta information
file_meta = FileMetaDataset()
file_meta.MediaStorageSOPClassUID = UID('1.2.840.10008.5.1.4.1.1.2')
file_meta.MediaStorageSOPInstanceUID = UID("1.2.3")
file_meta.ImplementationClassUID = UID("1.2.3.4")
file_meta.TransferSyntaxUID = '1.2.840.10008.1.2.1'

print("Setting dataset values...")
# Create the FileDataset instance (initially no data elements, but file_meta
# supplied)
ds = FileDataset(filename, {},
                 file_meta=file_meta, preamble=b"\0" * 128)

ds.PixelData = data_hu
ds.Rows, ds.Columns = data_hu.shape
ds.PixelSpacing = [pixel_width, pixel_width]
ds.SliceThickness = slice_width
ds.BitsAllocated = 16
ds.BitsStored = 12
ds.HighBit = 11
ds.PixelRepresentation = 0
ds.SmallestImagePixelValue = int(data_hu.min())
ds.LargestImagePixelValue = int(data_hu.max())
ds.SamplesPerPixel = 1
ds.PhotometricInterpretation = 'MONOCHROME2'
ds.WindowCenter = [40, -600]
ds.WindowWidth = [400, 1500]
ds.RescaleIntercept = '-1024.0'
ds.RescaleSlope = '1.0'
ds.RescaleType = 'HU'
# From: https://pydicom.github.io/pydicom/stable/auto_examples/image_processing/plot_downsize_image.html#sphx-glr-auto-examples-image-processing-plot-downsize-image-py
# Add the data elements -- not trying to set all required here. Check DICOM
# standard
ds.PatientName = "Test^Firstname"
ds.PatientID = "123456"

# Set the transfer syntax
ds.is_little_endian = True
ds.is_implicit_VR = True

# Set creation date/time
dt = datetime.datetime.now()
ds.ContentDate = dt.strftime('%Y%m%d')
timeStr = dt.strftime('%H%M%S.%f')  # long format with micro seconds
ds.ContentTime = timeStr

print("Writing test file", filename)
ds.save_as(filename)
print("File saved.")

# reopen the data just for checking
print('Load file {} ...'.format(filename))
ds = pydicom.dcmread(filename)
print(ds)

# %%
