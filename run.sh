export XCIST_UserPath=/home/brandon.nelson/Dev/DLIR_Ped_Generalizability/xcist

python raw_to_dicom.py

python DICOM_to_voxelized_phantom.py DICOM_to_voxelized_example_body.cfg 

python experiment_01_PedXCAT/experiment_pedxcat.py