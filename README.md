
See important notes from Paul FitzGerald on best usage of XCIST: [XCIST Readme](/home/brandon.nelson/Dev/DLIR_Ped_Generalizability/XCAT/CatSim/examples/evaluation/ReadMe.txt)


## Install

1. Create an environement (optional)

```bash
conda create -n XCIST
conda activate XCIST
```

git clone https://github.com/xcist/main.git XCIST

pip install XCIST/

conda install tqdm



# TODO

1. use [DICOM_to_voxelized_phantom.py](/gpfs_projects/brandon.nelson/XCIST/phantoms-voxelized/DICOM_to_voxelized/DICOM_to_voxelized_phantom.py) to turn .raw ground truth XCAT phantom files to XCIST format
    - *export to DICOM first*
    - see concents of </gpfs_projects/brandon.nelson/XCIST/phantoms-voxelized/> for examples