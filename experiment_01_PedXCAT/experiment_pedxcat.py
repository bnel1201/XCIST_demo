# %%
import catsim as xc
import sys
from pathlib import Path

# wdir = Path(__file__).parent

sys.path.extend(['/home/brandon.nelson/Dev/DLIR_Ped_Generalizability/XCAT/CatSim/examples/evaluation/pyfiles',
                 '/home/brandon.nelson/Dev/DLIR_Ped_Generalizability/XCAT/CatSim'])
# ^^^ From Paul FitzGerald Readme recommendation
# %%
# from my_commonTools import *
# %%
# userPath = getUserPath()
userPath = '/home/brandon.nelson/Dev/DLIR_Ped_Generalizability/xcist'
# xc.pyfiles.CommonTools.my_path.add_search_path(userPath)
# xc.pyfiles.CommonTools.my_path.add_search_path('/home/brandon.nelson/Dev/XCIST_demo/example_body')

ct = xc.CatSim('cfg/Phantom_male_infant_ref',
               'cfg/Physics_Default',
               'cfg/Protocol_Default',
               'cfg/Recon_Default',
               'cfg/Scanner_Default')

ct.resultsName = 'male_infant_test'
ct.cfg.experimentDirectory = 'male_infant_test'

ct.cfg.do_Sim = False
ct.cfg.do_Recon = True
ct.cfg.recon.saveImagePictureFiles = True
# ct.protocol.viewsPerRotation = 2
# ct.protocol.viewCount = ct.protocol.viewsPerRotation
# ct.protocol.stopViewId = ct.protocol.viewCount-1

# ct.physics.enableQuantumNoise = 1
# ct.physics.enableElectronicNoise = 1

# ct.physics.colSampleCount = 2
# ct.physics.rowSampleCount = 2
# ct.physics.srcXSampleCount = 2
# ct.physics.srcYSampleCount = 2
# ct.physics.viewSampleCount = 1

# %%
ct.run_all()
# %%
import numpy as np
import matplotlib.pyplot as plt

prep = xc.rawread(ct.resultsName+'.prep', [ct.protocol.viewCount, ct.scanner.detectorRowCount, ct.scanner.detectorColCount], 'float')
prep = prep[-1, :, :]
plt.plot(prep[1, :])
plt.show()
# %%
