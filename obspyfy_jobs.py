import sys
sys.path.append('/home/adrian/PhD/axisem3d_root/AxiSEM3D-master/cmb/output_handlers')
from obspyfy import obspyfy
import os

for dir in os.walk('simu_JOBS'):
    for subdir in dir[1]:
        obspyfy(os.path.abspath('simu_JOBS/' + subdir))
    break
