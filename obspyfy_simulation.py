import sys
sys.path.append('/home/users/scro4564/PhD/AxiSEM3D_Data_Handler')
from axisem3d_output.obspyfy import obspyfy

obspyfy('/disks/data/PhD/CMB/simu1D/FORWARD_CHECK', output_type='stations', 
        stations_paths='/disks/data/PhD/CMB/stations/STA_10DEG_CROSS.txt')                 
