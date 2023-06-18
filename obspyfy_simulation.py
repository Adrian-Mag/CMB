import sys
sys.path.append('/home/users/scro4564/PhD/AxiSEM3D_Data_Handler')
from AxiSEM3D_Data_Handler.obspyfy import obspyfy

obspyfy('/disks/data/PhD/CMB/simu3D_CMB/NORMAL_FAULT_100KM_CMB_50_0_30_30', output_type='stations', 
        stations_paths='/disks/data/PhD/CMB/stations/STA_10DEG_CROSS.txt')                 
