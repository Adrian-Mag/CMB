import sys
sys.path.append('/home/users/scro4564/PhD/AxiSEM3D_Data_Handler')
from obspyfy import obspyfy

obspyfy('/home/users/scro4564/PhD/CMB/simu1D_element/FORWARD', output_type='elements', 
        stations_paths='/home/users/scro4564/PhD/CMB/stations/STA_10DEG_CROSS.txt')                 