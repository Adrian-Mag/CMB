import sys
sys.path.append('/home/adrian/PhD/AxiSEM3D/Output_Handlers')
from AxiSEM3D_Data_Handler.obspyfy import obspyfy

obspyfy('/home/adrian/PhD/AxiSEM3D/CMB/simu1D/OBLIQUE_FAULT_100KM', output_type='stations', 
        stations_paths='/home/adrian/PhD/AxiSEM3D/CMB/stations/3D_UNIFORM_STA.txt')