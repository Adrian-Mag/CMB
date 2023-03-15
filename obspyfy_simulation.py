import sys
sys.path.append('/home/adrian/PhD/axisem3d_root/AxiSEM3D-master/cmb/output_handlers')
from obspyfy import obspyfy

obspyfy('/home/adrian/PhD/axisem3d_root/AxiSEM3D-master/cmb/simu1D_element/TEST_ELEMENT_AT2', output_type='elements', 
        stations_paths='/home/adrian/PhD/axisem3d_root/AxiSEM3D-master/cmb/stations/10DEG_GRID.txt')