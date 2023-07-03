from AxiSEM3D_Data_Handler.station_output import StationOutput
from AxiSEM3D_Data_Handler.element_output import ElementOutput

path = '/disks/data/PhD/CMB/simu1D_element/FORWARD/output/elements/entire_earth'
obj = ElementOutput(path_to_element_output=path)
obj.obspyfy('/disks/data/PhD/CMB/input3D_CMB/STA_10DEG_GRID.txt')
