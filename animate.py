from AxiSEM3D_Data_Handler.station_output import StationOutput
from AxiSEM3D_Data_Handler.element_output import ElementOutput

path = '/disks/data/PhD/CMB/simu1D_element/FORWARD/output/elements/entire_earth'
obj = ElementOutput(path_to_element_output=path)
#obj.stream([6371000, 0, 30], coord_in_deg=True, channels=['E']).plot()
obj.animation([0, 0, 0], [0, 0, 30], R_min=3400000, lower_range=0.3, upper_range=0.99, resolution=300, paralel_processing=True)