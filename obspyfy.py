from axisem3d_output.core.handlers.station_output import StationOutput
from axisem3d_output.core.handlers.element_output import ElementOutput

path = '/disks/data/PhD/CMB/simu3D_CMB/REAL_DATA/output/stations/Station_grid'
obj = StationOutput(path)
obj.obspyfy()
