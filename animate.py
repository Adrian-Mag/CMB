from axisem3d_output.core.handlers.element_output import ElementOutput

path = '/home/adrian/PhD/AxiSEM3D/CMB/recipe/output/elements/test'
obj = ElementOutput(path)
obj.animation(source_location=[0,0,0], station_location=[0,0,30],
              R_min=0, R_max=1000000, paralel_processing=False, channels=['UZ'])