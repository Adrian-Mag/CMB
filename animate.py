from axisem3d_output.core.handlers.element_output import ElementOutput

path = '/disks/data/PhD/CMB/simu1D_element/FORWARD/output/elements/forward_20s_1D'
obj = ElementOutput(path_to_element_output=path)
#obj.stream([6371000, 0, 30], coord_in_deg=True, channels=['E']).plot()
obj.animation([0, 0, 0], [0, 0, 40], R_min=3400000, lower_range=0.6, upper_range=0.99999, 
              resolution=200, frame_rate=10, video_duration=20,
              paralel_processing=True, channels=['UZ', 'UT', 'UR'])