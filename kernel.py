from element_interpolator import elements
import numpy as np

# Get the forward and backward data 
path_to_forward = '/home/adrian/PhD/axisem3d_root/AxiSEM3D-master/cmb/simu1D_element/ADJOINT_TEST'
path_to_backward = '/home/adrian/PhD/axisem3d_root/AxiSEM3D-master/cmb/simu1D_element/ADJOINT_TEST'

# Create element objects
forward_data = elements(path_to_forward, [2])
backward_data = elements(path_to_backward, [2])

