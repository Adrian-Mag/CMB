""" 
Code for creating vtk frames with wavefield in an Earth slice. 
The output can be dumped in paraview to make a video.

Code was already present in the example files of the AXISEM3D package
I've modified it slightly for myself. 
"""

import os
import pyvtk
import numpy as np
import xarray as xr

# The data structure in element-wise output is too complicated for xarray.open_mfdataset.
# Here we open the files as individual datasets and concatenate them on the variable level.
# This code is compatible with parallel netcdf build (single file output)

# load_wave_data=True:  read wave data and return numpy.ndarray
# load_wave_data=False: do not read wave data and return xarray.DataArray (use False if data is big)

def read_element_output(data_dir, load_wave_data=True):
    ################ open files ################
    # filenames
    nc_fnames = [f for f in os.listdir(data_dir) if 'axisem3d_synthetics.nc' in f]
    print('files to open: ', nc_fnames)

    # open files
    nc_files = []
    for nc_fname in nc_fnames:
        nc_files.append(xr.open_dataset(data_dir + '/' + nc_fname))
    
    ################ variables that are the same in the datasets ################
    # read Na grid (all azimuthal dimensions)
    na_grid = nc_files[0].data_vars['list_na_grid'].values.astype(int)

    # read time
    data_time = nc_files[0].data_vars['data_time'].values
    
    
    ################ variables to be concatenated over the datasets ################
    # define empty lists of xarray.DataArray objects
    xda_list_element_na = []
    xda_list_element_coords = []
    dict_xda_list_element = {}
    dict_xda_data_wave = {}
    for nag in na_grid:
        dict_xda_list_element[nag] = []
        dict_xda_data_wave[nag] = []
    
    # loop over nc files
    for nc_file in nc_files:
        # append DataArrays
        xda_list_element_na.append(nc_file.data_vars['list_element_na'])
        xda_list_element_coords.append(nc_file.data_vars['list_element_coords'])
        for nag in na_grid:
            dict_xda_list_element[nag].append(nc_file.data_vars['list_element__NaG=%d' % nag])
            dict_xda_data_wave[nag].append(nc_file.data_vars['data_wave__NaG=%d' % nag])
            
    # concat xarray.DataArray
    xda_list_element_na = xr.concat(xda_list_element_na, dim='dim_element')
    xda_list_element_coords = xr.concat(xda_list_element_coords, dim='dim_element')
    for nag in na_grid:
        dict_xda_list_element[nag] = xr.concat(dict_xda_list_element[nag], dim='dim_element__NaG=%d' % nag)
        dict_xda_data_wave[nag] = xr.concat(dict_xda_data_wave[nag], dim='dim_element__NaG=%d' % nag)
    # read data to numpy.ndarray
    list_element_na = xda_list_element_na.values.astype(int)
    list_element_coords = xda_list_element_coords.values
    dict_list_element = {}
    dict_data_wave = {}
    for nag in na_grid:
        dict_list_element[nag] = dict_xda_list_element[nag].values.astype(int)
        if load_wave_data:
            dict_data_wave[nag] = dict_xda_data_wave[nag].values
        
    ############### return ################
    if load_wave_data:
        return na_grid, data_time, list_element_na, list_element_coords, dict_list_element, dict_data_wave
    else:
        return na_grid, data_time, list_element_na, list_element_coords, dict_list_element, dict_xda_data_wave

# this function extracts waveforms with a given location (x, y)
def get_waveform_xy(xy, na_grid, list_element_na, list_element_coords, dict_data_wave,
                    channels=None, time_steps=None):
    # (x, y) to (s, phi)
    s = np.linalg.norm(xy)
    phi = np.arctan2(xy[1], xy[0])
    
    # point out of range
    if s > np.max(list_element_coords[:, :, 0]):
        return None
    
    # deal with default input
    if channels is None:
        channels = np.arange(dict_data_wave[na_grid[0]].shape[3])
    if time_steps is None:
        time_steps = np.arange(dict_data_wave[na_grid[0]].shape[4])
        
    
    ########## step 1: inplane interpolation ##########
    # find closest element using center coordinate
    s_center = list_element_coords[:, 2, 0]
    index_element = np.argmin(np.abs(s - s_center)) 
    # find the two GLL points, A and B, between which s is located
    s_element = list_element_coords[index_element, :, 0]
    index_A = np.argmin(np.abs(s - s_element))
    s_element_copy = s_element.copy()
    # set s of A to a crazy value to find the second closest point
    s_element_copy[index_A] = 1e100
    index_B = np.argmin(np.abs(s - s_element_copy))
    # interpolation factor at A
    factor_A = 1. / (s_element[index_B] - s_element[index_A]) * (s_element[index_B] - s) 
    factor_B = 1 - factor_A
    
    # na of closest elements
    # the FIVE columes are: 
    # 0 element tag in mesh
    # 1 actual nr
    # 2 stored nr (or na grid)
    # 3 element index in data (local)
    # 4 element index in data (global)
    element_na = list_element_na[index_element]
    
    # read waveforms at A and B and do inplane interpolation
    data_wave_A = dict_data_wave[element_na[2]][element_na[4], :, index_A][:, channels][:, :, time_steps]
    data_wave_B = dict_data_wave[element_na[2]][element_na[4], :, index_B][:, channels][:, :, time_steps]
    data_wave = data_wave_A * factor_A + data_wave_B * factor_B
    
    ########## step 2: Fourier interpolation ##########
    # complex type
    complex_type = np.complex32 if data_wave.dtype == np.complex64 else np.complex128
    
    # maximum Fourier order
    max_Fourier_order = element_na[1] // 2

    # initialize result with 0th order
    result = data_wave[0].copy()
    # add higher orders
    for order in np.arange(1, max_Fourier_order + 1):
        coeff = np.zeros(result.shape, dtype=complex_type)
        # real part
        coeff.real = data_wave[order * 2 - 1]
        # complex part of Fourier coefficients
        if order * 2 < len(data_wave): # check for Nyquist
            coeff.imag += data_wave[order * 2]
        result += (2. * np.exp(1j * order * phi) * coeff).real
    return result


# data dir
data_dir = '/home/adrian/PhD/axisem3d_root/AxiSEM3D-master/examples/10_S362ANI_element_global/simu1D/NORMALS/output/elements/orthogonal_azimuthal_slices'

# read
na_grid, data_time, list_element_na, list_element_coords, \
dict_list_element, dict_data_wave = read_element_output(data_dir)

# wave dimension to animation
wave_channel = 'RTZ'
wave_dim = 'RTZ'.index(wave_channel)

# time steps
ntime = len(data_time)

# phi of the slices
phi_slices = np.radians(np.arange(0, 360, 360))
nslice = len(phi_slices)

# GLL coords on elements
nelem = list_element_coords.shape[0]
ngll = list_element_coords.shape[1]
# flattened coords, (s, z)
element_coords_sz = list_element_coords.reshape((nelem * ngll), 2)

# connectivity list, shared by all slices
# with GLL_points_one_edge = [0, 2, 4] in the inparam.output.yaml,
# the GLL point layout on each element is
# 0--1--2
# |  |  |
# 3--4--5
# |  |  |
# 6--7--8
connectivity = []
for ielem in np.arange(nelem):
    start = ielem * 9
    connectivity.append([start + 0, start + 1, start + 4, start + 3])
    connectivity.append([start + 1, start + 2, start + 5, start + 4])
    connectivity.append([start + 3, start + 4, start + 7, start + 6])
    connectivity.append([start + 4, start + 5, start + 8, start + 7])

# loop over slices
for islice, phi in enumerate(phi_slices):
    # create vtk folder
    vtk_dir = data_dir + '/vtk/slice%d' % islice
    os.makedirs(vtk_dir, exist_ok=True)
    
    # vtk mesh
    xyz = np.ndarray((nelem * ngll, 3))
    xyz[:, 0] = element_coords_sz[:, 0] * np.cos(phi)
    xyz[:, 1] = element_coords_sz[:, 0] * np.sin(phi)
    xyz[:, 2] = element_coords_sz[:, 1]
    vtk_mesh = pyvtk.UnstructuredGrid(list(zip(xyz[:,0], xyz[:,1], xyz[:,2])), 
                                      quad=connectivity)

    # loop over elements to read wave data
    wave = np.ndarray((nelem * ngll, ntime))
    for ielem in np.arange(nelem):
        wave[(ielem * ngll):(ielem * ngll + ngll), :] = dict_data_wave[nslice][ielem, islice, :, wave_dim, :]
        
    # loop over time to write vtk
    for itime in np.arange(ntime):
        vtk = pyvtk.VtkData(vtk_mesh, pyvtk.PointData(pyvtk.Scalars(wave[:, itime], name='U' + wave_channel)))    
        vtk.tofile(vtk_dir + '/' + 'wave%d.vtk' % itime, 'binary')
        print('Done time step %d / %d' % (itime + 1, ntime), end='\r')
    print('\nDone slice %d / %d' % (islice + 1, len(phi_slices)))