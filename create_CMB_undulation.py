""" 
Code for creating undulation file for axisem 3D
Takes in some file that contains CMB data in lat lon format
and spews out .nc file with undulation ready to use with axisem3D
"""


import numpy as np
from scipy.interpolate import interp1d
from skimage.filters import gaussian
import matplotlib
matplotlib.use('tkagg')
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import cartopy.crs as ccrs
import os

# Number of grid points on latitude and longitude (should be the same
# with the one in the topography file)
nlat = 181
nlon = 361
nrow = nlat * nlon

# read raw data
txt_data = '/home/adrian/PhD/axisem3d_root/AxiSEM3D-master/cmb/cmb_models/cmb_topo_gaussian_7_CST.txt'
name = txt_data.split('/')[-1].split('.')[0]
bnd = np.loadtxt(txt_data)
bnd = bnd[:,2]

# Put data in matrix form and put average of the data at the poles
CMB = np.zeros((nlat, nlon))
# north pole
CMB[0, :] = bnd[0:nlon].sum() / nlon
# south pole
CMB[-1, :] = bnd[nrow - nlon:nrow].sum() / nlon
for i in range(1, nlat):
    CMB[i, :] = bnd[i * nlon:i * nlon + nlon]

# revert south to north
CMB = np.flip(CMB, axis=0)

# convert to SI
CMB *= 1e3

# structured grid
grid_lat = np.linspace(-90, 90, nlat)
grid_lon = np.linspace(-180, 180, nlon)

# Plot
ax = plt.axes(projection=ccrs.PlateCarree())
ax.set_global()
ax.coastlines()
ax.contourf(grid_lon, grid_lat, CMB)
ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True)
plt.show()

answer = input('Save it ?(y/n)')

if answer == "y":
    if os.path.isfile('cmb_models/' + name):
        print('file exists')
    else:
        # write to file
        nc = Dataset('cmb_models/' + name + '.nc', 'w')
        nc.createDimension('nlat', size=len(grid_lat))
        nc.createDimension('nlon', size=len(grid_lon))
        nc.createVariable('latitude', float, dimensions=('nlat'))
        nc['latitude'][:] = grid_lat
        nc.createVariable('longitude', float, dimensions=('nlon'))
        nc['longitude'][:] = grid_lon
        nc.createVariable('undulation_CMB', float, dimensions=('nlat', 'nlon'))
        nc['undulation_CMB'][:, :] = CMB
        nc.close()