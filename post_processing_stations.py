""" 
Code that creates animation of wave propagation across 
an array of stations. 

The code was already present in the example files of AXISEM3D
and I just modified it slightly to make it shorter 
"""

import os
import numpy as np
import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import animation
import cartopy.crs as ccrs
from matplotlib import colorbar, colors

data_dir_stations = '/home/adrian/PhD/axisem3d_root/AxiSEM3D-master/examples/10_S362ANI_element_global/simu1D/output/stations/Station_grid'    
dir_stations = '/home/adrian/PhD/axisem3d_root/AxiSEM3D-master/examples/10_S362ANI_element_global/simu1D/input/STATIONS.txt'
# since all of the station data is in a uniform format, we can easily use xarray's open_mfdataset to combine it all automatically
data = xr.open_mfdataset(data_dir_stations+"/axisem3d_synthetics.nc.*", engine="netcdf4", data_vars="different", concat_dim="dim_station", combine="nested")

# load station location information in AxiSEM coordinates
stations = pd.read_csv(dir_stations, delim_whitespace=True, header=0, names=["name","network","latitude","longitude","useless","depth"])

# find the permutation vector from the data output to the stations
stations["name_network"] = [net+'.'+nam for (nam, net) in zip(stations["name"], stations["network"])]
station_names_decoded = np.array([ls.decode("utf-8") for ls in data["list_station"].values])
permutation_vector = np.array([(nn==station_names_decoded).argmax() for nn in stations["name_network"]])
stations["permutation"] = permutation_vector

# put everything into simple numpy arrays for plotting, including sorting the wave data
lat_coords = stations["latitude"].values
lon_coords = stations["longitude"].values
t_coords = data["data_time"].values
displacement_data = data["data_wave"].values[stations["permutation"],:,:]

fig, axes = plt.subplots(2, 2, figsize=(8,9), dpi=300, constrained_layout = True, subplot_kw={'projection': ccrs.PlateCarree()})

i = 4000
animate = True

axes[0,0].set_title("Z", fontsize=16)
axes[0,1].set_title("R", fontsize=16)
axes[1,0].set_title("T", fontsize=16)
axes[1,1].set_title("PGV to Time", fontsize=16)

for col in range(2):
    for row in range(2):
        axes[row,col].set_extent((-180, 180, -80, 80), crs=ccrs.PlateCarree())
        axes[row,col].coastlines(resolution="10m")

# create fixed color bars using linear Normalize and dummy ScalarMapables
disp_norm = mpl.colors.Normalize(vmin=-1e-9, vmax=1e-9)
pgv_norm = mpl.colors.Normalize(vmin=0, vmax=0.1)
fig.colorbar(mpl.cm.ScalarMappable(norm=disp_norm, cmap='seismic'), ax=axes[1, 0], shrink=0.8, label="Velocity (m/s)", location="bottom")
fig.colorbar(mpl.cm.ScalarMappable(norm=pgv_norm, cmap='RdPu'), ax=axes[1, 1], shrink=0.8, label="PGV (m/s)", location="bottom")

if animate:
    # Use ArtistAnimation - flexible, but it is memory hungry as you have to make every plot element for every frame first
    frames = []
    props = dict(boxstyle='round', facecolor='wheat')
    for i in range(0,len(t_coords),100):
        timelabel = axes[0,0].text(0.025,0.95, f"{t_coords[i]:.2f} s", transform=axes[0,0].transAxes, ha="left", bbox=props)
        # tricontourf lets us use an unstructured set of points for contouring, and doesn't look too bad when we have a lot of stations
        tcf00 = axes[0,0].tricontourf(lat_coords, lon_coords, displacement_data[:,2,i], levels=100, cmap="seismic")
        tcf10 = axes[0,1].tricontourf(lat_coords, lon_coords, displacement_data[:,0,i], levels=100, cmap="seismic")
        tcf01 = axes[1,0].tricontourf(lat_coords, lon_coords, displacement_data[:,1,i], levels=100, cmap="seismic")
        tcf11 = axes[1,1].tricontourf(lat_coords, lon_coords, displacement_data[:,0,i], levels=100, cmap="RdPu")
        # collect each frame and append to frame list
        frames.append(tcf00.collections+tcf10.collections+tcf01.collections+tcf11.collections+[timelabel]) 
        print(str(i) + '/' + str(len(t_coords)))
    anim = animation.ArtistAnimation(fig, frames)
    Writer = animation.writers['ffmpeg']
    writer = Writer(fps=10, metadata=dict(artist='Me'), bitrate=1800)

    anim.save('finite_fault_movie.mp4', writer=writer)
else:    
    fig.suptitle(f"Time = {t_coords[i]:.2f} s", fontsize=16)
    axes[0,0].tricontourf(lat_coords, lon_coords, displacement_data[:,2,i], levels=100, cmap="seismic"),
    axes[0,1].tricontourf(lat_coords, lon_coords, displacement_data[:,0,i], levels=100, cmap="seismic"),
    axes[1,0].tricontourf(lat_coords, lon_coords, displacement_data[:,1,i], levels=100, cmap="seismic"),
    axes[1,1].tricontourf(lat_coords, lon_coords, displacement_data[:,0,i], levels=100, cmap="RdPu")
    fig.savefig("finite_fault_figure.png")