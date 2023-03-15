""" 
Code for creating a .txt file with stations data as a grid 
Ready to be used by axisem3D
"""

import numpy as np
import matplotlib
matplotlib.use('tkagg')
import matplotlib.pyplot as plt
import pandas as pd

# ARRAY
# Create Lat Lon grid grid
lons = np.arange(-180, 181, 10)
lats = np.arange(-90, 91, 10)
LAT, LON = np.meshgrid(lats, lons)

# Flatten in C style row-major 
lon_flat = LON.flatten('F')
lat_flat = LAT.flatten('F')

# Create station names with names less than 5 chars
names = []
nets = []
depth = []
useless = []
index = 0
for lat in lats:
    for lon in lons:
        index += 1
        names.append(str(index))
        nets.append('A')
        depth.append(0)
        useless.append(0)     


""" # CROSS
lons_1 = np.arange(-180, 181, 10)
lats_1 = np.zeros(len(lons_1))
lats_2 = np.arange(-90, 91, 10)
lons_2 = np.zeros(len(lats_2))
lon_flat = np.append(lons_1, lons_2)
lat_flat = np.append(lats_1, lats_2)
names = []
nets = []
depth = []
useless = []
index = 0
for lon, lat in zip(lon_flat, lat_flat):
        index += 1
        names.append(str(index))
        nets.append('A')
        depth.append(0)
        useless.append(0) 
 """

# Make and save data frame as .txt file   
data = {'#name': names, 'network': nets, 'latitude': lat_flat,
        'longitude': lon_flat, 'useless': useless, 'depth': depth}
df = pd.DataFrame(data=data)
df.to_csv('stations/10DEG_GRID.txt', sep=' ', index=False)
