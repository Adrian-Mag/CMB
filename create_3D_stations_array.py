"""
Creates evenly or randomly distributed points on a sphere, 
in a 3D ball, or in a 3D shell. 
"""

import math
import numpy as np
from mayavi import mlab
import pandas as pd


def rand_sphere(R_min, R_max, samples):
    points = np.random.uniform(low=-R_max, high=R_max, size=(samples, 3))
    
    # discard points outside the shell
    radii = np.sqrt((points**2).sum(axis=1))
    points = points[radii>R_min,:]
    radii = radii[radii>R_min]
    points = points[radii<R_max,:]
    radii = radii[radii<R_max]
    
    return points

def uniform_sphere(R_min, R_max, samples):
    points = []
    # create square mesh with sides 2R
    X = np.linspace(-R_max, R_max, samples)
    Y = np.linspace(-R_max, R_max, samples)
    Z = np.linspace(-R_max, R_max, samples)
    for x in X:
        for y in Y:
            for z in Z:
                r = np.sqrt(x**2 + y**2 + z**2)
                if r < R_max and r > R_min:
                    points.append([x,y,z])
    
    return np.asarray(points)

############
# Parameters
############
file_name = '3D_UNIFORM_STA.txt'
max_samples = 2000 # this is the no of samples on the biggest sphere
samples = 15 # this is for the uniform one
R_min = 3481e3 
R_max = 6371e3
distribution_type = 'uniform'

#######
# SOLVE
#######

if distribution_type == 'random':
    points = rand_sphere(R_min=R_min, R_max=R_max, samples=max_samples)
else:
    points = uniform_sphere(R_min=R_min, R_max=R_max, samples=samples)
##########
# PLOTTING
########## 

# Define the lat lon grid for the spheres
lat = np.arange(-90, 90.01, 1)*np.pi/180
lon = np.arange(0, 180.01, 1)*np.pi/180
LON, LAT = np.meshgrid(lon, lat)

 # plot
fig_earth = mlab.figure()

for R in [3471000, 6371000]:
    # Construct CMB and Surface matrices
    R_surface = R * np.ones(np.shape(LON))

    # Transform to cartesian coords
    X_surface = R_surface * np.cos(LAT) * np.cos(LON)
    Y_surface = R_surface * np.cos(LAT) * np.sin(LON)
    Z_surface = R_surface * np.sin(LAT)

    # Plot surface
    mlab.mesh(X_surface, Y_surface, Z_surface, color=(0.2,0.2,0.2), opacity=0.1)    

# Plot points
for point in points:
    mlab.points3d(point[0], point[1], point[2], color=(1,1,1), scale_factor=1e5, opacity=1)

mlab.show()

answer = input('Save it? (y/n): ')
if answer == 'y' or 'yes':
    # transform to spherical coordinates
    spherical_points = np.empty((len(points), 3))
    spherical_points[:,0] = np.round(6371000 - np.sqrt((points**2).sum(axis=1)), decimals=0) # depth in m
    spherical_points[:,1] = np.round(np.rad2deg(np.arctan(points[:,2] / np.sqrt(points[:,0]**2 + points[:,1]**2))),  decimals=0) # lat in deg
    spherical_points[:,2] = np.round(np.rad2deg(np.arctan2(points[:,1], points[:,0])),  decimals=0)
    names = np.arange(len(points))
    nets = np.asarray(['A' for _ in range(len(points))])
    useless = np.zeros(len(points))
    # Make and save data frame as .txt file   
    data = {'#name': names, 'network': nets, 'latitude': spherical_points[:,1],
            'longitude': spherical_points[:,2], 'useless': useless, 'depth': spherical_points[:,0]}
    df = pd.DataFrame(data=data)
    df.to_csv('stations/' + file_name, sep=' ', index=False)
    