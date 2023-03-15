""" Creates a txt space separated file of the form
    latitude[deg]    longitude[deg]   CMB[km]
        .                   .           .
    that contains the geometric description of the 
    CMB topography     
"""
import numpy as np
import matplotlib
matplotlib.use('tkagg')
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import os 

def analytical(LAT: np.ndarray, LON: np.ndarray) -> np.ndarray:
    """Analytical expression for CMB topography

    Args:
        LAT (np.ndarray): latitude matrix from meshgrid
        LON (np.ndarray): longitude matrix from meshgrid

    Returns:
        np.ndarray: CMB topography at each lat lon grid point 
    """
    sigma=0.057509
    a = 7
    print("Extent: ", 2*sigma*np.sqrt(np.log(10))*180/np.pi)
    return a*np.exp(-0*((LAT*np.pi/180)**2 + ((LON - 23)*np.pi/180)**2)/(2*sigma**2))

#give name
name = 'cmb_topo_gaussian_7_CST.txt'

# Define the lat lon grid 
lat = np.arange(-90, 90.01, 1)
lon = np.arange(-180, 180.01, 1)
LAT, LON = np.meshgrid(lat, lon)

# Create topography data from analytical function
cmb = analytical(LAT, LON)

# Set very small values to zero
#cmb[cmb<1e-6] = 0

# Put together data
a = np.zeros((len(LAT.flatten('F')), 3))
a[:,0] = np.transpose(LAT.flatten('F'))
a[:,1] = np.transpose(LON.flatten('F'))
a[:,2] = np.transpose(cmb.flatten('F'))
print(max(a[:,2]))
ax = plt.axes(projection=ccrs.PlateCarree())
ax.set_global()
ax.coastlines()
ax.contourf(lon, lat, np.transpose(cmb))
ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True)
plt.show()

answer = input('Save it ?(y/n)')

if answer == "y":
    # Save to text
    if os.path.isfile('cmb_models/' + name):
        print('file exists')
    else:
        np.savetxt('cmb_models/' + name, a, delimiter=' ', fmt="%s")