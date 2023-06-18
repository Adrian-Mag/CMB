import numpy as np
import matplotlib 
matplotlib.use('tkagg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def plot_gaussian_noise_on_earth():
    # Define the lat lon grid with higher resolution
    lat = np.linspace(-90, 90, 361)
    lon = np.linspace(-180, 180, 721)
    lon_grid, lat_grid = np.meshgrid(lon, lat)

    # Generate random Gaussian noise
    noise = np.random.normal(0, 1, size=lat_grid.shape)

    # Transform to cartesian coordinates
    R = 6371  # Earth radius in kilometers
    X = R * np.cos(np.radians(lat_grid)) * np.cos(np.radians(lon_grid))
    Y = R * np.cos(np.radians(lat_grid)) * np.sin(np.radians(lon_grid))
    Z = R * np.sin(np.radians(lat_grid))

    # Create a 3D plot
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    # Plot the surface
    ax.plot_surface(X, Y, Z, facecolors=plt.cm.viridis(noise), rstride=1, cstride=1)

    # Set labels and title
    ax.set_xlabel('X (km)')
    ax.set_ylabel('Y (km)')
    ax.set_zlabel('Z (km)')
    ax.set_title('Smooth Random Gaussian Noise on the Surface of the Earth')

    # Show the plot
    plt.show()
    
plot_gaussian_noise_on_earth()