import numpy as np
import matplotlib 
matplotlib.use('tkagg')
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import tqdm 

from AxiSEM3D_Data_Handler.element_output import element_output
from AxiSEM3D_Kernels.helper_functions import window_data, sph2cart, cart2sph


# Import element output file
element_output_path = 'output'
grid_format = [0, 2, 4]
source_loc = [6371000, 0, 0]
station_loc = [6371000, 0, 30]
R_max = 6371000
R_min = 3400000
N = 30
channel = 0 # RTZ

# Load data
element_obj = element_output(element_output_path, grid_format)
time = element_obj.data_time

# Form vectors for the two points (Earth frame)
point1 = sph2cart(source_loc[0], source_loc[1], source_loc[2])
point2 = sph2cart(station_loc[0], station_loc[1], station_loc[2])

# Do Gram-Schmidt orthogonalization to form slice basis (Earth frame)
base1 = point1 / np.linalg.norm(point1)
base2 = point2 - np.dot(point2, base1) * base1
base2 /= np.linalg.norm(base2)

# Generate inplane slice mesh (Slice frame)
inplane_dim1 = np.linspace(-R_max, R_max, N)
inplane_dim2 = np.linspace(-R_max, R_max, N)
inplane_DIM1, inplane_DIM2 = np.meshgrid(inplane_dim1, inplane_dim2)

# Initialize sensitivity values on the slice (Slice frame)
inplane_field = np.zeros((N, N))

# Load the data at all points
print('Loading data')
with tqdm(total=N**2) as pbar:
    for index1 in range(N):
        for index2 in range(N):
            [x, y, z] = inplane_dim1[index1] * base1 + inplane_dim2[index2] * base2  # Slice frame -> Earth frame
            rad, lat, lon = cart2sph(x, y, z)
            inplane_field[index1, index2] = element_obj(rad, lat, lon)
            pbar.update(1)
            
print('Create animation')
# Create a figure and axis
fig, ax = plt.subplots()
contour = ax.contourf(inplane_DIM1, inplane_DIM2, inplane_field[0])

def update(frame):
    ax.cla()
    contour = ax.contourf(inplane_DIM1, inplane_DIM2, inplane_field[frame])
    return contour

ani = animation.FuncAnimation(fig, update, frames=100, interval=10)

ani.save('animcation.mp4', writer='ffmpeg')

plt.show()