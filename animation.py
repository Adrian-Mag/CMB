import numpy as np
import matplotlib 
matplotlib.use('tkagg')
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from tqdm import tqdm 
import multiprocessing as mp


import sys
sys.path.append('/disks/data/PhD/AxiSEM3D-Kernels')
from AxiSEM3D_Data_Handler.element_output import ElementOutput
from AxiSEM3D_Kernels import sph2cart, cart2sph

############
# PARAMETERS
############

# Import element output file
element_output_path = '/disks/data/PhD/CMB/simu1D_element/BACKWARD_UNIT_DELAY'
name = 'backward_unit_delay'
grid_format = [0, 2, 4]
source_loc = [6371000, 0, 0]
station_loc = [6371000, 0, 30]
R_max = 6371000
R_min = 3400000
N = 200

video_duration = 20 # seconds
frame_rate = 10 # fps


################
# IMPLEMENTATION
################

def find_smallest_value(arr, percentage):
    # Flatten the array to a 1D array
    flattened = arr[np.isfinite(arr)].flatten()
    
    if len(flattened) == 0:
        return None

    # Sort the flattened array in ascending order
    sorted_arr = np.sort(flattened)
    
    # Compute the index that corresponds to 90% of the values
    percentile_index = int(len(sorted_arr) * (1 - percentage))
    
    # Get the value at the computed index
    smallest_value = sorted_arr[percentile_index]
    
    return smallest_value   

# Load data
element_obj = ElementOutput(element_output_path)
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
inplane_DIM1, inplane_DIM2 = np.meshgrid(inplane_dim1, inplane_dim2, indexing='xy')

# Initialize sensitivity values on the slice (Slice frame)
inplane_field = np.zeros((N, N, 3, len(time)))

# Load the data at all points
print('Loading data')     
pbar = tqdm(total=len(inplane_dim1) * len(inplane_dim2))
for index1 in range(len(inplane_dim1)):
    for index2 in range(len(inplane_dim2)):
        [x, y, z] = inplane_dim1[index1] * base1 + inplane_dim2[index2] * base2  # Slice frame -> Earth frame
        rad, lat, lon = cart2sph(x, y, z)
        if rad > R_min and rad < R_max:
            inplane_field[index2, index1, :, :] = element_obj.load_data_at_point([rad, np.rad2deg(lat), np.rad2deg(lon)])
        else:
            inplane_field[index2, index1, :, :] = np.full((3, len(time)), np.nan)
        pbar.update(1)
pbar.close()

print('Create animation')
# Create a figure and axis
cbar_min = find_smallest_value(np.log10(np.abs(inplane_field)), 0.5)
cbar_max = np.nanmax(np.log10(np.abs(inplane_field)))

fig, ax = plt.subplots()
contour = ax.contourf(inplane_DIM1, inplane_DIM2, np.nan_to_num(np.log10(np.abs(inplane_field[:, :, 0, 0]))), levels=np.linspace(cbar_min, cbar_max, 100), cmap='RdBu_r', extend='both')

def update(frame):
    ax.cla()
    contour = ax.contourf(inplane_DIM1, inplane_DIM2, np.nan_to_num(np.log10(np.abs(inplane_field[:, :, 0, frame * frame_step]))), levels=np.linspace(cbar_min, cbar_max, 100), cmap='RdBu_r', extend='both')
    plt.scatter(np.dot(point1, base1), np.dot(point1, base2))
    plt.scatter(np.dot(point2, base1), np.dot(point2, base2))
    print(100 * frame / (video_duration * frame_rate), '%')
    return contour


cbar = plt.colorbar(contour)

cbar_ticks = np.linspace(int(cbar_min), int(cbar_max), 5) # Example tick values
cbar_ticklabels = [str(cbar_tick) for cbar_tick in cbar_ticks] # Example tick labels
cbar.set_ticks(cbar_ticks)
cbar.set_ticklabels(cbar_ticklabels)
cbar.set_label('Intensity')

frame_step = int(len(time) / ( video_duration * frame_rate ))
ani = animation.FuncAnimation(fig, update, frames=video_duration * frame_rate, interval=1e3 / frame_rate)
ani.save(name + '_animation.mp4', writer='ffmpeg')