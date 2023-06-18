import numpy as np
import matplotlib 
matplotlib.use('tkagg')
import matplotlib.pyplot as plt
import matplotlib.animation as animation

from AxiSEM3D_Data_Handler.element_output import element_output

# Import element output file
element_output_path = 'output'
grid_format = [0, 2, 4]
element_obj = element_output(element_output_path, grid_format)







# Create a figure and axis
fig, ax = plt.subplots()

# Generate data for the contour plot
x = np.linspace(-2, 2, 100)
y = np.linspace(-2, 2, 100)
X, Y = np.meshgrid(x, y)
Z = np.sin(np.sqrt(X**2 + Y**2))

# Create the contour plot
contour = ax.contourf(X, Y, Z, cmap='viridis')

# Define the animation update function
def update(frame):
    ax.cla()  # Clear the axis
    # Generate updated data for the contour plot
    t = 0.1 * frame
    Z = np.sin(np.sqrt(X**2 + Y**2) + t)
    contour = ax.contourf(X, Y, Z, cmap='viridis')
    return contour

# Create the animation
ani = animation.FuncAnimation(fig, update, frames=100, interval=100)

# Save the animation as an MP4 file
ani.save('animation.mp4', writer='ffmpeg')

# Display the animation
plt.show()
