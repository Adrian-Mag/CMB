import sys
sys.path.append('/home/adrian/PhD/axisem3d_root/AxiSEM3D-master/cmb/output_handlers')
from element_output import element_output
import matplotlib
matplotlib.use('tkagg')
import matplotlib.pyplot as plt


path = '/home/adrian/PhD/axisem3d_root/AxiSEM3D-master/cmb/simu1D_element/ADJOINT_TEST'
grid_format = [2]
obj = element_output(path, grid_format)
data = obj.load_data_at_point([6370000,30,0])

plt.figure()
plt.plot(data[0,:])
plt.show()