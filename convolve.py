import numpy as np
from obspy import read
import matplotlib
matplotlib.use('tkagg')
import matplotlib.pyplot as plt

############
# DEFINE STF
############
def gaussian(dt, gamma):
    time = np.arange(0,)
    return np.exp(-gamma * (time)**2)


stream = read('/home/adrian/PhD/axisem3d_root/AxiSEM3D-master/cmb/simu1D/COLOMBI_GREEN/output/obspyfied/Station_grid.mseed')
stream = stream.select(station='10', channel='*Z')
fig_seismographs, ax = plt.subplots(len(stream))
if len(stream) > 2:
    for index, trace in enumerate(stream):
        id = str(trace.stats.network) + str(trace.stats.station) + str(trace.stats.location) + str(trace.stats.channel)
        time = np.arange(0, trace.data.size)*trace.stats.delta
        ax[index].plot(time, trace.data, lw=0.3, label=id)
else:
    trace = stream[0]
    time = np.arange(0, trace.data.size)*trace.stats.delta
    np.convolve(trace.data, gaussian(time[1] - time[0], gamma=0.5, T0=0))
    id = str(trace.stats.network) + str(trace.stats.station) + str(trace.stats.location) + str(trace.stats.channel)
    ax.plot(time, trace.data, lw=0.3, label=id)
plt.show()