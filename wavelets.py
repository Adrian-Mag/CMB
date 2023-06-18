import numpy as np
from obspy import read, read_events, read_inventory
import matplotlib
matplotlib.use('tkagg')
import matplotlib.pyplot as plt
import numpy as np
from obspy.taup import TauPyModel
from scipy import signal
from obspy.signal.tf_misfit import plot_tfr, plot_tf_misfits, tfpm


# GET DATA 
file1 = '/home/adrian/PhD/AxiSEM3D/CMB/simu1D/OBLIQUE_FAULT_100KM/output/obspyfied/OBLIQUE_FAULT_100KM.mseed'
file2 = '/home/adrian/PhD/AxiSEM3D/CMB/SIMU_JOBS/OBLIQUE_FAULT_100KM_CMB_5_0_30_30/output/obspyfied/OBLIQUE_FAULT_100KM_CMB_5_0_30_30.mseed'
cat_path = '/home/adrian/PhD/AxiSEM3D/CMB/SIMU_JOBS/OBLIQUE_FAULT_100KM_CMB_5_0_30_30/output/obspyfied/cat.xml'
inv = read_inventory('/home/adrian/PhD/AxiSEM3D/CMB/SIMU_JOBS/OBLIQUE_FAULT_100KM_CMB_5_0_30_30/output/obspyfied/STA_10DEG_GRID_inv.xml')


background_data = read(file1)
foreground_data = read(file2)
catalogue = read_events(cat_path)

network = 'A'
station = '354'
location = ''
channel = 'LXR'
model = TauPyModel('prem')
phase = "ScS"

background_trace = background_data.select(network, station, location, channel)[0]
foreground_trace = foreground_data.select(network, station, location, channel)[0]
# filter
#background_data.filter("bandpass", freqmin = 0.01, freqmax = 1)
#foreground_data.filter("bandpass", freqmin = 0.03333, freqmax = 0.05)
inventory = inv.select(network=network, station=station, location=location, channel=channel)

# Shift the 0 of the time axis to the event time
evt_time = catalogue[0].origins[0].time
start_time = background_trace.stats.starttime 
time_shift = start_time - evt_time
background_time = background_trace.times() + time_shift
foreground_time = foreground_trace.times() + time_shift
dt = background_time[1] - background_time[0]

# get backgorund wave data
background_data = background_trace.data
# Interpolate foregorund data to background time
foreground_data = np.interp(background_time, foreground_time, foreground_trace.data)

# find sampling freq and nyquist freq
fs = 1/dt
f_nyq = fs/2
print(f_nyq)


plot_tf_misfits(foreground_data, background_data, dt=dt, t0=min(background_time), fmin=1e-2, fmax=0.1, nf=100, w0=6,
                    norm='global', clim=0.)