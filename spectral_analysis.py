""" 
Code for spectral analysis on AXISEM3D NETCDF files
and on real data in obspy stream format
"""
import sys
sys.path.append('/home/adrian/PhD/axisem3d_root/AxiSEM3D-master/cmb/output_handlers')
from scipy.signal import zoom_fft
import numpy as np
from station_output import station_output
import matplotlib
matplotlib.use('tkagg')
import matplotlib.pyplot as plt
from obspy.imaging.spectrogram import spectrogram
from scipy.signal import argrelextrema 
import scipy
from obspy import read

# Choose between real or simulated data
data_type = 'simulated'

if data_type == 'simulated':
    # GET DATA FROM AXISEM SIMULATION (NETCDF)---------------------------------------------------
    # Path to simulation
    simu_paths = ['/home/adrian/PhD/axisem3d_root/AxiSEM3D-master/cmb/simu1D/60H_600S_NORMALS/output/stations/Station_grid/']

    # Stations to be analysed
    networks = ['A']
    station_names = ['STA_0_90']
    locations = ['']
    chn = 'R'
    chn_type = 'RTZ'

    # create netcdf station data object
    obj = station_output(simu_paths[0])
    # get time and waveforms as traces
    time = obj.time
    displacement_data = obj.stream(networks=networks, station_names=station_names, locations=locations, channels=chn)
    displacement_data = displacement_data[0]
else:
    # GET REAL DATA (NETCDF)-----------------------------------------------------------------

    path = 'real_seismic_data/'
    name = 'japan-5-70h_pre_processed'

    displacement_data = read(path + name + '.mseed')
    displacement_data = displacement_data[1]
    time = np.arange(0,len(displacement_data)) * displacement_data.stats.delta

# CODE ---------------------------------------------------------

# Compute Nyquist freq
f_nyq = 1/(2*displacement_data.stats.delta)

# Fourier transform 
spectrum_fft = np.fft.rfft(displacement_data.data)

# np fft freq
f_fft = np.linspace(0, f_nyq, len(spectrum_fft))    # frequency axis for plotting

# Find points of maximum in FT
locs = argrelextrema(abs(spectrum_fft), np.greater, order = 2)
spectrum_fft_extremum = abs(spectrum_fft[locs])
f_extremum = f_fft[locs]
f_extremum_filtered = f_extremum[spectrum_fft_extremum > max(spectrum_fft_extremum)/10]

# Info
print('Number of samples in time: ', displacement_data.stats.npts)
print('Nyquist freq: ', f_nyq)
print('FFT freq resolution: ', f_fft[1] - f_fft[0])
print('Periods present: ',1/f_extremum_filtered)

# PLOT---------------------------------------------------------------------------------------
fig1 = plt.figure('SPECTRUM')
plt.plot(f_fft * 1e3, abs(spectrum_fft)/max(abs(spectrum_fft)), color = 'g')
plt.xlim([0, 2]) # frequency in mHz
y_max = 1.1
plt.vlines(x = f_extremum_filtered * 1e3, ymin = 0, ymax = y_max, color = 'r', linewidth = 0.5)
plt.xlabel('Frequencies in mHz')
plt.ylabel('Amplitude spectrum, normalized to 1')

fig2 = plt.figure('displacement_data.data')
plt.plot(time, displacement_data.data)
plt.ylabel('Displacement')
plt.xlabel('Time in s')

plt.show()
