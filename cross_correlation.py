import numpy as np
from obspy.signal.cross_correlation import correlate_template, xcorr_max
from obspy import read, read_events, read_inventory
import matplotlib
matplotlib.use('tkagg')
import matplotlib.pyplot as plt
import numpy as np
from obspy.taup import TauPyModel


def _find_distance_in_degree(lat0: float, lon0: float, lat1: float, lon1:float) -> float:
    """Finds the distance in degrees between two geographical points

    Args:
        lat0 (float): point one latitude    
        lon0 (float): point one longitude
        lat1 (float): point two latitude
        lon1 (float): point two longitude

    Returns:
        float: distance in degerees
    """    
    return np.rad2deg(np.arccos(np.cos(np.deg2rad(lat0)) * np.cos(np.deg2rad(lat1)) * np.cos(np.deg2rad(lon0 - lon1)) \
                      + np.sin(np.deg2rad(lat0)) * np.sin(np.deg2rad(lat1))))

def cross_correlation(background_trace, foreground_trace, catalogue, inventory, 
                      network, station, location, channel, phase):
    
    
    # Shift the 0 of the time axis to the evet time
    evt_time = catalogue[0].origins[0].time
    start_time = background_trace.stats.starttime 
    time_shift = start_time - evt_time
    background_time = background_trace.times() + time_shift
    foreground_time = foreground_trace.times() + time_shift
    
    # find the positions of the statino and event
    evt_lat = np.deg2rad(catalogue[0].origins[0].latitude)
    evt_lon = np.deg2rad(catalogue[0].origins[0].longitude)
    evt_depth = catalogue[0].origins[0].depth
    
    sta_lat = inventory[0][0].latitude
    sta_lon = inventory[0][0].longitude
    
    # get backgorund wave data
    background_data = background_trace.data
    # Interpolate foregorund data to background time
    foreground_data = np.interp(background_time, foreground_time, foreground_trace.data)

    # Find the arrivals of the chosen phase
    distance_in_degree = _find_distance_in_degree(evt_lat, evt_lon, sta_lat, sta_lon)
    arrivals = model.get_travel_times(source_depth_in_km=evt_depth * 1e-3,
                                    distance_in_degree=distance_in_degree,
                                    phase_list=[phase])
    
    arrival_times = [arrival.time for arrival in arrivals]
    
    # get the time window based on arrival times
    window_left = min(arrival_times) - 500
    window_right = max(arrival_times) + 500
    left = np.argmin(np.abs(background_time - window_left))
    right = np.argmin(np.abs(background_time - window_right))
    # Find the difference between the foreground and background data 
    # compute correlation
    corr = correlate_template(background_data[left:right], 
                            foreground_data[left:right], 
                            mode = 'full',
                            normalize='full',
                            method = 'direct')
    delta_t_min = background_time[1] - background_time[0]
    print(xcorr_max(corr)[0]*delta_t_min)

    fig1 = plt.figure()
    plt.plot(background_time[left:right], background_data[left:right])
    plt.plot(background_time[left:right], foreground_data[left:right])

    fig2 = plt.figure()
    plt.plot(corr)
    plt.show()


file1 = '/home/adrian/PhD/axisem3d_root/AxiSEM3D-master/cmb/simu1D/XC_20S_1D_PRO/output/obspyfied/XC_20S_1D_PRO.mseed'
file2 = '/home/adrian/PhD/axisem3d_root/AxiSEM3D-master/cmb/simu1D_CMB/XC_20S_1D_CMB7_0_23_17/output/obspyfied/XC_20S_1D_CMB7_0_23_17.mseed'
cat_path = '/home/adrian/PhD/axisem3d_root/AxiSEM3D-master/cmb/simu1D/XC_20S_1D_PRO/output/obspyfied/cat.xml'
inv = read_inventory('/home/adrian/PhD/axisem3d_root/AxiSEM3D-master/cmb/simu1D_CMB/XC_20S_1D_CMB7_CST/output/obspyfied/10DEG_GRID_inv.xml')


background_data = read(file1)
foreground_data = read(file2)
catalogue = read_events(cat_path)

network = 'A'
station = '357'
location = ''
channel = 'LXR'
model = TauPyModel('prem')
phase = "ScS"

background_trace = background_data.select(network, station, location, channel)[0]
foreground_trace = foreground_data.select(network, station, location, channel)[0]
inventory = inv.select(network=network, station=station, location=location, channel=channel)

cross_correlation(background_trace, foreground_trace, catalogue, inventory, 
                      network, station, location, channel, phase)


