from obspy import read, read_events, read_inventory
import numpy as np
from obspy.taup import TauPyModel
from obspy.signal.cross_correlation import correlate_template, xcorr_max
import matplotlib
matplotlib.use('tkagg')
import matplotlib.pyplot as plt
import cartopy.crs as crs
import cartopy.feature as cfeature
from matplotlib import cm
from matplotlib.colors import ListedColormap

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


# read in the background model and event 
bkg_stream = read('/home/adrian/PhD/axisem3d_root/AxiSEM3D-master/cmb/simu1D/XC_20S_1D_PRO/output/obspyfied/XC_20S_1D_PRO.mseed')
fg_stream = read('/home/adrian/PhD/axisem3d_root/AxiSEM3D-master/cmb/simu1D_CMB/XC_20S_1D_CMB7_0_23_17/output/obspyfied/XC_20S_1D_CMB7_0_23_17.mseed')
cat = read_events('/home/adrian/PhD/axisem3d_root/AxiSEM3D-master/cmb/simu1D_CMB/XC_20S_1D_CMB7_CST/output/obspyfied/cat.xml')
inv = read_inventory('/home/adrian/PhD/axisem3d_root/AxiSEM3D-master/cmb/simu1D_CMB/XC_20S_1D_CMB7_CST/output/obspyfied/10DEG_GRID_inv.xml')
# The cross correlation will be done over the entire array, but must choose the component
chn = 'R'
bkg_stream = bkg_stream.select(channel='*Z')
fg_stream = fg_stream.select(channel='*Z')

# give a t_span [s]
t_span = 20
# give a model
model = TauPyModel('prem')
# choose a phase
phase = "PcP"


# get the event time
evt_time = cat[0].origins[0].time
evt_depth = cat[0].origins[0].depth
evt_lat = np.deg2rad(cat[0].origins[0].latitude)
evt_lon = np.deg2rad(cat[0].origins[0].longitude)

# initiate lists that will contain lat/lon/corr data
corr_lat = []
corr_lon = []
corr_shift = []

# Go over all stations
for bkg_trace in bkg_stream:
    # select only one trace
    bkg_id = bkg_trace.id
    print(bkg_id)
    fg_trace = fg_stream.select(id=bkg_id)[0]
    # get location of the station of this trace
    [network, station, location, channel] = bkg_id.split('.')
    inv_station = inv.select(network=network, station=station, location=location, channel=channel)
    sta_lat = inv_station[0][0].latitude
    sta_lon = inv_station[0][0].longitude

    # get the correct times
    start_time = bkg_trace.stats.starttime 
    time_shift = start_time - evt_time
    bkg_time = bkg_trace.times() + time_shift
    fg_time = fg_trace.times() + time_shift
    delta_t = bkg_time[1] - bkg_time[0]
    
    # get bkg wave data
    bkg_data = bkg_trace.data
    # Interpolate fg data to bkg time
    fg_data = np.interp(bkg_time, fg_time, fg_trace.data)
    
    # find t_center for a given phase
    distance_in_degree = _find_distance_in_degree(evt_lat, evt_lon, sta_lat, sta_lon)
    arrivals = model.get_travel_times(source_depth_in_km=evt_depth * 1e-3,
                                    distance_in_degree=distance_in_degree,
                                    phase_list=[phase])
    if arrivals:
        t_center = arrivals[0].time + t_span
        
        # find interval for cross correlation 
        t_right_bkg  = np.argmin(np.abs(bkg_time - t_center - 2*t_span))
        t_left_bkg  = np.argmin(np.abs(bkg_time - t_center + 2*t_span))
        t_right  = np.argmin(np.abs(bkg_time - t_center - t_span))
        t_left  = np.argmin(np.abs(bkg_time - t_center + t_span))
        
        # compute correlation
        corr = correlate_template(bkg_data[t_left_bkg:t_right_bkg], 
                                fg_data[t_left:t_right], 
                                mode = 'full',
                                normalize='full',
                                method = 'direct')
        # save data 
        corr_lat.append(sta_lat)
        corr_lon.append(sta_lon)
        corr_shift.append(xcorr_max(corr)[0]*delta_t)
    else:
        corr_lat.append(sta_lat)
        corr_lon.append(sta_lon)
        corr_shift.append(0)

print([min(corr_shift), max(corr_shift)])


rdbu = cm.get_cmap('RdBu', 512)
newcmp = ListedColormap(rdbu(np.linspace(0, 1, 256)))

fig = plt.figure(figsize=(10,8))

ax = fig.add_subplot(1,1,1, projection=crs.Robinson())
ax.set_global()
ax.add_feature(cfeature.COASTLINE)
ax.gridlines(crs=crs.PlateCarree(), draw_labels=True)

psm = plt.scatter(x=corr_lon, y=corr_lat,
            s=50,
            c=corr_shift, cmap=newcmp, vmin=-0.3, vmax=0.3,
            transform=crs.PlateCarree()) ## Important
fig.colorbar(psm, ax=ax)
plt.show()


    
    
    