""" 
Code template for pre-processing seismic data
"""
from obspy.core.inventory import Inventory, Network, Station, Channel
from obspy import read, read_inventory, read_events
from obspy.clients.fdsn import Client
import matplotlib
matplotlib.use('tkagg')
import matplotlib.pyplot as plt
from obspy.geodetics.base import gps2dist_azimuth


c = Client("IRIS")

###########
# READ DATA
###########
path = '/home/adrian/PhD/axisem3d_root/AxiSEM3D-master/cmb/real_seismic_data/SUMATRA_2009/'
name = 'SUMATRA_2009'

streams = read(path + name + '.mseed')
processed_streams = streams.copy()
processed_streams.clear()
invs = read_inventory(path + name + '_inv' + '.xml')
cat = read_events(path + name + '_cat' + '.xml')

##############
# GET METADATA
##############
evt_time = cat.events[0]['origins'][0]['time']
evt_lat = cat.events[0]['origins'][0]['latitude']
evt_lon = cat.events[0]['origins'][0]['longitude']

#############
# PRE PROCESS
#############
INSTRUMENT_REMOVAL = True
FILTERING = False
RESAMPLING = False
DETREND = True
TAPER = True
ROTATE = True
# INSTRUMENT REMOVAL--------------------------------------------------------------------------
if INSTRUMENT_REMOVAL == True:
    pre_filt = [1e-3, 1e-2, 10, 100]
    streams.remove_response(inventory=invs, pre_filt=pre_filt, output="DISP",
                    water_level=60, plot=True) 
    print('Instrument removal done')
else:
    print('No instrument removal done')
# FILTERING--------------------------------------------------------------------------
if FILTERING == True:
    streams.filter('lowpass', freq=1/2)
    print('Filtering done')
else:
    print('No filtering done')
# RESAMPLING-------------------------------------------------------------------------
if RESAMPLING == True:
    streams.decimate(factor=10, strict_length=False)
    print('Resampling done')
else:
    print('No resampling done')
# DETREND/DEMEAN--------------------------------------------------------------------------
if DETREND == True:
    streams.detrend()
    print('Detrending done')
else:
    print('No detrending done')
# TAPER----------------------------------------------------------------------------------------
if TAPER == True:
    streams.taper(max_percentage=0.001)
    print('Tapering done')
else:
    print('No tapering done')
# ROTATE----------------------------------------------------------------------------------
if ROTATE == True:
    for inv in invs:
        for sta in inv:
            sta_lat = sta.latitude
            sta_lon = sta.longitude
            back_azimuth = gps2dist_azimuth(lat1=evt_lat, lon1=evt_lon, lat2=sta_lat, lon2=sta_lon, 
                                    a=6378137.0, f=0.0033528106647474805)[2]
            # find traces of this station
            current_stream = streams.copy()
            current_stream.clear()
            # get code of the inventory (the station code), assuming just one station per inventory
            inv_station = sta.code
            for trace in streams:
                # get code of trace (the station code)
                st_station = trace.id.split(".")[1]
                if st_station == inv_station:
                    current_stream.append(trace)

            current_stream.rotate(method='->ZNE', inventory=invs)
            current_stream.rotate(method='NE->RT', inventory=invs, back_azimuth=back_azimuth)

            # Stitch back pre processed streams
            processed_streams.__iadd__(current_stream)

    print('Rotation into ZRT done')
else:
    print('No rotation into ZRT done')


#####################
# CORRECT INVENTORIES
#####################
correct_inv = Inventory(
        # We'll add networks later.
        networks=[],
        source="get_instaseis_data")
for trace in processed_streams:
    # Go trace by trace and get correct trace info
    [network, station, location, channel] = trace.id.split('.')
    # Find net + sta in inv corresponding to the trace and get geographical info
    inv_select = invs.select(network=network, station=station)
    [sta_lat, sta_lon, sta_elevation ]= [inv_select[0][0].latitude, inv_select[0][0].longitude, inv_select[0][0].elevation]

    # Create network if not already existent
    net_exists = False
    for existing_network in correct_inv:
        if existing_network.code == network:
            net_exists = True
            net = existing_network
    if net_exists is False:
        net = Network(
        code=network,
        stations=[])
        # add new network to inventory
        correct_inv.networks.append(net)

    # Create station if not already existent
    sta_exists = False
    for existing_station in net:
        if existing_station.code == station:
            sta_exists = True
            sta = existing_station
    if sta_exists is False:
        sta = Station(
            code=station,
            latitude=sta_lat,
            longitude=sta_lon,
            elevation=sta_elevation)
        # add station to network
        net.stations.append(sta)
    
    # Create channel
    cha = Channel(code=channel,
                location_code=location, # this is wrong, but is not needed
                latitude=sta_lat, # this is wrong, but is not needed
                longitude=sta_lon, # this is wrong, but is not needed
                elevation=sta_elevation, # this is wrong, but is not needed
                depth=0) # this is wrong, but is not needed)
    # add channel to station
    sta.channels.append(cha)

######
# PLOT
######
processed_streams.plot()
plt.show()


######
# SAVE
######
value = input('Save it? (y/n): ')
if value == 'y':
    processed_streams.write(path + name + '_preprocessed' + '.mseed')
    correct_inv.write(path + name + '_preprocessed_inv' + '.xml', format="STATIONXML")
