import matplotlib
matplotlib.use('tkagg')
import matplotlib.pyplot as plt
from obspy import read_inventory, UTCDateTime
from obspy.clients.fdsn import Client
plt.style.use('ggplot')
plt.rcParams['figure.figsize'] = (10, 8)
import instaseis
from obspy.core.inventory import Inventory, Network, Station, Channel
import os 
from obspy.core.event import Catalog, Event, Origin, FocalMechanism, MomentTensor, Tensor
from obspy.geodetics import FlinnEngdahl
import shutil
import obspy 

###############
# OPEN DATABASE
###############
model = 'prem_a_2s'
db = instaseis.open_db('syngine://' + model)

###################
# SOURCE PARAMETERS
###################
CREATE_NEW_SOURCE = False
SOURCE_TYPE = 'FINITE_SOURCE' # can be POINT_SOURCE, FINITE_SOURCE, GREENS

# MANUAL POINT SOURCE
latitude=0
longitude=0
depth_in_m=100e3
m_rr=1e20
m_tt=1e20
m_pp=1e20
m_rt=0
m_rp=0
m_tp=0

# QUAKEML POINT SOURCE
path_to_quakeml_cat = '/home/adrian/PhD/axisem3d_root/AxiSEM3D-master/cmb/simu1D/COLOMBI/output/obspyfied/cat.xml'

# USGS.param FINITE SOURCE
path_to_usgsparam = '/home/adrian/PhD/axisem3d_root/AxiSEM3D-master/cmb/real_sources/SUMATRA_2009.param'

#####################
# STATIONS PARAMETERS
#####################
CREATE_NEW_INVENTORY = False

# MANUAL INVENTORY
c = Client("IRIS")
net = "II,IU"
sta = "BFO,CTAO,MAJO"
loc = "00"
chn = "BH?"

# IMPORT INVENTORY
path_to_inv = '/home/adrian/PhD/axisem3d_root/AxiSEM3D-master/cmb/real_seismic_data/SUMATRA_2009/SUMATRA_2009_preprocessed_inv.xml'


###############
# CREATE SOURCE
###############
if CREATE_NEW_SOURCE is True:

    if SOURCE_TYPE == 'POINT_SOURCE' or SOURCE_TYPE == 'GREENS':
        # Point source manual~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Create the instaseis source 
        source = instaseis.Source(latitude=latitude,
                                  longitude=longitude,
                                  depth_in_m=depth_in_m,
                                  m_rr=m_rr,
                                  m_tt=m_tt,
                                  m_pp=m_pp,
                                  m_rt=m_rt,
                                  m_rp=m_rp,
                                  m_tp=m_tp)

        # Create the catalogue xml 
        cat = Catalog()
        event = Event()
        origin = Origin()
        
        origin.time = UTCDateTime("1970-01-01T00:00:00.0Z") # set as the default in obspy
        origin.latitude = latitude
        origin.longitude = longitude
        origin.depth = depth_in_m
        origin.depth_type = "operator assigned"
        origin.evaluation_mode = "manual"
        origin.evaluation_status = "preliminary"
        origin.region = FlinnEngdahl().get_region(origin.longitude, origin.latitude)
        
        focal_mechanisms = FocalMechanism()
        tensor = Tensor()
        moment_tensor = MomentTensor()

        tensor.m_rr = m_rr
        tensor.m_tt = m_tt
        tensor.m_pp = m_pp
        tensor.m_rt = m_rt
        tensor.m_rp = m_rp
        tensor.m_tp = m_tp
        
        moment_tensor.tensor = tensor
        focal_mechanisms.moment_tensor = moment_tensor
                            
        # make associations, put everything together
        cat.append(event)
        event.origins = [origin]
        event.focal_mechanisms = [focal_mechanisms]

    elif SOURCE_TYPE == 'FINITE_SOURCE':
        pass
        
else:

    if SOURCE_TYPE == 'POINT_SOURCE' or SOURCE_TYPE == 'GREENS':
        # Point Source from quakeml file~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        source = instaseis.Source.parse(path_to_quakeml_cat)

    elif SOURCE_TYPE == 'FINITE_SOURCE':
        # Finite source from USGS.param file~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        finite_source = instaseis.FiniteSource.from_usgs_param_file(path_to_usgsparam)
        nsamp = int(db.info.period / finite_source[0].dt) * 50
        finite_source.resample_sliprate(dt=finite_source[0].dt, nsamp=nsamp)
        finite_source.lp_sliprate(freq=1.0 / db.info.period)
        finite_source.resample_sliprate(dt=db.info.dt, nsamp=db.info.npts)
        finite_source.compute_centroid()
        source = finite_source

##################
# CREATE INVENTORY
##################
if CREATE_NEW_INVENTORY is True:
    # Create new FDSN inventory~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    duration = 10 * 60  # no importance in our case because we do not need the response
    event_time = UTCDateTime("2011-03-11T05:46:23.2")  # instaseis does not care about time
    inv = c.get_stations(network=net,
                            station=sta,
                            location=loc,
                            channel=chn,
                            starttime=event_time,
                            endtime=event_time + duration, level="channel")

    receivers = instaseis.Receiver.parse(inv)

else:
    # Parse existing inventory to instaseis receiver object~~~~~~~~~~~~~~~~~~~
    inv = read_inventory(path_to_inv)
    receivers = instaseis.Receiver.parse(inv)


####################
# COMPUTE SEISMOGRAM
####################
if SOURCE_TYPE == 'POINT_SOURCE':
    # Use point source~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    flag = False
    for receiver in receivers:
        print('Computing for receiver: ', receiver.station)
        st_point = db.get_seismograms(
            source, receiver, components=('RTZ'))
        if flag is False:
            st = st_point.copy()
            st.clear()
            st.__iadd__(st_point)
            flag = True
        else:
            st.__iadd__(st_point)
elif SOURCE_TYPE == 'FINITE_SOURCE':
    # Use finite source~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    flag = False
    for receiver in receivers:
        print('Computing for receiver: ', receiver.station)
        st_finite = db.get_seismograms_finite_source(
            source, receiver, components=('RTZ'), dt=1.0)
        if flag is False:
            st = st_finite.copy()
            st.clear()
            st.__iadd__(st_finite)
            flag = True
        else:
            st.__iadd__(st_finite)
elif SOURCE_TYPE == "GREENS":
    # Base URL and model choice.
    BASE_URL = 'http://service.iris.edu/irisws/syngine/1/query?format=miniseed&model=' + model +  '&'
    # A: Get it directly from Syngine.
    flag = False
    for receiver in receivers:
        print('Computing for receiver: ', receiver.station)
        MT_URL = BASE_URL + (
            "sourcelatitude=%g&sourcelongitude=%g&sourcedepthinmeters=%g&"
            "sourcemomenttensor=%g,%g,%g,%g,%g,%g&"
            "receiverlatitude=%g&receiverlongitude=%g&"
            "components=ZRT"% (source.latitude, source.longitude, source.depth_in_m, 
                            source.m_rr, source.m_tt, source.m_pp, source.m_rt, source.m_rp, source.m_tp,
                            receiver.latitude, receiver.longitude,))
        st_green= obspy.read(MT_URL.replace("+", ""))
        for trace in st_green:
            trace.id = receiver.network + '.' + receiver.station + '.' + receiver.location + '.' + trace.id.split('.')[-1]
        if flag is False:
            st = st_green.copy()
            st.clear()
            st.__iadd__(st_green)
            flag = True
        else:
            st.__iadd__(st_green)

#####################
# CORRECT INVENTORIES
#####################
correct_inv = Inventory(
        # We'll add networks later.
        networks=[],
        source="get_instaseis_data")
for trace in st:
    # Go trace by trace and get correct trace info
    [network, station, location, channel] = trace.id.split('.')
    # Find net + sta in inv corresponding to the trace and get geographical info
    inv_select = inv.select(network=network, station=station)
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


##################
# PLOT SEISMOGRAMS
##################
if len(st) > 10:
    ans = input('There are more than 10 traces. This may result in a memory error while plotting. \n \
                Do you still wish to plot the first 10 traces? (y/n): ')
    if ans == 'y':
        st[0:10].plot()
        plt.show()
else:
    st.plot()
    plt.show()  


##############
# SAVE STREAMS
##############
value = input('Save it? (y/n): ')
if value == 'y':
    name = input('Name: ')
    os.mkdir('./instaseis_data/' + name)
    st.write('./instaseis_data/' + name + '/' + name + '_INSTA' + '.mseed')
    correct_inv.write('./instaseis_data/' + name + '/' + name + '_INSTA_inv' + '.xml', format="STATIONXML")
    if CREATE_NEW_SOURCE is True:
        cat.write('./instaseis_data/' + name + '/' + name + '_INSTA_cat' + '.xml', format='QUAKEML')
    else:
        if SOURCE_TYPE == 'POINT_SOURCE':
            shutil.copy(path_to_quakeml_cat, './instaseis_data/' + name + '/' + name + '_INSTA_cat' + '.xml')
        else:
            shutil.copy(path_to_usgsparam, './instaseis_data/' + name + '/' + name + '_INSTA' + '.param')
