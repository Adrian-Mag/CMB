""" 
Simple template for downloading data from IRIS, including 
event, station, and wavefield data.
"""

from obspy import UTCDateTime, read
from obspy.clients.fdsn import Client
import matplotlib
matplotlib.use('tkagg')
import matplotlib.pyplot as plt
import os 

############
# EVENT INFO
############
c_event = Client("IRIS")
# Event time
event_time = UTCDateTime("2009-09-30T10:16:09")
# Get the event information. The temporal and magnitude constraints make it unique
cat = c_event.get_events(starttime=event_time - 10, endtime=event_time + 10,
                         minmagnitude=7.5, maxmagnitude=7.7, minlatitude=-1, 
                         maxlatitude=0, minlongitude=99, maxlongitude=100)


##############
# STATION INFO
##############
c = Client("IRIS")
net = "LB,CN,UW"
sta = "BMN,DRLN,BLOW"
loc = "*"
chn = "BH?"
pre_duration = 0
duration = 30*60
inv = c.get_stations(network=net, station=sta, location=loc, channel=chn,
                        starttime=event_time - pre_duration, endtime=event_time + duration,
                        level="response")
st = c.get_waveforms(network=net, station=sta, location=loc,
                        channel=chn, starttime=cat[0].origins[0].time - pre_duration,
                        endtime=event_time + duration)


######################
# PLOT AND INFO OUTPUT
######################
print(cat)
print(inv)
print(st)
st.plot()

plt.show()

######
# SAVE
######
value = input('Save it? : ')
if value=='y':
    name = input('Name: ')
    os.mkdir('./real_seismic_data/' + name)
    st.write('./real_seismic_data/' + name + '/' + name + '.mseed')
    cat.write('./real_seismic_data/' + name + '/' + name + '_cat' + '.xml', format="QUAKEML")
    inv.write('./real_seismic_data/' + name + '/' + name + '_inv' + '.xml', format="STATIONXML") 
