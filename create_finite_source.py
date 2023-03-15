""" 
Takes in either a manual list of point sources or a CMT solutions file
and parses it into a format that can be placed in axisem3D's inparam.source.yaml
"""
import numpy as np

# create function to format individual sources (finite rupture = collection of point sources approximation)

def format_input_point(i, input_data):
    input_string =     f"""
    - point{i}:
        location:
            latitude_longitude: [{input_data["lat"]}, {input_data["lon"]}]
            depth: {input_data["depth"]}
            ellipticity: false
            depth_below_solid_surface: true
            undulated_geometry: true
        mechanism:
            type: MOMENT_TENSOR
            data: [{input_data["MT"][0]:.3e}, {input_data["MT"][1]:.3e}, {input_data["MT"][2]:.3e}, {input_data["MT"][3]:.3e}, {input_data["MT"][4]:.3e}, {input_data["MT"][5]:.3e}]
            unit: {input_data["MT_unit"]:.3e}
        source_time_function:
            class_name: GaussianSTF
            half_duration: {input_data["hd"]:.3e}
            decay_factor: 1.628
            time_shift: {input_data["ts"]:.3e}
            use_derivative_integral: GAUSSIAN"""
    return input_string

TYPE = 'CMTSOLUTION'

if TYPE == 'MANUAL':
    # MANUAL------------------------------------------------------------------------------------------------------------------------
    # list of sources in the finite rupture; just an example shown here below but the true inparam.source.yaml shows an actual finite rupture 
    input_data_list = [{"lat":1,"lon":1, "depth":1, "MT":np.random.normal(0,1,6), "hd":2.0,"ts":-1, "MT_unit":1.0}, 
                    {"lat":1,"lon":1, "depth":1, "MT":np.random.normal(0,1,6), "hd":2.0,"ts":-1, "MT_unit":1.0}]
elif TYPE == 'CMTSOLUTION':
    # FROM USGS CMTSOLUTION FILE----------------------------------------------------------------------------------------------------------
    input_data_list = []
    index = 0
    f = open("/home/adrian/PhD/axisem3d_root/AxiSEM3D-master/cmb/real_sources/CMTSOLUTION_SUMATRA_2009", "r")
    for line in f:
        key = line.split(' ')[0]
        if key == 'event':
            event_name = line.split(' ')[-1][0:-1]
        elif key == 'latitude:':
            lat = float(line.split(' ')[-1][0:-1])
        elif key == 'longitude:':
            lon = float(line.split(' ')[-1][0:-1])
        elif key == 'depth:':
            depth = float(line.split(' ')[-1][0:-1])*1e3
        elif key == 'Mrr:':
            Mrr = float(line.split(' ')[-1][0:-1])
        elif key == 'Mtt:':
            Mtt = float(line.split(' ')[-1][0:-1])
        elif key == 'Mpp:':
            Mpp = float(line.split(' ')[-1][0:-1])
        elif key == 'Mrt:':
            Mrt = float(line.split(' ')[-1][0:-1])
        elif key == 'Mrp:':
            Mrp = float(line.split(' ')[-1][0:-1])
        elif key == 'Mtp:':
            Mtp = float(line.split(' ')[-1][0:-1])
        elif key == 'half':
            hd = float(line.split(' ')[-1][0:-1])
        elif key == 'time':
            ts = float(line.split(' ')[-1][0:-1])
        elif key == '\n':
            if index == 0:
                index += 1
            else:
                index += 1
                dictionary = {"lat":lat,"lon":lon, "depth":depth, "MT":[Mrr, Mtt, Mpp, Mrt, Mrp, Mtp], "hd":hd,"ts":ts, "MT_unit":0.0000001}
                input_data_list.append(dictionary)



# contents of the file must be manually copied into inparam.source.yaml under list_of_sources:

with open('SUMATRA_2009.txt', 'w') as file:
    for i, source in enumerate(input_data_list):
        file.write(format_input_point(i+1, source))
        