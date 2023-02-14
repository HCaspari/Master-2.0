
import math
import numpy as np
import netCDF4 as nc #read .nc files with weather data
from datetime import datetime, timedelta, date

#Input stats
mean_wind_speed = 10 #knots
mean_wind_direction = 90 #degrees
time_intervall =0.10 #hours
mesh_size = 10 #nautical miles
Start_east = 0
Start_north = 0
Start_position = (Start_east,Start_north) #Current Position
GlobalPositionVect = [(0,0)]
clock = 0
filename_AIS = "../env/input_files/ais_data_v4.csv"
travel_iteration            = 0
#vessel parameters:
vessel_length               = 101.26
vessel_draft                = 10.15
vessel_weight_disp          = 8856
#vessel_velocity             = 1*0.514444444
rho_air                     = 1.025
rho_water                   = 1025
rotor_amount = 4
h   = 35    #height flettner
d   = 5     #diameter of flettner
A   = h*d  #cross sectional area of flettner
Cl  = 12.5
Cd  = 0.2
Cm  = 0.2
alpha = 3.5

#datafiles

NW_file = "Weather_Data/NW_01_07_2020__01_07_2022.nc"  #Keys =['northward_wind', 'time', 'lat', 'lon']
EW_file = "Weather_Data/EW_01_07_2020__01_07_2022.nc"  #Keys =['eastward_wind', 'time', 'lat', 'lon']
dataset_NW = nc.Dataset(NW_file)
dataset_EW = nc.Dataset(EW_file)


##################
    #geospatial_lat_min: 58.9375
    #geospatial_lat_max: 70.1875
    #geospatial_lat_resolution: 0.125
    #geospatial_lat_units: degrees_north
    #geospatial_lon_min: 3.0625
    #geospatial_lon_max: 20.9375
    #geospatial_lon_resolution: 0.125
    #geospatial_lon_units: degrees_east

##################
#handling datafiles:

#print(dataset_EW)
#for var in dataset_NW.variables.values():
#    print(var)
#print("north",dataset_NW.variables.keys())
#print("East",dataset_EW.variables.keys())
northward_wind  = dataset_NW.variables["northward_wind"]
northward_time  = dataset_NW.variables["time"]
northward_lat   = dataset_NW.variables["lat"]
northward_lon   = dataset_NW.variables["lon"]

eastward_wind  = dataset_EW.variables["eastward_wind"]
eastward_time  = dataset_EW.variables["time"]
eastward_lat   = dataset_EW.variables["lat"]
eastward_lon   = dataset_EW.variables["lon"]
#print(northward_lon[:])
#print(northward_lat[:])


#print(f"The size of Eastward time is: {eastward_time.shape}")

#print(f"The first and last elements of Eastward time are {eastward_time[0]}, and {eastward_time[-1]}\n"
#      f"this gives {eastward_time[-1]-eastward_time[0]} seconds")
#print(f"given that start time is 01/01/1990 00:00:00, then adding 962 409 600 seconds \n"
#      f"gives the start date {date(1990,1,1)+ timedelta(days=11139)} "
#      f"and the end date {date(1990,1,1)+ timedelta(days=11869.9583333)}")
#print(f"The first and last elements of Eastward lat are {eastward_lat[0]}, and {eastward_lat[-1]}")
#print(f"The first and last elements of Eastward lon are {eastward_lon[0]}, and {eastward_lon[-1]}")
#print(f"The first and last elements of Eastward time are {eastward_time[0]}, and {eastward_time[-1]}"
#      f"it is {eastward_time[-1] - eastward_time[0]} elements")
#print(f"The second element of Eastward time is {eastward_time[1]}")
#print(f"The increment of Eastward time is {eastward_time[1]-eastward_time[0]} which is {(eastward_time[1]-eastward_time[0])/60}  minutes or "
#      f"{(eastward_time[1]-eastward_time[0])/3600} hours")
#
#print(f"The size of latitudes is: {eastward_lat.shape}")
#print(f"The size of longditudes  is: {eastward_lon.shape}")
#print(dataset_NW["northward_wind"][0,60,20])

#print(f"hour since 1990{date(1990,1,1)+ timedelta(days=42731850)}")




#dist_vect,travel_time,latitudes_vect,longditudes_vect,heading_vect = read_cols(filename_AIS)
#vector_of_positions(latitudes_vect,longditudes_vect)

#data = download_data() #Downloads data locally

#calculates weather at time=tid, latitude and londitude of trip
#def getweather(tid,latitude, longditude):
    #to get the correct index one needs to increase values by 0.125
    #for example, latitude 58.9375 = latitude[0]
#    lat_pos_north = round((latitude-dataset_NW["lat"][0])*8)/8
#    lon_pos_north = round((longditude-dataset_NW["lon"][0])*8)/8
#    lat_pos_north = dataset_NW.variables["lat"][0]              #Access correct position in vector of north wind
#    lon_pos_north = int((longditude-dataset_NW["northward_wind"]["lon"][0])*8)             #Access correct position in vector of east wind
#    lat_pos_east = int((longditude-dataset_EW["eastward_wind"]["lat"][0])*8)             #Access correct position in vector of east wind
#    lon_pos_east = int((longditude-dataset_EW["eastward_wind"]["lon"][0])*8)             #Access correct position in vector of east wind


    #tid         += 962409600 #measurements go from 01/07/2020 - to 01/07/2022
    #latitude    = round(latitude*8)/8
    #longditude  = round(longditude*8)/8
    #if tid < 0 or tid >= len(dataset_NW["northward_wind"][:,lat_pos,lon_pos]):
    #    print(f"time was out of bounds at {tid}, time must be within the span of one year (less than 1460)")
    #    return 1
#    if latitude < 58.9375 or latitude > 70.1875:
#        print("latitude out of bounds, latitude between 58.9375 and 70.1875")
#        return 1
#    elif longditude < 3.0625 or longditude > 20.9375: #været må være hentet på posisjonen longditude
#        print("longditude out of bounds, longditude between 3.0625 and 20.9375")
#        return 1
#    elif lat_pos_north >= len(dataset_NW["northward_wind"][tid,:,lon_pos_north]):
#        print(f"lat posistion is out of bound at {lat_pos_north} degrees")
#        return 1
#    elif lon_pos_north >= len(dataset_NW["northward_wind"][tid, lat_pos_north, :]):
#        print(f"lon posistion is out of bound at {lon_pos_north} degrees")
#        return 1#
#    WSN = dataset_NW["northward_wind"][tid,lat_pos_north,lon_pos_north]
#    WSE = dataset_EW["eastward_wind"][tid,lat_pos_east,lon_pos_east]


#    return WSN,WSE
#print(northward_lat[0])
#print(northward_lon[0])
#print(northward_lat[:])
#print(northward_lon[:])
#print(len(eastward_time[:]))

#translates radians to degrees for python
def r2d(rad):
    degree = rad*180/np.pi
    return degree
#translates degrees to rads for pythong
def d2r(degree):
    rad = np.pi*degree/180
    return rad

def getweather(tid,latitude, longditude):
    #to get the correct index one needs to increase values by 0.125
    #for example, latitude 58.9375 = latitude[0]
    lat_pos = int((latitude-eastward_lat[0])*8)              #Access correct position in vector of north wind
    lon_pos = int((longditude-eastward_lon[0])*8)             #Access correct position in vector of east wind
    #tid += 962409600 #measurements go from 01/07/2020 - to 01/07/2022
    #if tid < 0 or tid >= len(dataset_NW["northward_wind"][:,lat_pos,lon_pos]):
    #    print(f"time was out of bounds at {tid}, time must be within the span of one year (less than 1460)")
    #    return 1
    if latitude < eastward_lat[0] or latitude > eastward_lat[-1]:
        print("latitude out of bounds, latitude between 58.9375 and 70.1875")
        return 1
    elif longditude < eastward_lon [0] or longditude > eastward_lon[-1]: #været må være hentet på posisjonen longditude
        print("longditude out of bounds, longditude between 3.0625 and 20.9375")
        return 1
    elif lat_pos >= len(dataset_NW["northward_wind"][tid,:,lon_pos]):
        print(f"lat posistion is out of bound at {lat_pos} degrees")
        return 1
    elif lon_pos >= len(dataset_NW["northward_wind"][tid, lat_pos, :]):
        print(f"lon posistion is out of bound at {lon_pos} degrees")
        return 1
    elif tid >= eastward_time[-1]: #Denne fanger opp situasjoner nær slutten av året der vi mangler værdata for et nytt år. her vil
        #tiden vi henter værdatafra loope slik at istedenfor å PRØVE å hente været for 02/07/2021, så vil den loope tilbake og finne været for
        #den 02/07/2020.
        tid -= eastward_time[-1]

    WSN = dataset_NW["northward_wind"][tid,lat_pos,lon_pos]
    WSE = dataset_EW["eastward_wind"][tid,lat_pos,lon_pos]

    return WSN,WSE

def datetime_seconds(dtime):
    dt = dtime
    epoch = date(2020,7,1)
    delta = (dt-epoch)
    return delta.total_seconds()

def True_wind_speed(WSN,WSE):
    TWS = np.sqrt(WSN**2+WSE**2)
    return TWS

#Function that return true wind direction in degrees from wind speed north/east (WSN/WSE)
def True_wind_direction(vessel_heading,wind_speed_north,wind_speed_east):
    """
    :param vessel_heading: Vessel heading
    :param wind_speed_north: speed of wind in northward direction (negaitive means south)
    :param wind_speed_east: speed of wind in eastern direction (negative means west)
    :return: true wind direction in radians
    """
    wind_angle = math.atan2(wind_speed_north,wind_speed_east)
    true_wind_direction = wind_angle-vessel_heading

    return true_wind_direction

def Apparent_Wind_Speed(true_wind_speed, vessel_speed, true_wind_direction):
    """

    :param true_wind_speed: Wind speed in relation to vessel heading
    :param vessel_speed: speed vessel sails
    :param true_wind_direction: heading of wind in relation to vessel
    :return: apparend wind speed: speed of wind in relation to vessel speed and heading
    """

    #AWS_func = TWS_func - sailing_speed_func * np.sin(np.pi / 180 * sailing_direction_func)
    #from Seddiek et al. function 1 and 2
    AWS = np.sqrt(true_wind_speed ** 2 + vessel_speed**2-2*true_wind_speed*vessel_speed*np.cos(true_wind_direction))

    return AWS

#function that gives AWA by Seddiek
def Apparent_wind_angle(TWS, AWS, VS):
    """

    :param      TWS: Speed of wind in relation to global axis
    :param      AWS: speed of wind in relation to vessel speed and heading
    :param      VS: speed of vessel
    :return:    AWA: Apparent wind angle
    """
    AWA = math.acos((TWS ** 2 - AWS ** 2 - VS ** 2) / (-2 * AWS * VS))
    if AWA >= 360:
        AWA -= 360
    return r2d(AWA)

#Function that calculates AWA selfmade :D
def alpha(vessel_speed,vessel_heading, NWS,EWS):
    """
    :param vessel_speed: Speed of vessel
    :param vessel_heading: Heading of vessel
    :param NWS: Northern wind speed (decomposed)
    :param EWS: Eastern wind speed (decomposed)
    :return: Apparent wind angle in degrees
    """
    Vsx = np.cos(vessel_heading)*vessel_speed
    Vsy = np.sin(vessel_heading)*vessel_speed
    Wsx = EWS
    Wsy = NWS
    Vawx = Vsx + Wsx
    Vawy = Vsy + Wsy
    alpha = math.atan2(Vawy,Vawx)
    return r2d(alpha)
