
import math
import numpy as np
import netCDF4 as nc #read .nc files with weather data
from datetime import datetime, timedelta, date
from route_handling import mac_windows_file_handle


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

filename_AIS = mac_windows_file_handle("env/input_files/ais_data_v4.csv")
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


#NW_file = mac_windows_file_handle("Weather_Data/NW_01_07_2020__01_07_2022.nc") #Keys =['northward_wind', 'time', 'lat', 'lon']
#EW_file = mac_windows_file_handle("Weather_Data/EW_01_07_2020__01_07_2022.nc")  #Keys =['eastward_wind', 'time', 'lat', 'lon']
NW_file_20_21 = mac_windows_file_handle("Weather_Data/NW_01.07.2020_01.07.2021.nc") #Weather data file NW part 1
NW_file_21_22 = mac_windows_file_handle("Weather_Data/NW_01.07.2021_01.07.2022.nc") #Weather data file NW part 2
EW_file_20_21 = mac_windows_file_handle("Weather_Data/EW_01.07.2020_01.07.2021.nc") #Weather data file EW part 1
EW_file_21_22 = mac_windows_file_handle("Weather_Data/EW_01.07.2021_01.07.2022.nc") #Weather data file EW part 2

dataset_NW_1 = nc.Dataset(NW_file_20_21)
dataset_NW_2 = nc.Dataset(NW_file_21_22)
dataset_EW_1 = nc.Dataset(EW_file_20_21)
dataset_EW_2 = nc.Dataset(EW_file_21_22)

## Test weather

file_Sleipnir_A  = mac_windows_file_handle("Weather_Data/Platform Check/Sleipner-A_2021_03.nc") #WEather Data taken from Sleipnir A
file_Troll_A     = mac_windows_file_handle("Weather_Data/Platform Check/Troll-A_2020_07.nc") #Weather data taken from Troll A

dataset_Sleipnir_A  = nc.Dataset(file_Sleipnir_A)
dataset_Troll_A     = nc.Dataset(file_Troll_A)

##
#print(dataset_Sleipnir_A)
#print(dataset_Troll_A)

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
#print(dataset_test)
###for var in dataset_test.variables.values():
#    print(var)
#print(dataset_EW)
#for var in dataset_NW.variables.values():
#    print(var)#
#    print(dataset_test_1.variables.keys())
#print(dataset_test_1.variables["time"])
#print(dataset_test_2.variables["time"])
#print(dataset_test_3.variables["time"])
#print("time first ", dataset_test_1.variables["time"][0])
#print("time last ", dataset_test_1.variables["time"][-1])

#print("time first of 2 ", dataset_test_2.variables["time"][0])
#print("time last of 2", dataset_test_2.variables["time"][-1])
#print("lon 1", dataset_test_1.variables["lon"][:])
#print("lat 1", dataset_test_1.variables["lat"][:])
#print(len(dataset_test_1.variables))
#print(len(dataset_test_1.dimensions))
#print(dataset_test_1.dimensions)
#print(dataset_test_1["time"])
#print(dataset_test_1["x"])
#print(dataset_test_1["y"])
#print(dataset_test_1["height"])

# datasets have 4 dimentions: time, x, y, height.
#time is days since 2018-12-01 00:00:00
#increases by 0.0416 days per value, meaning increases by hour
#x is latitude
#y is longitude
#height is either 10,20,40--- up to 600, we only need 10 at position 0, so write 0 (the 0th position)


#dataset_NW["northward_wind"][:,lat_pos,lon_pos])
#print(dataset_test.variables["y"][:])
#print(dataset_test.variables["x"][:])

#print("north",dataset_NW_1.variables.keys())
#print("")
#print("East",dataset_EW.variables.keys())
northward_wind_1  = dataset_NW_1.variables["northward_wind"]
northward_time_1  = dataset_NW_1.variables["time"]
northward_lat_1   = dataset_NW_1.variables["lat"]
northward_lon_1   = dataset_NW_1.variables["lon"]
eastward_wind_1   = dataset_EW_1.variables["eastward_wind"]
eastward_time_1   = dataset_EW_1.variables["time"]
eastward_lat_1    = dataset_EW_1.variables["lat"]
eastward_lon_1    = dataset_EW_1.variables["lon"]

northward_wind_2  = dataset_NW_2.variables["northward_wind"]
northward_time_2  = dataset_NW_2.variables["time"]
northward_lon_2   = dataset_NW_2.variables["lon"]
northward_lat_2   = dataset_NW_2.variables["lat"]
eastward_wind_2   = dataset_EW_2.variables["eastward_wind"]
eastward_time_2   = dataset_EW_2.variables["time"]
eastward_lat_2    = dataset_EW_2.variables["lat"]
eastward_lon_2    = dataset_EW_2.variables["lon"]





#print("first 1 lat",northward_lat_1[0])
#print("last 1 lat",northward_lat_1[-1])
#print("first 1 lon",northward_lon_1[0])
#print("last 1 lon",northward_lon_1[-1])
#print("first 2 lat",northward_lat_2[0])
#print("last 2 lat",northward_lat_2[-1])
#print("first 2 lon",northward_lon_2[0])
#print("last 2 lon",northward_lon_2[-1])



#print(f"The size of Eastward time_1 is: {eastward_time_1.shape}")
#print(f"The size of Eastward time_2 is: {eastward_time_2.shape}")

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
#print(dataset_NW_1["northward_wind"][0,60,20])


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
    """
    :param rad: input angle radians
    :return: angle in degrees
    """
    degree = rad*180/np.pi
    return degree
#translates degrees to rads for pythong
def d2r(degree):
    """
    :param degree: input angle degrees
    :return: angle in radians
    """
    rad = np.pi*degree/180
    return rad


def getweather(tid,latitude, longditude):
    """
    :param tid: Time we want to access weather data
    :param latitude: Latitude where we want to access weather data
    :param longditude: Longdtitude where we want to access weather data
    :return: WSN and WSE in m/s at time = tid, and position (lat,lon)
    """

    if tid >= 17520:  #if end of year 2 is reached while sailing, weatherdata from the beginning of year one is used
        tid -= 17520

    if tid <= 8742:

        lat_pos = int((latitude-eastward_lat_1[0])*8)              #Access correct position in vector of north wind
        lon_pos = int((longditude-eastward_lon_1[0])*8)             #Access correct position in vector of east wind

        if latitude < eastward_lat_1[0] or latitude > eastward_lat_1[-1]:
            print("latitude out of bounds, latitude between 58.9375 and 70.1875")
            return 1
        elif longditude < eastward_lon_1[0] or longditude > eastward_lon_1[-1]: #været må være hentet på posisjonen longditude
            print("longditude out of bounds, longditude between 3.0625 and 20.9375")
            return 1
        elif lat_pos >= len(dataset_NW_1["northward_wind"][tid,:]):
            print("yay")
        elif lat_pos >= len(dataset_NW_1["northward_wind"][tid,:,lon_pos]):
            print(f"lat posistion is out of bound at {lat_pos} degrees")
            return 1
        elif lon_pos >= len(dataset_NW_1["northward_wind"][tid, lat_pos, :]):
            print(f"lon posistion is out of bound at {lon_pos} degrees")
            return 1

        WSN = dataset_NW_1["northward_wind"][tid,lat_pos,lon_pos]
        WSE = dataset_EW_1["eastward_wind"][tid,lat_pos,lon_pos]

    #if tid ==
    else:

        tid -=  8743        #indekserer tiden i andre filen fra start igjen

        lat_pos = int((latitude - eastward_lat_2[0]) * 8)  # Access correct position in vector of north wind
        lon_pos = int((longditude - eastward_lon_2[0]) * 8)  # Access correct position in vector of east wind

        WSN = dataset_NW_2["northward_wind"][tid, lat_pos, lon_pos]
        WSE = dataset_EW_2["eastward_wind"][tid, lat_pos, lon_pos]

    return WSN,WSE




def datetime_seconds(dtime):
    """
    :param dtime: inputs date during over 2 year period (for this dataset)
    :return: point of time in seconds from 01/07/2020. Gives posibility of accessing specific sailing speed at given date and time
    """
    dt = dtime
    epoch = date(1990,1,1)
    delta = (dt-epoch)
    return delta.total_seconds()


def True_wind_speed(WSN,WSE): 
    """
    :param WSN: Wind Speed North in m/s
    :param WSE: Wind Speed East in m/s
    :return: True wind speed in m/s
    """
    TWS = np.sqrt(WSN**2+WSE**2)

    return TWS

#Function that return true wind direction in degrees from wind speed north/east (WSN/WSE)
def True_wind_direction(vessel_heading,wind_speed_north,wind_speed_east):
    """
    :param vessel_heading: Vessel heading in degrees
    :param wind_speed_north: speed of wind in northward direction (negaitive means south)
    :param wind_speed_east: speed of wind in eastern direction (negative means west)
    :return: true wind direction [degrees]
    """
    wind_angle_rads = math.atan2(wind_speed_north,wind_speed_east)   # gives direction, 0 degrees equals east, 90 degrees = north ....
    wind_angle_degs = r2d(wind_angle_rads)
    wind_angle_degs = (wind_angle_degs + 360) % 360             #normalizes degrees

    true_wind_direction = wind_angle_degs-vessel_heading        #true wind direction in degrees

    true_wind_direction = (true_wind_direction + 360) % 360

    return true_wind_direction #in degrees

def Apparent_Wind_Speed(true_wind_speed, vessel_speed, true_wind_direction):
    """

    :param true_wind_speed: Wind speed in relation to vessel heading
    :param vessel_speed: speed vessel sails
    :param true_wind_direction: heading of wind in relation to vessel
    :return: apparend wind speed m/s [float]: speed of wind in relation to vessel speed and heading
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
    :return:    AWA: Apparent wind angle [degree]
    """
    AWA = math.acos((TWS ** 2 - AWS ** 2 - VS ** 2) / (-2 * AWS * VS))
    AWA = (AWA + 360)%360
    return r2d(AWA)


#Function that calculates AWA selfmade
def alpha(vessel_speed,vessel_heading, NWS,EWS):
    """
    :param vessel_speed: Speed of vessel
    :param vessel_heading: Heading of vessel
    :param NWS: Northern wind speed (decomposed)
    :param EWS: Eastern wind speed (decomposed)
    :return: Apparent wind angle in degrees
    """
    Vsx = np.cos(vessel_heading)*vessel_speed   #Decompose sailing speed
    Vsy = np.sin(vessel_heading)*vessel_speed   #Decompose sailing speed
    Wsx = EWS                                   #Decompose Wind speed
    Wsy = NWS                                   #Decompose Wind speed
    Vawx = Vsx + Wsx                            #Recombine speeds
    Vawy = Vsy + Wsy                            #Recombine speeds
    alpha_temp = r2d(math.atan2(Vawy,Vawx))          #Change to degrees, calculate angle
    alpha      = (alpha_temp+360) % 360                   #Normalize angle between 0 and 360 degrees

    return alpha


def add_hours_to_date(date,hours):
    """

    :param date: start date
    :param hours: 6 hour intervalls
    :return: datetime
    """

    date_start = date
    hours = hours*6
    # add the specified number of hours to the input date
    new_datetime = date_start + timedelta(hours=hours)

    return new_datetime
#print(add_hours_to_date(date(2020,7,2),4))

def add_seconds_to_date(date,seconds):
    date_start = date
    # add the specified number of hours to the input date
    new_datetime = date_start + timedelta(seconds=seconds)
    return new_datetime

def find_index(date_now):
    """

    :param date_now:
    :return: correct index of csv file from copernicus in situ
    """
    date_start      = datetime(1950,1,1)
    second_delta    = (date_now-date_start).total_seconds()
    #indeksert per 10. sekund, ergo dele på 360, ikke 3600
    day_delta       = second_delta/3600/24
    return day_delta

#starttime = datetime(2020,7,1,00,00,00)
#tid = 8785.57
#latitude = 46
#longitude = 126
#print(len(northward_lat_2))
#print(len(eastward_lat_2))
#print(len(northward_lon_2))
#print(len(eastward_lon_2))
#print(northward_time_2[:])
#print(add_hours_to_date(starttime,tid))
#WSN = dataset_NW_2["northward_wind"][tid, lat_pos, lon_pos]
#WSE = dataset_EW_2["eastward_wind"][tid, lat_pos, lon_pos]

#getweather(tid,latitude,longitude)
#for i in range(100000):
#   WSN = dataset_NW_2["northward_wind"][i, lat_pos, lon_pos]
#    WSE = dataset_EW_2["eastward_wind"][i, lat_pos, lon_pos]
#    print(f"works with tid =  {i}")

#WSN = dataset_NW_2["northward_wind"][41.5699, lat_pos, lon_pos]
#WSE = dataset_EW_2["eastward_wind"][41.5699, lat_pos, lon_pos]

#print("time dif time 1", northward_time_1[1]-northward_time_1[0])
#print("time dif time 2", northward_time_2[1]-northward_time_2[0])
#print("first time element 1",northward_time_1[0])
#print("last time element 1",northward_time_1[-1])
#print("first element time 2", northward_time_2[0])
#print("time dif end 1 start 2", northward_time_2[0]-northward_time_1[-1])
#print("XXXXXXXXXXXX")


#starttime = datetime(1990,1,1,00,00,00)

#starttime_data_2 = add_hours_to_date(starttime,(northward_time_1[-1]/3600))#
#print("start2",starttime_data_2)


#print("XXXXXXXXXXXX")
#starttime = northward_time_1[0]/3600
#starttime = datetime(1990,1,1,00,00,00)

#print("end of time 1 in hours",northward_time_1[-1]/3600)
#print("start of time 2 in hours",northward_time_2[0]/3600)
#print("start of time 1,",add_hours_to_date(starttime,northward_time_1[0]/3600))
#print(f"add {(northward_time_1[-1]-northward_time_1[0])/3600} hours")
#print("end of time 1,",add_hours_to_date(starttime,northward_time_1[-1]/3600))
#print("start of time 2,",add_hours_to_date(starttime,northward_time_2[0]/3600))
#print(f"add {(northward_time_2[-1]-northward_time_2[0])/3600} hours")
#print("end of time 2,",add_hours_to_date(starttime,northward_time_2[-1]/3600))

#print("add 8748 timer to start time 1,", add_hours_to_date(datetime(2020,7,1,18),8748))
#print("start of time 2,",add_hours_to_date(starttime,northward_time_2[0]/3600))
#print("end of time 2,",add_hours_to_date(starttime,northward_time_2[-1]/3600))

#starttime = datetime(2020,7,1,00,00,00)
#print(add_hours_to_date(starttime,8760))
def find_position(ten_mins,TIME):
    # create a masked array from the input TIME array
    time = np.ma.masked_invalid(TIME)

    # remove masked values from the time array
    time = time[~time.mask]

    position = np.where(time == ten_mins)
    try:
        index    = position[0][0]
    except (IndexError):
        #print("Could not find index for", ten_mins)
        index = np.abs(time-ten_mins).argmin()
        return index

    #print("The index of ",ten_mins, "ten minute intervalls in the weather matrix is:", index)

    return index
#print(dataset_test_2["x"][:])

#print(dataset_test_2["y"][:])

def get_Weather_date_index(datenow):
    start_date = datetime(2020,1,1)

    #finding number of hours between start date and current date
    index = datenow-start_date
    index = int(index.total_seconds()/3600)
    return index

get_Weather_date_index(datetime(2020,7,2,12))
#print("\n")
#print(f"the date 2.7.2020  is the {datetime_seconds(date(2020,7,2))} second")
#print(f"the date 2.7.2020  is the 579 element in the datasets")

#print(dataset_Troll_A.variables)


def find_average_weather_Troll():
    True_wind_speed_vector_Troll = []
    for i in range(1,31):
        for j in range(1,24):
            if i == 1 and j <4 :
                break
            index       = find_index(datetime(2020,7,i,j))
            position    = find_position(index,Time_T)
            WSN,WSE     = getweather(position,Troll_lat,Troll_lon)
            True_wind_speed_vector_Troll.append(True_wind_speed(WSN,WSE))

    return True_wind_speed_vector_Troll

def find_average_weather_Sleipnir():
    True_wind_speed_vector_Sleipnir = []
    for i in range(1,31):
        for j in range(1,24):
            if i == 1 and j <4 :
                break
            index       = find_index(datetime(2021,3,i,j))
            position    = find_position(index,Time_S)
            WSN,WSE     = getweather(position,Troll_lat,Troll_lon)
            True_wind_speed_vector_Sleipnir.append(True_wind_speed(WSN,WSE))
    return True_wind_speed_vector_Sleipnir



start_date_T    = datetime(1950,1,1,00,00,00)
Time_T          = dataset_Troll_A.variables["TIME"]
WSPD_T          = dataset_Troll_A.variables["WSPD"]
WSPD_T_QC       = dataset_Troll_A.variables["WSPD_QC"]
VDIR_T          = dataset_Troll_A.variables["WDIR"]
VDIR_T_QC       = dataset_Troll_A.variables["WDIR_QC"]

start_date_S    = datetime(1950,1,1,00,00,00)

Time_S          = dataset_Sleipnir_A.variables["TIME"]
WSPD_S          = dataset_Sleipnir_A.variables["WSPD"]
WSPD_S_QC       = dataset_Sleipnir_A.variables["WSPD_QC"]
VDIR_S          = dataset_Sleipnir_A.variables["WDIR"]
VDIR_S_QC       = dataset_Sleipnir_A.variables["WDIR_QC"]


#days_second_july    = find_index(datetime(2020, 7, 1, 12))
#days_tenth_march    = find_index(datetime(2021, 3, 10,12))
#second_july = find_position(days_second_july,Time_T)
#tenth_march = find_position(days_tenth_march,Time_S)
#print(len(Time_S))

Troll_lat = 60.64350
Troll_lon = 3.71930
Sleipnir_lat = 58.37110
Sleipnir_lon = 1.90910
#WSN_Troll,WSE_Troll = getweather(second_july,Troll_lat,Troll_lon)
#WSN_Sleipnir,WSE_Sleipnir = getweather(tenth_march,Sleipnir_lat,Sleipnir_lon)


#ind = find_index(datetime(2021,3,10,12))

#WSN2,WSE2 = getweather(find_index(datetime(2020,7,1,12)),Sleipnir_lat,Sleipnir_lon)
#WSN1,WSE1 = getweather(find_index(datetime(2021,3,10,12)),Troll_lat,Troll_lon)
#print("True wind speed is ", True_wind_speed(WSN1,WSE1),"at Troll")
#print("True wind speed is ", True_wind_speed(WSN2,WSE2),"at Sleipner")
#print("True wind direction is ",True_wind_direction(0,WSN1,WSE1),"at troll")#
#print("True wind direction is ",True_wind_direction(0,WSN2,WSE2),"at sleipner")#

def get_old_weather_vect_Sleipner():
    Vect_Sleipnir_Old = []
    for i in range(1, 31):
        for j in range(1, 24):
            index = get_Weather_date_index(datetime(2020,7,i,j))
            WSN, WSE = getweather(index, Sleipnir_lat, Sleipnir_lon)
            Vect_Sleipnir_Old.append(True_wind_speed(WSN, WSE))

    return Vect_Sleipnir_Old


def get_old_weather_vect_Troll():
    Vect_Troll_Old = []
    for i in range(1, 31):
        for j in range(1, 24):
            index = get_Weather_date_index(datetime(2021,3,i,j))
            WSN, WSE = getweather(index, Troll_lat, Troll_lon)
            Vect_Troll_Old.append(True_wind_speed(WSN, WSE))
    return Vect_Troll_Old

Vect_Troll = find_average_weather_Troll()
Vect_Troll_Old = get_old_weather_vect_Troll()
Vect_Sleipner = find_average_weather_Sleipnir()
Vect_Sleipner_Old = get_old_weather_vect_Sleipner()



print("Average wind Troll new june of 2020 ", np.average(Vect_Troll))
print("Average wind Troll old june of 2020", np.average(Vect_Troll_Old))
print("Average wind Sleipnir new march of 2021", np.average(Vect_Sleipner))
print("Average wind Sleipnir old march of 2021", np.average(Vect_Sleipner_Old))