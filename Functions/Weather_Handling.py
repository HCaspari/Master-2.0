
import math
import numpy as np
import netCDF4 as nc #read .nc files with weather data
from datetime import datetime, timedelta, date
from route_handling import mac_windows_file_handle
from plot_functions import plot_something, plot_histogram, plot_Vect_Daily,plot_Vect_Weekly,plot_Vect_hourly, plot_Vect_hourly_single
from file_handling import write_to_file, combine_troll,combine_Sleipnir
import matplotlib.pyplot as plt

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
Troll_lat = 60.64350
Troll_lon = 3.71930
Sleipnir_lat = 58.37110
Sleipnir_lon = 1.90910

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


#Clean Data
start_date_T    = datetime(1950,1,1,00,00,00)
Time_T          = dataset_Troll_A.variables["TIME"]
Time_T_masked   = np.ma.masked_invalid(Time_T[:])
Time_T_clean    = Time_T_masked[~Time_T_masked.mask]
Time_T_clean    = np.insert(Time_T_clean,0, 25749.0 )

WSPD_T          = dataset_Troll_A.variables["WSPD"]
WSPD_T_masked   = np.ma.masked_invalid(WSPD_T[:])
WSPD_T_clean    = WSPD_T_masked[~WSPD_T_masked.mask]
WSPD_T_clean    = np.insert(WSPD_T_clean, 0,WSPD_T_clean[0] )
#insert copy of first value into beginning of dataset for precision

WDIR_T        = dataset_Troll_A.variables["WDIR"]#from true north
WDIR_T_masked = np.ma.masked_invalid(WDIR_T[:])
WDIR_T_clean  = WDIR_T_masked[~WDIR_T_masked.mask]
WDIR_T_clean = np.insert(WDIR_T_clean, 0,WDIR_T_clean[0] )

start_date_S    = datetime(1950, 1, 1, 00, 00, 00)
Time_S          = dataset_Sleipnir_A.variables["TIME"]
Time_S_masked   = np.ma.masked_invalid(Time_S[:])
Time_S_clean    = Time_S_masked[~Time_S_masked.mask]

WSPD_S = dataset_Sleipnir_A.variables["WSPD"]
WSPD_S_masked = np.ma.masked_invalid(WSPD_S[:])
WSPD_S_clean = WSPD_S_masked[~WSPD_S_masked.mask]
#insert copy of first value into beginning of dataset for precision

# WSPD_S_QC       = dataset_Sleipnir_A.variables["WSPD_QC"]
WDIR_S        = dataset_Sleipnir_A.variables["WDIR"]  # from true north
WDIR_S_masked = np.ma.masked_invalid(WDIR_S[:])
WDIR_S_clean  = WDIR_S_masked[~WDIR_S_masked.mask]







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

def add_days_to_date(days):
    date = datetime(1950,1,1)
    new_datetime = date + timedelta(days=days)
    return new_datetime

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
    day_delta       = second_delta/3600/24

    return day_delta

def find_position(ten_mins,TIME):
    # create a masked array from the input TIME array
    time = np.ma.masked_invalid(TIME)

    # remove masked values from the time array
    time = time[~time.mask]

    position = np.where(time == ten_mins)
    try:
        index    = position[0][0]
    except (IndexError):
        index = np.abs(time-ten_mins).argmin()
        return index

    return index

def get_Weather_date_index(datenow):
    start_date = datetime(2020,1,1)

    #finding number of hours between start date and current date
    index = datenow-start_date
    index = int(index.total_seconds()/3600)
    return index

def find_average_weather_Troll_Measured():
    True_wind_speed_vector_Troll     = []
    True_wind_direction_vector_Troll = []
    for i in range(1,31):
        if i == 1:
            for j in range(1, 24):
                index = find_index(datetime(2020, 7, i, j))
                position = find_position(index, Time_T_clean)
                True_wind_speed_vector_Troll.append(WSPD_T_clean[position])
                True_wind_direction_vector_Troll.append(WDIR_T_clean[position])
        if i > 1:
            for j in range(1,24):
                index = find_index(datetime(2020,7,i,j))
                position = find_position(index, Time_T_clean)
                True_wind_speed_vector_Troll.append(WSPD_T_clean[position])
                True_wind_direction_vector_Troll.append(WDIR_T_clean[position])
    return True_wind_speed_vector_Troll, True_wind_direction_vector_Troll

def find_average_weather_Sleipnir_Measured():
    True_wind_speed_vector_Sleipnir = []
    True_wind_direction_vector_Sleipnir = []
    for i in range(1,31):
        if i == 1:
            for j in range(0, 24):
                index       = find_index(datetime(2021,3,i,j))
                position    = find_position(index,Time_S_clean)
                True_wind_speed_vector_Sleipnir.append(round(WSPD_S_clean[position],3))
                True_wind_direction_vector_Sleipnir.append(round(WDIR_S_clean[position],3))
        if i >= 1:
            for j in range(1, 24):
                index       = find_index(datetime(2021, 3, i, j))
                position    = find_position(index, Time_S_clean)
                True_wind_speed_vector_Sleipnir.append(round(WSPD_S_clean[position],3))
                True_wind_direction_vector_Sleipnir.append(round(WDIR_S_clean[position],3))
    return True_wind_speed_vector_Sleipnir, True_wind_direction_vector_Sleipnir

def get_calculated_weather_vect_Sleipner():
    WSPD_Vect_Sleipnir_Old = []
    WDir_Vect_Sleipner_Old = []
    for i in range(1, 31):
        for j in range(1, 24):
            index = get_Weather_date_index(datetime(2020,7,i,j))
            WSN, WSE = getweather(index, Sleipnir_lat, Sleipnir_lon)
            WSPD_Vect_Sleipnir_Old.append(True_wind_speed(WSN, WSE))
            WDir_Vect_Sleipner_Old.append(True_wind_direction(0, WSN, WSE))
    return WSPD_Vect_Sleipnir_Old,WDir_Vect_Sleipner_Old

def get_calculated_weather_vect_Troll():
    WSPD_Vect_Troll_Old = []
    WDir_Vect_Troll_Old = []
    for i in range(1, 31):
        for j in range(1, 24):
            index = get_Weather_date_index(datetime(2021,3,i,j))
            WSN, WSE = getweather(index, Troll_lat, Troll_lon)
            WSPD_Vect_Troll_Old.append(True_wind_speed(WSN, WSE))
            WDir_Vect_Troll_Old.append(True_wind_direction(0,WSN,WSE))
    return WSPD_Vect_Troll_Old,WDir_Vect_Troll_Old

def get_calculated_weather_vect_Troll_year(year):
    WSPD_Vect_Troll_Old_year = []
    WDir_Vect_Troll_Old_year = []
    feb = 28
    if year % 4 == 0 or year == 2000:
        feb = 29
    for k in range(1,13):
        if k == 2:
            for i in range(1, feb):
                for j in range(1, 24):
                    index = get_Weather_date_index(datetime(year, k, i, j))
                    WSN, WSE = getweather(index, Troll_lat, Troll_lon)
                    WSPD_Vect_Troll_Old_year.append(True_wind_speed(WSN, WSE))
                    WDir_Vect_Troll_Old_year.append(True_wind_direction(0, WSN, WSE))
        elif k%2 == 0:
            for i in range(1, 31):
                for j in range(1, 24):
                    index = get_Weather_date_index(datetime(year,k,i,j))
                    WSN, WSE = getweather(index, Troll_lat, Troll_lon)
                    WSPD_Vect_Troll_Old_year.append(True_wind_speed(WSN, WSE))
                    WDir_Vect_Troll_Old_year.append(True_wind_direction(0,WSN,WSE))
        elif k%2 == 1:
            for i in range(1, 30):
                for j in range(1, 24):
                    index = get_Weather_date_index(datetime(year,k,i,j))
                    WSN, WSE = getweather(index, Troll_lat, Troll_lon)
                    WSPD_Vect_Troll_Old_year.append(True_wind_speed(WSN, WSE))
                    WDir_Vect_Troll_Old_year.append(True_wind_direction(0,WSN,WSE))
    return WSPD_Vect_Troll_Old_year,WDir_Vect_Troll_Old_year


def Hourly_Concat(ten_mins):
    # Create a list of indices corresponding to each day
    indices = [i for i in range(0, len(ten_mins), 6)]

    # Calculate the average for each day
    hourly_avg = [sum(ten_mins[i:i + 6]) / 6 for i in indices]

    return hourly_avg

def Time_Concat(ten_mins):
    # Create a list of indices corresponding to each day
    indices = [i for i in range(0, len(ten_mins), 6)]

    # Calculate the average for each day
    hourly_avg = [ten_mins[i] for i in indices]

    return hourly_avg

WSPD_Vect_Sleipner_Measured, WDIR_Vect_Sleipner_Measured = find_average_weather_Sleipnir_Measured()
WSPD_Vect_Sleipner_Calculated, WDIR_Vect_Sleipner_Calculated = get_calculated_weather_vect_Sleipner()

WSPD_Vect_Troll_Measured, WDIR_Vect_Troll_Measured = find_average_weather_Troll_Measured()
WSPD_Vect_Troll_Calculated, WDIR_Vect_Troll_Calculated = get_calculated_weather_vect_Troll()

a = []
for i in range(31):
    a.append(i)
#plot_something("Sleipner_Calculated", WSPD_Vect_Sleipner_Calculated,"Day","WSPD")
#plot_something("Sleipner_Measured", WSPD_Vect_Sleipner_Measured,"Day","WSPD")
plot_something("Troll Measured", WSPD_Vect_Troll_Measured,"Hour","WSPD")


#print("Average wind Troll new july of 2020 ", np.average(WSPD_Vect_Troll_Measured))
#print("Average wind Troll old july of 2020", np.average(WSPD_Vect_Troll_Calculated))
#print("Average wind Sleipnir new march of 2021", np.average(WSPD_Vect_Sleipner_Measured))
#print("Average wind Sleipnir old march of 2021", np.average(WSPD_Vect_Sleipner_Calculated))

x_axis = []
y_axis = []
for i in range(31):
    x_axis.append(i)


#print(add_days_to_date(Time_S_clean[0]))
#print(add_days_to_date(Time_S_clean[-1]))

#print(add_days_to_date(Time_T_clean[0]))
#print(add_days_to_date(Time_T_clean[-1]))

#plot_Vect_Daily(WSPD_Vect_Sleipner_Calculated,WSPD_Vect_Sleipner_Measured,"Day","WSPD", "WSPD Calculated", " WSPD Measured","Sleipner March 2021")
#plot_Vect_Daily(WSPD_Vect_Troll_Calculated,WSPD_Vect_Troll_Measured,"Day","WSPD", "WSPD Calculated", "WSPD Measured", "Troll July 2020")




#WSPD_Troll_Hourly_Calculated_2019, WDir_Vect_Troll_Calculated_2019 = get_calculated_weather_vect_Troll_year(2019)
WSPD_Troll_Hourly_Calculated_2020, WDir_Vect_Troll_Calculated_2020 = get_calculated_weather_vect_Troll_year(2020)
#WSPD_Troll_Hourly_Calculated_2021, WDir_Vect_Troll_Calculated_2021 = get_calculated_weather_vect_Troll_year(2021)

#WSPD_T_conc_Measured_2021, WDIR_T_conc_Measured_2021, TIME_T_conc_Measured_2021 = combine_troll("2021 Troll", 2021)
WSPD_T_conc_Measured_2020, WDIR_T_conc_Measured_2020, TIME_T_conc_Measured_2020 = combine_troll("2020 Troll", 2020)
#WSPD_T_conc_measured_2019, WDIR_T_conc_Measured_2019, TIME_T_conc_Measuredd_2019 = combine_troll("2019 Troll", 2019)

#2019
#WSPD_Troll_Hourly_measured_2019 = Hourly_Concat(WSPD_T_conc_measured_2019)
#WDIR_Troll_Hourly_Measured_2019 = Hourly_Concat(WDIR_T_conc_Measured_2019)
#TIME_Troll_Hourly_measured_2019 = Time_Concat(TIME_T_conc_Measured_2019)

#2020
WSPD_Troll_Hourly_measured_2020 = Hourly_Concat(WSPD_T_conc_Measured_2020)
#WDIR_Troll_Hourly_Measured_2020 = Hourly_Concat(WDIR_T_conc_Measured_2020)
#TIME_Troll_Hourly_measured_2020 = Time_Concat(TIME_T_conc_Measured_2020)

#2021
#WSPD_Troll_Hourly_measured_2021 = Hourly_Concat(WSPD_T_conc_Measured_2021)
#WDIR_Troll_Hourly_Measured_2021 = Hourly_Concat(WDIR_T_conc_Measured_2021)
#TIME_Troll_Hourly_measured_2021 = Time_Concat(TIME_T_conc_Measured_2021)


#plot_Vect_Weekly(WSPD_Troll_Hourly_Calculated_2020,WSPD_Troll_Hourly_measured_2020,  "Week",
#                "WSPD in m/s","Calculated","Measured","Troll 2020 by week")

#plot_histogram(WSPD_Troll_Hourly_measured_2020,WSPD_Troll_Hourly_Calculated_2020, "Distribution of Wind Speeds Measured vs Calculated")

#plot_Vect_Daily(WSPD_Vect_Sleipner_Calculated,WSPD_Vect_Sleipner_Measured,"Hour","Average WSPD", "Average WSPD Calculated", "Average WSPD Measured","Sleipner March 2021")
#plot_Vect_Daily(WSPD_Vect_Troll_Calculated,WSPD_Vect_Troll_Measured,"Hour","Average WSPD", "Average WSPD Calculated", "Average WSPD Measured", "Troll July 2020")

#plot_Vect_hourly(WSPD_Vect_Sleipner_Calculated,WSPD_Vect_Sleipner_Measured,"Hour","WSPD", "WSPD Calculated", " WSPD Measured","Sleipner March 2021")
plot_Vect_hourly(WSPD_Vect_Troll_Calculated,WSPD_Vect_Troll_Measured,"Hour","WSPD", "WSPD Calculated", "WSPD Measured", "Troll July 2020")

#plot_Vect_hourly_single(WSPD_Vect_Troll_Measured,"Hour","WSPD", "WSPD measured", "Troll July 2020")