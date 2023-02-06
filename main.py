#yay 2.0

import math

import pandas as pd
#from pandas import DataFrame
#import random
import numpy as np
import matplotlib.pyplot as plt
import geopy.distance #package to calculate distance between two lat/lon points
#from env.old_code.MetoceanDownloader import MetoceanDownloader
import netCDF4 as nc #read .nc files with weather data
from numpy.ma.core import MaskedConstant
from scipy.interpolate import interp1d
from bisect import bisect_left
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
filename_AIS = "env/input_files/ais_data_v4.csv"
travel_iteration            = 0
#vessel parameters:
vessel_length               = 101.26
vessel_draft                = 10.15
vessel_weight_disp          = 8856
#vessel_velocity             = 1*0.514444444
rho_air                     = 1.025
rho_water                   = 1025

#comparisson_vessel_power_by_speed = vessel_velocity**3*0.8*0.55

#datafiles

NW_file = "Weather_Data/NW_01_07_2020__01_07_2022.nc" #Keys =['northward_wind', 'time', 'lat', 'lon']
EW_file = "Weather_Data/EW_01_07_2020__01_07_2022.nc" #Keys =['eastward_wind', 'time', 'lat', 'lon']
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
#print(dataset_NW.variables.keys())
#print(dataset_EW.variables.keys())
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
#print(f"The size of latitudes is: {eastward_lat.shape}")
#print(f"The size of longditudes  is: {eastward_lon.shape}")
#print(dataset_NW["northward_wind"][0,60,20])

#print(f"hour since 1990{date(1990,1,1)+ timedelta(days=42731850)}")

#write data found to file called filename_func
def write_to_file(data,filename_func):
    pd.DataFrame(data).to_csv(filename_func)
    return 0

#reads columns of dataframe to find distances, latitudes, longditudes of trip
def read_cols(filename_ais_data):
    df= pd.read_csv(filename_ais_data)
    #a = df[['dt',"distance","travel_time"]]
    distances  = df[["distance"]]
    travel_time = df[["travel_time"]]
    latitudes   = df[["lat"]]
    longditudes = df[["lon"]]
    heading     = df[["heading"]]
    dist_vect_func = []
    travel_time_vect_func   = []
    latitudes_vect_func     = []
    longditudes_vect_func   = []
    heading_vect_func       = []
    for i in range (len(distances)):
        dist_vect_func.append(distances.loc[i].iat[0])
        travel_time_vect_func.append(travel_time.loc[i].iat[0])
        latitudes_vect_func.append(round(latitudes.loc[i].iat[0],4))
        longditudes_vect_func.append(round(longditudes.loc[i].iat[0],4))
        heading_vect_func.append(heading.loc[i].iat[0])
    heading_file = "env/input_files/heading.csv"
    write_to_file(heading, heading_file)
    return dist_vect_func,travel_time_vect_func,latitudes_vect_func,longditudes_vect_func,heading_vect_func

#create vector of positions over trip, and removes duplicate positions, saves positions to file
def vector_of_positions(lats,lons):
    position_array_func = []
    for i in range(len(lats)):
        position_array_func.append((round(lats[i],3),round(lons[i],3)))
    #remove duplicate positions
    #print(len(position_array_func))
    j= 0
    while j < len(position_array_func)-1:
        lat1 = position_array_func[j][0]
        lat2 = position_array_func[j+1][0]
        lon1 = position_array_func[j][1]
        lon2 = position_array_func[j+1][1]
        if lat1 == lat2 and lon1 == lon2:
            del position_array_func[j]
        else:
            j += 1
    #print(len(position_array_func))
    position_array_file  = "env/position_array"
    write_to_file(position_array_func,position_array_file)
    return position_array_func

#dist_vect,travel_time,latitudes_vect,longditudes_vect,heading_vect = read_cols(filename_AIS)
#vector_of_positions(latitudes_vect,longditudes_vect)

#data = download_data() #Downloads data locally

#calculates weather at time=tid, latitude and londitude of trip
#def getweather(tid,latitude, longditude):
    #to get the correct index one needs to increase values by 0.125
    #for example, latitude 58.9375 = latitude[0]
#    lat_pos_north = round((latitude-dataset_NW["lat"][0])*8)/8
    lon_pos_north = round((longditude-dataset_NW["lon"][0])*8)/8
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
def getweather(tid,latitude, longditude):
    #to get the correct index one needs to increase values by 0.125
    #for example, latitude 58.9375 = latitude[0]
    lat_pos = int((latitude-eastward_lat[0])*8)              #Access correct position in vector of north wind
    lon_pos = int((longditude-eastward_lon[0])*8)             #Access correct position in vector of east wind


    #tid += 962409600 #measurements go from 01/07/2020 - to 01/07/2022
    if tid < 0 or tid >= len(dataset_NW["northward_wind"][:,lat_pos,lon_pos]):
        print(f"time was out of bounds at {tid}, time must be within the span of one year (less than 1460)")
        return 1
        print("latitude out of bounds, latitude between 58.9375 and 70.1875")
        return 1
    elif longditude < 3.0625 or longditude > 20.9375: #været må være hentet på posisjonen longditude
        print("longditude out of bounds, longditude between 3.0625 and 20.9375")
        return 1
    elif lat_pos >= len(dataset_NW["northward_wind"][tid,:,lon_pos]):
        print(f"lat posistion is out of bound at {lat_pos} degrees")
        return 1
    elif lon_pos >= len(dataset_NW["northward_wind"][tid, lat_pos, :]):
        print(f"lon posistion is out of bound at {lon_pos} degrees")
        return 1
    WSN = dataset_NW["northward_wind"][tid,lat_pos,lon_pos]
    WSE = dataset_EW["eastward_wind"][tid,lat_pos,lon_pos]

    return WSN,WSE

#a,b = getweather(0,58.938,4)

#getweather(datetime)
#%date%hour%minute%second%

def datetime_seconds(dtime):
    dt = dtime
    epoch = date(2020,7,1)
    delta = (dt-epoch)
    return delta.total_seconds()
    


#calculates direction between two points, (20.323,23.243),(34.235, 43.345)
def calc_bearing(pointA, pointB):
    deg2rad = math.pi / 180
    latA = pointA[0] * deg2rad
    latB = pointB[0] * deg2rad
    lonA = pointA[1] * deg2rad
    lonB = pointB[1] * deg2rad

    delta_ratio = math.log(math.tan(latB/ 2 + math.pi / 4) / math.tan(latA/ 2 + math.pi / 4))
    delta_lon = abs(lonA - lonB)

    delta_lon %= math.pi
    bearing = math.atan2(delta_lon, delta_ratio)/deg2rad
    return bearing
def read_position_vect_from_file(filename_func):
    readdata = pd.read_csv(filename_func, usecols=["0", "1"])
    readdata = readdata[readdata != 0]

    position_array =  []
    #print(f"im reading {filename_func} now")
    for i in range(len(readdata)):
        position_array.append((round(readdata.loc[i].iat[0], 3),round(readdata.loc[i].iat[1], 3)))
    return position_array
#print(calc_bearing(pos1, pos2))

#translates radians to degrees for python
def r2d(rad):
    degree = rad*180/np.pi
    return degree
#translates degrees to rads for pythong
def d2r(degree):
    rad = np.pi*degree/180
    return rad
#Finds windspeeds at position in position array and time in tid in weater data
def Find_WSV_over_trip_at_time_TID(position_array_func, tid):
    Wind_North_vector_func = [] #gives wind_speed_north as position (position_Array [tid])
    Wind_East_vector_func  = []
    Wind_tot_vector_func = []
    for j in range(len(position_array_func)):
        lat_func = position_array_func[j][0]
        lon_func = position_array_func[j][1]
        WSN_func, WSE_func= getweather(tid, lat_func,lon_func) #gives wind at time tid for all lat_vectors and lon vectors
        Wind_North_vector_func.append(WSN_func)
        Wind_East_vector_func.append(WSE_func)
        Wind_tot_vector_func.append(np.sqrt(WSN_func**2 + WSE_func**2)) #returns total windspeed at position index j
    return Wind_North_vector_func,Wind_East_vector_func, Wind_tot_vector_func

#function that gives true windspeed in m/s for wind north,east
def AWS(TWS_func, sailing_speed_func,sailing_direction_func):
    AWS_func = TWS_func - sailing_speed_func * np.sin(np.pi / 180 * sailing_direction_func)
    return AWS_func
#Function that return true wind direction in degrees from wind speed north/east (WSN/WSE)
def AWD(WSNA, WSEA,sailing_dir_func):
    AWD_func = round((180/np.pi)*(math.atan2(WSNA,WSEA))+sailing_dir_func,2)
    if AWD_func >= 360:
        AWD_func -= 360
    return AWD_func
#function that returns true wind direction and speed using earlier functions
def Apparent_wind_and_direction_arrays_over_trip(position_array_func, tid):#, vessel_velocity):
    Wind_North_vector_func, Wind_East_vector_func, Wind_tot_vector_func = Find_WSV_over_trip_at_time_TID(position_array_func, tid)          # windspeeds at (lat,lon) for t € 1050930 through 1051824
    AWD_Trip  = AWD(Wind_North_vector_func, Wind_East_vector_func)  #true winddirection at (lat,lon) for t € 1050930 through 1051824
    AWS_Trip  = AWS(Wind_tot_vector_func, AWD_Trip)#, vessel_velocity)          #true windspeed at (lat,lon) for t € 1050930 through 1051824
    return AWS_Trip,AWD_Trip
#find force at each position of route
def Force_at_position(AWS, AWD):
    lift = 0.5 * rho_air * A * AWS ** 2 * Cl                                 #lift force from traut
    drag = 0.5 * rho_air * A * AWS ** 2 * Cd                                 #drag force from traut
    #P_input_fletner = 0.5 * rho_air * A * TWS ** 3 * Cm * alpha              #Input to flettner is energy to spin rotors
    if 0 <= AWD <= 90 or 270 <= AWD <= 360:                                  #drag is set to negative if the wind is coming ahead, and positive if not
        drag *= -1
    Force_4_flettners = (lift + drag)*4/1000 #kN                             #KiloNewton
    #The force generated by the flettner is the total output force minus the force needed to spin the rotor.

    perp_force      = Force_4_flettners*abs(np.cos(AWD))                     #perpendicular force flettner in KN
    if type(perp_force) != MaskedConstant:
        perp_force = round(min(perp_force, 1400),2)  # Force is maximum 1400 kN, so as not to break flettners (cite: Norsepower technical)
    forward_force   = Force_4_flettners*abs(np.sin(AWD))            #directional force flettner in KN

    if type(forward_force) != MaskedConstant:
        forward_force   = round(min(forward_force,1400),2)              #Force is maximum 1400 kN, so as not to break flettners

    return forward_force, perp_force
#returns avg forces over each point of a trip
def Force_over_trip(position_vector, time):
    Forward_force_array   = np.zeros(len(position_vector))
    Perp_force_array      = np.zeros(len(position_vector))
    Speed_sailed_array    = np.zeros(len(position_vector))
    TWS, AWD_Array        = Apparent_wind_and_direction_arrays_over_trip(position_vector, time)
    for Pos_num in range(len( position_vector)): #Finds total force over each trip
        if Pos_num >= 1:
            Forward_force_position, Perp_force_position = Force_at_position(TWS, AWD_Array, Pos_num)
            if Forward_force_position == Forward_force_array[Pos_num-1]:
                Speed_sailed_array[Pos_num] = Speed_sailed_array[Pos_num-1]  # KNOTS
                Forward_force_array[Pos_num] = Forward_force_array[Pos_num-1]
                Perp_force_array[Pos_num] = Perp_force_array[Pos_num-1]
            else:
                speed_sailed = Speed_sailed(Perp_force_position,Forward_force_position) #KNOTS
                distance = speed_sailed*time
                #print(distance)
                Speed_sailed_array[Pos_num]  = speed_sailed #KNOTS
                Forward_force_array[Pos_num] = Forward_force_position
                Perp_force_array[Pos_num]    = Perp_force_position
        else:
            Forward_force_position, Perp_force_position = Force_at_position(TWS, AWD_Array, Pos_num)
            speed_sailed = Speed_sailed(Perp_force_position, Forward_force_position)  # KNOTS
            Speed_sailed_array[Pos_num] = speed_sailed  # KNOTS
            Forward_force_array[Pos_num] = Forward_force_position
            Perp_force_array[Pos_num] = Perp_force_position
    Forward_force_total = np.sum(Forward_force_array)
    Perp_force_total    = np.sum(Perp_force_array)
    Speed_sailed_total  = np.sum(Speed_sailed_array) #KNOTS
    Avg_perp_force_trip = Perp_force_total/len(Perp_force_array)
    Avg_forw_force_trip = Forward_force_total/len(Forward_force_array)
    Avg_WS              = sum(TWS)/len(TWS)
    Avg_Speed_sailed    = Speed_sailed_total/len(Speed_sailed_array) #KNOTS
    return Avg_forw_force_trip, Avg_perp_force_trip, Avg_WS, Avg_Speed_sailed#KNOTS , Speed_sailed_array
#function to calculate forward force, average force, perpendicular force over time array, and add to files
def Force_over_year(six_hour_periods):#1460 for a year, 7307 for whole period
    forw_force_array_time           = np.zeros(len(six_hour_periods))
    perp_force_array_time           = np.zeros(len(six_hour_periods))
    avg_AWS_array_over_time         = np.zeros(len(six_hour_periods))
    speed_sailed_array_over_time    = np.zeros(len(six_hour_periods))
    for time in range(0,len(six_hour_periods)):#skal være 7307 Finds force at each position for each day of year
        avg_perp_force_trip, avg_forward_force_trip , avg_WS, avg_speed_sailed_over_trip = Force_over_trip(position_array, time)
        forw_force_array_time[time]             = avg_forward_force_trip
        perp_force_array_time[time]             = avg_perp_force_trip
        avg_AWS_array_over_time[time]           = avg_WS
        speed_sailed_array_over_time[time]      = avg_speed_sailed_over_trip
        if time%10 == 0:
            #print(f"speed at this iteration is {vessel_velocity*5.44} knots")
            print(f"avg. force at day {time/4} is equal to {forw_force_array_time[time]}")
            print(f"perp. force at day {time/4} is equal to {perp_force_array_time[time]}")
            print(f"average speed sailed at day {time/4} is equal to {speed_sailed_array_over_time[time]} knots")
            print(f"Current time= {datetime.now().time()}")
    # writing to pandas, uncomment to rewrite file
    filename_avg_forward_force  = f"env/output_files5/avg_forward_force.txt"
    filename_perp_force         = f"env/output_files5/avg_perp_force.txt"
    filename_speed_sailed       = f"env/output_files5/speed_sailed_over_time.txt"
    write_to_file(forw_force_array_time, filename_avg_forward_force)
    write_to_file(perp_force_array_time, filename_perp_force)     #winddir 31.12
    write_to_file(speed_sailed_array_over_time,filename_speed_sailed)
    return forw_force_array_time, perp_force_array_time, speed_sailed_array_over_time
#plot power from Force_over_trip with respect to time
def plot_power(title, y_axis, x_label, y_label):
    x_axis = []
    for i in range(1,10):
        x_axis.append(730*i)
    plt.plot(y_axis)
    plt.xticks(x_axis)
    plt.title(title)
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.show()
    return 0
#plot average power over the 4 year period
def plot_avg_power(tid,power):
    plt.plot(tid,power)
    plt.title(" avg Poweroutput from flettnerrotor")
    plt.xlabel("triptime of year")
    plt.ylabel("kN produced by flettner on average")
    plt.show()
#plot weekly and daily power over 4 years
def plot_weekly_and_daily_avg_power(power):
    #create weekly average power vector
    monthly_power_average = []
    weekly_power_average = []
    daily_power_average = []
    for i in range(0,len(power)-4,4):
        daily_power_sum = 0
        for j in range(4):
            daily_power_sum += power[i+j]
        daily_power_average.append(daily_power_sum/4)

    for i in range(0,len(daily_power_average)-7,7):
        weekpowersum = 0
        for j in range(7):
            weekpowersum += daily_power_average[i+j]
        weekly_power_average.append(weekpowersum/7)

    for i in range(0,len(weekly_power_average)-4,4):
        monthly_power_sum = 0
        for j in range(4):
            monthly_power_sum += weekly_power_average[i+j]
        monthly_power_average.append(monthly_power_sum/4)
    plt.plot(daily_power_average)
    plt.title(" avg daily poweroutput from flettnerrotor")
    plt.xlabel("day of dataset")
    plt.ylabel("kW produced by flettner on average that day")
    plt.show()
    plt.plot(weekly_power_average)
    plt.title(" avg Poweroutput from flettnerrotor per week")
    plt.xlabel("week of dataset")
    plt.ylabel("kW produced by flettner on average that week")
    plt.show()
    plt.plot(monthly_power_average)
    plt.title(" avg Poweroutput from flettnerrotor per month")
    plt.xlabel("month of dataset")
    plt.ylabel("kW produced by flettner on average that month")
    plt.show()
#clears entries with nav status 2 or 5 for AIS data
def find_vals():
    df = pd.read_csv("env/input_files/ais_data_v3.csv")
    del df["Unnamed: 0"]
    df = df[df.nav_status !=2] #removes every entry with nav_Status = 2
    df = df[df.nav_status !=5] #removes every entry with nav_status = 5
    df.to_csv("ais_data_v4.csv")
    return 0
#function to read data from file (save time after running program through)
def read_array_from_file(filename_func):
    readdata = pd.read_csv(filename_func)
    readdata = readdata[readdata != 0]
    data_array = []
    print(f"im reading {filename_func} now")
    for i in range(len(readdata)):
        data_array.append(round(readdata.loc[i].iat[1],3))
    #print(data_array)
    return data_array
#dist_vect,travel_time_vect,latitudes_vect,longditudes_vect, heading_vect = read_cols(filename_AIS)
#position_array = vector_of_positions(latitudes_vect,longditudes_vect)[:1460]
#print("position array = ", position_array)
#print("heading vector = ", heading_vect)x|



#units time since 1.1.1900 00:00:00. Start date: 11.03.2013 to 26.06.2014 12:00:00 in 6 hour intervalls
#lowerbound = 999354.0, Upper bound = 1003242, Length: 3888
#time_Vector = ds2["time"][:]
#for i in range(len(time_Vector)):
#    time_Vector[i] -= 981774

#lat_Vector = ds2["lat"][:]
#lon_Vector = ds2["lon"][:]

#returns propulsion force at every position along vector of positions from the true wind speed and direction
# here one should input rawdata from flettner efficiency, for now using approximated values from paper
# using source DOI: 10.1016/j.apenergy.2013.07.026, "Propulsive power contribution of a kite and a Flettner rotor on selected shipping routes"
rotor_amount = 4
rho_air = 1.025 #density of air
h   = 35    #height flettner
d   = 5     #diameter of flettner
A   = h*d  #cross sectional area of flettner
Cl  = 12.5
Cd  = 0.2
Cm  = 0.2
alpha = 3.5

#finds resistance by speed using admirality formula
def resistance_by_speed(sailing_speed):
    admirality_coeff = 385.46
    power = (vessel_weight_disp**(2/3)*sailing_speed**3)/admirality_coeff
    force = round(power/(sailing_speed/5.44444444),2)
    #print(f"vessel power required at {sailing_speed} knots equals:{power} Watt")
    #print(f"resistance force required at {sailing_speed} knots equals:{force} kN")
    return force
#plots resistance of speed found by admirality formula
def plot_resistance():
    speeds = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16]
    resistance = [0,6.05,24.18,54.41,96.74,151.15,217.65,299.25,386.98,489.72,604.6,731.56,870.62,1021.77,1185.01,1360.34,1547.77]
    plt.plot(speeds,resistance)
    plt.title("Resistance in kN")
    plt.xlabel("Vessel speed")
    plt.ylabel("kN")
    plt.show()
    return 0
#function that finds driftangle beta that is large enough to counteract drifting force
def solve_beta_by_perp_force(perp_force_func, vessel_velocity):
    Fy = perp_force_func*1000 #Newton
    RHS = Fy/(0.5*vessel_length*vessel_draft*vessel_velocity**2)
    # Solving second degree function to find beta
    a = 0.461
    b = 1.091
    c = -RHS
    Beta_1 = (-b + np.sqrt(b ** 2 - 4 * a * c)) / (2 * a)
    Beta_2 = (-b - np.sqrt(b ** 2 - 4 * a * c)) / (2 * a)
    Beta_func = round(max(Beta_1, Beta_2),3)  # the angle will always be the positive solution
    return Beta_func
#function that finds xfold resistance by beta (Drift angle)
def xfold(beta):
    if beta < 8:
        modeltest = [0,1, 1.2, 1.25, 1.45, 1.9]
        modelX = [0, 2, 4, 6, 8, 10]

        skogman = [0, 1, 1.18, 1.3, 1.59]
        skogmanX = [0, 2, 4, 6, 8]

        #wagner_mariner = [0, 1, 1.18, 1.3, 1.7]
        #wagnerX = [0, 2, 4, 6, 8]

        #wagner_arbitrary = [0, 1, 1.21, 1.599, 2.1]
        #wagnerarbX = [0, 2, 4, 6, 8]

        xfoldmodel = interp1d(modelX, modeltest)
        xfoldskogman = interp1d(skogmanX, skogman)
        xfold_final = (xfoldmodel(beta) + xfoldskogman(beta)) / 2
    else:
        #if driftangle exceeds 8 degrees, then resistance is doubled
        xfold_final = 2
    return xfold_final
# finds speed sailed given perp force and forward force IN KNOTS
def Speed_sailed(perp_force, forward_force):
    # Empty vectors to store values
    sailing_resistance_vector = []
    total_resistance_vector = []

    # Set vessel speed [knots] in intervall from 0,20 with stepsize 0.1
    vessel_velocity = np.linspace(0.1, 20, 200)
    # ratio hydrodynamic resistance to total res is approximately 0.85
    ratio_hydrodyn_to_tot_res = 0.85

    # loop through all vessel velocitites [knots] and find correpsonding drifting and sailing resistance
    # add these together: total resistance
    for velocity in vessel_velocity:
        sailing_resistance = resistance_by_speed(velocity) * ratio_hydrodyn_to_tot_res
        sailing_resistance_vector.append(sailing_resistance)
        drift_angle = solve_beta_by_perp_force(perp_force, velocity)
        resistance_multiplier = xfold(drift_angle)
        if resistance_multiplier < 1:
            resistance_multiplier = 1
        total_resistance = resistance_multiplier * sailing_resistance
        total_resistance_vector.append(total_resistance)

        #stop iteration when total_resistance > forward force
        if total_resistance > forward_force:
            speed_achieved = velocity
            #speed_achieved = take_closest(total_resistance_vector,forward_force) #IN KNOTS
            return speed_achieved
    speed_achieved = take_closest(total_resistance_vector,forward_force) #IN KNOTS
    return speed_achieved #IN KNOTS
def Speed_sailed_point(perp_force, forward_force, sailing_speed):

    # ratio hydrodynamic resistance to total res is approximately 0.85
    ratio_hydrodyn_to_tot_res = 0.85

    # Empty vectors to store values
    sailing_resistance_vector = []
    total_resistance_vector = []

    # Set vessel speed [knots] in intervall from 0,20 with stepsize 0.1
    vessel_velocity = np.linspace(0.1, 20, 200)
    #
    for velocity in vessel_velocity:
        sailing_resistance = resistance_by_speed(velocity) * ratio_hydrodyn_to_tot_res
        sailing_resistance_vector.append(sailing_resistance)
        drift_angle = solve_beta_by_perp_force(perp_force, velocity)
        resistance_multiplier = xfold(drift_angle)
        if resistance_multiplier < 1:
            resistance_multiplier = 1
        total_resistance = resistance_multiplier * sailing_resistance
        total_resistance_vector.append(total_resistance)

        #stop iteration when total_resistance > forward force
        if total_resistance > forward_force:
            speed_achieved = velocity
            #speed_achieved = take_closest(total_resistance_vector,forward_force) #IN KNOTS
            return speed_achieved
    speed_achieved = take_closest(total_resistance_vector,forward_force) #IN KNOTS


    return speed_achieved #IN KNOTS

#function interpolates over y values in list and finds closest value in table, returns x value (basically inverse func)
def take_closest(myList, myNumber):

    pos = bisect_left(myList, myNumber)
    if pos == 0:
        return myList[0]
    if pos == len(myList):
        return myList[-1]
    before = myList[pos - 1]
    after = myList[pos]
    if after - myNumber < myNumber - before:
        return myList.index(after)/10
    else:
        return myList.index(before)/10

def start_from_files():
    avg_speed_sailed        = "env/output_files5/speed_sailed_over_time.txt"
    avg_force_data          = "env/output_files5/avg_forward_force.txt"
    avg_perp_force_data     = "env/output_files5/avg_perp_force.txt"

    avg_speed_sailed_array  = read_array_from_file(avg_speed_sailed)
    avg_force_array         = read_array_from_file(avg_force_data)
    avg_perp_force_array    = read_array_from_file(avg_perp_force_data)

    print(f"speed sailed array          {avg_speed_sailed_array}")
    print(f"propulsive force array      {avg_force_array}")
    print(f"perpendicular force array   {avg_perp_force_array}")

    print(f"Average speed sailed over 5 years:                  {np.average(avg_speed_sailed_array)} knots")
    print(f"Average force produced over 5 years:                {np.average(avg_force_array)} kN")
    print(f"Averge Perpendicular force experienced over 5 years {np.average(avg_perp_force_array)} kN")

    print(f"length sailing speed array      {len(avg_speed_sailed_array)}")
    print(f"length force array              {len(avg_force_array)}")
    print(f"length perpendicular force array{len(avg_perp_force_array)}")


    plot_power("Average speed sailed over route", avg_speed_sailed_array,"Month of period", "knots")
    plot_power("Average propulsive force over route",avg_force_array, "Month of period","kN produced by flettner")
    plot_power("Average perpendicular force over route", avg_perp_force_array,"Month of period", "kN produced by flettner")

    return 0
time_Vector = [0,1]
def startfromscratch():
    avg_forward_force_over_time, avg_perp_force_over_time, speed_sailed_over_time = Force_over_year(time_Vector)
    return avg_forward_force_over_time, avg_perp_force_over_time

def plot_power_2(title, y_axis, x_label, y_label):
    speed = [1,3,6,9,12]
    plt.plot(speed,y_axis)
    plt.title(title)
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.show()
    return 0

def plot_percent(title, y_axis, x_label, y_label):
    speed = [3,6,9,12]
    plt.plot(speed,y_axis)
    plt.title(title)
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.show()
    return 0

def iterate_drift_angle(vessel_velocity):
    iterations = 80
    Fy_vector = []
    for beta in range(0,iterations):
        beta = beta/10
        Cl_beta = 0.461*beta*+1.091*beta*abs(beta)
        Fy      = Cl_beta*0.5*vessel_length*vessel_draft*vessel_velocity**2
        Fy_vector.append((beta,Fy))
    for j in range(0,iterations,4):
        print(f"A drift angle of {j/10} degrees provides "
              f"{round(Fy_vector[j][1]/1000,3)} kN of force righting force" )
    filename_avg_righting_force = "/Weather/env/old_code/output_files5/righting_force.txt"
    write_to_file(Fy_vector,filename_avg_righting_force)
    return Fy_vector

#################################################

#Run to start from files one speed
#start_from_files()

#file =  "env/old_code/output_files5/speed_sailed_over_time.txt
#file2 = "env/old_code/output_files5/avg_forward_force.txt"
#file3 = "env/old_code/output_files5/avg_perp_force.txt"
#average_Speed = read_array_from_file(file)
#average_force = read_array_from_file(file2)
#average_perp  = read_array_from_file(file3)

#Run to start from files with different speeds()
#start_from_files_dif_speed()

#Run to start program over
#startfromscratch()

#Run to iterate over different speeds
#run_dif_speed()
#################################################

#Calculate vesel speeds achieved:
#calculate_speeds_achieved()

#plot_resistance()




#Windspeed_North, Windspeed_East = getweather(0,60,0)



#Function that calculates time spent sailing from port a to port b
def main(route, time):
    tot_sailing_dist        = 0
    bad_weather_positions   = []
    poor_sailing_time       = 0
    poor_sailing_distance   = 0
    sailing_speed           = 4 #initializing sailing speed of 10 knots (will change after one iteration)
    for i in range(len(route)-1): #create itteration through route
        position_first      = route[i]
        position_next       = route[i+1]
        sailing_distance    = round(geopy.distance.geodesic(position_first,position_next).nautical,3)    #In Nautical Miles
        sailing_direction   = calc_bearing(position_first,position_next)                                 #In Degrees (North is 0)
        WSE_func,WSN_func   = getweather(time,position_first[0],position_first[1])                       #Gives Wind speed East and North
        TWS                 = np.sqrt(WSE_func**2 + WSN_func**2)                                         #Finds True Windspeed (pythagoras)
        AWA                 = np.arctan2(TWS,sailing_speed)
        forward_force_func,perpendicular_force_func = Force_at_position(TWS,AWA)                         #Forward and Perpendicular force from Flettners
        if type(forward_force_func) != MaskedConstant or type(perpendicular_force_func) != MaskedConstant:
            sailing_speed    = Speed_sailed_point(perpendicular_force_func,forward_force_func, initial_speed)    #Sailing Speed obtained in KNOTS
            sailing_time     = sailing_distance/sailing_speed                                                    #time used to sail trip added
            time             += sailing_time
        tot_sailing_dist     += sailing_distance
        #check for extreme time usage
        if sailing_speed < 1:
            #bad_weather_positions.append(f"position is: {(position_first,position_next)},sailing speed is {sailing_speed},"
            #                             f" sailing time is {sailing_time}, and sailing distance is {sailing_distance}. \n")
            poor_sailing_time += sailing_time
            poor_sailing_distance += sailing_distance

        #print(*bad_weather_positions, sep = "\n")
        #print(f"total time sailed at less than 1 knot is {poor_sailing_time}\n"
               #f"this time is used to sail {poor_sailing_distance} nautical miles")

        #print(f"total sailing time on this route is {time_of_trip} and distance sailed is [{tot_sailing_dist}")
    return time,tot_sailing_dist, poor_sailing_time, poor_sailing_distance, sailing_speed


filename = "env/position_array"  # file containing position array of route
position_array = read_position_vect_from_file(filename) # reads said file

#Wind_North_vector_func,Wind_East_vector_func, Wind_tot_vector_func = Find_WSV_over_trip_at_time_TID(position_array,0)


#time_of_trip,total_sailing_distance = main(position_array)
#print(f"Time of trip is {time_of_trip},\n"
#      f" giving a total trip distance is {total_sailing_distance}\n"
#      f"with an average sailing speed of {total_sailing_distance/time_of_trip} knots")


def read_route(csv):
    df = pd.read_csv(csv)
    latitudes   = df["latitude"]
    longditudes = df["longditude"]
    positions = []
    for i in range(len(latitudes)):
        positions.append((latitudes[i],longditudes[i]))
    return positions



#main = time_of_trip,tot_sailing_dist, poor_sailing_time, poor_sailing_distance

#change to simulation(route) so that you dont need to hardcode
def simulation(csv):
    route_travel            = read_route(csv)
    time_of_simulation      = 2*365*24*60*60 #two years in seconds
    time_of_trip            = np.zeros(time_of_simulation)
    tot_sailing_dist        = np.zeros(time_of_simulation)
    poor_sailing_time       = np.zeros(time_of_simulation)
    poor_sailing_distance   = np.zeros(time_of_simulation)
    sailing_speed_vect      = np.zeros(time_of_simulation)
    for time in range(0,time_of_simulation,3600): #calculating every hour

        time_of_trip_1,tot_sailing_dist_1, poor_sailing_time_1, poor_sailing_distance_1, sailing_speed = main(route_travel,time)
        time_of_trip[time]             = time_of_trip_1
        tot_sailing_dist[time]         = tot_sailing_dist_1
        poor_sailing_time[time]        = poor_sailing_time_1
        poor_sailing_distance[time]    = poor_sailing_distance_1
        sailing_speed_vect[time]       = sailing_speed
        if time%10 == 0 and time > 0:
            poor_sailing_speed = sum(poor_sailing_distance) / sum(poor_sailing_time)
            print(f"speed sailing {csv}, distance of {tot_sailing_dist[time]} is {np.average(sailing_speed_vect[0:time])} knots")
            print(f"total time sailed at less than 1 knot is {poor_sailing_time[time]}\n"
                  f"this time is used to sail {poor_sailing_distance[time]} nautical miles\n"
                  f"at an average speed of {poor_sailing_speed} knots")
            print(datetime.now())
    poor_sailing_speed = sum(poor_sailing_distance)/sum(poor_sailing_time)
    print(f"speed sailing {csv} is {np.average(sailing_speed_vect)}")
    print(f"total time sailed at less than 1 knot is {sum(poor_sailing_time)}\n"
          f"this time is used to sail {sum(poor_sailing_distance)} nautical miles\n"
          f"at an average speed of {poor_sailing_speed} knots")
    return time_of_trip,tot_sailing_dist,sailing_speed_vect

start_position      = (0,60)    #Cooridnates of port a
end_position        = (10,56)   #Coordinates of port b
initial_speed       = 4
Trond_aalesund      = "Rute_Trondheim_Aalesund.csv"
Aalesund_Floro      = "Rute_Aalesund_floro.csv"
Floro_Bergen        = "Rute_Floro_Bergen.csv"
Bergen_Stavanger    = "Rute_Bergen_Stavanger.csv"


Trip_time_vector, Tot_sailing_distance_vector, sailing_speed_vector = simulation(Trond_aalesund)
#simulation(Aalesund_Floro)
#simulation(Floro_Bergen)
#simulation(Bergen_Stavanger)


print(f"Trip time for each repetition {Trip_time_vector}")
print(f"Total Sailing dist for each repetition {Tot_sailing_distance_vector[:10]}, should be equal")
print(f"Sailing speed for each repetition {sailing_speed_vector} in knots")


print(r2d(np.arctan2(0.75,0.50)),r2d(np.arctan2(0.75,-0.50)),r2d(np.arctan2(-0.75,-0.50)),r2d(np.arctan2(-0.75,0.50)))
print(r2d(np.arctan2(10,4)))

def prøve_å_forstå(Vessel_heading,Vessel_Speed,Wind_heading,Wind_speed):
    route_travel_ex     = read_route(Trond_aalesund)
    Vessel_speed_x      = np.cos(Vessel_heading)*Vessel_Speed
    Vessel_speed_y      = np.sin(Vessel_heading)*Vessel_Speed
    Wind_speed_x        = np.cos(Wind_heading)*Wind_speed
    Wind_speed_y        = np.sin(Wind_heading)*Wind_speed
    x                   = Wind_speed_x + Vessel_speed_x
    y                   = Wind_speed_y + Vessel_speed_y
    AWA                 = r2d(np.arctan2(y,x))
    print(f"Apparent wind angle is {AWA} degrees")
    return 0
#prøve_å_forstå(-45,4,60,10)















print("Finished <3<3")
