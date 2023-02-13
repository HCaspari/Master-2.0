
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
import searoute as sr
import chardet
import folium as folium
import webbrowser
import gmplot as gmp


from file_handling import write_to_file
from Resistance_functions import resistance_by_speed,solve_beta_by_perp_force, xfold
from plot_functions import plot_power, plot_power_2, plot_percent, plot_avg_power, plot_weekly_and_daily_avg_power, plot_resistance
from Weather_Handling import getweather, r2d, d2r, True_wind_direction, True_wind_speed, Apparent_Wind_Speed, Apparent_wind_angle, alpha

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
rotor_amount = 4
h   = 35    #height flettner
d   = 5     #diameter of flettner
A   = h*d  #cross sectional area of flettner
Cl  = 12.5
Cd  = 0.2
Cm  = 0.2
alpha_const = 3.5

#comparisson_vessel_power_by_speed = vessel_velocity**3*0.8*0.55

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

#print(calc_bearing(pos1, pos2))






#find force at each position of route

#clears entries with nav status 2 or 5 for AIS data


#returns propulsion force at every position along vector of positions from the true wind speed and direction
# here one should input rawdata from flettner efficiency, for now using approximated values from paper
# using source DOI: 10.1016/j.apenergy.2013.07.026, "Propulsive power contribution of a kite and a Flettner rotor on selected shipping routes"



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


def Force_at_position(AWS, AWD):
    lift = 0.5 * rho_air * A * AWS ** 2 * Cl                                 #lift force from traut
    drag = 0.5 * rho_air * A * AWS ** 2 * Cd                                 #drag force from traut
    #P_input_fletner = 0.5 * rho_air * A * TWS ** 3 * Cm * alpha_const              #Input to flettner is energy to spin rotors
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


#Function that calculates time spent sailing from port a to port b, through predetermined route
def main(route, time):
    tot_sailing_dist        = 0
    poor_sailing_time       = 0
    poor_sailing_distance   = 0
    sailing_speed_vector    = []
    route_sailing_time      = []
    vessel_speed           = 4 #initializing sailing speed of 4 knots (will change after one iteration)

    for i in range(len(route)-1): #create itteration through route
        position_first      = route[i]
        position_next       = route[i+1]
        sailing_distance    = round(geopy.distance.geodesic(position_first,position_next).nautical,3)    #In Nautical Miles
        vessel_heading   = calc_bearing(position_first,position_next)                                 #In Degrees (North is 0)
        WSE_func,WSN_func   = getweather(time,position_first[0],position_first[1])                       #Gives Wind speed East and North
        TWS                 = np.sqrt(WSE_func**2 + WSN_func**2)                                         #Finds True Windspeed (pythagoras)
        AWA                 = alpha(vessel_speed,vessel_heading,WSN_func,WSE_func )                                              #Finds Apparent wind angle
        forward_force_func,perpendicular_force_func = Force_at_position(TWS,AWA)                         #Forward and Perpendicular force from Flettners
        if type(forward_force_func) != MaskedConstant or type(perpendicular_force_func) != MaskedConstant:
            vessel_speed    = Speed_sailed_point(perpendicular_force_func,forward_force_func, initial_speed)    #Sailing Speed obtained in KNOTS
            sailing_speed_vector.append(vessel_speed)
            sailing_time     = sailing_distance/vessel_speed                                                    #time used to sail trip added
            route_sailing_time.append(sailing_time)
        tot_sailing_dist     += sailing_distance
        #check for extreme time usage
        if vessel_speed < 1:
            poor_sailing_time += sailing_time
            poor_sailing_distance += sailing_distance
    total_time_sailed_route = sum(route_sailing_time)
    return total_time_sailed_route,tot_sailing_dist, poor_sailing_time, poor_sailing_distance, sailing_speed_vector



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
    positions = np.asarray(positions)
    return positions



#main = time_of_trip,tot_sailing_dist, poor_sailing_time, poor_sailing_distance
def simulation(csv):
    """

    :param csv: A csv
    :return:    time_of_trip
                tot_sailing_dist
                sailing_speed_simulation_vector

    """
    hour_intervall                  = 1                                                #at what hourly interval should we simulate?
    route_travel                    = read_route(csv)
    time_of_simulation              = 17520                                             #two years in hours
    time_of_trip                    = np.zeros(int(time_of_simulation/hour_intervall))
    tot_sailing_dist                = np.zeros(int(time_of_simulation/hour_intervall))
    poor_sailing_time               = np.zeros(int(time_of_simulation/hour_intervall))
    poor_sailing_distance           = np.zeros(int(time_of_simulation/hour_intervall))
    sailing_speed_simulation_vector = np.zeros(int(time_of_simulation/hour_intervall))
    for time in range(0,int(time_of_simulation/hour_intervall)) : #calculating every 12 hours

        time_of_trip_1,tot_sailing_dist_1, poor_sailing_time_1, poor_sailing_distance_1, sailing_speed_vector = main(route_travel,time)
        time_of_trip[int(time)]              = time_of_trip_1
        tot_sailing_dist[int(time)]          = tot_sailing_dist_1
        poor_sailing_time[int(time)]         = poor_sailing_time_1
        poor_sailing_distance[int(time)]     = poor_sailing_distance_1
        sailing_speed_simulation_vector[time] = np.average(sailing_speed_vector)
        if time%1000 == 0 and time > 0:
            poor_sailing_speed = sum(poor_sailing_distance) / (sum(poor_sailing_time))
            print(f"speed sailing {csv}, distance of {tot_sailing_dist[time]}\n"
                  f" is {np.average(sailing_speed_vector[0:time])} knots")
            print(f"total time sailed at less than 1 knot is {poor_sailing_time[time]}\n"
                  f"this time is used to sail {poor_sailing_distance[time]} nautical miles\n"
                  f"at an average speed of {poor_sailing_speed} knots")
            print(f"{datetime.now()},time is {time}")
    poor_sailing_speed = sum(poor_sailing_distance)/(sum(poor_sailing_time))
    print(f"speed sailing {csv} is {np.average(sailing_speed_simulation_vector)}")
    print(f"total time sailed at less than 1 knot is {sum(poor_sailing_time)}\n"
          f"this time is used to sail {sum(poor_sailing_distance)} nautical miles\n"
          f"at an average speed of {poor_sailing_speed} knots")
    return time_of_trip,tot_sailing_dist,sailing_speed_simulation_vector



#change to simulation(route) so that you dont need to hardcode
end_position        = (10,56)
start_position      = (0,60)    #Cooridnates of port a
#Coordinates of port b
initial_speed       = 4
Trond_aalesund      = "Route_Trondheim_Aalesund.csv"
Aalesund_Floro      = "Route_Aalesund_floro.csv"
Floro_Bergen        = "Route_Floro_Bergen.csv"
Bergen_Stavanger    = "Route_Bergen_Stavanger.csv"

route_Trond_Aal         = read_route(Trond_aalesund)
route_AalFloro          = read_route(Aalesund_Floro)
route_floro_bergen      = read_route(Floro_Bergen)
route_bergen_stavanger  = read_route(Bergen_Stavanger)
#print(route_bergen_stavanger)

def runsimulation():
    Trip_time_vector_TA, Tot_sailing_distance_vector_TA, sailing_speed_simulation_vector_TA = simulation(Trond_aalesund)
    Trip_time_vector_AF, Tot_sailing_distance_vector_AF, sailing_speed_simulation_vector_AF = simulation(Aalesund_Floro)
    Trip_time_vector_FB, Tot_sailing_distance_vector_FB, sailing_speed_simulation_vector_FB = simulation(Floro_Bergen)
    Trip_time_vector_BS, Tot_sailing_distance_vector_BS, sailing_speed_simulation_vector_BS = simulation(Bergen_Stavanger)

    #Trip_time_vector_TA2, Tot_sailing_distance_vector_TA2, sailing_speed_simulation_vector_TA2 = simulation(Trond_aalesund)

    #print(f"Trip time for each repetition {Trip_time_vector_TA2}")
    #print(f"Total Sailing dist for each repetition {Tot_sailing_distance_vector_TA2[:10]}, should be equal")
    #print(f"Sailing speed for each repetition {sailing_speed_simulation_vector_TA2} in knots")


    #sailing_speed_Trondheim_Aalesund_fil_2 = "Output_files/Trondheim_Aalesund_reise_2"


    sailing_speed_trond_aalesund_fil    = "Output_files/Trondheim_Aalesund_reise"
    sailing_speed_Aalesund_Floro_fil    = "Output_files/Aalesund_Floro_reise"
    sailing_speed_Floro_Bergen_fil      = "Output_files/Floro_Bergen_reise"
    sailing_speed_Bergen_Stavanger_fil  = "Output_files/Bergen_Stavanger_reise"

    #write_to_file(sailing_speed_simulation_vector_TA2,sailing_speed_Trondheim_Aalesund_fil_2)

    write_to_file(sailing_speed_simulation_vector_TA,sailing_speed_trond_aalesund_fil)
    write_to_file(sailing_speed_simulation_vector_AF,sailing_speed_Aalesund_Floro_fil)
    write_to_file(sailing_speed_simulation_vector_FB,sailing_speed_Floro_Bergen_fil)
    write_to_file(sailing_speed_simulation_vector_BS,sailing_speed_Bergen_Stavanger_fil)


    #sailing_Speed = read_array_from_file(sailing_speed_trond_aalesund_fil)
    #print(np.average(sailing_Speed))
    return 0



#AWS (wind speed, vessel speed, wind direction)
WSN_temp            = 5
WSE_temp            = -5
vessel_speed_temp   = 0
vessel_heading_temp = 0
twd_temp            = True_wind_direction(vessel_heading_temp,WSN_temp,WSE_temp) #Heading, WSN, WSE
degrees             = r2d(twd_temp)
tws_temp            = True_wind_speed(WSN_temp,WSE_temp)
AWS_temp            = Apparent_Wind_Speed(tws_temp,vessel_speed_temp,twd_temp) #TWS, VS, Twd
print("AWS", AWS_temp)
print("TWS",tws_temp)
print("TWD",r2d(twd_temp))
sediek              = Apparent_wind_angle(tws_temp,AWS_temp,vessel_speed_temp)
egendefinert        = alpha(vessel_speed_temp,vessel_heading_temp,WSN_temp,WSE_temp)
print("apparent wind angle using alpha, egendefinert",egendefinert)
print("apparent wind angle using function from sediek",sediek)



print("Finished <3<3")

