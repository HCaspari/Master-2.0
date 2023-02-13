
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


from file_handling import write_to_file, read_route
from route_handling import calc_bearing
from Force_functions import Sailing_resistance,Beta_solver, Drift_resistance_multiplier, Force_produced, Speed_achieved
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
alpha_const = 3.5

#comparisson_vessel_power_by_speed = vessel_velocity**3*0.8*0.55





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
        vessel_heading      = calc_bearing(position_first,position_next)                                 #In Degrees (North is 0)
        WSE,WSN             = getweather(time,position_first[0],position_first[1])                       #Gives Wind speed East and North
        TWS                 = True_wind_speed(WSN,WSE)                              #Finds True Windspeed (pythagoras)
        TWD                 = True_wind_direction(vessel_heading,WSN,WSE)
        AWS                 = Apparent_Wind_Speed(TWS,vessel_speed,TWD)
        AWA                 = alpha(vessel_speed,vessel_heading,WSN,WSE )                                              #Finds Apparent wind angle
        forward_force_func,perpendicular_force_func = Force_produced(AWS, AWA)                         #Forward and Perpendicular force from Flettners
        if type(forward_force_func) != MaskedConstant or type(perpendicular_force_func) != MaskedConstant:
            vessel_speed    = Speed_achieved(perpendicular_force_func, forward_force_func)    #Sailing Speed obtained in KNOTS
            sailing_speed_vector.append(vessel_speed)
            sailing_time     = sailing_distance/vessel_speed                                                    #time used to sail trip added
            route_sailing_time.append(sailing_time)
        if type(forward_force_func) == MaskedConstant or type(perpendicular_force_func) == MaskedConstant:
            print("ouchie, we have a mask", i)
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





#main = time_of_trip,tot_sailing_dist, poor_sailing_time, poor_sailing_distance
def simulation(csv):
    """

    :param csv: A csv
    :return:    time_of_trip
                tot_sailing_dist
                sailing_speed_simulation_vector

    """
    initial_speed       = 4

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


def runsimulation(route):
    """
        :param route:   0 runs all routes
                        1 runs Trondheim Ålesund
                        2 runs Ålesund Florø
                        3 runs Florø Bergen
                        4 runs Bergen Stavanger
        :return: saved files with simulation results
        """

    Trond_aalesund      = "Route_data/Route_Trondheim_Aalesund.csv"
    Aalesund_Floro      = "Route_data/Route_Aalesund_floro.csv"
    Floro_Bergen        = "Route_data/Route_Floro_Bergen.csv"
    Bergen_Stavanger    = "Route_data/Route_Bergen_Stavanger.csv"

    if route == 1 or route == 0:
        Trip_time_vector_TA, Tot_sailing_distance_vector_TA, sailing_speed_simulation_vector_TA = simulation(Trond_aalesund)
        sailing_speed_trond_aalesund_fil    = "Output_files/Trondheim_Aalesund_reise_3"
        write_to_file(sailing_speed_simulation_vector_TA, sailing_speed_trond_aalesund_fil)

    if route == 2 or route == 0:
        Trip_time_vector_AF, Tot_sailing_distance_vector_AF, sailing_speed_simulation_vector_AF = simulation(Aalesund_Floro)
        sailing_speed_Aalesund_Floro_fil    = "Output_files/Aalesund_Floro_reise"
        write_to_file(sailing_speed_simulation_vector_AF, sailing_speed_Aalesund_Floro_fil)

    if route == 3 or route == 0:
        Trip_time_vector_FB, Tot_sailing_distance_vector_FB, sailing_speed_simulation_vector_FB = simulation(Floro_Bergen)
        sailing_speed_Floro_Bergen_fil      = "Output_files/Floro_Bergen_reise"
        write_to_file(sailing_speed_simulation_vector_FB, sailing_speed_Floro_Bergen_fil)

    if route == 4 or route == 0:
        Trip_time_vector_BS, Tot_sailing_distance_vector_BS, sailing_speed_simulation_vector_BS = simulation(Bergen_Stavanger)
        sailing_speed_Bergen_Stavanger_fil  = "Output_files/Bergen_Stavanger_reise"
        write_to_file(sailing_speed_simulation_vector_BS, sailing_speed_Bergen_Stavanger_fil)

    return 0



def test_func():
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
    return 0

runsimulation(1)
print("Finished <3<3")

