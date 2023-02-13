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

from plot_functions import plot_power,plot_avg_power,plot_resistance,plot_weekly_and_daily_avg_power
from file_handling import write_to_file, read_cols, read_position_vect_from_file
from Force_functions import resistance_by_speed,solve_beta_by_perp_force, xfold, iterate_drift_angle
from Old_route_calc_funcs import Force_over_year,Force_over_trip

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
Trondheim_location = 63.437686821303096, 10.402184694640052
Aalesund_location  = 62.93245830958637, 6.3481997169859055


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
    position_array_file  = "../env/position_array"
    write_to_file(position_array_func,position_array_file)
    return position_array_func

def generate_intermediate_points(start_point, end_point, num_points):
    """
    Generates intermediate points between two given coordinate points using linear interpolation.

    Parameters:
        start_point (tuple): The starting coordinate point, in the form (latitude, longitude)
        end_point (tuple): The ending coordinate point, in the form (latitude, longitude)
        num_points (int): The number of intermediate points to generate

    Returns:
        list: A list of intermediate points, in the form [(latitude, longitude), ...]
    """
    intermediate_points = []
    latitude_delta = (end_point[0] - start_point[0]) / (num_points + 1)
    longitude_delta = (end_point[1] - start_point[1]) / (num_points + 1)

    for i in range(1, num_points + 1):
        latitude = start_point[0] + i * latitude_delta
        longitude = start_point[1] + i * longitude_delta
        intermediate_points.append((latitude, longitude))

    return intermediate_points

def generate_intricate_route(route,points):
    """
    :param route: vector of positions over route
    :param points: number of points between separate positions
    :return: detailed route as vector of coordinates
    """
    newroute = [route[0][:]]
    for position in range(8):
        start_position  = route[position]
        end_position    = route[position+1]
        intermediate_points = generate_intermediate_points(start_position,end_position, points)
        newroute.append(intermediate_points)
    newroute.append(end_position)
    return newroute
#long_route_trond_aal = generate_intricate_route(route_Trond_Aal, 10)#
#print(long_route_trond_aal)
#T_A_Oppgradert = generate_intricate_route(route_Trond_Aal,5)

def createmap(trip_vector):
    """
    :param trip_vector: A vector of coordinates over route
    :return: a map that shows given route
    """

    foliummap = folium.Map(location=Trondheim_location, tiles="Stamen Terrain", zoom_start=9)
    #foliummap.show_in_browser()

    route = trip_vector
    for coordinates in route:
        foliummap.add_child(
                folium.Marker(
                    location=coordinates,
                    icon=folium.Icon(color="%s" % "blue"),
                )
        )
        foliummap.add_child(
            folium.Marker(
                location=Aalesund_location,
                icon=folium.Icon(color="%s" % "blue"),
            )
        )

    foliummap.show_in_browser()
    return 0
#createmap(T_A_Oppgradert)
