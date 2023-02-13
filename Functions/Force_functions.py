
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


#finds resistance by speed using admirality formula
def resistance_by_speed(sailing_speed):
    admirality_coeff = 385.46
    power = (vessel_weight_disp**(2/3)*sailing_speed**3)/admirality_coeff
    force = round(power/(sailing_speed/5.44444444),2)
    #print(f"vessel power required at {sailing_speed} knots equals:{power} Watt")
    #print(f"resistance force required at {sailing_speed} knots equals:{force} kN")
    return force

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
