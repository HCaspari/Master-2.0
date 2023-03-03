
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
#import searoute as sr
import chardet
import folium as folium
import webbrowser
#import gmplot as gmp
from file_handling import write_to_file
from pathlib import Path, PureWindowsPath

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
filename_AIS = PureWindowsPath(Path("../env/input_files/ais_data_v4.csv"))
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
A_rotor   = h * d  #cross sectional area of flettner
Cl_flettner  = 12.5
Cd_flettner  = 0.2
Cm_flettner  = 0.2
alpha_flettner = 3.5


#finds resistance by speed using admirality formula
def Sailing_resistance(sailing_speed):
    """
    :param sailing_speed: vessel sailing speed
    :return: sailing resistance force found from admirality formula in kN
    """
    admirality_coeff = 385.46
    power = (vessel_weight_disp**(2/3)*sailing_speed**3)/admirality_coeff
    force = round(power/(sailing_speed/5.44444444),2)
    #print(f"vessel power required at {sailing_speed} knots equals:{power} Watt")
    #print(f"resistance force required at {sailing_speed} knots equals:{force} kN")
    return force

#function that finds driftangle beta that is large enough to counteract drifting force
def Beta_solver(perp_force, vessel_velocity):
    """
    :param perp_force: perpendicular force observed by vessel
    :param vessel_velocity: vessel sailing speed
    :return: drift angle beta
    """
    Fy = perp_force * 1000 #Newton
    RHS = Fy/(0.5*vessel_length*vessel_draft*vessel_velocity**2)
    # Solving second degree function to find beta
    a = 0.461
    b = 1.091
    c = -RHS
    Beta_1 = (-b + np.sqrt(b ** 2 - 4 * a * c)) / (2 * a)
    Beta_2 = (-b - np.sqrt(b ** 2 - 4 * a * c)) / (2 * a)
    Beta = round(max(Beta_1, Beta_2),3)  # the angle will always be the positive solution
    return Beta

#function that finds Drift_resistance_multiplier resistance by beta (Drift angle)
def Drift_resistance_multiplier(beta):
    """
    :param beta: drift angle
    :return: drag resistance multiplying factor due to drifting
    """
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
    """
    :param vessel_velocity: speed of vessel
    :return: vector of forces that are produced when sailing with applied sideforces
    """
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
    filename_avg_righting_force = PureWindowsPath(Path("/Weather/env/old_code/output_files5/righting_force.txt"))
    write_to_file(Fy_vector,filename_avg_righting_force)
    return Fy_vector

def Force_produced(AWS, AWD):
    """
    :param AWS: Apparent wind speed in m/s
    :param AWD: Apparent wind direction (in degrees)
    :return: forward and perpendicular force produced by flettners (in kN)
    """

    lift = 0.5 * rho_air * A_rotor * AWS ** 2 * Cl_flettner                                 #lift force from traut
    drag = 0.5 * rho_air * A_rotor * AWS ** 2 * Cd_flettner                                 #drag force from traut
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


def Force_produced_new_CL(AWS, AWD):
    """
    :param AWS: Apparent wind speed in m/s
    :param AWD: Apparent wind direction (in degrees)
    :return: forward, perpendicular force produced by flettners (in kN), power in kW
    """
    SR = 3*np.pi*5/(2*AWS)

    #Found in Tillig Design operation and analysis of wind assisted carrier
    Cl = -0.0046*SR**5 + 0.1145*SR**4 - 0.9817*SR**3 + 3.1309*SR**2 - 0.1039*SR
    Cd = -0.0017*SR**5 + 0.0464*SR**4 - 0.4424*SR**3 + 1.7243*SR**2 - 1.641*SR + 0.6375
    Cp =  0.0001*SR**5 - 0.0004*SR**4 + 0.0143*SR**3 - 0.0168*SR**2 + 0.0234*SR
    lift = 0.5 * rho_air * A_rotor * AWS ** 2 * Cl  # lift force from traut
    drag = 0.5 * rho_air * A_rotor * AWS ** 2 * Cd  # drag force from traut

    power_consumption = Cp * (rho_air/2) * A_rotor * AWS**3       #Input to flettner is energy to spin rotors


    if 0 <= AWD <= 90 or 270 <= AWD <= 360:  # drag is set to negative if the wind is coming ahead, and positive if not
        drag *= -1
    Force_4_flettners = (lift + drag) * 4 / 1000  # kN                             #KiloNewton
    # The force generated by the flettner is the total output force minus the force needed to spin the rotor.

    perp_force = Force_4_flettners * abs(np.cos(AWD))  # perpendicular force flettner in KN
    if type(perp_force) != MaskedConstant:
        perp_force = round(min(perp_force, 1400),
                           2)  # Force is maximum 1400 kN, so as not to break flettners (cite: Norsepower technical)
    forward_force = Force_4_flettners * abs(np.sin(AWD))  # directional force flettner in KN

    if type(forward_force) != MaskedConstant:
        forward_force = round(min(forward_force, 1400), 2)  # Force is maximum 1400 kN, so as not to break flettners



    return forward_force, perp_force, power_consumption

#function interpolates over y values in list and finds closest value in table, returns x value (basically inverse func)
def take_closest(myList, myNumber):
    """
    :param myList: list of values
    :param myNumber: number to be found in said list
    :return: index of value in list closes to myNumber
    """

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

def Speed_achieved_old(perp_force, forward_force):
    """
    :param perp_force: sideforces observed by the vessel from flettners (in kN)
    :param forward_force:  propulsive force observed by vessel from flettners (in kN)
    :return: speed achieved by vessel when observing these forces
    """

    # ratio hydrodynamic resistance to total res is approximately 0.85
    ratio_hydrodyn_to_tot_res = 0.85

    # Empty vectors to store values
    sailing_resistance_vector = []
    total_resistance_vector = []

    # Set vessel speed [knots] in intervall from 0,20 with stepsize 0.1
    vessel_velocity = np.linspace(0.1, 20, 200)
    total_resistance = 0
    #
    for velocity in vessel_velocity:
        total_sailing_resistance = Sailing_resistance(velocity) * ratio_hydrodyn_to_tot_res
        sailing_resistance_vector.append(total_sailing_resistance)
        drift_angle = Beta_solver(perp_force, velocity)
        resistance_multiplier = Drift_resistance_multiplier(drift_angle)
        if resistance_multiplier < 1:
            resistance_multiplier = 1
            total_resistance        = resistance_multiplier * total_sailing_resistance
            total_resistance_vector.append(total_resistance)


        # stop iteration when total_resistance > forward force
        if total_resistance > forward_force:
            speed_achieved = velocity
            return speed_achieved

    speed_achieved = take_closest(total_resistance_vector, forward_force)  # IN KNOTS


    return speed_achieved #IN KNOTS
