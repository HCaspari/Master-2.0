import numpy as np
import geopy.distance #package to calculate distance between two lat/lon points
import pandas as pd
from numpy.ma.core import MaskedConstant
from datetime import datetime
from file_handling import write_to_file, read_route, test_read_files, write_to_file_2, add_timestamp_to_dataframe
from Force_functions import Force_produced, Speed_achieved_Valid
from Weather_Handling import getweather, r2d, True_wind_direction, True_wind_speed, Apparent_Wind_Speed, Apparent_wind_angle, Apparent_Wind_Angle, add_hours_to_date

from route_handling import calc_vessel_heading, mac_windows_file_handle, calc_vessel_heading_2

timestart = datetime.now()

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
alpha_const = 3.5

mac_windows_file_handle(filename_AIS)


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


#Function that calculates time spent sailing from port a to port b, through route
def main(route, iteration, date_of_simulation,routenumber):
    """
    :param route: Vector of route coordinates [(x1,y1),(x2,y2),...,(xn,yn)]
    :param iteration: Number of repetitions of simulation

    :return:For each iteration of simulation:
            total_time_sailed_route: float
            tot_sailing_dist: vector of sailing distances over route
            poor_sailing_time: float, time sailed less than 1 knot
            poor_sailing_distance: float, distance sailed at less than 1 knot
            sailing_speed_vector: vector of speed sailed at each point along route
            TWS: Vector of true wind speed
            TWD: Vector of true wind direction
            Route sailing time: Time used to sail each route iteration
            Coordinate sailing time: Time used to sail between two points on route
            datestamp_vector: Vector containing times vessel is at each point over route
            battery_use_vector: Total batterypower neeed for sailing route this iteration

    """
    # initializing sailing speed of 4 knots (will change after one iteration)
    vessel_speed                    = 4
    route_sailing_time              = iteration
    tot_sailing_dist                = 0
    poor_sailing_time               = 0
    poor_sailing_distance           = 0
    sailing_speed_vector            = []
    coordinate_sailing_time         = []
    apparent_wind_speed_observed    = []
    true_wind_speed_vector          = []
    true_wind_direction_vector      = []
    datestamp_vector                = [date_of_simulation]
    battery_use_vector              = []
    forward_force_vect              = []
    total_time_sailed_route         = 0

    for i in range(len(route)-1): #create itteration through route
        position_first      = route[i]
        position_next       = route[i+1]
        sailing_distance    = geopy.distance.geodesic(position_first,position_next).nautical#In Nm
        vessel_heading      = calc_vessel_heading_2(position_first, position_next)
        WSE,WSN             = getweather(route_sailing_time,position_first[0],position_first[1])
        TWS                 = True_wind_speed(WSN,WSE)
        TWD                 = True_wind_direction(vessel_heading,WSN,WSE)
        AWS                 = Apparent_Wind_Speed(TWS,vessel_speed,TWD)
        AWA                 = Apparent_Wind_Angle(vessel_speed, vessel_heading, WSN, WSE)
        Forward_Force,Perp_Force = Force_produced(AWS, AWA)



        #With Kite added aswell as Flettner, assume 1.4 times more force generated
        if routenumber == 14:
            Forward_Force       = 1.4*Forward_Force
            if iteration % 1000 == 0 and i == 1:
                print("Running kite (1.4 force)")

        #With utilization of modern vessel design.
        if routenumber == 16 or routenumber == 151 or routenumber == 152:
            Forward_Force        = 2.5*Forward_Force


        #With Kite and Slimmer vessel aswell as Flettner,
        # assume 3.5 (2.5*1.4) times more force generated
        if routenumber == 17:
            if iteration % 1000 == 0 and i == 1:
                print("running Kite, and newbuild (3.5 force)")
            Forward_Force        = 3.5*Forward_Force


        vessel_speed,total_resistance,battery_need_power = \
            Speed_achieved_Valid(Perp_Force, Forward_Force)


        if type(Forward_Force) == MaskedConstant or type(Perp_Force) == MaskedConstant:
            print("ouchie, we have a mask", i)
            return 1

        #If wind observed equals zero, sailing time is set to one
        # and calculations are reated with new experienced wind at next intervall

        # If ship experiences no wind, it waits an hour recalculating forces
        if vessel_speed == 0:
            sailing_time = 1.00
        else:
            sailing_time     = round(sailing_distance/vessel_speed,3)

        # check for extreme time usage, whenever sailed speed is less than 1 knot.
        if vessel_speed < 1:
            poor_sailing_time += sailing_time
            poor_sailing_distance += sailing_distance

        #When batteries are not utilized, the batterypower calculated is set to 0
        if routenumber != 15 and routenumber != 152:
            battery_need_power = 0

        #When batteries are utilized, speeds below 2 knots are st to 2 knots
        if routenumber == 15 or routenumber == 152:
            if vessel_speed < 2:
                vessel_speed = 2


        #Append data to vectors
        sailing_speed_vector.append(vessel_speed)
        coordinate_sailing_time.append(sailing_time)
        true_wind_speed_vector.append(TWS)
        true_wind_direction_vector.append(TWD)
        apparent_wind_speed_observed.append(AWS)
        tot_sailing_dist += sailing_distance
        battery_use_vector.append(round(battery_need_power*(sailing_distance/2),3))
        forward_force_vect.append(round(Forward_Force,3))


        route_sailing_time += sailing_time  #Sailing time of total route
        if len(coordinate_sailing_time) != 0:
            total_time_sailed_route = sum(coordinate_sailing_time)

        # Adding timestamp to first column
        if i > 0:
            time_beginning  = datestamp_vector[i-1]
            time_new        = add_hours_to_date(time_beginning,sailing_time)
            datestamp_vector.append(time_new)

    return total_time_sailed_route,tot_sailing_dist, poor_sailing_time, \
           poor_sailing_distance, sailing_speed_vector,true_wind_speed_vector,\
           true_wind_direction_vector, route_sailing_time, coordinate_sailing_time, \
           datestamp_vector, battery_use_vector





#main = time_of_trip,tot_sailing_dist, poor_sailing_time, poor_sailing_distance
def simulation(csv,routenumber,interval):
    """
    :param csv: A csv
    :param routenumber: 1 runs Trondheim Ålesund
                        2 runs Ålesund Florø
                        3 runs Florø Bergen
                        4 runs Bergen Stavanger
                        5 runs Aberdeen Færøyene
                        6 runs Amsterdam Newcastle
                        7 runs Danmark AMsterdam
                        8 runs Færøyene Ålesund
                        81 runs Fær - åles retur
                        9 runs Newcastle Aberdeen
                        10 runs Ålesund Danmark
    :param interval: What intervall to write files
    :return:    time_of_trip
                tot_sailing_dist
                VS_simulation_vector
    """

    starttime = datetime(2020,7,1,00,00,00)
    iteration                       = 0
    hour_intervall                  = interval                                                 #at what hourly interval should we simulate?
    route_travel                    = read_route(csv)
    time_of_simulation              = 17520                                             #two years in hours
    time_of_trip                    = []
    tot_sailing_dist                = []
    poor_sailing_time               = []
    poor_sailing_distance           = []
    #VS_simulation_vector            = np.zeros(int(time_of_simulation/hour_intervall), dtype=object)
    datestamp_simulation_vector     = []
    VS_simulation_vector            = []
    TWS_simulation_vector           = []
    TWD_simulation_vector           = []
    battery_use_iteration_vect      = []
    #date_of_simulation              = starttime
    poor_sailing_speed = 0
    Batteryneed_simulation_Vector            = []

    for iteration in range(0,time_of_simulation,hour_intervall) : #repeating simulation for each hour_intervall through a year
        date_of_simulation = add_hours_to_date(starttime,iteration)
        time_of_trip_1,tot_sailing_dist_1, poor_sailing_time_1, poor_sailing_distance_1, sailing_speed_vector, TWS_vector, TWD_vector,\
        route_sailing_time, coordinate_sailing_time, datestamp_vector, battery_use_vector = main(route_travel,iteration, date_of_simulation,routenumber)
        time_of_trip.append(time_of_trip_1)
        tot_sailing_dist.append(tot_sailing_dist_1)
        poor_sailing_time.append(poor_sailing_time_1)
        poor_sailing_distance.append(poor_sailing_distance_1)
        VS_simulation_vector.extend(sailing_speed_vector)
        TWS_simulation_vector.extend(TWS_vector)
        TWD_simulation_vector.extend(TWD_vector)
        datestamp_simulation_vector.extend(datestamp_vector)
        Batteryneed_simulation_Vector.extend(battery_use_vector) #dette er totalt batteri brukt for denne iterasjonen

        battery_use_iteration_vect.extend(battery_use_vector)               #Dette er vectoren med batteribruk per punkt over ruten

        if iteration%100 == 0:
            print("progress is made, iteration:", iteration, "/", time_of_simulation," on route ", routenumber, "at time",  datetime.now())
        date_of_simulation = add_hours_to_date(date_of_simulation,hour_intervall)



    print(f"Average speed sailing {csv} over {iteration} iterations is {np.average(VS_simulation_vector[0])}")
    print(f"Throughout all iterations, the vessel sails less than one knot for an average of {np.average(poor_sailing_time)} hours\n"
          f"this iteration is used to sail on average {np.average(poor_sailing_distance)} nautical miles\n"
          f"at an average speed of {poor_sailing_speed} knots")
    print(f"the average batterypower needed for each iteration of the simulation is {np.average(Batteryneed_simulation_Vector)}"
          f"with a max used batterypower over one simulation of {max(Batteryneed_simulation_Vector)}"
          f"and a min used batterypower over one simulation of {min(Batteryneed_simulation_Vector)}")

#Run each function
    if routenumber == 151:
        print("here")
 # Trond_Ålesund

    if routenumber == 1:
        file_speed_Trond_Aalesund = mac_windows_file_handle("Output_files/Trondheim_Ålesund_iteration_1/savespeed_TrondAales.csv")
        file_TWS_Trond_Aalesund = mac_windows_file_handle("Output_files/Trondheim_Ålesund_iteration_1/saveTWS_TrondAales.csv")
        file_TWD_Trond_Aalesund = mac_windows_file_handle("Output_files/Trondheim_Ålesund_iteration_1/saveTWD_TrondAales.csv")
        datestamp_file_1 = mac_windows_file_handle("Output_files/Datestamps/datestamp1.csv")
        write_to_file(VS_simulation_vector, file_speed_Trond_Aalesund)  #VS vector file
        write_to_file(TWS_simulation_vector,file_TWS_Trond_Aalesund)    #TWS vector file
        write_to_file(TWD_simulation_vector,file_TWD_Trond_Aalesund)    #TWD vector file
        write_to_file(datestamp_simulation_vector, datestamp_file_1)                 #Datestamp file
        add_timestamp_to_dataframe(file_speed_Trond_Aalesund, datestamp_file_1)
        add_timestamp_to_dataframe(file_TWS_Trond_Aalesund, datestamp_file_1)
        add_timestamp_to_dataframe(file_TWD_Trond_Aalesund, datestamp_file_1)

    if routenumber == 110:
        file_speed_Trond_Aalesund = mac_windows_file_handle(
            "Output_files/Trondheim_Ålesund_iteration_10/savespeed_TrondAales.csv")
        file_TWS_Trond_Aalesund = mac_windows_file_handle(
            "Output_files/Trondheim_Ålesund_iteration_10/saveTWS_TrondAales.csv")
        file_TWD_Trond_Aalesund = mac_windows_file_handle(
            "Output_files/Trondheim_Ålesund_iteration_10/saveTWD_TrondAales.csv")
        datestamp_file_21 = mac_windows_file_handle("Output_files/Datestamps/datestamp21.csv")
        write_to_file(VS_simulation_vector, file_speed_Trond_Aalesund)  # VS vector file
        write_to_file(TWS_simulation_vector, file_TWS_Trond_Aalesund)  # TWS vector file
        write_to_file(TWD_simulation_vector, file_TWD_Trond_Aalesund)  # TWD vector file
        write_to_file(datestamp_simulation_vector, datestamp_file_21)  # Datestamp file
        add_timestamp_to_dataframe(file_speed_Trond_Aalesund, datestamp_file_21)
        add_timestamp_to_dataframe(file_TWS_Trond_Aalesund, datestamp_file_21)
        add_timestamp_to_dataframe(file_TWD_Trond_Aalesund, datestamp_file_21)

    if routenumber == 1100:
        file_speed_Trond_Aalesund = mac_windows_file_handle(
            "Output_files/Trondheim_Ålesund_iteration_100/savespeed_TrondAales.csv")
        file_TWS_Trond_Aalesund = mac_windows_file_handle(
            "Output_files/Trondheim_Ålesund_iteration_100/saveTWS_TrondAales.csv")
        file_TWD_Trond_Aalesund = mac_windows_file_handle(
            "Output_files/Trondheim_Ålesund_iteration_100/saveTWD_TrondAales.csv")
        datestamp_file_22 = mac_windows_file_handle("Output_files/Datestamps/datestamp22.csv")
        write_to_file(VS_simulation_vector, file_speed_Trond_Aalesund)  # VS vector file
        write_to_file(TWS_simulation_vector, file_TWS_Trond_Aalesund)  # TWS vector file
        write_to_file(TWD_simulation_vector, file_TWD_Trond_Aalesund)  # TWD vector file
        write_to_file(datestamp_simulation_vector, datestamp_file_22)  # Datestamp file
        add_timestamp_to_dataframe(file_speed_Trond_Aalesund, datestamp_file_22)
        add_timestamp_to_dataframe(file_TWS_Trond_Aalesund, datestamp_file_22)
        add_timestamp_to_dataframe(file_TWD_Trond_Aalesund, datestamp_file_22)

    if routenumber == 11000:
        file_speed_Trond_Aalesund = mac_windows_file_handle(
            "Output_files/Trondheim_Ålesund_iteration_1000/savespeed_TrondAales.csv")
        file_TWS_Trond_Aalesund = mac_windows_file_handle(
            "Output_files/Trondheim_Ålesund_iteration_1000/saveTWS_TrondAales.csv")
        file_TWD_Trond_Aalesund = mac_windows_file_handle(
            "Output_files/Trondheim_Ålesund_iteration_1000/saveTWD_TrondAales.csv")
        datestamp_file_23 = mac_windows_file_handle("Output_files/Datestamps/datestamp23.csv")
        write_to_file(VS_simulation_vector, file_speed_Trond_Aalesund)  # VS vector file
        write_to_file(TWS_simulation_vector, file_TWS_Trond_Aalesund)  # TWS vector file
        write_to_file(TWD_simulation_vector, file_TWD_Trond_Aalesund)  # TWD vector file
        write_to_file(datestamp_simulation_vector, datestamp_file_23)  # Datestamp file
        add_timestamp_to_dataframe(file_speed_Trond_Aalesund, datestamp_file_23)
        add_timestamp_to_dataframe(file_TWS_Trond_Aalesund, datestamp_file_23)
        add_timestamp_to_dataframe(file_TWD_Trond_Aalesund, datestamp_file_23)

# Ålesund_floro
    elif routenumber == 2:
        file_speed_Aalesund_Floro = mac_windows_file_handle("Output_files/Ålesund_Florø/savespeed_AalesFloro.csv")
        file_TWS_Aalesund_Floro = mac_windows_file_handle("Output_files/Ålesund_Florø/saveTWS_AalesFloro.csv")
        file_TWD_Aalesund_Floro = mac_windows_file_handle("Output_files/Ålesund_Florø/saveTWD_AalesFloro.csv")
        datestamp_file_2 = mac_windows_file_handle("Output_files/Datestamps/datestamp2.csv")

        write_to_file(VS_simulation_vector, file_speed_Aalesund_Floro)  #VS vector file
        write_to_file(TWS_simulation_vector,file_TWS_Aalesund_Floro)    #TWS vector file
        write_to_file(TWD_simulation_vector,file_TWD_Aalesund_Floro)    #TWD vector file
        write_to_file(datestamp_simulation_vector, datestamp_file_2)                 #Datestamp file
        add_timestamp_to_dataframe(file_speed_Aalesund_Floro, datestamp_file_2)
        add_timestamp_to_dataframe(file_TWS_Aalesund_Floro, datestamp_file_2)
        add_timestamp_to_dataframe(file_TWD_Aalesund_Floro, datestamp_file_2)
# Floro_bergen
    elif routenumber == 3:
        file_speed_Floro_Bergen = mac_windows_file_handle("Output_files/Florø_Bergen/savespeed_FloroBergen.csv")
        file_TWS_Floro_Bergen = mac_windows_file_handle("Output_files/Florø_Bergen/saveTWS_FloroBergen.csv")
        file_TWD_Floro_Bergen = mac_windows_file_handle("Output_files/Florø_Bergen/saveTWD_FloroBergen.csv")
        datestamp_file_3 = mac_windows_file_handle("Output_files/Datestamps/datestamp3.csv")
        write_to_file(VS_simulation_vector, file_speed_Floro_Bergen)    #VS vector file
        write_to_file(TWS_simulation_vector,file_TWS_Floro_Bergen)      #TWS vector file
        write_to_file(TWD_simulation_vector,file_TWD_Floro_Bergen)      #TWD vector file
        write_to_file(datestamp_simulation_vector, datestamp_file_3)                 #Datestamp file
        add_timestamp_to_dataframe(file_speed_Floro_Bergen, datestamp_file_3)
        add_timestamp_to_dataframe(file_TWS_Floro_Bergen, datestamp_file_3)
        add_timestamp_to_dataframe(file_TWD_Floro_Bergen, datestamp_file_3)
# Bergen_stavanger
    elif routenumber == 4:
        file_speed_Bergen_Stavanger = mac_windows_file_handle("Output_files/Bergen_Stavanger/savespeed_BrgStvg.csv")
        file_TWS_Bergen_Stavanger = mac_windows_file_handle("Output_files/Bergen_Stavanger/saveTWS_BergenStavanger.csv")
        file_TWD_Bergen_Stavanger = mac_windows_file_handle("Output_files/Bergen_Stavanger/saveTWD_BergenStavanger.csv")
        datestamp_file_4 = mac_windows_file_handle("Output_files/Datestamps/datestamp4.csv")

        write_to_file(VS_simulation_vector,file_speed_Bergen_Stavanger) #VS vector file
        write_to_file(TWS_simulation_vector,file_TWS_Bergen_Stavanger)  #TWS vector file
        write_to_file(TWD_simulation_vector,file_TWD_Bergen_Stavanger)  #TWD vector file
        write_to_file(datestamp_simulation_vector, datestamp_file_4)       #Datestamp file
        add_timestamp_to_dataframe(file_speed_Bergen_Stavanger, datestamp_file_4)
        add_timestamp_to_dataframe(file_TWS_Bergen_Stavanger, datestamp_file_4)
        add_timestamp_to_dataframe(file_TWD_Bergen_Stavanger, datestamp_file_4)
# Aberdeen_Færøyene
    elif routenumber == 5:
        file_speed_Aber_Faer = mac_windows_file_handle("Output_files/Aberdeen_Færøyene/savespeed_Aberdeen_Færøyene.csv")
        file_TWS_Aber_Faer = mac_windows_file_handle("Output_files/Aberdeen_Færøyene/saveTWD_Aberdeen_Færøyene.csv")
        file_TWD_Aber_Faer = mac_windows_file_handle("Output_files/Aberdeen_Færøyene/saveTWS_Aberdeen_Færøyene.csv")
        datestamp_file_5 = mac_windows_file_handle("Output_files/Datestamps/datestamp5.csv")

        write_to_file(VS_simulation_vector,file_speed_Aber_Faer) #VS vector file
        write_to_file(TWS_simulation_vector,file_TWS_Aber_Faer)  #TWS vector file
        write_to_file(TWD_simulation_vector,file_TWD_Aber_Faer)  #TWD vector file
        write_to_file(datestamp_simulation_vector, datestamp_file_5)          #Datestamp file
        add_timestamp_to_dataframe(file_speed_Aber_Faer, datestamp_file_5)
        add_timestamp_to_dataframe(file_TWS_Aber_Faer, datestamp_file_5)
        add_timestamp_to_dataframe(file_TWD_Aber_Faer, datestamp_file_5)
# Amsterdam Newcastle
    elif routenumber == 6:
        file_speed_Amst_New = mac_windows_file_handle(
            "Output_files/Amsterdam_Newcastle/savespeed_Amsterdam_Newcastle.csv")
        file_TWS_Amst_New = mac_windows_file_handle("Output_files/Amsterdam_Newcastle/saveTWD_Amsterdam_Newcastle.csv")
        file_TWD_Amst_New = mac_windows_file_handle("Output_files/Amsterdam_Newcastle/saveTWS_Amsterdam_Newcastle.csv")
        datestamp_file_6 = mac_windows_file_handle("Output_files/Datestamps/datestamp6.csv")

        write_to_file(VS_simulation_vector,file_speed_Amst_New) #VS vector file
        write_to_file(TWS_simulation_vector,file_TWS_Amst_New)  #TWS vector file
        write_to_file(TWD_simulation_vector,file_TWD_Amst_New)  #TWD vector file
        write_to_file(datestamp_simulation_vector, datestamp_file_6)         #Datestamp file
        add_timestamp_to_dataframe(file_speed_Amst_New, datestamp_file_6)
        add_timestamp_to_dataframe(file_TWS_Amst_New, datestamp_file_6)
        add_timestamp_to_dataframe(file_TWD_Amst_New, datestamp_file_6)
# Danmark Amsterdam
    elif routenumber == 7:
        file_speed_Dk_Amst = mac_windows_file_handle("Output_files/Danmark_Amsterdam/savespeed_Danmark_Amsterdam.csv")
        file_TWS_Dk_Amst = mac_windows_file_handle("Output_files/Danmark_Amsterdam/saveTWD_Danmark_Amsterdam.csv")
        file_TWD_Dk_Amst = mac_windows_file_handle("Output_files/Danmark_Amsterdam/saveTWS_Danmark_Amsterdam.csv")
        datestamp_file_7 = mac_windows_file_handle("Output_files/Datestamps/datestamp7.csv")

        write_to_file(VS_simulation_vector,file_speed_Dk_Amst) #VS vector file
        write_to_file(TWS_simulation_vector,file_TWS_Dk_Amst)  #TWS vector file
        write_to_file(TWD_simulation_vector,file_TWD_Dk_Amst)  #TWD vector file
        write_to_file(datestamp_simulation_vector, datestamp_file_7)        #Datestamp file
        add_timestamp_to_dataframe(file_speed_Dk_Amst, datestamp_file_7)
        add_timestamp_to_dataframe(file_TWS_Dk_Amst, datestamp_file_7)
        add_timestamp_to_dataframe(file_TWD_Dk_Amst, datestamp_file_7)
# Færøyene Ålesund with battery
    elif routenumber == 8:

        file_speed_Faer_Aal = mac_windows_file_handle(
            "Output_files/Færøyene_Ålesund_2x_with_BATTERY/savespeed_Færøyene_Ålesund.csv")
        file_TWS_Faer_Aal = mac_windows_file_handle(
            "Output_files/Færøyene_Ålesund_2x_with_BATTERY/saveTWD_Færøyene_Ålesund.csv")
        file_TWD_Faer_Aal = mac_windows_file_handle(
            "Output_files/Færøyene_Ålesund_2x_with_BATTERY/saveTWS_Færøyene_Ålesund.csv")
        datestamp_file_8 = mac_windows_file_handle("Output_files/Datestamps/datestamp8.csv")

        write_to_file(VS_simulation_vector,file_speed_Faer_Aal) #VS vector file
        write_to_file(TWS_simulation_vector,file_TWS_Faer_Aal)  #TWS vector file
        write_to_file(TWD_simulation_vector,file_TWD_Faer_Aal)  #TWD vector file
        write_to_file(datestamp_simulation_vector, datestamp_file_8)         #Datestamp file
        add_timestamp_to_dataframe(file_speed_Faer_Aal,datestamp_file_8)
        add_timestamp_to_dataframe(file_TWS_Faer_Aal,datestamp_file_8)
        add_timestamp_to_dataframe(file_TWD_Faer_Aal,datestamp_file_8)
#Faraoe Ålesund Return
    elif routenumber == 81:

        file_speed_Faer_Aal_return = mac_windows_file_handle("Output_files/Færøyene_Ålesund_return/savespeed_Færøyene_Ålesund_return.csv")
        #file_TWS_Faer_Aal_return = mac_windows_file_handle("Output_files/Færøyene_Ålesund_return/saveTWS_Færøyene_Ålesund_return.csv")
        #file_TWD_Faer_Aal_return = mac_windows_file_handle("Output_files/Færøyene_Ålesund_return/saveTWD_Færøyene_Ålesund_return.csv")
        datestamp_file_8 = mac_windows_file_handle("Output_files/Datestamps/datestamp8.csv")

        write_to_file(VS_simulation_vector,file_speed_Faer_Aal_return) #VS vector file
        #write_to_file(TWS_simulation_vector,file_TWS_Faer_Aal_return)  #TWS vector file
        #write_to_file(TWD_simulation_vector,file_TWD_Faer_Aal_return)  #TWD vector file
        write_to_file(datestamp_simulation_vector, datestamp_file_8)         #Datestamp file
        add_timestamp_to_dataframe(file_speed_Faer_Aal_return,datestamp_file_8)
        #add_timestamp_to_dataframe(file_TWS_Faer_Aal_return,datestamp_file_8)
        #add_timestamp_to_dataframe(file_TWD_Faer_Aal_return,datestamp_file_8)
# Newcastle Aberdeen
    elif routenumber == 9:

        file_speed_New_Aber = mac_windows_file_handle("Output_files/Newcastle_Aberdeen/SS_Newcastle_Aberdeen.csv")
        file_TWS_New_Aber = mac_windows_file_handle("Output_files/Newcastle_Aberdeen/TWD_Newcastle_Aberdeen.csv")
        file_TWD_New_Aber = mac_windows_file_handle("Output_files/Newcastle_Aberdeen/TWS_Newcastle_Aberdeen.csv")
        datestamp_file_9 = mac_windows_file_handle("Output_files/Datestamps/datestamp9.csv")

        write_to_file(VS_simulation_vector,file_speed_New_Aber) #VS vector file
        write_to_file(TWS_simulation_vector,file_TWS_New_Aber)  #TWS vector file
        write_to_file(TWD_simulation_vector,file_TWD_New_Aber)  #TWD vector file
        write_to_file(datestamp_simulation_vector, datestamp_file_9)         #Datestamp file
        add_timestamp_to_dataframe(file_speed_New_Aber, datestamp_file_9)
        add_timestamp_to_dataframe(file_TWS_New_Aber, datestamp_file_9)
        add_timestamp_to_dataframe(file_TWD_New_Aber, datestamp_file_9)
# Ålesund Danmark
    elif routenumber == 10:


        file_speed_Aal_Dk = mac_windows_file_handle("Output_files/Ålesund_Danmark/savespeed_Ålesund_Danmark.csv")
        file_TWS_Aal_Dk = mac_windows_file_handle("Output_files/Ålesund_Danmark/saveTWD_Ålesund_Danmark.csv")
        file_TWD_Aal_Dk = mac_windows_file_handle("Output_files/Ålesund_Danmark/saveTWS_Ålesund_Danmark.csv")
        datestamp_file_10 = mac_windows_file_handle("Output_files/Datestamps/datestamp10.csv")

        write_to_file(VS_simulation_vector,file_speed_Aal_Dk) #VS vector file
        write_to_file(TWS_simulation_vector,file_TWS_Aal_Dk)  #TWS vector file
        write_to_file(TWD_simulation_vector,file_TWD_Aal_Dk)  #TWD vector file
        write_to_file(datestamp_simulation_vector, datestamp_file_10)       #Datestamp file
        add_timestamp_to_dataframe(file_speed_Aal_Dk, datestamp_file_10)
        add_timestamp_to_dataframe(file_TWS_Aal_Dk, datestamp_file_10)
        add_timestamp_to_dataframe(file_TWD_Aal_Dk, datestamp_file_10)
# Floro port
    elif routenumber == 11:


        file_speed_Floro_port = mac_windows_file_handle("Output_files/Floro_port/savespeed_port.csv")
        file_TWS_Floro_port = mac_windows_file_handle("Output_files/Floro_port/TWD_Floro_port.csv")
        file_TWD_Floro_port = mac_windows_file_handle("Output_files/Floro_port/TWS_Floro_port.csv")
        datestamp_file_11 = mac_windows_file_handle("Output_files/Datestamps/datestamp11.csv")

        write_to_file(VS_simulation_vector,file_speed_Floro_port) #VS vector file
        write_to_file(TWS_simulation_vector,file_TWS_Floro_port)  #TWS vector file
        write_to_file(TWD_simulation_vector,file_TWD_Floro_port)  #TWD vector file
        write_to_file(datestamp_simulation_vector, datestamp_file_11)          #Datestamp file
        add_timestamp_to_dataframe(file_speed_Floro_port, datestamp_file_11)
        add_timestamp_to_dataframe(file_TWS_Floro_port, datestamp_file_11)
        add_timestamp_to_dataframe(file_TWD_Floro_port, datestamp_file_11)
# Poor route
    elif routenumber == 12:

        file_speed_Poor = mac_windows_file_handle("Output_files/Poor_Route/savespeed_Poor_Route.csv")
        file_TWS_Poor = mac_windows_file_handle("Output_files/Poor_Route/saveTWS_Poor_Route.csv")
        file_TWD_Poor = mac_windows_file_handle("Output_files/Poor_Route/saveTWD_Poor_Route.csv")
        file_battery  = mac_windows_file_handle("Output_files/Poor_Route/Battery_use.csv")
        datestamp_file_12 = mac_windows_file_handle("Output_files/Datestamps/datestamp12.csv")

        write_to_file(VS_simulation_vector, file_speed_Poor)  # VS vector file
        write_to_file(TWS_simulation_vector, file_TWS_Poor)  # TWS vector file
        write_to_file(TWD_simulation_vector, file_TWD_Poor)  # TWD vector file
        write_to_file(Batteryneed_simulation_Vector,file_battery)
        write_to_file(datestamp_simulation_vector, datestamp_file_12)  # Datestamp file
        add_timestamp_to_dataframe(file_speed_Poor, datestamp_file_12)
        add_timestamp_to_dataframe(file_TWS_Poor, datestamp_file_12)
        add_timestamp_to_dataframe(file_TWD_Poor, datestamp_file_12)
        add_timestamp_to_dataframe(file_battery,datestamp_file_12)
# Good route
    elif routenumber == 13:
        file_speed_Fav = mac_windows_file_handle("Output_files/Good_Route/savespeed_Good_Route.csv")
        file_TWS_Fav = mac_windows_file_handle("Output_files/Good_Route/saveTWS_Good_Route.csv")
        file_TWD_Fav = mac_windows_file_handle("Output_files/Good_Route/saveTWD_Good_Route.csv")
        file_battery = mac_windows_file_handle("Output_files/Good_Route/Battery_use.csv")
        datestamp_file_13 = mac_windows_file_handle("Output_files/Datestamps/datestamp13.csv")
        write_to_file(VS_simulation_vector, file_speed_Fav)  # VS vector file
        write_to_file(TWS_simulation_vector, file_TWS_Fav)  # TWS vector file
        write_to_file(TWD_simulation_vector, file_TWD_Fav)  # TWD vector file
        write_to_file(datestamp_simulation_vector, datestamp_file_13)  # Datestamp file
        write_to_file(Batteryneed_simulation_Vector,file_battery)
        add_timestamp_to_dataframe(file_speed_Fav, datestamp_file_13)
        add_timestamp_to_dataframe(file_TWS_Fav, datestamp_file_13)
        add_timestamp_to_dataframe(file_TWD_Fav, datestamp_file_13)
        add_timestamp_to_dataframe(file_battery,datestamp_file_13)
# Good route with Kite
# Set forcefunction to 1.4
    elif routenumber == 14:
        file_speed_Fav = mac_windows_file_handle("Output_files/Færøyene_Ålesund_With_Kite/savespeed_Good_Route.csv")
        file_TWS_Fav = mac_windows_file_handle("Output_files/Færøyene_Ålesund_With_Kite/saveTWS_Good_Route.csv")
        file_TWD_Fav = mac_windows_file_handle("Output_files/Færøyene_Ålesund_With_Kite/saveTWD_Good_Route.csv")
        file_battery = mac_windows_file_handle("Output_files/Færøyene_Ålesund_With_Kite/Battery_use.csv")
        datestamp_file_14 = mac_windows_file_handle("Output_files/Datestamps/datestamp14.csv")
        write_to_file(VS_simulation_vector, file_speed_Fav)  # VS vector file
        write_to_file(TWS_simulation_vector, file_TWS_Fav)  # TWS vector file
        write_to_file(TWD_simulation_vector, file_TWD_Fav)  # TWD vector file
        write_to_file(datestamp_simulation_vector, datestamp_file_14)  # Datestamp file
        write_to_file(Batteryneed_simulation_Vector,file_battery)
        add_timestamp_to_dataframe(file_speed_Fav, datestamp_file_14)
        add_timestamp_to_dataframe(file_TWS_Fav, datestamp_file_14)
        add_timestamp_to_dataframe(file_TWD_Fav, datestamp_file_14)
        add_timestamp_to_dataframe(file_battery,datestamp_file_14)
# Færøyene with Battery
# Set batteryfunction on
    elif routenumber == 15:
        file_speed_Fav = mac_windows_file_handle("Output_files/Færøyene_Ålesund_With_Battery/savespeed_Good_Route.csv")
        file_TWS_Fav = mac_windows_file_handle("Output_files/Færøyene_Ålesund_With_Battery/saveTWS_Good_Route.csv")
        file_TWD_Fav = mac_windows_file_handle("Output_files/Færøyene_Ålesund_With_Battery/saveTWD_Good_Route.csv")
        file_battery = mac_windows_file_handle("Output_files/Færøyene_Ålesund_With_Battery/Battery_use.csv")
        datestamp_file_15 = mac_windows_file_handle("Output_files/Datestamps/datestamp15.csv")
        write_to_file(VS_simulation_vector, file_speed_Fav)  # VS vector file
        write_to_file(TWS_simulation_vector, file_TWS_Fav)  # TWS vector file
        write_to_file(TWD_simulation_vector, file_TWD_Fav)  # TWD vector file
        write_to_file(datestamp_simulation_vector, datestamp_file_15)  # Datestamp file
        write_to_file(battery_use_iteration_vect,file_battery)
        add_timestamp_to_dataframe(file_speed_Fav, datestamp_file_15)
        add_timestamp_to_dataframe(file_TWS_Fav, datestamp_file_15)
        add_timestamp_to_dataframe(file_TWD_Fav, datestamp_file_15)
        add_timestamp_to_dataframe(file_battery,datestamp_file_15)

    elif routenumber == 151:
        file_speed_Fav = mac_windows_file_handle("Output_files/Færøyene_Ålesund_Newbuild/savespeed_Færøyene_Ålesund.csv")
        file_battery = mac_windows_file_handle("Output_files/Færøyene_Ålesund_With_Battery/Battery_use.csv")
        datestamp_file_18 = mac_windows_file_handle("Output_files/Datestamps/datestamp18.csv")
        write_to_file(VS_simulation_vector, file_speed_Fav)  # VS vector file
        write_to_file(datestamp_simulation_vector, datestamp_file_18)  # Datestamp file
        write_to_file(Batteryneed_simulation_Vector, file_battery)
        add_timestamp_to_dataframe(file_speed_Fav, datestamp_file_18)
        add_timestamp_to_dataframe(file_battery, datestamp_file_18)

# Færøyene with Battery and newbuild
# Set batteryfunction and newbuild
    elif routenumber == 152:
        file_speed_Fav = mac_windows_file_handle("Output_files/Færøyene_Ålesund_Newbuild_Battery/savespeed_Good_Route.csv")
        file_battery = mac_windows_file_handle("Output_files/Færøyene_Ålesund_Newbuild_Battery/Battery_use.csv")
        datestamp_file_20 = mac_windows_file_handle("Output_files/Datestamps/datestamp20.csv")

        write_to_file(datestamp_simulation_vector, datestamp_file_20)  # Datestamp file
        write_to_file(VS_simulation_vector, file_speed_Fav)  # VS vector file
        write_to_file(Batteryneed_simulation_Vector, file_battery)

        add_timestamp_to_dataframe(file_speed_Fav, datestamp_file_20)
        add_timestamp_to_dataframe(file_battery, datestamp_file_20)

    # Good Route With Newbuild:
# Set Forcefunction to 2.5
    elif routenumber == 16:
        file_speed_Fav = mac_windows_file_handle("Output_files/Good_Route_Newbuild/savespeed_Good_Route.csv")
        file_TWS_Fav = mac_windows_file_handle("Output_files/Good_Route_Newbuild/saveTWS_Good_Route.csv")
        file_TWD_Fav = mac_windows_file_handle("Output_files/Good_Route_Newbuild/saveTWD_Good_Route.csv")
        file_battery = mac_windows_file_handle("Output_files/Good_Route_Newbuild/Battery_use.csv")
        datestamp_file_16 = mac_windows_file_handle("Output_files/Datestamps/datestamp16.csv")
        write_to_file(VS_simulation_vector, file_speed_Fav)  # VS vector file
        write_to_file(TWS_simulation_vector, file_TWS_Fav)  # TWS vector file
        write_to_file(TWD_simulation_vector, file_TWD_Fav)  # TWD vector file
        write_to_file(datestamp_simulation_vector, datestamp_file_16)  # Datestamp file
        write_to_file(Batteryneed_simulation_Vector, file_battery)
        add_timestamp_to_dataframe(file_speed_Fav, datestamp_file_16)
        add_timestamp_to_dataframe(file_TWS_Fav, datestamp_file_16)
        add_timestamp_to_dataframe(file_TWD_Fav, datestamp_file_16)
        add_timestamp_to_dataframe(file_battery, datestamp_file_16)

# Good Route With Flettner, Kite and Battery and Newbuild
 #Set forcefunction to 2.5*1.4, add battery option
    elif routenumber == 17:
        file_speed_Fav = mac_windows_file_handle("Output_files/Good_Route_With_Kite_Battery_Newbuild/savespeed_Good_Route.csv")
        file_TWS_Fav   = mac_windows_file_handle("Output_files/Good_Route_With_Kite_Battery_Newbuild/saveTWS_Good_Route.csv")
        file_TWD_Fav   = mac_windows_file_handle("Output_files/Good_Route_With_Kite_Battery_Newbuild/saveTWD_Good_Route.csv")
        file_battery   = mac_windows_file_handle("Output_files/Good_Route_With_Kite_Battery_Newbuild/Battery_use.csv")
        datestamp_file_17 = mac_windows_file_handle("Output_files/Datestamps/datestamp17.csv")
        write_to_file(VS_simulation_vector, file_speed_Fav)  # VS vector file
        write_to_file(TWS_simulation_vector, file_TWS_Fav)  # TWS vector file
        write_to_file(TWD_simulation_vector, file_TWD_Fav)  # TWD vector file
        write_to_file(datestamp_simulation_vector, datestamp_file_17)  # Datestamp file
        write_to_file(Batteryneed_simulation_Vector, file_battery)
        add_timestamp_to_dataframe(file_speed_Fav, datestamp_file_17)
        add_timestamp_to_dataframe(file_TWS_Fav, datestamp_file_17)
        add_timestamp_to_dataframe(file_TWD_Fav, datestamp_file_17)
        add_timestamp_to_dataframe(file_battery, datestamp_file_17)
#Newcastle Ålesund
    elif routenumber == 18:
        file_speed_Fav = mac_windows_file_handle("Output_files/Newcastle_Ålesund/SS_Newcastle_Ålesund.csv")
        datestamp_file_17 = mac_windows_file_handle("Output_files/Datestamps/datestampNew.csv")
        write_to_file(VS_simulation_vector, file_speed_Fav)  # VS vector file
        write_to_file(datestamp_simulation_vector, datestamp_file_17)  # Datestamp file
        add_timestamp_to_dataframe(file_speed_Fav, datestamp_file_17)

    return 0

#change to simulation(route) so that you dont need to hardcode

def runsimulation(route, interval):
    """
        :param route:   1 runs Trondheim Ålesund
                        2 runs Ålesund Florø
                        3 runs Florø Bergen
                        4 runs Bergen Stavanger
                        5 runs Aberdeen Færøyene
                        6 runs Amsterdam Newcastle
                        7 runs Danmark AMsterdam
                        8 runs Færøyene Ålesund
                        9 runs Newcastle Aberdeen
                        10 runs Ålesund Danmark
                        11 runds floro port
        :param interval: At what time interval do we run simulation
        :return: saved files with simulation results
        """

    if route == 1:
        Trond_aalesund = mac_windows_file_handle("Route_data/Trondheim_Ålesund_Route/Route_Trond_Aales.csv")
        print("Running simulation for route Trondheim Aalesund now")
        simulation(Trond_aalesund,route,interval)
        print("Simulation Trondheim to Ålesund is now complete.\n")
    if route == 110:
        Trond_aalesund = mac_windows_file_handle("Route_data/Trondheim_Ålesund_Route/Route_Trond_Aales.csv")
        print("Running simulation for route Trondheim Aalesund now 10")
        simulation(Trond_aalesund,route,interval)
        print("Simulation Trondheim to Ålesund is now complete. 10 \n")
    if route == 1100:
        Trond_aalesund = mac_windows_file_handle("Route_data/Trondheim_Ålesund_Route/Route_Trond_Aales.csv")
        print("Running simulation for route Trondheim Aalesund now 100")
        simulation(Trond_aalesund,route,interval)
        print("Simulation Trondheim to Ålesund is now complete. 100\n")
    if route == 11000:
        Trond_aalesund = mac_windows_file_handle("Route_data/Trondheim_Ålesund_Route/Route_Trond_Aales.csv")
        print("Running simulation for route Trondheim Aalesund now 1000")
        simulation(Trond_aalesund,route,interval)
        print("Simulation Trondheim to Ålesund is now complete. 1000 \n")
    if route == 2:
        Aalesund_Floro      = mac_windows_file_handle("Route_data/Ålsund_Florø_Rute/Route_Ålesund_Florø.csv")
        print("Running simulation for route Ålesund Florø now")
        simulation(Aalesund_Floro,route,interval)
        print("Simulation Ålesund to Florø is now complete.\n")
    if route == 3:
        Floro_Bergen        = mac_windows_file_handle("Route_data/Florø_Bergen_Rute/Route_Floro_Bergen.csv")
        print("Running simulation for route Florø Bergen now")
        simulation(Floro_Bergen,route,interval)
        print("Simulation Florø to Bergen d is now complete.\n")
    if route == 4:
        Bergen_Stavanger    = mac_windows_file_handle("Route_data/Bergen_Stavanger_Route/Route_Bergen_Stavanger.csv")
        print("Running simulation for route Bergen Stavanger now")
        simulation(Bergen_Stavanger,route,interval)
        print("Simulation Bergen to Stavanger is now complete.\n")
    if route == 5:
        Aberdeen_Faer       = mac_windows_file_handle("Route_data/Aberdeen_Færøyene_Route/Route_Aberdeen_Færøyene.csv")
        print("Running simulation for route Aberdeen_Faerøyene now")
        simulation(Aberdeen_Faer, route,interval)
        print("Simulation Aberdeen to Færøyene is now complete.\n")
    if route == 6:
        Amst_New            = mac_windows_file_handle("Route_data/Amsterdam_Newcastle_Route/Route_Amsterdam_Newcastle.csv")
        print("Running simulation for route Amsterdam_Newcastle now")
        simulation(Amst_New, route,interval)
        print("Simulation Amsterdam to Newcastle is now complete.\n")
    if route == 7:
        DK_Amst             = mac_windows_file_handle("Route_data/Danmark_Amsterdam_Route/Route_Danmark_Amsterdam.csv")
        print("Running simulation for route Danmark_Amsterdam now")
        simulation(DK_Amst, route,interval)
        print("Simulation Danmark to Amsterdam is now complete.\n")
    if route == 8:
        Faer_Aale        = mac_windows_file_handle("Route_data/Færøyene_Ålesund_Route/Route_Færøyene_Ålesund.csv")
        print("Running simulation for route Færøyene_Ålesund now")
        simulation(Faer_Aale, route,interval)
        print("Simulation Færøyene to Ålesund is now complete.\n")
    if route == 81:
        Faer_Aale_return = mac_windows_file_handle("Route_data/Færøyene_Ålesund_Return_Route/Route_Færøyene_Ålesund_retur.csv")
        print("Running simulation for route Færøyene_Ålesund Return now")
        simulation(Faer_Aale_return, route, interval)
        print("Simulation Færøyene to Ålesund REturn is now complete.\n")
    if route == 9:
        New_Aber            = mac_windows_file_handle("Route_data/Newcastle_Aberdeen_Route/Route_Newcastle_Ålesund.csv")
        print("Running simulation for route Newcastle_Aberdeen now")
        simulation(New_Aber, route,interval)
        print("Simulation Newcastle to Aberdeen is now complete.\n")
    if route == 10:
        Aale_DK             = mac_windows_file_handle("Route_data/Ålesund_Danmark_Route/Route_Ålesund_DK.csv")
        print("Running simulation for route Ålesund_Danmark now")
        simulation(Aale_DK, route,interval)
        print("Simulation Ålesund to Danmark is now complete.\n")
    if route == 11:
        Floro_port          = mac_windows_file_handle("Route_data/Floro_port/Floro_port.csv")
        print("Running simulation for route in florø port now")
        simulation(Floro_port, route, interval)
        print("Simulation for route in florø port is now complete")
# Poor Route
    if route == 12:
        Poor_route          = mac_windows_file_handle("Route_data/Poor_Route/Poor_route.csv")
        print("Running simulation for poor route now")
        simulation(Poor_route, route, interval)
        print("Simulation for Poor Route  is now complete")
# Good route
    if route == 13:
        Good_route    = mac_windows_file_handle("Route_data/Good_Route/Good_route.csv")
        print("Running simulation for Good route now")
        simulation(Good_route, route, interval)
        print("Simulation for Good route is now complete")
# Færøyene Ålesund  with Kite
    if route == 14:
        Kite   = mac_windows_file_handle("Route_data/Færøyene_Ålesund_Route/Route_Færøyene_Ålesund.csv")
        print("Running simulation Ålesund - Færøyene with Kite now")
        simulation(Kite, route, interval)
        print("Simulation for Kite is now complete")
# Færøyene Ålesund Route With Battery
    if route == 15:
        Battery_route    = mac_windows_file_handle("Route_data/Færøyene_Ålesund_Route/Route_Færøyene_Ålesund.csv")
        print("Running simulation for Battery now")
        simulation(Battery_route, route, interval)
        print("Simulation for Battery is now complete")
# Færøyene Ålesund Route With Newbuild
    if route == 151:
        Fær_åles_New = mac_windows_file_handle("Route_data/Færøyene_Ålesund_Route/Route_Færøyene_Ålesund.csv")
        print("Running simulation for Battery now")
        simulation(Fær_åles_New, route, interval)
        print("Simulation for Battery is now complete")
# Færøyene Ålesund Route With Newbuild
    if route == 152:
        Fær_Åles_with_Bat = mac_windows_file_handle("Route_data/Færøyene_Ålesund_Route/Route_Færøyene_Ålesund.csv")
        print("Running simulation for Battery now")
        simulation(Fær_Åles_with_Bat, route, interval)
        print("Simulation for Battery is now complete")
# Good Route With Newbuild
    if route == 16:
        Newbuild = mac_windows_file_handle("Route_data/Good_Route/Good_route.csv")
        print("Running simulation for Newbuild now")
        simulation(Newbuild, route, interval)
        print("Simulation for Newbuild is now complete")
# Good Route With Flettner, Kite and Battery and Newbuild
    if route == 17:
        Flettner_kite_battery_Newbuild    = mac_windows_file_handle("Route_data/Good_Route/Good_route.csv")
        print("Running simulation for Everything now")
        simulation(Flettner_kite_battery_Newbuild, route, interval)
        print("Simulation for Everythin is now complete")
    if route == 18:
        Newcastle_Åles = mac_windows_file_handle("Route_data/Newcastle_Ålesund_Route/Route_Newcastle_Ålesund.csv")
        print("Running simulation for Newcastle to ålesund now")
        simulation(Newcastle_Åles, route, interval)
        print("Simulation for Newcastle to Ålesund is now complete")
    return 0









#reset_index()
steps = 10
print(f"Everything is working 1241, 16/03 running simulation with steps {steps}")
#runsimulation(1,1)
#runsimulation(110,10)
#runsimulation(1100,100)
#runsimulation(11000,1000)
#runsimulation(2,steps)
#runsimulation(3,steps)
#runsimulation(4,steps)
#runsimulation(5,steps)
#runsimulation(6,steps)
#runsimulation(7,steps)
#runsimulation(8,steps) # Færøyene to Ålesund
#runsimulation(81,steps) # Færøyene to Ålesund Return
#runsimulation(9,steps)
#runsimulation(10,steps)
#runsimulation(11,steps)

#runsimulation(12,steps) #normal bad route planning

#runsimulation(15,steps) # Færøyene Ålesund battery
#runsimulation(16,steps) #Newbuild, change to 2.5 force
#runsimulation(13,steps) #normal good route planning
#runsimulation(151,steps) # Fær åles newbuild
#runsimulation(152,steps) #Fær åles newbuild with battery
#runsimulation(17,10)
#runsimulation(18,10)


timeend = datetime.now()
print("Finished <3<3")
print("Bonus print: Mathias er digg, kl. 1501, den 28/04/2023")
print(f"time used is {timeend-timestart}")


