import numpy as np
import geopy.distance #package to calculate distance between two lat/lon points
import pandas as pd
from numpy.ma.core import MaskedConstant
from datetime import datetime
from file_handling import write_to_file, read_route, test_read_files, write_to_file_2, add_timestamp_to_dataframe
from Force_functions import Force_produced, Speed_achieved_old
from Weather_Handling import getweather, r2d, True_wind_direction, True_wind_speed, Apparent_Wind_Speed, Apparent_wind_angle, alpha, add_hours_to_date

from route_handling import calc_vessel_heading, mac_windows_file_handle, calc_vessel_heading_2



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


#Function that calculates time spent sailing from port a to port b, through predetermined route
def main(route, iteration, date_of__return):
    """
    :param route: Vector of route coordinates [(x1,y1),(x2,y2),...,(xn,yn)]
    :param iteration: Number of repetitions of _return

    :return:    For each iteration of _return:
                total_time_sailed_route: float
                tot_sailing_dist: vector of sailing distances over route
                poor_sailing_time: float, time sailed less than 1 knot
                poor_sailing_distance: float, distance sailed at less than 1 knot
                sailing_speed_vector: vector of speed sailed at each point along route
                TWS: Vector of true wind speed
                TWD: Vector of true wind direction

    """


    vessel_speed                    = 4 #initializing sailing speed of 4 knots (will change after one iteration)
    route_sailing_time              = iteration
    tot_sailing_dist                = 0
    poor_sailing_time               = 0
    poor_sailing_distance           = 0
    sailing_speed_vector            = []
    coordinate_sailing_time         = []
    apparent_wind_speed_observed    = []
    true_wind_speed_vector          = []
    true_wind_direction_vector      = []
    datestamp_vector                = [date_of__return]

    #np.append(tot_sailing_dist,[1])

    for i in range(len(route)-1): #create itteration through route

        position_first      = route[i]
        position_next       = route[i+1]
        sailing_distance    = round(geopy.distance.geodesic(position_first,position_next).nautical,3)    #In Nautical Miles
        vessel_heading      = calc_vessel_heading_2(position_first, position_next)                       #In Degrees (North is 90)
        WSE,WSN             = getweather(route_sailing_time,position_first[0],position_first[1])         #Wind speed East and North
        TWS                 = True_wind_speed(WSN,WSE)                                                   #True Windspeed (pythagoras)
        TWD                 = True_wind_direction(vessel_heading,WSN,WSE)                                #True Wind Direction
        AWS                 = Apparent_Wind_Speed(TWS,vessel_speed,TWD)                                  #Apparent wind speed
        AWA                 = alpha(vessel_speed,vessel_heading,WSN,WSE )                                #Apparent wind angle
        Forward_Force,Perp_Force = Force_produced(AWS, AWA)                           #Forward and Perpendicular force from Flettners
        vessel_speed    = round(Speed_achieved_old(Perp_Force, Forward_Force),3)    #Sailing Speed obtained in KNOTS


        if type(Forward_Force) == MaskedConstant or type(Perp_Force) == MaskedConstant:
            print("ouchie, we have a mask", i)
            return 1

        #Speed and time calculations:

        if vessel_speed == 0:
            sailing_time = 1.00                                           #If ship experiences no wind, it waits an hour before attempting to sail with new wind
            print(f"Slowsteam :( at location {position_first}, at time {route_sailing_time}")
        else:
            sailing_time     = round(sailing_distance/vessel_speed,3)           #time used to sail trip added

        # check for extreme time usage
        if vessel_speed < 1:
            poor_sailing_time += sailing_time
            poor_sailing_distance += sailing_distance

        #Append data to vectors
        sailing_speed_vector.append(vessel_speed)
        coordinate_sailing_time.append(sailing_time)
        true_wind_speed_vector.append(TWS)
        true_wind_direction_vector.append(TWD)
        apparent_wind_speed_observed.append(AWS)
        tot_sailing_dist += sailing_distance

        route_sailing_time += sailing_time                          #Sailing time of total route
        if len(coordinate_sailing_time) != 0:
            total_time_sailed_route = sum(coordinate_sailing_time)

        # Adding timestamp to first column
        if i > 0:
            time_beginning  = datestamp_vector[i-1]
            time_new        = add_hours_to_date(time_beginning,sailing_time)
            datestamp_vector.append(time_new)

    #Create vector of lists

    #tot_sailing_dist                = np.array(tot_sailing_dist)
    #sailing_speed_vector            = np.array(sailing_speed_vector)
    #true_wind_speed_vector          = np.array(true_wind_speed_vector)
    #true_wind_direction_vector      = np.array(true_wind_direction_vector)
    #route_sailing_time              = np.array(route_sailing_time)
    #coordinate_sailing_time         = np.array(coordinate_sailing_time)

    return total_time_sailed_route,tot_sailing_dist, poor_sailing_time, poor_sailing_distance, sailing_speed_vector,\
           true_wind_speed_vector, true_wind_direction_vector, route_sailing_time, coordinate_sailing_time, datestamp_vector




#Wind_North_vector_func,Wind_East_vector_func, Wind_tot_vector_func = Find_WSV_over_trip_at_time_TID(position_array,0)
#time_of_trip,total_sailing_distance = main(position_array)
#print(f"Time of trip is {time_of_trip},\n"
#      f" giving a total trip distance is {total_sailing_distance}\n"
#      f"with an average sailing speed of {total_sailing_distance/time_of_trip} knots")





#main = time_of_trip,tot_sailing_dist, poor_sailing_time, poor_sailing_distance
def _return(csv,routenumber,interval):
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
                        9 runs Newcastle Aberdeen
                        10 runs Ålesund Danmark
    :param interval: What intervall to write files
    :return:    time_of_trip
                tot_sailing_dist
                VS__return_vector
    """

    starttime = datetime(2020,7,1,00,00,00)
    iteration                       = 0
    hour_intervall                  = interval                                                 #at what hourly interval should we simulate?
    route_travel                    = read_route(csv)
    time_of__return              = 17520                                             #two years in hours
    time_of_trip                    = np.zeros(int(time_of__return/hour_intervall))
    tot_sailing_dist                = np.zeros(int(time_of__return/hour_intervall))
    poor_sailing_time               = np.zeros(int(time_of__return/hour_intervall))
    poor_sailing_distance           = np.zeros(int(time_of__return/hour_intervall))
    #VS__return_vector            = np.zeros(int(time_of__return/hour_intervall), dtype=object)
    datestamp__return_vector     = []
    VS__return_vector            = []
    TWS__return_vector           = []
    TWD__return_vector           = []
    date_of__return = starttime
    poor_sailing_speed = 0


    for iteration in range(0,int(time_of__return/hour_intervall)) : #repeating _return for each hour_intervall through a year
    #for iteration in range(8400,8600) : #repeating _return for each hour_intervall through a year

        time_of_trip_1,tot_sailing_dist_1, poor_sailing_time_1, poor_sailing_distance_1, sailing_speed_vector, TWS_vector, TWD_vector, route_sailing_time, coordinate_sailing_time, datestamp = main(route_travel,iteration, date_of__return)
        time_of_trip[int(iteration)]            = time_of_trip_1
        tot_sailing_dist[int(iteration)]        = tot_sailing_dist_1
        poor_sailing_time[int(iteration)]       = poor_sailing_time_1
        poor_sailing_distance[int(iteration)]   = poor_sailing_distance_1
        VS__return_vector.extend(sailing_speed_vector)
        TWS__return_vector.extend(TWS_vector)
        TWD__return_vector.extend(TWD_vector)
        datestamp__return_vector.extend(datestamp)

#Komenterer ut dette når jeg simultankjører
        #if iteration%1000 == 0 and iteration != 0 or iteration == 1:
        #    poor_sailing_speed = 0
        #    if poor_sailing_time_1 > 0:
        #        poor_sailing_speed = poor_sailing_distance_1 / poor_sailing_time_1
        #    print(f"[1] speed sailing the distance of {tot_sailing_dist[iteration]} nm\n"
        #          f" is {np.average(sailing_speed_vector[0:iteration])} knots")
        #    print(f"[2] total iteration sailed at less than 1 knot is {poor_sailing_time[iteration]}\n"
        #          f"[3] this iteration is used to sail {poor_sailing_distance[iteration]} nm\n"
        #          f"[4] at an average speed of {poor_sailing_speed} knots")
        #    #print(f"[5] The wind at this point in time was measured to be {TWS__return_vector[iteration]} with an AWA of {TWD__return_vector[iteration]}")
        #    print(f"[6] {datetime.now()},iteration is {iteration}")

        if iteration%100 == 0:
            print("progress is made, iteration:", iteration, "on route", routenumber)
        date_of__return = add_hours_to_date(date_of__return,hour_intervall)



    print(f"Average speed sailing {csv} over {iteration} iterations is {np.average(VS__return_vector[0])}")
    print(f"Throughout all iterations, the vessel sails less than one knot for an average of {np.average(poor_sailing_time)} hours\n"
          f"this iteration is used to sail on average {np.average(poor_sailing_distance)} nautical miles\n"
          f"at an average speed of {poor_sailing_speed} knots")


    #Read files:
    #Trond_Ålesund return
    file_speed_Trond_Aalesund   = mac_windows_file_handle("Output_files/Trondheim_Ålesund_return/savespeed_TrondAales_return.csv")
    file_TWS_Trond_Aalesund   = mac_windows_file_handle("Output_files/Trondheim_Ålesund_return/saveTWS_TrondAales_return.csv")
    file_TWD_Trond_Aalesund     = mac_windows_file_handle("Output_files/Trondheim_Ålesund_return/saveTWD_TrondAales_return.csv")
    datestamp_file_1            = mac_windows_file_handle("Output_files/Datestamps_return/datestamp1_return.csv")


    #Ålesund_floro return
    file_speed_Aalesund_Floro   = mac_windows_file_handle("Output_files/Ålesund_Florø_return/savespeed_AalesFloro_return.csv")
    file_TWS_Aalesund_Floro     = mac_windows_file_handle("Output_files/Ålesund_Florø_return/saveTWS_AalesFloro_return.csv")
    file_TWD_Aalesund_Floro     = mac_windows_file_handle("Output_files/Ålesund_Florø_return/saveTWD_AalesFloro_return.csv")
    datestamp_file_2            = mac_windows_file_handle("Output_files/Datestamps_return/datestamp2_return.csv")


    #Floro_bergen return
    file_speed_Floro_Bergen     = mac_windows_file_handle("Output_files/Florø_Bergen_return/savespeed_FloroBergen_return.csv")
    file_TWS_Floro_Bergen       = mac_windows_file_handle("Output_files/Florø_Bergen_return/saveTWS_FloroBergen_return.csv")
    file_TWD_Floro_Bergen       = mac_windows_file_handle("Output_files/Florø_Bergen_return/saveTWD_FloroBergen_return.csv")
    datestamp_file_3            = mac_windows_file_handle("Output_files/Datestamps_return/datestamp3_return.csv")


    #Bergen_stavanger return
    file_speed_Bergen_Stavanger = mac_windows_file_handle("Output_files/Bergen_Stavanger_return/savespeed_BrgStvg_return.csv")
    file_TWS_Bergen_Stavanger   = mac_windows_file_handle("Output_files/Bergen_Stavanger_return/saveTWS_BergenStavanger_return.csv")
    file_TWD_Bergen_Stavanger   = mac_windows_file_handle("Output_files/Bergen_Stavanger_return/saveTWD_BergenStavanger_return.csv")
    datestamp_file_4            = mac_windows_file_handle("Output_files/Datestamps_return/datestamp4_return.csv")


    #Aberdeen_Færøyene return
    file_speed_Aber_Faer    = mac_windows_file_handle("Output_files/Aberdeen_Færøyene_return/savespeed_Aberdeen_Færøyene_return.csv")
    file_TWS_Aber_Faer      = mac_windows_file_handle("Output_files/Aberdeen_Færøyene_return/saveTWD_Aberdeen_Færøyene_return.csv")
    file_TWD_Aber_Faer      = mac_windows_file_handle("Output_files/Aberdeen_Færøyene_return/saveTWS_Aberdeen_Færøyene_return.csv")
    datestamp_file_5        = mac_windows_file_handle("Output_files/Datestamps_return/datestamp5_return.csv")


    #Amsterdam Newcastle return
    file_speed_Amst_New = mac_windows_file_handle("Output_files/Amsterdam_Newcastle_return/savespeed_Amsterdam_Newcastle_return.csv")
    file_TWS_Amst_New   = mac_windows_file_handle("Output_files/Amsterdam_Newcastle_return/saveTWD_Amsterdam_Newcastle_return.csv")
    file_TWD_Amst_New    = mac_windows_file_handle("Output_files/Amsterdam_Newcastle_return/saveTWS_Amsterdam_Newcastle_return.csv")
    datestamp_file_6    = mac_windows_file_handle("Output_files/Datestamps_return/datestamp6_return.csv")


    #Danmark Amsterdam return 
    file_speed_Dk_Amst  = mac_windows_file_handle("Output_files/Danmark_Amsterdam_return/savespeed_Danmark_Amsterdam_return.csv")
    file_TWS_Dk_Amst    = mac_windows_file_handle("Output_files/Danmark_Amsterdam_return/saveTWD_Danmark_Amsterdam_return.csv")
    file_TWD_Dk_Amst    = mac_windows_file_handle("Output_files/Danmark_Amsterdam_return/saveTWS_Danmark_Amsterdam_return.csv")
    datestamp_file_7    = mac_windows_file_handle("Output_files/Datestamps_return/datestamp7_return.csv")


    #Færøyene Ålesund return

    file_speed_Faer_Aal = mac_windows_file_handle("Output_files/Færøyene_Ålesund_return/savespeed_Færøyene_Ålesund_return.csv")
    file_TWS_Faer_Aal   = mac_windows_file_handle("Output_files/Færøyene_Ålesund_return/saveTWD_Færøyene_Ålesund_return.csv")
    file_TWD_Faer_Aal   = mac_windows_file_handle("Output_files/Færøyene_Ålesund_return/saveTWS_Færøyene_Ålesund_return.csv")
    datestamp_file_8    = mac_windows_file_handle("Output_files/Datestamps_return/datestamp8_return.csv")


    #Newcastle Aberdeen return

    file_speed_New_Aber = mac_windows_file_handle("Output_files/Newcastle_Aberdeen_return/SS_Newcastle_Aberdeen_return.csv")
    file_TWS_New_Aber   = mac_windows_file_handle("Output_files/Newcastle_Aberdeen_return/TWD_Newcastle_Aberdeen_return.csv")
    file_TWD_New_Aber   = mac_windows_file_handle("Output_files/Newcastle_Aberdeen_return/TWS_Newcastle_Aberdeen_return.csv")
    datestamp_file_9    = mac_windows_file_handle("Output_files/Datestamps_return/datestamp9_return.csv")


    #Ålesund Danmark return 

    file_speed_Aal_Dk   = mac_windows_file_handle("Output_files/Ålesund_Danmark_return/savespeed_Ålesund_Danmark_return.csv")
    file_TWS_Aal_Dk     = mac_windows_file_handle("Output_files/Ålesund_Danmark_return/saveTWD_Ålesund_Danmark_return.csv")
    file_TWD_Aal_Dk     = mac_windows_file_handle("Output_files/Ålesund_Danmark_return/saveTWS_Ålesund_Danmark_return.csv")
    datestamp_file_10   = mac_windows_file_handle("Output_files/Datestamps_return/datestamp10_return.csv")

    #Floro port return

    file_speed_Floro_port   = mac_windows_file_handle("Output_files/Floro_port/savespeed_port.csv")
    file_TWS_Floro_port     = mac_windows_file_handle("Output_files/Floro_port/TWD_Floro_port.csv")
    file_TWD_Floro_port     = mac_windows_file_handle("Output_files/Floro_port/TWS_Floro_port.csv")
    datestamp_file_11       = mac_windows_file_handle("Output_files/Datestamps/datestamp11.csv")


    #Write to files

    if routenumber == 1:

        write_to_file(VS__return_vector, file_speed_Trond_Aalesund)  #VS vector file
        write_to_file(TWS__return_vector,file_TWS_Trond_Aalesund)    #TWS vector file
        write_to_file(TWD__return_vector,file_TWD_Trond_Aalesund)    #TWD vector file
        write_to_file(datestamp__return_vector, datestamp_file_1)                 #Datestamp file
        add_timestamp_to_dataframe(file_speed_Trond_Aalesund, datestamp_file_1)
        add_timestamp_to_dataframe(file_TWS_Trond_Aalesund, datestamp_file_1)
        add_timestamp_to_dataframe(file_TWD_Trond_Aalesund, datestamp_file_1)

    elif routenumber == 2:
        write_to_file(VS__return_vector, file_speed_Aalesund_Floro)  #VS vector file
        write_to_file(TWS__return_vector,file_TWS_Aalesund_Floro)    #TWS vector file
        write_to_file(TWD__return_vector,file_TWD_Aalesund_Floro)    #TWD vector file
        write_to_file(datestamp__return_vector, datestamp_file_2)                 #Datestamp file
        add_timestamp_to_dataframe(file_speed_Aalesund_Floro, datestamp_file_2)
        add_timestamp_to_dataframe(file_TWS_Aalesund_Floro, datestamp_file_2)
        add_timestamp_to_dataframe(file_TWD_Aalesund_Floro, datestamp_file_2)

    elif routenumber == 3:
        write_to_file(VS__return_vector, file_speed_Floro_Bergen)    #VS vector file
        write_to_file(TWS__return_vector,file_TWS_Floro_Bergen)      #TWS vector file
        write_to_file(TWD__return_vector,file_TWD_Floro_Bergen)      #TWD vector file
        write_to_file(datestamp__return_vector, datestamp_file_3)                 #Datestamp file
        add_timestamp_to_dataframe(file_speed_Floro_Bergen, datestamp_file_3)
        add_timestamp_to_dataframe(file_TWS_Floro_Bergen, datestamp_file_3)
        add_timestamp_to_dataframe(file_TWD_Floro_Bergen, datestamp_file_3)

    elif routenumber == 4:
        write_to_file(VS__return_vector,file_speed_Bergen_Stavanger) #VS vector file
        write_to_file(TWS__return_vector,file_TWS_Bergen_Stavanger)  #TWS vector file
        write_to_file(TWD__return_vector,file_TWD_Bergen_Stavanger)  #TWD vector file
        write_to_file(datestamp__return_vector, datestamp_file_4)       #Datestamp file
        add_timestamp_to_dataframe(file_speed_Bergen_Stavanger, datestamp_file_4)
        add_timestamp_to_dataframe(file_TWS_Bergen_Stavanger, datestamp_file_4)
        add_timestamp_to_dataframe(file_TWD_Bergen_Stavanger, datestamp_file_4)

    elif routenumber == 5:
        write_to_file(VS__return_vector,file_speed_Aber_Faer) #VS vector file
        write_to_file(TWS__return_vector,file_TWS_Aber_Faer)  #TWS vector file
        write_to_file(TWD__return_vector,file_TWD_Aber_Faer)  #TWD vector file
        write_to_file(datestamp__return_vector, datestamp_file_5)          #Datestamp file
        add_timestamp_to_dataframe(file_speed_Aber_Faer, datestamp_file_5)
        add_timestamp_to_dataframe(file_TWS_Aber_Faer, datestamp_file_5)
        add_timestamp_to_dataframe(file_TWD_Aber_Faer, datestamp_file_5)

    elif routenumber == 6:
        write_to_file(VS__return_vector,file_speed_Amst_New) #VS vector file
        write_to_file(TWS__return_vector,file_TWS_Amst_New)  #TWS vector file
        write_to_file(TWD__return_vector,file_TWD_Amst_New)  #TWD vector file
        write_to_file(datestamp__return_vector, datestamp_file_6)         #Datestamp file
        add_timestamp_to_dataframe(file_speed_Amst_New, datestamp_file_6)
        add_timestamp_to_dataframe(file_TWS_Amst_New, datestamp_file_6)
        add_timestamp_to_dataframe(file_TWD_Amst_New, datestamp_file_6)

    elif routenumber == 7:
        write_to_file(VS__return_vector,file_speed_Dk_Amst) #VS vector file
        write_to_file(TWS__return_vector,file_TWS_Dk_Amst)  #TWS vector file
        write_to_file(TWD__return_vector,file_TWD_Dk_Amst)  #TWD vector file
        write_to_file(datestamp__return_vector, datestamp_file_7)        #Datestamp file
        add_timestamp_to_dataframe(file_speed_Dk_Amst, datestamp_file_7)
        add_timestamp_to_dataframe(file_TWS_Dk_Amst, datestamp_file_7)
        add_timestamp_to_dataframe(file_TWD_Dk_Amst, datestamp_file_7)

    elif routenumber == 8:
        write_to_file(VS__return_vector,file_speed_Faer_Aal) #VS vector file
        write_to_file(TWS__return_vector,file_TWS_Faer_Aal)  #TWS vector file
        write_to_file(TWD__return_vector,file_TWD_Faer_Aal)  #TWD vector file
        write_to_file(datestamp__return_vector, datestamp_file_8)         #Datestamp file
        add_timestamp_to_dataframe(file_speed_Faer_Aal,datestamp_file_8)
        add_timestamp_to_dataframe(file_TWS_Faer_Aal,datestamp_file_8)
        add_timestamp_to_dataframe(file_TWD_Faer_Aal,datestamp_file_8)

    elif routenumber == 9:
        write_to_file(VS__return_vector,file_speed_New_Aber) #VS vector file
        write_to_file(TWS__return_vector,file_TWS_New_Aber)  #TWS vector file
        write_to_file(TWD__return_vector,file_TWD_New_Aber)  #TWD vector file
        write_to_file(datestamp__return_vector, datestamp_file_9)         #Datestamp file
        add_timestamp_to_dataframe(file_speed_New_Aber, datestamp_file_9)
        add_timestamp_to_dataframe(file_TWS_New_Aber, datestamp_file_9)
        add_timestamp_to_dataframe(file_TWD_New_Aber, datestamp_file_9)

    elif routenumber == 10:
        write_to_file(VS__return_vector,file_speed_Aal_Dk) #VS vector file
        write_to_file(TWS__return_vector,file_TWS_Aal_Dk)  #TWS vector file
        write_to_file(TWD__return_vector,file_TWD_Aal_Dk)  #TWD vector file
        write_to_file(datestamp__return_vector, datestamp_file_10)       #Datestamp file
        add_timestamp_to_dataframe(file_speed_Aal_Dk, datestamp_file_10)
        add_timestamp_to_dataframe(file_TWS_Aal_Dk, datestamp_file_10)
        add_timestamp_to_dataframe(file_TWD_Aal_Dk, datestamp_file_10)

    elif routenumber == 11:
        write_to_file(VS__return_vector,file_speed_Floro_port) #VS vector file
        write_to_file(TWS__return_vector,file_TWS_Floro_port)  #TWS vector file
        write_to_file(TWD__return_vector,file_TWD_Floro_port)  #TWD vector file
        write_to_file(datestamp__return_vector, datestamp_file_11)          #Datestamp file
        add_timestamp_to_dataframe(file_speed_Floro_port, datestamp_file_11)
        add_timestamp_to_dataframe(file_TWS_Floro_port, datestamp_file_11)
        add_timestamp_to_dataframe(file_TWD_Floro_port, datestamp_file_11)

    return 0



#change to _return(route) so that you dont need to hardcode


def run_return(route, interval):
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
        :param interval: At what time interval do we run _return
        :return: saved files with _return results
        """

    Trond_aalesund      = mac_windows_file_handle("Route_data/Trondheim_Ålesund_Rute/Route_Trond_Aales.csv")
    Aalesund_Floro      = mac_windows_file_handle("Route_data/Ålsund_Florø_Rute/Route_Ålesund_Florø.csv")
    Floro_Bergen        = mac_windows_file_handle("Route_data/Florø_Bergen_Rute/Route_Floro_Bergen.csv")
    Bergen_Stavanger    = mac_windows_file_handle("Route_data/Bergen_Stavanger_Rute/Route_Bergen_Stavanger.csv")
    Aberdeen_Faer       = mac_windows_file_handle("Route_data/Aberdeen_Færøyene_Route/Route_Aberdeen_Færøyene.csv")
    Amst_New            = mac_windows_file_handle("Route_data/Amsterdam_Newcastle_Route/Route_Amsterdam_Newcastle.csv")
    DK_Amst             = mac_windows_file_handle("Route_data/Danmark_Amsterdam_Route/Route_Danmark_Amsterdam.csv")
    Faer_Aale           = mac_windows_file_handle("Route_data/Færøyene_Ålesund_Route/Route_Færøyene_Ålesund.csv")
    New_Aber            = mac_windows_file_handle("Route_data/Newcastle_Aberdeen_Route/Route_Newcastle_Aberdeen.csv")
    Aale_DK             = mac_windows_file_handle("Route_data/Ålesund_Danmark_Route/Route_Ålesund_DK.csv")
    Floro_port          = mac_windows_file_handle("Route_data/Floro_port/Floro_port.csv")


    if route == 1:
        print("Running _return for route Trondheim Aalesund now")
        _return(Trond_aalesund,route,interval)
        print("Simulation Trondheim to Ålesund is now complete.\n")

    if route == 2:
        print("Running _return for route Ålesund Florø now")
        _return(Aalesund_Floro,route,interval)
        print("Simulation Ålesund to Florø is now complete.\n")


    if route == 3:
        print("Running _return for route Florø Bergen now")
        _return(Floro_Bergen,route,interval)
        print("Simulation Florø to Bergen d is now complete.\n")

    if route == 4:
        print("Running _return for route Bergen Stavanger now")
        _return(Bergen_Stavanger,route,interval)
        print("Simulation Bergen to Stavanger is now complete.\n")


    if route == 5:
        print("Running _return for route Aberdeen_Faerøyene now")
        _return(Aberdeen_Faer, route,interval)
        print("Simulation Aberdeen to Færøyene is now complete.\n")

    if route == 6:
        print("Running _return for route Amsterdam_Newcastle now")
        _return(Amst_New, route,interval)
        print("Simulation Amsterdam to Newcastle is now complete.\n")

    if route == 7:
        print("Running _return for route Danmark_Amsterdam now")
        _return(DK_Amst, route,interval)
        print("Simulation Danmark to Amsterdam is now complete.\n")

    if route == 8:
        print("Running _return for route Færøyene_Ålesund now")
        _return(Faer_Aale, route,interval)
        print("Simulation Færøyene to Ålesund is now complete.\n")

    if route == 9:
        print("Running _return for route Newcastle_Aberdeen now")
        _return(New_Aber, route,interval)
        print("Simulation Newcastle to Aberdeen is now complete.\n")

    if route == 10:
        print("Running _return for route Ålesund_Danmark now")
        _return(Aale_DK, route,interval)
        print("Simulation Ålesund to Danmark is now complete.\n")

    if route == 11:
        print("Running _return for route in florø port now")
        _return(Floro_port, route, interval)
        print("Simulation for route in florø port is now complete")


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





#reset_index()
print("Everything is working 1241, 16/03")
steps = 1000
run_return(1,steps)
run_return(2,steps)
run_return(3,steps)
run_return(4,steps)
run_return(5,steps)
run_return(6,steps)
run_return(7,steps)
run_return(8,steps)
run_return(9,steps)
run_return(10,steps)
#run_return(11,steps)


import threading

def run_return_multiple(a, b):
    run_return(a,b)
    print(f"Simulation of route{a} every {b} iteration is completed")
def run_multiple():
    j = 100
    inputs = [(1,j),(2,j), (3,j), (4,j), (5,j), (6,j), (7,j), (8,j), (9,j),(10,j), (11,j)]

    threads = []

    for input in inputs:
        thread = threading.Thread(target=run_return_multiple, args=input)
        thread.start()
        threads.append(thread)

    for thread in threads:
        thread.join()
    return 0
#run_multiple()





print("Finished <3<3")
print("Bonus print: Mathias er digg, kl. 1112, den 20032023")


