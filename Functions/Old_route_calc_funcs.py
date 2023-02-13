
import numpy as np
from numpy.ma.core import MaskedConstant
from datetime import datetime
from plot_functions import plot_power
from file_handling import write_to_file, read_array_from_file, read_position_vect_from_file
from Weather_Handling import getweather, Apparent_wind_angle, Apparent_Wind_Speed
#from main import Speed_sailed

filename = "../env/position_array"  # file containing position array of route
position_array = read_position_vect_from_file(filename) # reads said file

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
    filename_avg_forward_force  = f"../env/output_files5/avg_forward_force.txt"
    filename_perp_force         = f"../env/output_files5/avg_perp_force.txt"
    filename_speed_sailed       = f"../env/output_files5/speed_sailed_over_time.txt"
    write_to_file(forw_force_array_time, filename_avg_forward_force)
    write_to_file(perp_force_array_time, filename_perp_force)     #winddir 31.12
    write_to_file(speed_sailed_array_over_time,filename_speed_sailed)
    return forw_force_array_time, perp_force_array_time, speed_sailed_array_over_time

#Finds windspeeds at position in position array and time in tid in weater data
#function that gives true windspeed in m/s for wind north,east
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


def Apparent_wind_and_direction_arrays_over_trip(position_array_func, tid):#, vessel_velocity):
    Wind_North_vector_func, Wind_East_vector_func, Wind_tot_vector_func = Find_WSV_over_trip_at_time_TID(position_array_func, tid)          # windspeeds at (lat,lon) for t € 1050930 through 1051824
    AWD_Trip  = Apparent_wind_angle(Wind_North_vector_func, Wind_East_vector_func)  #true winddirection at (lat,lon) for t € 1050930 through 1051824
    AWS_Trip  = Apparent_Wind_Speed(Wind_tot_vector_func, AWD_Trip)#, vessel_velocity)          #true windspeed at (lat,lon) for t € 1050930 through 1051824
    return AWS_Trip,AWD_Trip

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