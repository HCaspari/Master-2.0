
import numpy as np
from plot_functions import plot_power
from file_handling import write_to_file, read_array_from_file, read_position_vect_from_file
from Weather_Handling import getweather, Apparent_wind_angle, Apparent_Wind_Speed
from pathlib import PureWindowsPath, Path



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
    avg_speed_sailed        = PureWindowsPath(Path("../env/output_files5/speed_sailed_over_time.txt"))
    avg_force_data          = PureWindowsPath(Path("../env/output_files5/avg_forward_force.txt"))
    avg_perp_force_data     = PureWindowsPath(Path("../env/output_files5/avg_perp_force.txt"))

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


def create_routes():

    route_DK_Stav   = auto_create_Route(Danmark,Stavanger,Mid_Danmark_stvg)
    route_Amst_New  = auto_create_Route(Amsterdam,Newcastle)
    route_Dk_Amst   = auto_create_Route(Danmark,Amsterdam)
    route_New_Aber  = auto_create_Route(Newcastle,Aberdeen)
    route_Aber_Faer  = auto_create_Route(Aberdeen,Færøyene,Midpoint_Fa_Ab)
    route_Faer_Aales  = auto_create_Route(Færøyene,Aalesund)
    route_Aales_Dk   = auto_create_Route(Aalesund,Danmark, Midpoint_2_Dan_Aal, Midpoint_1_Dan_Aal)

    filename_DK_Stav    = mac_windows_file_handle("Route_data/Danmark_Stavanger_Route/Route_Danmark_Stavanger.csv")
    filename_Amst_New   = mac_windows_file_handle("Route_data/Amsterdam_Newcastle_Route/Route_Amsterdam_Newcastle.csv")
    filename_Dk_Amst    = mac_windows_file_handle("Route_data/Danmark_Amsterdam_Route/Route_Danmark_Amsterdam.csv")
    filename_Faer_Aal     = mac_windows_file_handle("Route_data/Færøyene_Ålesund_Route/Route_Færøyene_Ålesund.csv")
    filename_New_Aber   = mac_windows_file_handle("Route_data/Newcastle_Aberdeen_Route/Route_Newcastle_Aberdeen.csv")
    filename_Aal_Dk      = mac_windows_file_handle("Route_data/Ålesund_Danmark_Route/Route_Ålesund_Florø.csv")
    filename_Aber_Faer   = mac_windows_file_handle("Route_data/Aberdeen_Færøyene_Route/Route_Aberdeen_Færøyene.csv")

    write_to_file(route_Aber_Faer,filename_Aber_Faer)
    write_to_file(route_DK_Stav,filename_DK_Stav)
    write_to_file(route_Amst_New,filename_Amst_New)
    write_to_file(route_Dk_Amst,filename_Dk_Amst)
    write_to_file(route_New_Aber,filename_New_Aber)
    write_to_file(route_Faer_Aales,filename_Faer_Aal)
    write_to_file(route_Aales_Dk,filename_Aal_Dk)

    return route_Aber_Faer,route_Faer_Aales,route_Aales_Dk,route_Dk_Amst,route_DK_Stav,route_New_Aber,route_Amst_New



