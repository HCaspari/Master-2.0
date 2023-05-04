
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


def extract_data():
    # Load the first NetCDF file
    ds1 = xr.open_dataset('../Weather Check/2020 Troll/NO_TS_MO_Troll-A_202001.nc')
    ds2 = xr.open_dataset('../Weather Check/2020 Troll/NO_TS_MO_Troll-A_202002.nc')
    ds3 = xr.open_dataset('../Weather Check/2020 Troll/NO_TS_MO_Troll-A_202003.nc')
    ds4 = xr.open_dataset('../Weather Check/2020 Troll/NO_TS_MO_Troll-A_202004.nc')
    ds5 = xr.open_dataset('../Weather Check/2020 Troll/NO_TS_MO_Troll-A_202005.nc')
    ds6 = xr.open_dataset('../Weather Check/2020 Troll/NO_TS_MO_Troll-A_202006.nc')
    ds7 = xr.open_dataset('../Weather Check/2020 Troll/NO_TS_MO_Troll-A_202007.nc')
    ds8 = xr.open_dataset('../Weather Check/2020 Troll/NO_TS_MO_Troll-A_202008.nc')
    ds9 = xr.open_dataset('../Weather Check/2020 Troll/NO_TS_MO_Troll-A_202009.nc')
    ds10 = xr.open_dataset('../Weather Check/2020 Troll/NO_TS_MO_Troll-A_202010.nc')
    ds11 = xr.open_dataset('../Weather Check/2020 Troll/NO_TS_MO_Troll-A_202011.nc')
    ds12 = xr.open_dataset('../Weather Check/2020 Troll/NO_TS_MO_Troll-A_202012.nc')

    variable_names = list(ds1.variables.keys())
    #northward_wind_1  = dataset_NW_1.variables["northward_wind"]
    # Print the names of all variables in the NetCDF file
    print(variable_names)


    WSPD_1 = ds1.variables["WSPD"][:,0].values
    WSPD_2 = ds2.variables["WSPD"][:,0].values
    WSPD_3 = ds3.variables["WSPD"][:,0].values
    WSPD_4 = ds4.variables["WSPD"][:,0].values
    WSPD_5 = ds5.variables["WSPD"][:,0].values
    WSPD_6 = ds6.variables["WSPD"][:,0].values
    WSPD_7 = ds7.variables["WSPD"][:,0].values
    WSPD_8 = ds8.variables["WSPD"][:,0].values
    WSPD_9 = ds9.variables["WSPD"][:,0].values
    WSPD_10 = ds10.variables["WSPD"][:,0].values
    WSPD_11 = ds11.variables["WSPD"][:,0].values
    WSPD_12 = ds12.variables["WSPD"][:,0].values

    combined = np.concatenate((WSPD_1,WSPD_2,WSPD_3,WSPD_4,WSPD_5,WSPD_6,WSPD_7,WSPD_8,WSPD_9,WSPD_10,WSPD_11,WSPD_12),0)

    # Load the datasets into a list
    datasets = [xr.open_dataset(fp) for fp in file_paths]

    # Extract the WSPD variable data for each dataset using a loop
    WSPD_list = []
    for ds in datasets:
        WSPD = ds['WSPD'][:, 0].values
        WSPD_list.append(WSPD)
    return 0

def søppel():
    # Read the data from the two CSV files into Pandas dataframes
    df1 = pd.read_csv('Output_files/Floro_port/FloroPort_windspeed_winddirection.csv', sep=';', decimal=',')
    df2 = pd.read_csv('Output_files/Floro_port/TWD_Floro_port.csv')

    # Convert the datetime column to a datetime object
    df1['DateTime'] = pd.to_datetime(df1['DateTime'])
    df2['DateTime'] = pd.to_datetime(df2['DateTime'])

    # Define the datetime period to plot (e.g., from January 1, 2022 to March 31, 2022)
    start_date = pd.to_datetime('2021-10-01')
    end_date = pd.to_datetime('2021-11-01')

    # Filter the dataframes to include only data within the datetime period
    df1 = df1.loc[(df1['DateTime'] >= start_date) & (df1['DateTime'] <= end_date)]
    df2 = df2.loc[(df2['DateTime'] >= start_date) & (df2['DateTime'] <= end_date)]

    # Plot the wind speed data from file 1 in blue
    plt.plot(df1['DateTime'], df1['Wind direction'], color='blue', label='File 1')

    # Plot the wind speed data from file 2 in red
    plt.plot(df2['DateTime'], df2['Wind direction'], color='red', label='File 2')

    # Set the axis labels and legend
    plt.xlabel('Datetime')
    plt.ylabel('Wind direction')
    plt.legend()

    # Show the plot
    plt.show()
    return 0

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

    power_consumption = Cp * (rho_air/2) * A_rotor * AWS**3       #Input to flettner is energy to spin rotors (traut et al)


    if 0 <= AWD <= 90 or 270 <= AWD <= 360:  # drag is set to negative if the wind is coming ahead, and positive if not
        drag *= -1
    Force_4_flettners = (lift + drag) * rotor_amount / 1000  # kN                             #KiloNewton
    # The force generated by the flettner is the total output force minus the force needed to spin the rotor.

    perp_force = Force_4_flettners * abs(np.cos(AWD))  # perpendicular force flettner in KN
    if type(perp_force) != MaskedConstant:
        perp_force = round(min(perp_force, 1400),
                           2)  # Force is maximum 1400 kN, so as not to break flettners (cite: Norsepower technical)
    forward_force = Force_4_flettners * abs(np.sin(AWD))  # directional force flettner in KN

    if type(forward_force) != MaskedConstant:
        forward_force = round(min(forward_force, 1400), 2)  # Force is maximum 1400 kN, so as not to break flettners



    return forward_force, perp_force, power_consumption

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