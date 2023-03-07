import pandas as pd
import numpy as np
from route_handling import mac_windows_file_handle



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
alpha = 3.5


#write data found to file called filename_func
def write_to_file(data,filename_func):
    """
    :param data: data that is to be written to a file
    :param filename_func:  name of file that data is to be written to
    :return: 0
    """
    pd.DataFrame(data).to_csv(filename_func)
    return 0

#reads columns of dataframe to find distances, latitudes, longditudes of trip
def read_cols(filename_ais_data):
    """
    this function is specific to reading ais data
    :param filename_ais_data: filename of data, containing distance, traveltime, latitudes, longditudes,  heading of vessel
    :return: vectors showing distance and travel time between lats and longs, the latitudes and longditudes of those points and the heading that the vesesl has at these points.
    """
    df= pd.read_csv(filename_ais_data)
    #a = df[['dt',"distance","travel_time"]]
    distances  = df[["distance"]]
    travel_time = df[["travel_time"]]
    latitudes   = df[["lat"]]
    longditudes = df[["lon"]]
    heading     = df[["heading"]]
    dist_vect_func = []
    travel_time_vect_func   = []
    latitudes_vect_func     = []
    longditudes_vect_func   = []
    heading_vect_func       = []
    for i in range (len(distances)):
        dist_vect_func.append(distances.loc[i].iat[0])
        travel_time_vect_func.append(travel_time.loc[i].iat[0])
        latitudes_vect_func.append(round(latitudes.loc[i].iat[0],4))
        longditudes_vect_func.append(round(longditudes.loc[i].iat[0],4))
        heading_vect_func.append(heading.loc[i].iat[0])
    heading_file = mac_windows_file_handle("env/input_files/heading.csv")
    write_to_file(heading, heading_file)
    return dist_vect_func,travel_time_vect_func,latitudes_vect_func,longditudes_vect_func,heading_vect_func

def read_position_vect_from_file(filename_func):
    """
    :param filename_func: reads position vessel has in ais data
    :return: vector showing position data of vessel
    """
    readdata = pd.read_csv(filename_func)#, usecols=["1","0"])
    readdata = readdata[readdata != 0]

    position_array =  []
    for i in range(len(readdata)):
        position_array.append((round(readdata.loc[i].iat[0], 3),round(readdata.loc[i].iat[1], 3)))
    return position_array

def find_vals():
    """
    cleans ais data for specific values
    :return: 0
    """
    df = pd.read_csv(mac_windows_file_handle("env/input_files/ais_data_v3.csv"))
    del df["Unnamed: 0"]
    df = df[df.nav_status !=2] #removes every entry with nav_Status = 2
    df = df[df.nav_status !=5] #removes every entry with nav_status = 5
    df.to_csv(mac_windows_file_handle("ais_data_v4.csv"))
    return 0

def read_array_from_file(filename_func):
    """
    function to read data from file (save time after running program through)
    :param filename_func: name of file we want to read data from
    :return: array of that data
    """
    readdata = pd.read_csv(filename_func)
    readdata = readdata[readdata != 0]
    data_array = []
    #print(f"im reading {filename_func} now")
    for i in range(len(readdata)):
        data_array.append(round(readdata.loc[i].iat[1],3))
    #print(data_array)
    return data_array

def read_route(csv):
    """
    reads route data from a csv file
    :param csv: file containing route data
    :return: vector of positions thorugh route
    """
    df = pd.read_csv(csv)
    latitudes   = df["latitude"]
    longditudes = df["longditude"]
    positions = []
    for i in range(len(latitudes)):
        positions.append((latitudes[i],longditudes[i]))
    positions = np.asarray(positions)
    return positions

filname = mac_windows_file_handle("Output_files/Bergen_Stavanger_reise")
ARRAY = read_array_from_file(filname)

def test_read_files(routenumber):

    file_speed_Trond_Aalesund = mac_windows_file_handle("Output_files/savespeed_TrondAales.csv")
    file_TWS_Trond_Aalesund = mac_windows_file_handle("Output_files/saveTWS_TrondAales.csv")
    file_TWD_Trond_Aalesund = mac_windows_file_handle("Output_files/saveTWD_TrondAales.csv")
    file_speed_Aalesund_Floro = mac_windows_file_handle("Output_files/savespeed_AalesFloro.csv")
    file_TWS_Aalesund_Floro = mac_windows_file_handle("Output_files/saveTWS_AalesFloro.csv")
    file_TWD_Aalesund_Floro = mac_windows_file_handle("Output_files/saveTWD_AalesFloro.csv")
    file_speed_Floro_Bergen = mac_windows_file_handle("Output_files/savespeed_FloroBergen.csv")
    file_TWS_Floro_Bergen = mac_windows_file_handle("Output_files/saveTWS_FloroBergen.csv")
    file_TWD_Floro_Bergen = mac_windows_file_handle("Output_files/saveTWD_FloroBergen.csv")
    file_speed_Bergen_Stavanger = mac_windows_file_handle("Output_files/savespeed_BrgStvg.csv")
    file_TWS_Bergen_Stavanger = mac_windows_file_handle("Output_files/saveTWS_BergenStavanger.csv")
    file_TWD_Bergen_Stavanger = mac_windows_file_handle("Output_files/saveTWD_BergenStavanger.csv")

    if routenumber == 1:
        a = read_array_from_file(file_speed_Trond_Aalesund)
        b = read_array_from_file(file_TWS_Trond_Aalesund)
        c = read_array_from_file(file_TWD_Trond_Aalesund)
        print("read_test complete (speed, TWS, TWD)", a,"\n", b,"\n",c)
    elif routenumber == 2:
        d = read_array_from_file(file_speed_Aalesund_Floro)
        e = read_array_from_file(file_TWS_Aalesund_Floro)
        f =read_array_from_file(file_TWD_Aalesund_Floro)
        print("read_test complete, (speed, TWS, TWD)", d,"\n", e,"\n",f)
    elif routenumber == 3:
        g = read_array_from_file(file_speed_Floro_Bergen)
        h = read_array_from_file(file_TWS_Floro_Bergen)
        i = read_array_from_file(file_TWD_Floro_Bergen)
        print("read_test complete, (speed, TWS, TWD)", g,"\n", h,"\n", i)
    elif routenumber == 4:
        j = read_array_from_file(file_speed_Bergen_Stavanger)
        k = read_array_from_file(file_TWS_Bergen_Stavanger)
        l = read_array_from_file(file_TWD_Bergen_Stavanger)
        print("read_test complete, (speed, TWS, TWD)", j,"\n", k,"\n", l)