import pandas as pd
import numpy as np
import datetime

import csv
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
    df = pd.DataFrame(data)
    pd.DataFrame(data).to_csv(filename_func)
    return 0


def write_to_file_2(datatype, data, filename_func):
    """
    :param data: data that is to be written to a file
    :param filename_func:  name of file that data is to be written to
    :return: 0
    """

    df = pd.DataFrame({datatype: data})
    df.to_csv(filename_func)
    return 0


#reads columns of dataframe to find distances, latitudes, longitudes of trip
def read_cols(filename_ais_data):
    """
    this function is specific to reading ais data
    :param filename_ais_data: filename of data, containing distance, traveltime, latitudes, longitudes,  heading of vessel
    :return: vectors showing distance and travel time between lats and longs, the latitudes and longitudes of those points and the heading that the vesesl has at these points.
    """
    df= pd.read_csv(filename_ais_data)
    #a = df[['dt',"distance","travel_time"]]
    distances  = df[["distance"]]
    travel_time = df[["travel_time"]]
    latitudes   = df[["lat"]]
    longitudes = df[["lon"]]
    heading     = df[["heading"]]
    dist_vect_func = []
    travel_time_vect_func   = []
    latitudes_vect_func     = []
    longitudes_vect_func   = []
    heading_vect_func       = []
    for i in range (len(distances)):
        dist_vect_func.append(distances.loc[i].iat[0])
        travel_time_vect_func.append(travel_time.loc[i].iat[0])
        latitudes_vect_func.append(round(latitudes.loc[i].iat[0],4))
        longitudes_vect_func.append(round(longitudes.loc[i].iat[0],4))
        heading_vect_func.append(heading.loc[i].iat[0])
    heading_file = mac_windows_file_handle("env/input_files/heading.csv")
    write_to_file(heading, heading_file)
    return dist_vect_func,travel_time_vect_func,latitudes_vect_func,longitudes_vect_func,heading_vect_func

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

def read_datestamp_from_file(datestamp_file):
    datestamp = pd.read_csv(datestamp_file)
    datestamp = datestamp[datestamp != 0]
    datestamp_Array = []
    for i in range(len(datestamp)):
        datestamp_Array.append(datestamp.loc[i].iat[1])

    return datestamp_Array

def add_timestamp_to_dataframe(datafile, timestamp_file):
    """
    Adds column at beginning of csv file with timestamp
    :param datafile: csv containing data
    :param timestamp_file: timestamp file with time stamps
    :return: 0
    """

    timestamp_array = read_datestamp_from_file(timestamp_file)

    df = pd.read_csv(datafile)
    #print(df)
    df = df.drop(df.columns[0],axis=1)              #Drop first column of nothing
    #print(df)
    len_datafile = len(df.iloc[:,0])
    len_timestamp = len(timestamp_array)
    while len_datafile > len_timestamp:
        df = df.drop(df.index[-1])
        len_datafile = len(df.iloc[:, 0])
    #print(df)
    timestamp_array = timestamp_array[:len_datafile]
    df.insert(0, "timestamp", timestamp_array)       #insert timestamp vector
    #print(df)
    df.to_csv(datafile)
    return 0

#datafile = mac_windows_file_handle("Output_files/Færøyene_Ålesund/savespeed_Færøyene_Ålesund.csv")
#datestamp_file = mac_windows_file_handle("Output_files/datestamp9.csv")
#add_timestamp_to_dataframe(datafile, datestamp_file)

def read_route(csv):
    """
    reads route data from a csv file
    :param csv: file containing route data
    :return: vector of positions thorugh route
    """
    df = pd.read_csv(csv)
    latitudes   = df["latitude"]
    longitudes = df["longitude"]
    positions = []
    for i in range(len(latitudes)):
        positions.append((latitudes[i],longitudes[i]))
    positions = np.asarray(positions)
    return positions

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

    return 0

def create_array_with_datetime():
    # define a numpy array
    data = np.array([[1, 2], [3, 4], [5, 6]])

    # create a datetime object for today's date and time
    today = datetime.datetime.now()

    # create a new array with an extra column for the datetime
    dt_data = np.zeros((data.shape[0], data.shape[1] + 1), dtype=object)

    # set the original data in the new array
    dt_data[:, :-1] = data

    # set the datetime in the last column of the new array
    dt_data[:, -1] = today

    # create a vector of data
    vector = np.array([7, 8, 9])

    # add the vector as a column to the right of the datetime array
    final_data = np.concatenate((dt_data, np.expand_dims(vector, axis=1)), axis=1)

    # print the final array
    print(final_data)

    return 0



def drop_duplicates():
    """
    # Read in the CSV file

    :return: 0
    """
    df0 = pd.read_csv(mac_windows_file_handle('Output_files/Floro_port/TWS_Floro_port.csv'))
    df1 = pd.read_csv(mac_windows_file_handle('Output_files/Floro_port/TWD_Floro_port.csv'))


    # Drop every third row
    df0.drop_duplicates(subset=[df0.columns[-1]], keep='last', inplace=True)
    df1.drop_duplicates(subset=[df1.columns[-1]], keep='last', inplace=True)

    # Write the updated dataframe back to a new CSV file
    df0.to_csv(mac_windows_file_handle('Output_files/Floro_port/TWS_Floro_port.csv'), index=False)
    df1.to_csv(mac_windows_file_handle('Output_files/Floro_port/TWD_Floro_port.csv'), index=False)


def reset_index():
    """
    resets index of csv to 1
    :return: 0
    """
    # Read in the CSV file
    df = pd.read_csv(mac_windows_file_handle('Output_files/Floro_port/TWS_Floro_port.csv'))
    df1 = pd.read_csv(mac_windows_file_handle('Output_files/Floro_port/TWS_Floro_port.csv'))

    # Reset the index and rename the first column to 'index'
    index = []
    for i in range(len(df)):
        index.append(i)

    df.insert(loc=0, column='new_column_name', value=index)
    df1.insert(loc=0, column='new_column_name', value=index)
    # Write the updated dataframe back to a new CSV file
    df.to_csv(mac_windows_file_handle('Output_files/Floro_port/TWS_Floro_port.csv'), index=False)
    df1.to_csv(mac_windows_file_handle('Output_files/Floro_port/TWS_Floro_port.csv'), index=False)

def flip_route(csv_file_old_route,csv_file_new_route):
    """

    :param csv_file_old_route: Filepath of old route file
    :param csv_file_new_route: Filepath to new, empty route file
    :return: flips old route file to create return route file.
    """
    # Set the file path of the CSV file
    file_path_old = csv_file_old_route
    file_path_new = csv_file_new_route
    dfnew = pd.DataFrame(columns=["latitude", "longitude"])
    # Get the total number of rows in the CSV file
    num_rows = sum(1 for line in open(file_path_old))
    latitude = []
    longitude = []
    # Read the CSV file from the last row to the first row

    dfold = pd.read_csv(file_path_old)
    #dfnew = pd.read_csv(file_path_new)
    for i in range(0, num_rows-1, 1):
        latitude.append(dfold.loc[i,"latitude"])
        longitude.append(dfold.loc[i,"longitude"])
    latitude.reverse()
    longitude.reverse()
    latitude_df = {"latitude": [latitude]}
    longitude_df = {"latitude": [longitude]}

    dfnew = pd.DataFrame(latitude)
    dfnew["longitude"] = longitude

    dfnew.to_csv(file_path_new)



#In this example, we first define a numpy array data with shape (3,2).
# Then, we create a datetime object for today's date and time using the
# datetime.datetime.now() function.

#Next, we create a new numpy array dt_data with an extra column to hold
# the datetime. We set the original data in the first columns of the new
# array and set the datetime in the last column of the new array.

#Then, we create a vector of data vector with shape (3,). We add this
# vector as a column to the right of the datetime array using the
# np.concatenate() function. Finally, we print the final array with
# shape (3,4) that has the original data in the first two columns,
# the datetime in the third column, and the vector data in the fourth column.
