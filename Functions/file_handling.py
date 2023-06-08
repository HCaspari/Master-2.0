import pandas as pd
import numpy as np
from datetime import datetime
import xarray as xr
import re
import csv
from route_handling import mac_windows_file_handle
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import threading
import folium

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
    filename_func = "../" + filename_func
    try:
        readdata = pd.read_csv(filename_func, usecols=['speed'])
    except ValueError:
        try:
            readdata = pd.read_csv(filename_func, usecols=["0"])
        except ValueError:
            readdata = pd.read_csv(filename_func, usecols=["average sailing speed"])
    data_array = [readdata]

    return data_array

#speed_optimal = read_array_from_file("../Output_files/Good_Route/savespeed_Good_Route.csv")
#print(np.mean(speed_optimal))
#speed_poor = read_array_from_file("../Output_files/Poor_Route/savespeed_Poor_Route.csv")
#print(np.mean(speed_poor))
#newbuild_speed = read_array_from_file("../Output_files/Good_Route_Newbuild/savespeed_Good_Route.csv")
#print(np.mean(newbuild_speed))

def find_average(name,filename):
    df = pd.read_csv("../"+filename)
    try:
        print("mean is:", df["speed"].mean(), "on", name)
        return df["speed"].mean()
    except KeyError:
        try:
            print(f"mean is:", df["0"].mean(), "on", name)
            return df["0"].mean()
        except KeyError:
            return 0

    return 0

#good = find_average("Good route", "Output_files/Good_Route/savespeed_Good_Route.csv")
#bad = find_average("poor route","Output_files/Poor_Route/savespeed_Poor_Route.csv")
#newbuild = find_average("Newbuild","Output_files/Good_Route_Newbuild/savespeed_Good_Route.csv")
#battery = find_average("Battery","Output_files/Færøyene_Ålesund_With_Battery/savespeed_Good_Route.csv")
##kite = find_average("With kite", "Output_files/Good_Route_With_Kite/savespeed_Good_Route.csv")

#aber_Fær = find_average("Aberdeen-Færøyene", "Output_files/Aberdeen_Færøyene/savespeed_Aberdeen_Færøyene.csv")
#amst_new = find_average("Amsterdam-Newcastle", "Output_files/Amsterdam_Newcastle/savespeed_Amsterdam_Newcastle.csv")
#berg_stvg = find_average("Bergen_Stavanger","Output_files/Bergen_Stavanger/savespeed_BrgStvg.csv")
#dk_amst = find_average("Danmark_Amsterdam","Output_files/Danmark_Amsterdam/savespeed_Danmark_Amsterdam.csv")
#Florø_Brg = find_average("Florø-Bergen", "Output_files/Florø_Bergen/savespeed_FloroBergen.csv")
#new_aber = find_average("Newcastle to Aberdeen", "Output_files/Newcastle_Aberdeen/SS_Newcastle_Aberdeen.csv")
#Trond_aales = find_average("Trondheim_Ålesund_iteration_100","Output_files/Trondheim_Ålesund_iteration_100/savespeed_TrondAales.csv")
#aales_Floro = find_average("Ålesund til Florø","Output_files/Ålesund_Florø/savespeed_AalesFloro.csv")


#Fær_åles = find_average("Færøyene Ålesund", "Output_files/Færøyene_Ålesund/savespeed_Færøyene_Ålesund.csv")
#åles_Dan = find_average("Ålesund Danmark", "Output_files/Ålesund_Danmark/savespeed_Ålesund_Danmark.csv")

#print("Average of all", (good+bad+newbuild+aber_Fær+amst_new+berg_stvg+dk_amst+Florø_Brg+new_aber+Trond_aales+aales_Floro+Fær_åles+åles_Dan)/13)

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




#dictionary = {date of start: datetime, sailing speeds observed: np.array, average speed route iteration: float}

#01012020, [0.4,3.1,3.4,....] , 5 knop
#01012020, [0.4,3.1,3.4,....] , 1 knop
#01012020, [0.4,3.1,3.4,....] , 5.3 knop
#01012020, [0.4,3.1,3.4,....] , 2.3 knop
#01012020, [0.4,3.1,3.4,....] , 3 knop
#01012020, [0.4,3.1,3.4,....] , 4 knop

#reliability given as percentage of time vessel sails route over x knots on average.
#reliability for route x, speed 5 knots is 30%
#reliability for route x, speed 4 knots is 55%

#compare this is corresponding for normal shipping.
#https://atlas-network.com/container-carriers-sailing-reliability-over-and-beyond/

def write_to_dictionary(filename):
    name = filename[16:]
    end_index = name.find("/")
    name = name[:end_index]
    print("starting write_to_dictionary of", name)
    routename = {}
    df = pd.read_csv(filename)
    length_dataframe = len(df)
    VS_total   = 0
    k = 0
    iterations_this_loop = 0
    for i in range(1,length_dataframe):     #loop over
        timestamp_previous    = df.loc[i-1]["timestamp"]
        timestamp_next        = df.loc[i]["timestamp"]
        sailing_speed = df.loc[i-1]["0"]
        VS_total += sailing_speed
        iterations_this_loop += 1
        if timestamp_next < timestamp_previous: #if simulation resets, calculate avg speed
            index = df.loc[k]["timestamp"]
            routename[index] = {}
            routename[index]["average sailing speed"] = VS_total/iterations_this_loop
            #reset indices
            k = i
            VS_total = 0
            iterations_this_loop = 0

        if i % 50000 == 0:
            print("progress", i, "of", length_dataframe, datetime.now())

    df = pd.DataFrame.from_dict(routename, orient="index")

    df.to_csv(mac_windows_file_handle("Output_files/"+name+"/avgspeed.csv"))
    return 0

def run_write_to_multiple(filename):
    write_to_dictionary(filename)
    print(f"Simulation of route{filename}")

def run_multiple():
    a = "Output_files/Færøyene_Ålesund_return/savespeed_Færøyene_Ålesund_return.csv"
    b = "Output_files/Færøyene_Ålesund/savespeed_Færøyene_Ålesund.csv"
    inputs = [(a),(b)]
               #(mac_windows_file_handle("Output_files/Newcastle_Aberdeen/SS_Newcastle_Aberdeen.csv")),
               #(mac_windows_file_handle("Output_files/Newcastle_Aberdeen_return/SS_Newcastle_Aberdeen_return.csv")),
               #(mac_windows_file_handle("Output_files/Trondheim_Ålesund_iteration_100/savespeed_TrondAales.csv")),
               #(mac_windows_file_handle("Output_files/Trondheim_Ålesund_return/savespeed_TrondAales_return.csv")),
               #(mac_windows_file_handle("Output_files/Ålesund_Danmark/savespeed_Ålesund_Danmark.csv")),
               #(mac_windows_file_handle("Output_files/Ålesund_Danmark_return/savespeed_Ålesund_Danmark_return.csv")),
               #(mac_windows_file_handle("Output_files/Ålesund_Florø/savespeed_AalesFloro.csv")),
               #(mac_windows_file_handle("Output_files/Ålesund_Florø_return/savespeed_AalesFloro_return.csv"))]

    threads = []

    for input in inputs:
        a = input
        thread = threading.Thread(target=run_write_to_multiple, args=a)
        thread.start()
        threads.append(thread)

    for thread in threads:
        thread.join()
    return 0

def write_all_to_dict():
    #print(datetime.now())
    ##write_to_dictionary(mac_windows_file_handle("Output_files/Aberdeen_Færøyene/savespeed_Aberdeen_Færøyene.csv"))
    ##print(datetime.now())
    ##write_to_dictionary(mac_windows_file_handle("Output_files/Aberdeen_Færøyene_return/savespeed_Aberdeen_Færøyene_return.csv"))
    ##print(datetime.now())
    #write_to_dictionary(mac_windows_file_handle("Output_files/Amsterdam_Newcastle/savespeed_Amsterdam_Newcastle.csv"))
    #print(datetime.now())
    #write_to_dictionary(mac_windows_file_handle("Output_files/Amsterdam_Newcastle_return/savespeed_Amsterdam_Newcastle_return.csv"))
    #print(datetime.now())
    #write_to_dictionary(mac_windows_file_handle("Output_files/Bergen_Stavanger/savespeed_BrgStvg.csv"))
    #print(datetime.now())
    #write_to_dictionary(mac_windows_file_handle("Output_files/Bergen_Stavanger_return/savespeed_BrgStvg_return.csv"))
    #print(datetime.now())
    #write_to_dictionary(mac_windows_file_handle("Output_files/Danmark_Amsterdam/savespeed_Danmark_Amsterdam.csv"))
    #print(datetime.now())
    #write_to_dictionary(mac_windows_file_handle("Output_files/Danmark_Amsterdam_return/savespeed_Danmark_Amsterdam_return.csv"))
    #print(datetime.now())
    #write_to_dictionary(mac_windows_file_handle("Output_files/Floro_port/savespeed_port.csv"))
    #print(datetime.now())
    #write_to_dictionary(mac_windows_file_handle("Output_files/Florø_Bergen/savespeed_FloroBergen.csv"))
    #print(datetime.now())

    #run This next

    #write_to_dictionary(mac_windows_file_handle("Output_files/Florø_Bergen_return/savespeed_FloroBergen_return.csv"))
    #print(datetime.now())
    #write_to_dictionary(mac_windows_file_handle("Output_files/Færøyene_Ålesund/savespeed_Færøyene_Ålesund.csv"))
    #print(datetime.now())
    #write_to_dictionary(mac_windows_file_handle("Output_files/Færøyene_Ålesund_return/savespeed_Færøyene_Ålesund_return.csv"))
    #print(datetime.now())
    #write_to_dictionary(mac_windows_file_handle("Output_files/Newcastle_Aberdeen/SS_Newcastle_Aberdeen.csv"))
   # print(datetime.now())
   # write_to_dictionary(mac_windows_file_handle("Output_files/Newcastle_Aberdeen_return/SS_Newcastle_Aberdeen_return.csv"))
   # print(datetime.now())
   # write_to_dictionary(mac_windows_file_handle("Output_files/Trondheim_Ålesund_iteration_100/savespeed_TrondAales.csv"))
   # print(datetime.now())
   # write_to_dictionary(mac_windows_file_handle("Output_files/Trondheim_Ålesund_return/savespeed_TrondAales_return.csv"))
   # print(datetime.now())
    #write_to_dictionary(mac_windows_file_handle("Output_files/Ålesund_Danmark/savespeed_Ålesund_Danmark.csv"))
    #print(datetime.now())
    #write_to_dictionary(mac_windows_file_handle("Output_files/Ålesund_Danmark_return/savespeed_Ålesund_Danmark_return.csv"))
    #print(datetime.now())
    write_to_dictionary(mac_windows_file_handle("Output_files/Ålesund_Florø/savespeed_AalesFloro.csv"))
    print(datetime.now())
    write_to_dictionary(mac_windows_file_handle("Output_files/Ålesund_Florø_return/savespeed_AalesFloro_return.csv"))
    print(datetime.now())
    return 0

def combine_troll(year_Place, year):
    # Define the file paths
    file_paths = [f'../Weather Check/{year_Place}/NO_TS_MO_Troll-A_{year}01.nc',
                  f'../Weather Check/{year_Place}/NO_TS_MO_Troll-A_{year}02.nc',
                  f'../Weather Check/{year_Place}/NO_TS_MO_Troll-A_{year}03.nc',
                  f'../Weather Check/{year_Place}/NO_TS_MO_Troll-A_{year}04.nc',
                  f'../Weather Check/{year_Place}/NO_TS_MO_Troll-A_{year}05.nc',
                  f'../Weather Check/{year_Place}/NO_TS_MO_Troll-A_{year}06.nc',
                  f'../Weather Check/{year_Place}/NO_TS_MO_Troll-A_{year}07.nc',
                  f'../Weather Check/{year_Place}/NO_TS_MO_Troll-A_{year}08.nc',
                  f'../Weather Check/{year_Place}/NO_TS_MO_Troll-A_{year}09.nc',
                  f'../Weather Check/{year_Place}/NO_TS_MO_Troll-A_{year}10.nc',
                  f'../Weather Check/{year_Place}/NO_TS_MO_Troll-A_{year}11.nc',
                  f'../Weather Check/{year_Place}/NO_TS_MO_Troll-A_{year}12.nc']

    # Load the datasets into a list
    datasets = [xr.open_dataset(fp) for fp in file_paths]

    # Extract the WSPD variable data for each dataset using a loop
    TIME_list = []
    WSPD_list = []
    WDIR_list = []


    for ds in datasets:
        WSPD = ds['WSPD'][:, 0].values
        WDIR = ds['WDIR'][:, 0].values
        TIME = ds['TIME'][:].values
        WSPD_list.append(WSPD)
        WDIR_list.append(WDIR)
        TIME_list.append(TIME)

    WSPD_concat = np.concatenate(WSPD_list)
    WDIR_concat = np.concatenate(WDIR_list)
    TIME_concat = np.concatenate(TIME_list)

    return WSPD_concat,WDIR_concat,TIME_concat

def combine_Sleipnir():
    # Define the file paths
    file_paths = ['../Weather Check/2021 Sleipner/NO_TS_MO_Sleipner-A_202101.nc',
                  '../Weather Check/2021 Sleipner/NO_TS_MO_Sleipner-A_202102.nc',
                  '../Weather Check/2021 Sleipner/NO_TS_MO_Sleipner-A_202103.nc',
                  '../Weather Check/2021 Sleipner/NO_TS_MO_Sleipner-A_202104.nc',
                  '../Weather Check/2021 Sleipner/NO_TS_MO_Sleipner-A_202105.nc',
                  '../Weather Check/2021 Sleipner/NO_TS_MO_Sleipner-A_202106.nc',
                  '../Weather Check/2021 Sleipner/NO_TS_MO_Sleipner-A_202107.nc',
                  '../Weather Check/2021 Sleipner/NO_TS_MO_Sleipner-A_202108.nc',
                  '../Weather Check/2021 Sleipner/NO_TS_MO_Sleipner-A_202109.nc',
                  '../Weather Check/2021 Sleipner/NO_TS_MO_Sleipner-A_202110.nc',
                  '../Weather Check/2021 Sleipner/NO_TS_MO_Sleipner-A_202112.nc']


    # Load the datasets into a list
    datasets = [xr.open_dataset(fp) for fp in file_paths]

    # Extract the WSPD variable data for each dataset using a loop
    TIME_list = []
    WSPD_list = []
    WDIR_list = []


    for ds in datasets:
        # Get the names of all variables in the NetCDF file
        variable_names = list(ds.variables.keys())

        # Print the names of all variables in the NetCDF file
        print(variable_names)
        WSPD = ds['WSPD'][:, 0].values
        WDIR = ds['WDIR'][:, 0].values
        TIME = ds['TIME'][:].values
        try:
            WSPD_list.append(WSPD)
            WDIR_list.append(WDIR)
            TIME_list.append(TIME)
        except KeyError:
            print("nah")




    WSPD_S_concatenated = np.concatenate(WSPD_list)
    WDIR_S_concatenated = np.concatenate(WDIR_list)
    TIME_S_concatenated = np.concatenate(TIME_list)
    return WSPD_S_concatenated,WDIR_S_concatenated,TIME_S_concatenated

#WSPD_concat_2021,WDIR_concat_2021,TIME_concat_2021 = combine_troll("2021 Troll",2021)

# Read the CSV file into a DataFrame
#df = pd.read_csv('../Route_data/Poor_Route/Route_Poor_Route.csv')

# Drop the column you want to remove
#df = df.drop(df.columns[0], axis=1)

# Write the updated DataFrame to a new CSV file
#df.to_csv('../Route_data/Poor_Route/Good_route.csv', index=False)

#a = read_array_from_file("../Output_files/Færøyene_Ålesund_With_Battery/Battery_use.csv")
#print(np.average(a))

#print("Trond Åles", np.average(read_array_from_file("Output_files/Trondheim_Ålesund_iteration_100/avgspeed.csv")))

#Trond_åles = read_array_from_file("Output_files/Trondheim_Ålesund_iteration_100/savespeed_TrondAales.csv")
    #print(type(Fær_Aales))
#print(type(Trond_åles))
#print(len(Trond_åles[0][:]))

#Trond_åles = read_array_from_file("Output_files/Trondheim_Ålesund_iteration_100/savespeed_TrondAales.csv")
#print(type(Trond_åles))
#print(Trond_åles.shape)
#print(Trond_åles[0][:])
#Trond_åles_t = Trond_åles.T
#print("here",Trond_åles_t[0][:])
#Trond_åles_list = [item for sublist in Trond_åles for sublist_2 in sublist for item in sublist_2]
#print(Trond_åles_list[:10])
def read_array_As_list(file,rownumber):
    data = []

    with open("../"+ file, 'r') as file1:
        reader = csv.reader(file1)
        next(reader)  # Skip the header row

        for row in reader:
            value = row[rownumber]  # Access the third column (index 2)
            if type(value) != str:
                value = float(value)
            data.append(value)
    return data

def plot_box():


    Trond_åles = read_array_As_list("Output_files/Trondheim_Ålesund_iteration_100/savespeed_TrondAales.csv",2)
    Åles_Floro = read_array_As_list("Output_files/Ålesund_Florø/savespeed_AalesFloro.csv",2)
    Floro_Bergen = read_array_As_list("Output_files/Florø_Bergen/savespeed_FloroBergen.csv",2)
    Bergen_stvg = read_array_As_list("Output_files/Bergen_Stavanger/savespeed_BrgStvg.csv",2)
    Aber_Fær = read_array_As_list("Output_files/Aberdeen_Færøyene/savespeed_Aberdeen_Færøyene.csv",2)
    Amster_New = read_array_As_list("Output_files/Amsterdam_Newcastle/savespeed_Amsterdam_Newcastle.csv",2)
    Dan_Amster = read_array_As_list("Output_files/Danmark_Amsterdam/savespeed_Danmark_Amsterdam.csv",2)
    New_Aber = read_array_As_list("Output_files/Newcastle_Aberdeen/SS_Newcastle_Aberdeen.csv",2)
    Aales_Danm = read_array_As_list("Output_files/Ålesund_Danmark/savespeed_Ålesund_Danmark.csv",2)
    Fær_Aales = read_array_As_list("Output_files/Færøyene_Ålesund/savespeed_Færøyene_Ålesund.csv",2)




    #Fær_Aales_newbuild = read_array_As_list("Output_files/Færøyene_Ålesund_Newbuild/savespeed_Færøyene_Ålesund.csv")
    #Fær_Aales_Battery   = read_array_As_list("Output_files/Færøyene_Ålesund_With_Battery/savespeed_Good_Route.csv")
    #Fær_Aales_Kite      = read_array_As_list("Output_files/Færøyene_Ålesund_With_Kite/savespeed_Good_Route.csv")

    # Combine all route data into a list
    data = [Trond_åles,Åles_Floro,Floro_Bergen,Bergen_stvg,Aber_Fær,Amster_New,Dan_Amster
        ,New_Aber,Aales_Danm,Fær_Aales]#, Fær_Aales_Battery, Fær_Aales_Kite, Fær_Aales_newbuild]

    mean = [np.array(np.mean(route)) for route in data]
    print(mean)
    data = [np.array(route) for route in data]

    # Transpose the data
    #data = np.array(data).T

    # Create a figure and axes
    fig, ax = plt.subplots(figsize=(10, 8))

    # Add the mean as a red point on the plot for each data series
    #for i, m in enumerate(mean):
    #    ax.plot(i + 1, m, 'ro', label='Mean')

    # Create a boxplot
    ax.boxplot(data, sym='k.')

    # Set labels and title
    ax.set_xlabel('Routes')
    ax.set_ylabel('Sailing Speed (knots)')
    ax.set_title('Sailing Speed')

    # Set x-axis tick labels
    x_labels = ["Trondheim to Ålesund", "Ålesund to Floro", "Floro to Bergen", "Bergen  to Stavanger", "Aberdeen to Faraoe", "Amsterdam to Newcastle",
                        "Denmark to Amsterdam", "Newcastle to Aberdeen", "Ålesund to Denmark", "Faraoe to Ålesund"]#, "Faraoe-Ålesund_Battery","Faraoe-Ålesund_Kite", "Faraoe-Ålesund_newbuild"]
    ax.set_xticklabels(x_labels)
    ax.set_xticklabels(x_labels, rotation=30, ha='right')

    # Adjust the bottom margin to make room for the table
    plt.subplots_adjust(bottom=0.2)

    # Display the plot
    plt.show()


def plot_box_return():
    Trondheim_Åles_Retur = read_array_As_list("Output_files/Trondheim_Ålesund_return/savespeed_TrondAales_return.csv",2)
    Åles_Floro_Retur = read_array_As_list("Output_files/Ålesund_Florø_return/savespeed_AalesFloro_return.csv",2)
    Floro_Bergen_Retur = read_array_As_list("Output_files/Florø_Bergen_return/savespeed_FloroBergen_return.csv",2)
    Bergen_Stavanger_Retur = read_array_As_list("Output_files/Bergen_Stavanger_return/savespeed_BrgStvg_return.csv",2)
    Aber_Fær_retur = read_array_As_list("Output_files/Aberdeen_Færøyene_return/savespeed_Aberdeen_Færøyene_return.csv",2)
    Dan_Amster_retur = read_array_As_list("Output_files/Danmark_Amsterdam_return/savespeed_Danmark_Amsterdam_return.csv",2)
    Amster_New_retur = read_array_As_list("Output_files/Amsterdam_Newcastle_return/savespeed_Amsterdam_Newcastle_return.csv",2)
    New_Aber_retur = read_array_As_list("Output_files/Amsterdam_Newcastle_return/savespeed_Amsterdam_Newcastle_return.csv",2)
    Aales_Danm_retur = read_array_As_list("Output_files/Ålesund_Danmark_return/savespeed_Ålesund_Danmark_return.csv",2)
    Fær_Aales_retur = read_array_As_list("Output_files/Færøyene_Ålesund_return/savespeed_Færøyene_Ålesund_return.csv",2)

    # Combine all route data into a list
    data = [Trondheim_Åles_Retur, Åles_Floro_Retur, Floro_Bergen_Retur, Bergen_Stavanger_Retur, Aber_Fær_retur,
            Amster_New_retur, Dan_Amster_retur, New_Aber_retur, Aales_Danm_retur, Fær_Aales_retur]


    data = [np.array(route) for route in data]
    average = [np.mean(route) for route in data]
    print("average return", average)



    # Create a figure and axes
    fig, ax = plt.subplots(figsize=(10, 8))


    # Create a boxplot
    ax.boxplot(data, sym='k.')

    # Set labels and title
    ax.set_xlabel('Return Routes')
    ax.set_ylabel('Sailing Speed (knots)')
    ax.set_title('Sailing Speed')

    # Set x-axis tick labels
    x_labels = ["Ålesund to Trondheim", "Floro to Ålesund", "Bergen to Floro", "Stavanger to Bergen",
                "Faraoe to Aberdeen ", "Newcastle to Amsterdam ",
                "Amsterdam to Denmark", "Aberdeen to Newcastle", "Denmark to Ålesund",
                "Ålesund to Faraoe"]
    ax.set_xticklabels(x_labels)
    ax.set_xticklabels(x_labels, rotation=30, ha='right')

    # Adjust the bottom margin to make room for the table
    plt.subplots_adjust(bottom=0.2)

    # Display the plot
    plt.show()



#plot_box()
#plot_box_return()

