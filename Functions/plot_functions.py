import matplotlib.pyplot as plt
import pandas as pd
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


#plot power from Force_over_trip with respect to time
def plot_power(title, y_axis, x_label, y_label):
    """
    :param title: Title of graph
    :param y_axis: What goes on the y-axis
    :param x_label: x_label
    :param y_label: y_label
    :return: Shows graph
    """
    x_axis = []
    for i in range(1,10):
        x_axis.append(730*i)
    plt.plot(y_axis)
    plt.xticks(x_axis)
    plt.title(title)
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.show()
    return 0


#This function will output a histogram that show the speed distribution for a given route
def histogram(filename, title):
    """
    :param filename: name of output file containing data to plot
    :param title: title of histogram
    :return: plots histogram of data
    """
    columns = ['number', 'speed']
    df1 = pd.read_csv(filename,header = None, names=columns)
    #speeds_observed = []
    #speeds_observed.append(df1["speed"].value_counts())
    #print(speeds_observed)

    plt.hist(df1['speed'], edgecolor='black', range=[0,8], bins=8)
    plt.title(title)
    plt.xlabel('Speed [Knots]')
    plt.ylabel('Occurences')
    plt.show()
    return 0
Trond_Aalesund      = "Trondheim Ålesund"
Aalesund_Floro      = "Ålesund Florø"
Floro_Bergen        = "Florø Bergen"
Bergen_Stavanger    = "Bergen Stavanger"
vector = []
Test = "Test"

#histogram(mac_windows_file_handle('Output_files/Trondheim_Aalesund_reise', Trond_Aalesund))
#histogram(mac_windows_file_handle('Output_files/Aalesund_Floro_reise', Aalesund_Floro))
#histogram(mac_windows_file_handle('Output_files/Florø_Bergen/savespeed_FloroBergen.csv'), Floro_Bergen)
#histogram(mac_windows_file_handle('Output_files/savespeed_TrondAales.csv'), Trond_Aalesund)
#histogram(mac_windows_file_handle("Output_files/savespeed_TrondAales.csv"),Test)