import datetime
import datetime as datetime
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from route_handling import mac_windows_file_handle
import matplotlib as mpl
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

def plot_something(title, x_axis, y_axis, x_label, y_label):
    """
    :param title: Title of graph
    :param x_axis: What goes on the x-axis
    :param y_axis: What goes on the y-axis
    :param x_label: x_label
    :param y_label: y_label
    :return: Shows graph
    """

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

    plt.hist(df1['speed'], edgecolor='black', range=[0,10], bins=10)
    plt.title(title)
    plt.xlabel('Speed [Knots]')
    plt.ylabel('Occurences')
    plt.show()
    return 0

def plot_histogram(data1,data2,title):

    # Define the range of values to bin the data into
    bins = [0,2,4,6,8,10,12,14,16,18,20,22,24]

    # Calculate the histogram of the data
    hist_1, bin_edges_1 = np.histogram(data1, bins=bins)
    hist_2, bin_edges_2 = np.histogram(data2, bins=bins)

    # Calculate the percentage of values in each bin
    bin_percentages_1 = hist_1 / np.sum(hist_1) * 100
    bin_percentages_2 = hist_2 / np.sum(hist_2) * 100


    # Create a bar plot of the percentage of values in each bin
    #plt.bar(bins[:-1], bin_percentages_1,edgecolor='black', width=np.diff(bins), align='edge')
    #plt.bar(bins[:-1], bin_percentages_2,edgecolor='black', width=np.diff(bins), align='edge')

    # Create a bar plot of the percentage of values in each bin
    plt.bar(bins[:-1], bin_percentages_1, color = "red" , edgecolor='black', width=np.diff(bins), align='edge', label = "Measured")
    plt.bar(bins[:-1], bin_percentages_2, edgecolor='black', width=np.diff(bins), align='edge', alpha=0.5, hatch=':', label = "Calculated",
            lw=0.5)

    #plt.hist(data,edgecolor='black', range=[0,10], bins=bins)
    # Add axis labels and a title
    plt.legend(loc='upper right')
    plt.xlabel('Speed in knots')
    plt.ylabel('Percentage of Values')
    plt.title(title)

    #Show plot
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


#histogram(mac_windows_file_handle("Output_files/Floro_port/TWS_Floro_port.csv"),Test)




def plot_wind_data(filename):
    # Read the CSV file with semicolon separator and comma as decimal separator
    df = pd.read_csv(filename, sep=';', decimal=',')

    # Convert wind speed values to knots
    df['TWS'] = df['TWS'] * 1.94384

    # Extract wind speed and direction columns and remove missing values
    wind_speed = df['TWS'].dropna()
    #wind_dir = df['Wind direction'].dropna()

    # Plot wind speed histogram
    plt.hist(wind_speed, range=[0,20], bins=20)
    plt.title('Wind Speed')
    plt.xlabel('Speed (knots)')
    plt.ylabel('Frequency')
    plt.show()

    # Plot wind direction histogram
    plt.hist(wind_dir, range=[0,360], bins=36)
    plt.title('Wind Direction')
    plt.xlabel('Direction (degrees)')
    plt.ylabel('Frequency')
    plt.show()

def histogramangle(filename, title):
    """
    :param filename: name of output file containing data to plot
    :param title: title of histogram
    :return: plots histogram of data
    """

    columns = ['number', 'Angle']
    df1 = pd.read_csv(filename,header = None, names=columns)
    #speeds_observed = []
    #speeds_observed.append(df1["speed"].value_counts())
    #print(speeds_observed)

    plt.hist(df1['Angle'], edgecolor='black', range=[0,360], bins=36)
    plt.title(title)
    plt.xlabel('Angles')
    plt.ylabel('Occurences')
    plt.show()
    return 0


#plot_wind_data("../Output_files/Aberdeen_Færøyene/saveTWS_Aberdeen_Færøyene.csv")
#histogramangle("Output_files/Floro_port/TWD_Floro_port.csv", Test)
#plot_wind_data('Output_files/Floro_port/FloroPort_windspeed_winddirection.csv')

#filnavn = "Output_files/Newcastle_Aberdeen/SS_Newcastle_Aberdeen.csv"
#filnavn2 = "Output_files/Færøyene_Ålesund/savespeed_Færøyene_Ålesund.csv"
#navn    = "newcastle_aberdeen"
#navn2   = "Færøyene Ålesund"

#histogram(mac_windows_file_handle(filnavn),navn)
#histogram(mac_windows_file_handle(filnavn2),navn2)

