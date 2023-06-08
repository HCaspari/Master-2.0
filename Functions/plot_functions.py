
#import datetime as datetime
from datetime import datetime
from datetime import timedelta
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import numpy
import numpy as np
import pandas as pd
from route_handling import mac_windows_file_handle
import matplotlib as mpl
import math
import csv
from file_handling import read_array_from_file, read_array_As_list
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

def plot_something(title, y_axis, x_label, y_label):
    """
    :param title: Title of graph
    :param x_axis: What goes on the x-axis
    :param y_axis: What goes on the y-axis
    :param x_label: x_label
    :param y_label: y_label
    :return: Shows graph
    """
    plt.ylim(0,20)
    plt.plot(y_axis)
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

    plt.hist(df1['speed'], edgecolor='black', range=[0,12], bins=12)
    plt.title(title)
    plt.xlabel('Speed [knots]')
    plt.ylabel('Occurences')
    plt.show()
    return 0




#Trond_Aalesund      = "Trondheim Ålesund"
#Aalesund_Floro      = "Ålesund Florø"
#Floro_Bergen        = "Florø Bergen"
#Bergen_Stavanger    = "Bergen Stavanger"
#vector = []

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

    plt.hist(df1['Angle'], edgecolor='black', range=[0,360], bins=36)
    plt.title(title)
    plt.xlabel('Angles')
    plt.ylabel('Occurences')
    plt.show()
    return 0

# assume your data is stored in a list called WSPD_Vect_Troll_Old
def plot_Vect_Weekly(datavector_one, datavector_two, xlabel, ylabel, title1, title2, title):

    weekly_avg_one = []
    weekly_avg_two = []
    if len(datavector_one) > len(datavector_two):
        datavector_one = datavector_one[:len(datavector_two)]
    if len(datavector_two) > len(datavector_one):
        datavector_two = datavector_two[:len(datavector_one)]

    # Create a list of indices corresponding to each day
    indices = [i for i in range(0, len(datavector_one), 168)]

    # Calculate the average for each week in first list (168 hours)
    for i in range(len(datavector_one)):
        #  Checks vector for nan values, if discovered, uses prior value
        if math.isnan(datavector_one[i]):
            datavector_one[i] = datavector_one[i-1]
        if math.isnan(datavector_two[i]):
            datavector_two[i] = datavector_two[i-1]

    #check for remaining nan-values
    for i in range(len(datavector_one)):
        #  Checks vector for nan values, if discovered, uses prior value
        if math.isnan(datavector_one[i]):
            print("nan value at index",i)
        if math.isnan(datavector_two[i]):
            print("nan value at index",i)

    for i in indices:
        weekly_avg_one.append(sum(datavector_one[i:i + 168]) / 168)
        weekly_avg_two.append(sum(datavector_two[i:i + 168]) / 168)

    # Create the figure and axes objects
    fig, ax = plt.subplots()

    # Plot the first graph
    ax.plot(weekly_avg_one, label=title1)

    # Plot the second graph
    ax.plot(weekly_avg_two, label=title2, color='red')

    # Add a legend
    ax.legend()

    # Calculate the correlation coefficient

    corr = np.corrcoef(datavector_one, datavector_two)[0, 1]

    # Display the correlation coefficient
    print(f"The correlation coefficient between y1 and y2 is {corr}")

    # Show the plot

    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.show()

    return 0

def plot_Vect_Daily(datavector_one, datavector_two, xlabel, ylabel, title1, title2, title):

    if len(datavector_one) > len(datavector_two):
        datavector_one = datavector_one[:len(datavector_two)]
    if len(datavector_two) > len(datavector_one):
        datavector_two = datavector_two[:len(datavector_one)]

    # Create a list of indices corresponding to each day
    indices = [i for i in range(0, len(datavector_one), 24)]

    # Calculate the average for each day
    daily_avg = [sum(datavector_one[i:i + 24]) / 24 for i in indices]

    # Calculate the average for each day in the second list
    second_daily_avg = [sum(datavector_two[i:i + 24]) / 24 for i in indices]

    # Create the figure and axes objects
    fig, ax = plt.subplots(figsize=(9, 8))

    # Manually set the dates for the x-axis
    #dates = ['2020-07-01', '2020-08-20', '2020-10-9', '2020-11-28', '2020-01-17', '2021-02-05', '2021-03-08', "2021-04-27", "2021-06-16"]

    # Generate the corresponding dates for the x-axis
    start_date = datetime.date(2020, 7, 1)
    dates = [start_date + datetime.timedelta(days=i) for i in range(len(daily_avg))]

    # Format the x-axis labels as dates
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))

    # Rotate the x-axis tick labels if needed
    plt.xticks(rotation=45)

    # Plot the first graph
    ax.plot(dates, daily_avg, label=title1)

    # Plot the second graph
    ax.plot(dates, second_daily_avg, label=title2, color='red')

    # Add a legend
    ax.legend()

    # Calculate the correlation coefficient

    corr = np.corrcoef(datavector_one, datavector_two)[0, 1]

    # Display the correlation coefficient
    print(f"The correlation coefficient between y1 and y2 is {corr}")

    # Show the plot

    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)

    # Adjust the spacing to push the plot up
    plt.subplots_adjust(bottom=0.15)

    #Show graph
    plt.show()

    return 0

def plot_Vect_hourly(datavector_one, datavector_two, xlabel, ylabel, title1, title2, title):
    if len(datavector_one) > len(datavector_two):
        datavector_one = datavector_one[:len(datavector_two)]
    if len(datavector_two) > len(datavector_one):
        datavector_two = datavector_two[:len(datavector_one)]

    # Generate corresponding dates
    initial_date = datetime.datetime(2020, 7, 1)
    days = list(range(len(datavector_one)))
    dates = [initial_date + datetime.timedelta(days=d) for d in days]

    # Create the figure and axes objects
    fig, ax = plt.subplots(figsize=(10,12))

    # Plot the first graph with corresponding dates
    ax.plot(dates, datavector_one, label=title1)

    # Plot the second graph with corresponding dates
    ax.plot(dates, datavector_two, label=title2, color='red')

    # Add a legend
    ax.legend()

    # Format the x-axis as dates
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))
    ax.xaxis.set_major_locator(mdates.AutoDateLocator())

    # Show the plot
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)

    # Calculate the correlation coefficient
    corr = np.corrcoef(datavector_one, datavector_two)[0, 1]

    # Display the correlation coefficient
    print(f"The correlation coefficient between {title1} and {title2} is {corr}")

    plt.title(title + " hourly")
    plt.show()

    return 0


def plot_Vect_hourly_single(datavector_one, xlabel, ylabel, title1, title):


    # Create the figure and axes objects
    fig, ax = plt.subplots()

    # Plot the first graph
    ax.plot(datavector_one, label=title1, color = "red")

    # Add a legend
    ax.legend()

    # Show the plot
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)

    plt.ylim(0, 5)
    plt.title(title + " hourly")
    plt.show()

    return 0
#Define the data vectors
#months = range(31)
#VLSFO = [0,228257.6625,456515.325,684772.9875,913030.65,1141288.313,1369545.975,1597803.638,1826061.3,2054318.963,2282576.625,2510834.2882739091,952967349.613,3195607.275,3423864.938,3652122.6,3880380.263,4108637.925,4336895.588,4565153.25,4793410.913,5021668.575,5249926.238,5478183.9,5706441.563,5934699.225,6162956.888,6391214.55,6619472.213,6847729.875,7075987.538,7304245.2,7532502.863,7760760.525,7989018.188,8217275.85,8445533.513,8673791.175,8902048.838,9130306.5]
#Battery = [4500000,4554225.6,4608451.2,4662676.8,4716902.4,4771128,4825353.6,4879579.2,4933804.8,4988030.4,5042256,5096481.6,5150707.2,5204932.8,5259158.4,5313384,5367609.6,5421835.2,5476060.8,5530286.4,5584512,5638737.6,5692963.2,5747188.8,5801414.4,5855640,5909865.6,5964091.2,6018316.8,6072542.4,6126768,6180993.6,6235219.2,6289444.8,6343670.4,6397896,6452121.6,6506347.2,6560572.8,6614798.4,6669024]
#e_Ammonia = [4500000,4521473.338,4542946.675,4564420.013,4585893.35,4607366.688,4628840.026,4650313.363,4671786.701,4693260.038,4714733.376,4736206.714,4757680.051,4779153.389,4800626.726,4822100.064,4843573.402,4865046.739,4886520.077,4907993.415,4929466.752,4950940.09,4972413.427,4993886.765,5015360.103,5036833.44,5058306.778,5079780.115,5101253.453,5122726.791,5144200.128,5165673.466,5187146.803,5208620.141,5230093.479,5251566.816,5273040.154,5294513.491,5315986.829,5337460.167,5358933.504]
#E_LH_2 = [4500000,4521473.338,4542946.675,4564420.013,4585893.35,4607366.688,4628840.026,4650313.363,4671786.701,4693260.038,4714733.376,4736206.714,4757680.051,4779153.389,4800626.726,4822100.064,4843573.402,4865046.739,4886520.077,4907993.415,4929466.752,4950940.09,4972413.427,4993886.765,5015360.103,5036833.44,5058306.778,5079780.115,5101253.453,5122726.791,5144200.128,5165673.466,5187146.803,5208620.141,5230093.479,5251566.816,5273040.154,5294513.491,5315986.829,5337460.167,5358933.504]
#e_methanol = [4500000,4532210.006,4564420.013,4596630.019,4628840.026,4661050.032,4693260.038,4725470.045,4757680.051,4789890.058,4822100.064,4854310.07,4886520.077,4918730.083,4950940.09,4983150.096,5015360.103,5047570.109,5079780.115,5111990.122,5144200.128,5176410.135,5208620.141,5240830.147,5273040.154,5305250.16,5337460.167,5369670.173,5401880.179,5434090.186,5466300.192,5498510.199,5530720.205,5562930.211,5595140.218,5627350.224,5659560.231,5691770.237,5723980.244,5756190.25,5788400.256]

def plot_fuelprices():
    # Define the data
    years = ['2020-2030', '2030-2040', '2040-2050']
    vlsfo = [0.62, 1.23, 1.77]
    vlsfo_no_tax  = [0.48,0.81,1.08]
    vlsfo_without_flettner = [1.64,3.24,4.66]
    fossillng = [0.42, 0.81, 1.19]
    fossillng_no_tax = [0.32,0.51,0.70]
    fossillng_without_flettner = [1.10, 2.13, 3.12]
    emethanol = [2.75, 2.39, 2.02]
    elng = [2.26, 1.95, 1.63]
    eammoniakk = [1.79, 1.51, 1.21]
    batteries = [0.72, 0.83, 0.94]
    ngammoniakk = [0.98, 0.97, 0.97]


    ##vlsfo = [10,10,10]
    #vlsfo_no_tax  = [10,10,10]
    #vlsfo_without_flettner = [0,0,0]
    #fossillng = [10,10,10]
    #fossillng_no_tax = [10,10,10]
    #fossillng_without_flettner = [0,0,0]
    #emethanol = [10,10,10]
    #elng = [10,10,10]
    #eammoniakk = [10,10,10]
    #batteries = [10,10,10]
    #ngammoniakk = [0,0,0]

    # Plot the data
    fig, ax = plt.subplots(figsize=(10, 6))

    ax.plot(years, vlsfo, color="red", label="VLSFO + carbon tax")
    ax.plot(years, vlsfo_no_tax, color="red", linestyle="--", label="VLSFO no tax")
    ax.plot(years, vlsfo_without_flettner, color="red",linestyle="-.", label="VLSFO no flettner")
    ax.plot(years, fossillng, color="peachpuff", label="LNG no tax")
    ax.plot(years, fossillng_no_tax, color="peachpuff",linestyle="--", label="LNG no tax")
    ax.plot(years, fossillng_without_flettner, color="peachpuff",linestyle="-.", label="LNG + carbon tax")
    ax.plot(years, emethanol, color="limegreen", label="E-methanol")
    ax.plot(years, elng, color="seagreen", label="E-LNG")
    ax.plot(years, eammoniakk, color="darkgreen", label="E-ammoniakk")
    ax.plot(years, ngammoniakk, color="blue", label="NG-ammoniakk (CCS)")
    ax.plot(years, batteries, color="cyan", label="Batteries")

    # Add labels and titles
    ax.set_ylabel("MUSD")
    ax.set_title("VOYEX by period")

    #ax.legend(loc='center left', bbox_to_anchor=(0.5, 0.5))
    plt.show()

    return 0

def plot_Vect_daily_single(datavector_one, xlabel, ylabel, title):
    # Create the figure and axes objects
    fig, ax = plt.subplots()

    # Create a list of indices corresponding to each day
    indices = [i for i in range(0, len(datavector_one), 24)]

    # Calculate the average for each day
    daily_avg = [sum(datavector_one[i:i + 24]) / 24 for i in indices]

    # Plot the first graph
    ax.plot(datavector_one, label=title, color="red")

    # Add a legend
    ax.legend()

    # Show the plot
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)

    plt.ylim(0, 25)
    plt.title(title + " daily")
    plt.show()

    return 0

def plot_histogram(data1,data2,title,data1label,data2label,max,interval):

    # Define the range of values to bin the data into
    bins = []
    for i in range(0,max,interval):
        bins.append(i)

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
    plt.bar(bins[:-1], bin_percentages_1, color = "red" , edgecolor='black', width=np.diff(bins), align='edge', label = data1label)
    plt.bar(bins[:-1], bin_percentages_2, edgecolor='black', width=np.diff(bins), align='edge', alpha=0.5, label = data2label,
            lw=0.5)

    #plt.hist(data,edgecolor='black', range=[0,10], bins=bins)
    # Add axis labels and a title
    plt.legend(loc='upper right')
    plt.xlabel('Sailing speed in knots')
    plt.ylabel('Percentage of values')
    plt.title(title)

    #Show plot
    plt.show()
    return 0

def plot_histogram_single(data1,title, label):
    # Define the range of values to bin the data into
    bins = []
    for i in range(10,40,2):
        bins.append(i/10)
    # Calculate the histogram of the data
    hist_1, bin_edges_1 = np.histogram(data1, bins=bins)


    # Calculate the percentage of values in each bin
    bin_percentages_1 = hist_1 / np.sum(hist_1) * 100

    # Create a bar plot of the percentage of values in each bin
    plt.bar(bins[:-1], bin_percentages_1, edgecolor='black', width=np.diff(bins), align='edge', label =  label)


    #plt.hist(data,edgecolor='black', range=[0,10], bins=bins)
    # Add axis labels and a title
    plt.legend(loc='upper right')
    plt.xlabel('Sailed speed in knots')
    plt.ylabel('Percentage of Values')
    plt.title(title)

    #Show plot
    plt.show()
    return 0

def print_averages():
    print("Trond Åles", np.average(Trond_Aales))
    print("Trond Åles retur", np.average(Trond_Aales_retur))
    print("Åles Floro", np.average(Aales_Floro))
    print("Åles Floro retur", np.average(Aales_Floro_retur))
    print("Floro bergen ", np.average(Floro_Bergen))
    print("Floro bergen retur ", np.average(Floro_Bergen_retur))
    print("Berg stav", np.average(Bergen_Stvg))
    print("Berg stav retur", np.average(Bergen_Stvg_retur))

    print("Aber - Fær", np.average(Aber_Faer))
    print("Aber - Fær_retur", np.average(Aber_Faer_retur))
    print("Amster - New",np.average(Amster_New))
    print("Amster - New return",np.average(Amster_New_retur))
    print("Danmark - Amster",np.average(Danmark_Amster))
    print("Danmark - Amster retur",np.average(Danmark_Amster_retur))
    print("Fær - Åles", np.average(Faer_Åles))
    print("Fær - Åles return", np.average(Faer_Åles_retur))
    print("New - Aber", np.average(Newcast_Aber))
    print("New - Aber retur", np.average(Newcast_Aber))
    print("Ålesund - Dk", np.average(Åles_Danmark))
    print("Ålesund - Dk retur", np.average(Åles_Danmark_retur))
    return 0


def plot_battery_old(filebattery,filespeed):

    #manage data

    sailingcombine_each_iteration(filebattery)
    combine_each_iteration(filespeed)
    # Create the figure and axes objects
    fig, ax1 = plt.subplots()

    # Plot the speed data on the first y-axis
    ax1.plot(dates, sailing_speed, color='blue')
    ax1.set_xlabel('Date', rotation=45)
    ax1.set_ylabel('Sailing Speed (m/s)', color='blue')
    ax1.tick_params('y', colors='blue')

    # Create a second y-axis
    ax2 = ax1.twinx()

    # Plot the battery power data on the second y-axis
    ax2.plot(dates, battery_power, color='red')
    ax2.set_ylabel('Battery Power (kWh)', color='red')
    ax2.tick_params('y', colors='red')

    # Add a title
    plt.title('Speed and Battery Power')

    # Display the plot
    plt.show()

def plot_battery_power(battery_array, speed_array):

    # read data from files
    dates            = []  #Dates
    #battery_power   = read_array_As_list(file_battery,2)  # Battery power in kWh
    #speed           = read_array_As_list(file_speed,2)    # Speed in m/s

    #Calculate average for each day of simulation:
    startdate = datetime.datetime(2020,7,1)

    for i in range(len(battery_array)):
        dates.append(startdate + datetime.timedelta(hours = i))


    #Checks if length of vectors is equal
    if len(battery_array) > len(speed_array):
        battery_array   = battery_array[:len(speed_array)]
        dates           = dates[:len(speed_array)]
    if len(speed_array) > len(battery_array):
        speed_array = speed_array[:len(battery_array)]

    #for i in range(len(battery_array)):
    #    battery_array[i] = battery_array[i]
    #    speed_array[i] = speed_array[i]

    # Create the figure and axes objects
    fig, ax1 = plt.subplots(figsize = (6.5,6.5))

    # Plot the speed data on the first y-axis
    ax1.plot(dates, speed_array, color='blue')
    ax1.set_xlabel('Date',)
    ax1.set_ylabel('Sailing Speed (m/s)', color='blue')
    ax1.tick_params('y', colors='blue')

    # Create a second y-axis
    ax2 = ax1.twinx()

    # Plot the battery power data on the second y-axis
    ax2.plot(dates, battery_array, color='red')
    ax2.set_ylabel('Battery Power (kWh)', color='red')
    ax2.tick_params('y', colors='red')

    # Rotate the x-axis tick labels if needed
    plt.xticks(rotation=60)
    for tick in ax1.get_xticklabels():
        tick.set_rotation(60)

    # Add a title
    plt.title('Speed and Battery Power')



    # Display the plot
    plt.show()




def moving_average(data, window_size):
    cumsum = np.cumsum(data, dtype=float)
    cumsum[window_size:] = cumsum[window_size:] - cumsum[:-window_size]
    moving_avg = cumsum[window_size - 1:] / window_size
    return moving_avg

def plot_Two_things(array1,array2, title):

    # Read data from files
    dates = []  # Dates

    interval = 10
    # Calculate average for each day of simulation
    startdate = datetime.datetime(2020, 7, 1)
    enddate = datetime.datetime(2022, 7, 1)
    stepsize = datetime.timedelta(hours=(interval*10))

    current_date = startdate
    while current_date <= enddate:
        dates.append(current_date)
        current_date += stepsize


    #create indices for every 10 entry
    indices = [i for i in range(0, len(array1), interval)]

    print(f"Battery average is {np.average(array1)}, and speed average is {np.average(array2)}")


    #find average of each row of averages
    array1_ten = [sum(array1[i:i + interval]) / interval for i in indices]
    array2_ten = [sum(array2[i:i + interval]) / interval for i in indices]


    window_size = 20
    smooth_array1 = moving_average(array1_ten, window_size)
    smooth_array2 = moving_average(array2_ten, window_size)


    smooth_array1 = smooth_array1[:-1]
    smooth_array2 = smooth_array2[:-1]
    dates1 = dates[:-1]
    dates2 = dates[:-1]


    fig, ax1 = plt.subplots(figsize=(8, 8))
    lim = 1700
    # Plot the speed data on the first y-axis
    ax1.plot(dates1[window_size - 1:], smooth_array1, color='blue')
    ax1.set_xlabel('Date', )
    ax1.set_ylabel('Battery power (kWh)', color='blue')
    ax1.tick_params('y', colors='blue')

    # Create a second y-axis
    ax2 = ax1.twinx()



    # Plot the battery power data on the second y-axis
    ax2.plot(dates2[window_size - 1:], smooth_array2, color='red')
    ax2.set_ylabel('Sailing speed (knots)', color='red')
    ax2.tick_params('y', colors='red')

    # Rotate the x-axis tick labels if needed
    plt.xticks(rotation=60)
    for tick in ax1.get_xticklabels():
        tick.set_rotation(60)

    #find average of each value:
    avg_array1 = np.average(array1)
    avg_array2 = np.average(array2)

    # Plot the average speed as a green dashed line
    ax1.axhline(y=avg_array1, color='blue', linestyle='dashed', label='Average speed')
    ax1.yaxis.set_tick_params(length=0)  # Remove ticks on the y-axis for the average speed

    # Plot the average battery usage as a blue solid line
    ax2.axhline(y=avg_array2, color='red', linestyle='dashed', label='Average battery power')
    ax2.yaxis.set_tick_params(length=0)  # Remove ticks on the y-axis for the average battery usage

    # Combine the legends for both y-axes
    handles1, labels1 = ax1.get_legend_handles_labels()
    handles2, labels2 = ax2.get_legend_handles_labels()
    handles = handles1 + handles2
    labels = labels1 + labels2

    # Manually add labels for the speed and battery function values
    labels.extend(['Speed', 'Battery power'])
    handles.extend([plt.Line2D([], [], color='red', linestyle='-'), plt.Line2D([], [], color='blue', linestyle='-')])

    # Create the combined legend
    ax1.legend(handles, labels, loc='upper left')

    # Add a title
    plt.title(title)

    # Display the plot
    plt.show()

    return 0





# National Routes
Aales_Floro         = read_array_from_file("Output_files/Ålesund_Florø/savespeed_AalesFloro.csv")
Aales_Floro_retur   = read_array_from_file("Output_files/Ålesund_Florø_return/savespeed_AalesFloro_return.csv")
Trond_Aales         = read_array_from_file("Output_files/Trondheim_Ålesund_iteration_100/savespeed_TrondAales.csv")
Trond_Aales_retur   = read_array_from_file("Output_files/Trondheim_Ålesund_return/savespeed_TrondAales_return.csv")
Floro_Bergen        = read_array_from_file("Output_files/Florø_Bergen/savespeed_FloroBergen.csv")
Floro_Bergen_retur  = read_array_from_file("Output_files/Florø_Bergen_return/savespeed_FloroBergen_return.csv")
Bergen_Stvg         = read_array_from_file("Output_files/Bergen_Stavanger/savespeed_BrgStvg.csv")
Bergen_Stvg_retur   = read_array_from_file("Output_files/Bergen_Stavanger_return/savespeed_BrgStvg_return.csv")
Fær_Åles            = read_array_from_file("Output_files/Færøyene_Ålesund/savespeed_Færøyene_Ålesund.csv")
Fær_Åles_retur      = read_array_from_file("Output_files/Færøyene_Ålesund_return/savespeed_Færøyene_Ålesund_return.csv")
Fær_Åles_bat        = read_array_from_file("Output_files/Færøyene_Ålesund_With_Battery/Battery_use.csv")



#International Routes
Aber_Faer           = read_array_from_file("Output_files/Aberdeen_Færøyene/savespeed_Aberdeen_Færøyene.csv")
Aber_Faer_retur     = read_array_from_file("Output_files/Aberdeen_Færøyene_return/savespeed_Aberdeen_Færøyene_return.csv")
Åles_Danmark        = read_array_from_file("Output_files/Ålesund_Danmark/savespeed_Ålesund_Danmark.csv")
Åles_Danmark_retur  = read_array_from_file("Output_files/Ålesund_Danmark_return/savespeed_Ålesund_Danmark_return.csv")
Trond_Åles_1_iter     = read_array_from_file("Output_files/Trondheim_Ålesund_iteration_1/savespeed_TrondAales.csv")
Trond_Åles_10_iter    = read_array_from_file("Output_files/Trondheim_Ålesund_iteration_10/savespeed_TrondAales.csv")
Trond_Åles_100_iter   = read_array_from_file("Output_files/Trondheim_Ålesund_iteration_100/savespeed_TrondAales.csv")
Trond_Åles_1000_iter  = read_array_from_file("Output_files/Trondheim_Ålesund_iteration_1000/savespeed_TrondAales.csv")

Fær_Åles_retur      = read_array_from_file("Output_files/Færøyene_Ålesund_return/savespeed_Færøyene_Ålesund_return.csv")
Newcast_Aber        = read_array_from_file("Output_files/Newcastle_Aberdeen/SS_Newcastle_Aberdeen.csv")
Newcast_Aber_retur  = read_array_from_file("Output_files/Newcastle_Aberdeen_return/SS_Newcastle_Aberdeen_return.csv")
Amster_New          = read_array_from_file("Output_files/Amsterdam_Newcastle/savespeed_Amsterdam_Newcastle.csv")
Amster_New_retur    = read_array_from_file("Output_files/Amsterdam_Newcastle_return/savespeed_Amsterdam_Newcastle_return.csv")
Danmark_Amster      = read_array_from_file("Output_files/Danmark_Amsterdam/savespeed_Danmark_Amsterdam.csv")
Danmark_Amster_retur= read_array_from_file("Output_files/Danmark_Amsterdam_return/savespeed_Danmark_Amsterdam_return.csv")
Aber_Åles        = read_array_from_file("Output_files/Newcastle_Ålesund/SS_Newcastle_Ålesund.csv")


#special Routes
Poor_route          = read_array_from_file("Output_files/Poor_Route/savespeed_Poor_Route.csv")
Good_route          = read_array_from_file("Output_files/Good_Route/savespeed_Good_Route.csv")
Good_route_kite     = read_array_from_file("Output_files/Good_Route_With_Kite/savespeed_Good_Route.csv")
Good_route_newbuild = read_array_from_file("Output_files/Good_Route_Newbuild/savespeed_Good_Route.csv")
Fær_Åles_battery   = read_array_from_file("Output_files/Færøyene_Ålesund_With_Battery/savespeed_Good_Route.csv")
Fær_Åles_newbuild  = read_array_from_file("Output_files/Færøyene_Ålesund_Newbuild/savespeed_Færøyene_Ålesund.csv")
Fær_Åles_new_bat   = read_array_from_file("Output_files/Færøyene_Ålesund_Newbuild_Battery/savespeed_Good_Route.csv")




def combine_each_iteration(filename):

    """
    :param filename: CSV file with datetime in column 1, information in column 2
    :return: array with average of each iteration
    get first value to find batterypower used for each route iteration
    get second value to find average speed sailed for each iteration
    """

    # Open the CSV file
    with open(filename, 'r') as file:
        reader = csv.reader(file)
        next(reader)  # Skip the header row

        # Initialize variables
        current_iteration = 1
        iteration_value = []
        average_values = []
        Battery_sum   = []
        Speed_average = []
        Speed_all       = []
        # Read the first row

        time_previous   = datetime.strptime(next(reader)[1], '%Y-%m-%d %H:%M:%S.%f')
        time_zero       = time_previous
        # Iterate over the remaining rows
        for row in reader:
            speed = float(row[2])
            #find time difference
            time_now = datetime.strptime(row[1], '%Y-%m-%d %H:%M:%S.%f')
            if type(time_previous) == str:
                time_previous = datetime.strptime(time_previous, '%Y-%m-%d %H:%M:%S.%f')
            change = time_now-time_zero

            if  time_previous > time_now or change == timedelta(hours=100) or change == timedelta(hours=1000):
                # Calculate average speed for the last iteration
                route_sum = sum(iteration_value)
                Battery_sum.append(route_sum)
                Speed_average.append(np.average(iteration_value))
                Speed_all.append(iteration_value)
                # reset vector of iterations
                iteration_value = []
                time_previous = time_now
                if change == 100:
                    time_zero += timedelta(hours=100)
                if change == 1000:
                    time_zero += timedelta(hours=1000)

            else:
                # Add speed to the current iteration
                iteration_value.append(speed)
                time_previous = time_now



    # Print the average speeds for each iteration
    #print(average_value)
    print(len(Battery_sum))
    return Battery_sum,Speed_average

