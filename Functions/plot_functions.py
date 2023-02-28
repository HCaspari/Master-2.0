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
#plot average power over the 4 year period
def plot_avg_power(tid,power):
    plt.plot(tid,power)
    plt.title(" avg Poweroutput from flettnerrotor")
    plt.xlabel("triptime of year")
    plt.ylabel("kN produced by flettner on average")
    plt.show()
#plot weekly and daily power over 4 years
def plot_weekly_and_daily_avg_power(power):
    #create weekly average power vector
    monthly_power_average = []
    weekly_power_average = []
    daily_power_average = []
    for i in range(0,len(power)-4,4):
        daily_power_sum = 0
        for j in range(4):
            daily_power_sum += power[i+j]
        daily_power_average.append(daily_power_sum/4)

    for i in range(0,len(daily_power_average)-7,7):
        weekpowersum = 0
        for j in range(7):
            weekpowersum += daily_power_average[i+j]
        weekly_power_average.append(weekpowersum/7)

    for i in range(0,len(weekly_power_average)-4,4):
        monthly_power_sum = 0
        for j in range(4):
            monthly_power_sum += weekly_power_average[i+j]
        monthly_power_average.append(monthly_power_sum/4)
    plt.plot(daily_power_average)
    plt.title(" avg daily poweroutput from flettnerrotor")
    plt.xlabel("day of dataset")
    plt.ylabel("kW produced by flettner on average that day")
    plt.show()
    plt.plot(weekly_power_average)
    plt.title(" avg Poweroutput from flettnerrotor per week")
    plt.xlabel("week of dataset")
    plt.ylabel("kW produced by flettner on average that week")
    plt.show()
    plt.plot(monthly_power_average)
    plt.title(" avg Poweroutput from flettnerrotor per month")
    plt.xlabel("month of dataset")
    plt.ylabel("kW produced by flettner on average that month")
    plt.show()
#plots resistance of speed found by admirality formula
def plot_resistance():
    speeds = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16]
    resistance = [0,6.05,24.18,54.41,96.74,151.15,217.65,299.25,386.98,489.72,604.6,731.56,870.62,1021.77,1185.01,1360.34,1547.77]
    plt.plot(speeds,resistance)
    plt.title("Resistance in kN")
    plt.xlabel("Vessel speed")
    plt.ylabel("kN")
    plt.show()
    return 0

def plot_power_2(title, y_axis, x_label, y_label):
    speed = [1,3,6,9,12]
    plt.plot(speed,y_axis)
    plt.title(title)
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.show()
    return 0

def plot_percent(title, y_axis, x_label, y_label):
    speed = [3,6,9,12]
    plt.plot(speed,y_axis)
    plt.title(title)
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.show()
    return 0


#This function will output a histogram that show the speed distribution for a given route
def historgram(filename, title):
    columns = ['number', 'speed']
    df1 = pd.read_csv(filename,header = None, names=columns)
    plt.hist(df1['speed'], edgecolor='black', range=[0,15], bins=15)
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

#historgram(mac_windows_file_handle('Output_files/Trondheim_Aalesund_reise', Trond_Aalesund))
#historgram(mac_windows_file_handle('Output_files/Aalesund_Floro_reise', Aalesund_Floro))
#historgram(mac_windows_file_handle('Output_files/Floro_Bergen_reise', Floro_Bergen))
historgram(mac_windows_file_handle('Output_files/savespeed_TrondAales.csv', Trond_Aalesund))
