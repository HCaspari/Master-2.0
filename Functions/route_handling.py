import math
import pandas as pd
import csv
import os
import numpy as np

import folium as folium
from geopy.distance import geodesic
import platform



#write data found to file called filename_func
def write_to_file(data,filename_func):
    pd.DataFrame(data).to_csv(filename_func)
    return 0

def read_position_vect_from_file(filename_func):
    readdata = pd.read_csv(filename_func)#, usecols=["1","0"])
    readdata = readdata[readdata != 0]

    position_array =  []
    for i in range(len(readdata)):
        position_array.append((round(readdata.loc[i].iat[0], 3),round(readdata.loc[i].iat[1], 3)))
    return position_array

def mac_windows_file_handle(path):
    if platform.system() == "Windows":
        return "../"+path
    elif platform.system() == "Darwin":
        return path
    else:
        return 1

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
filename_AIS =mac_windows_file_handle("env/input_files/ais_data_v4.csv")
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
Trondheim = 63.437686821303096, 10.402184694640052
Aalesund_location  = 62.93245830958637, 6.3481997169859055
Trond_Aalesund      = mac_windows_file_handle("Route_data/Trondheim_Ålesund_Route/Route_Trond_Aales.csv")
Aalesund_Floro      = mac_windows_file_handle("Route_data/Ålesund_Florø_Route/Route_Ålesund_Florø.csv")
Floro_Bergen        = mac_windows_file_handle("Route_data/Florø_Bergen_Route/Route_Floro_Bergen.csv")
Bergen_Stavanger    = mac_windows_file_handle("Route_data/Bergen_Stavanger_Route/Route_Bergen_Stavanger.csv")

Route_Trond_Aal     = read_position_vect_from_file(Trond_Aalesund)
Route_Aal_Floro     = read_position_vect_from_file(Aalesund_Floro)
Route_Floro_Bergen  = read_position_vect_from_file(Floro_Bergen)
Route_Bergen_Stvg   = read_position_vect_from_file(Bergen_Stavanger)

#create vector of positions over trip, and removes duplicate positions, saves positions to file
def vector_of_positions(lats,lons):
    position_array_func = []
    for i in range(len(lats)):
        position_array_func.append((round(lats[i],3),round(lons[i],3)))
    #remove duplicate positions
    #print(len(position_array_func))
    j= 0
    while j < len(position_array_func)-1:
        lat1 = position_array_func[j][0]
        lat2 = position_array_func[j+1][0]
        lon1 = position_array_func[j][1]
        lon2 = position_array_func[j+1][1]
        if lat1 == lat2 and lon1 == lon2:
            del position_array_func[j]
        else:
            j += 1
    #print(len(position_array_func))
    position_array_file  = mac_windows_file_handle("env/position_array")
    write_to_file(position_array_func,position_array_file)
    return position_array_func

#calculates direction between two points, (20.323,23.243),(34.235, 43.345)
def calc_vessel_heading(pointA, pointB):
    """

    :param pointA: Coordinates of a
    :param pointB: coordinates of b
    :return: value in degreees
    """
    deg2rad = math.pi / 180
    rad2deg = 180 / math.pi
    latA = pointA[0] * deg2rad
    lonA = pointA[1] * deg2rad

    latB = pointB[0] * deg2rad
    lonB = pointB[1] * deg2rad

    delta_ratio = math.log(math.tan(latB/ 2 + math.pi / 4) / math.tan(latA/ 2 + math.pi / 4))
    delta_lon = abs(lonA - lonB)

    delta_lon %= math.pi
    vessel_heading = math.atan2(delta_lon, delta_ratio) * rad2deg
    return vessel_heading #IN DEGREES

def calc_vessel_heading_2(pointA, pointB):

    """
    Created by https://gist.github.com/jeromer/2005586
    Calculates the bearing between two points.
    The formulae used is the following:
        θ = atan2(sin(Δlong).cos(lat2) cos(lat1).sin(lat2) − sin(lat1).cos(lat2).cos(Δlong))
    :Parameters:
      - `pointA: The tuple representing the latitude/longitude for the
        first point. Latitude and longitude must be in decimal degrees
      - `pointB: The tuple representing the latitude/longitude for the
        second point. Latitude and longitude must be in decimal degrees
    :Returns:
      The bearing in degrees
    :Returns Type:
      float
    """
    pointA = tuple(pointA)
    pointB = tuple(pointB)
    if (type(pointA) != tuple) or (type(pointB) != tuple):
        raise TypeError("Only tuples are supported as arguments")

    lat1 = math.radians(pointA[0])
    lat2 = math.radians(pointB[0])

    diffLong = math.radians(pointB[1] - pointA[1])

    x = math.sin(diffLong) * math.cos(lat2)
    y = math.cos(lat1) * math.sin(lat2) - (math.sin(lat1) * math.cos(lat2) * math.cos(diffLong))

    initial_bearing = math.atan2(x, y)

    # Now we have the initial bearing but math.atan2 return values
    # from -180° to + 180° which is not what we want for a compass bearing
    # The solution is to normalize the initial bearing as shown below
    initial_bearing = math.degrees(initial_bearing)
    compass_bearing = (initial_bearing + 360) % 360

    return compass_bearing

def generate_intermediate_points(start_point, end_point, num_points):
    """
    Generates intermediate points between two given coordinate points using linear interpolation.

    Parameters:
        start_point (tuple): The starting coordinate point, in the form (latitude, longitude)
        end_point (tuple): The ending coordinate point, in the form (latitude, longitude)
        num_points (int): The number of intermediate points to generate

    Returns:
        list: A list of intermediate points, in the form [(latitude, longitude), ...]
    """
    intermediate_points = []
    latitude_delta = (end_point[0] - start_point[0]) / (num_points + 1)
    longitude_delta = (end_point[1] - start_point[1]) / (num_points + 1)

    for i in range(1, num_points + 1):
        latitude = start_point[0] + i * latitude_delta
        longitude = start_point[1] + i * longitude_delta
        intermediate_points.append((round(latitude,3), round(longitude,3)))

    return intermediate_points
#test_route = [(0,0),(10,10),(20,20)]

def equal_dist():
    start_point = Route_Trond_Aal[0]
    end_point = Route_Trond_Aal[-1]
    distance = geodesic(start_point, end_point).nautical
    num_points = 10
    interval = distance / (num_points - 1)
    points = []
    for i in range(num_points):
        distance_from_start = i * interval
        fraction = distance_from_start / distance
        lat = start_point[0] + fraction * (end_point[0] - start_point[0])
        lon = start_point[1] + fraction * (end_point[1] - start_point[1])
        points.append((lat, lon))

    maptest = folium.Map(location=start_point, zoom_start=10)
    folium.PolyLine(points, color="blue", weight=2.5, opacity=1).add_to(maptest)
    maptest.show_in_browser()
    return 0
#equal_dist()


def generate_intricate_route(route,points):
    """
    :param route: vector of positions over route
    :param points: number of points between separate positions
    :return: detailed route as vector of coordinates
    """
    newroute = [route[0][:]]
    for position in range(len(route)-1):
        start_position  = route[position]
        end_position    = route[position+1]
        intermediate_points = generate_intermediate_points(start_position,end_position, points)
        newroute.extend(intermediate_points)
        newroute.append((round(end_position[0],3),round(end_position[1],3)))
    return newroute

#generate_intricate_route(Route_Trond_Aal,4)

def testthingy():
    print(f"Test Short route start {test_route[0]}")
    print(f"Test Short route end {test_route[-1]}")
    test_route_long = generate_intricate_route(test_route, 3)#
    print(f"Test Long route vector: {test_route_long}")
    print(f"Test lenght of short route: {len(test_route)}")
    print(f"Test lenght of long route: {len(test_route_long)}")

    print(f"Short route start {Route_Trond_Aal[0]}")
    print(f"Short route end {Route_Trond_Aal[-1]}")
    Route_Trond_Aal_Long = generate_intricate_route(Route_Trond_Aal, 3)#
    print(f"Long route vector: {Route_Trond_Aal_Long}")
    print(f"lenght of short route: {len(Route_Trond_Aal)}")
    print(f"lenght of long route: {len(Route_Trond_Aal_Long)}")
    return 0
#testthingy()
#Route_Trond_Aal_Long = generate_intricate_route(Route_Trond_Aal, 3)#
def createmap(trip_vector):
    """
    :param trip_vector: A vector of coordinates over route
    :return: a map that shows given route
    """

    foliummap = folium.Map(location=[59.6131, 5.0301], tiles="Stamen Terrain", zoom_start=5)
    #foliummap.show_in_browser()

    route = trip_vector
    for coordinates in route:
        foliummap.add_child(
                folium.Marker(
                    location=coordinates,
                    icon=folium.Icon(color="%s" % "blue"),
                )
        )

    foliummap.show_in_browser()
    return 0
#createmap(Route_Trond_Aal)
#createmap(Route_Trond_Aal_Long)

#createmap(Route_Trond_Aal)
#Specific_route_TA = generate_intricate_route(Route_Trond_Aal,15)
#createmap(Specific_route_TA)

def distance_calc(route):
    totaldist = 0
    for i in range(len(route)-1):
        start_coord = route[i]
        end_coord = route[i + 1]
        totaldist += geodesic(start_coord,end_coord).kilometers

    print(f"Total distance = {totaldist}")
#distance_between_points(Route_Trond_Aal)
#distance_between_points(Route_Aal_Floro)
#distance_between_points(Route_Floro_Bergen)
##distance_between_points(Route_Bergen_Stvg)
def distance_between_points(route):
    dist_vect = []
    for i in range(len(route)-1):
        start_coord = route[i]
        end_coord = route[i + 1]
        dist_vect.append(geodesic(start_coord,end_coord).kilometers)
    print(dist_vect)

#distance_calc(Route_Trond_Aal)
#distance_calc(Route_Aal_Floro)
#distance_calc(Route_Floro_Bergen)
#distance_calc(Route_Bergen_Stvg)

intricate_Trond_aal = generate_intricate_route(Route_Trond_Aal,15)
intricate_Aal_Floro = generate_intricate_route(Route_Aal_Floro,15)
intricate_Floro_Brg = generate_intricate_route(Route_Floro_Bergen,15)
intricate_Brg_Stvg  = generate_intricate_route(Route_Bergen_Stvg,10)
route_Trond_Aals_intricate  = mac_windows_file_handle("Route_data/route_Trond_Aales_Intricate")
route_Aals_Floro_intricate  = mac_windows_file_handle("Route_data/route_Aales_Floro_Intricate")
route_Floro_Brg_intricate   = mac_windows_file_handle("Route_data/route_Floro_Brg_Intricate")
route_Brg_Stvg_intricate    = mac_windows_file_handle("Route_data/route_Brg_Stv_Intricate")
route_Stvg_Dk_intricate     = mac_windows_file_handle("Route_data/Danmark_Stavanger_Route/Route_Danmark_Stavanger.csv")

#write_to_file(intricate_Trond_aal,route_Trond_Aals_intricate)
#write_to_file(intricate_Aal_Floro,route_Aals_Floro_intricate)
#write_to_file(intricate_Floro_Brg,route_Floro_Brg_intricate)
#write_to_file(intricate_Brg_Stvg,route_Brg_Stvg_intricate)

#createmap(Route_Bergen_Stvg)



#####

#Route creation : Aalesund -- Færøyene -- Aberdeen -- Newcastle -- Amsterdam -- Esbjerg (Danmark) -- Aalesund
#Route Coordinates:
def auto_create_Route(start_point_coordinates, end_point_coordinates, midpoint1 = (), midpoint2 = ()):
    Route = [0]
    #if midpoints do not exist, create direct route
    if len(midpoint1) == 0:
        route = [start_point_coordinates, end_point_coordinates]
        distance = int(np.floor(geodesic(start_point_coordinates, end_point_coordinates).kilometers // 7))
        Route = generate_intricate_route(route, distance)
        #createmap(Route)

    # if midpoint 1 exists, and midpoint2 does not, create indirect route
    elif len(midpoint1) == 2 and len(midpoint2) != 2:

        route_step_1    = [start_point_coordinates, midpoint1] #first segment of 1 step route
        route_step_2    = [midpoint1,end_point_coordinates]    #second segment of 1 step route

        distance_step_1 = geodesic(start_point_coordinates, midpoint1).kilometers #distance of first step
        distance_step_2 = geodesic(midpoint1, end_point_coordinates).kilometers   #distance of second step

        distance_tot    = int(np.floor((distance_step_1+distance_step_2)//7))          #distance div by 7

        Route1 = generate_intricate_route(route_step_1, distance_tot)             #route part 1
        Route2 = generate_intricate_route(route_step_2, distance_tot)             #route part 2

        Route  = Route1 + Route2
        #createmap(Route)

    #if midpoint 1 and 2 exists, create indirect, 3 point route
    elif len(midpoint1) == 2 and len(midpoint2) == 2:
        route_step_1 = [start_point_coordinates, midpoint1]  # first segment of 1 step route
        route_step_2 = [midpoint1, midpoint2]  # second segment of 1 step route
        route_step_3 = [midpoint2, end_point_coordinates]  # second segment of 1 step route

        distance_step_1 = geodesic(start_point_coordinates, midpoint1).kilometers  # distance of first step
        distance_step_2 = geodesic(midpoint1, midpoint2).kilometers  # distance of second step
        distance_step_3 = geodesic(midpoint2, end_point_coordinates).kilometers  # distance of second step

        distance_tot = int(np.floor((distance_step_1 + distance_step_2 + distance_step_3) // 7))  # distance div by 7

        Route1 = generate_intricate_route(route_step_1, distance_tot)  # route part 1
        Route2 = generate_intricate_route(route_step_2, distance_tot)  # route part 2
        Route3 = generate_intricate_route(route_step_3, distance_tot)  # route part 3

        Route = Route1 + Route2 + Route3
        #createmap(Route)

    #Creates a tripvector of points from start_point to end_point with stepdistance of max 7 km.
    # (meaning weather will always be recalculated when entering new 0.125 degree lat/lon which is weather granularity
    #print("no route found")
    return Route

#Route point coordinates:

Aalesund        = (62.48342017643765, 5.922191001412515)
Stavanger       = (58.92502271580142, 5.580428183923716)
Bergen          = (60.134509001387705, 4.958929708820197)
Mid_Danmark_stvg= (58.534289780314545, 4.595293578632795)
Færøyene        = (62.01263234289485, -6.774748315943504)
Midpoint_Fa_Ab  = (58.80320072664799, -1.1708712520663145)
Aberdeen        = (57.15470822505185, -2.1185733842766115)
Newcastle       = (54.97876964024249, -1.3553315768007936)
Amsterdam       = (52.47169283545967, 4.537125962881187)
Danmark         = (55.47341900177149, 8.301921187910168)
Midpoint_1_Dan_Aal  = (60.44513450388027, 3.418535867789226)
Midpoint_2_Dan_Aal  = (62.324854695720006, 5.0124722561975235)
Florø           = (61.476094547417446, 4.716076439675587)

#createmap([Færøyene,Stavanger])
#createmap([Færøyene,(61.242101, -0.405119),Stavanger])
#Fær_Stvg = auto_create_Route(Færøyene,Stavanger,(61.242101, -0.405119))
#createmap(Fær_Stvg)
#to create route:
#1. Run create_folder_with_files(Route_name,directory_path)
#   Directory path is usually = "C:\\Users\\haako\\Documents\\10. Semester\\MASTAH\\Master-2.0\\Output_files"
#3. See if intermediate points are needed
#4. run route_start_end = auto_create_route(Startcoordinate, endcoordinate, optional midcoordinate 1, optional midcoordinate 2)
#5. run write_to_file(route_start_end), this writes route coordinates to file
#6. run main(routename, simulation count)



#route_stv_Dk        = auto_create_Route(Stavanger,Danmark,Mid_Danmark_stvg)
#filename_Stv_Dk     = mac_windows_file_handle("Route_data/Stavanger_Danmark_Route/Route_Stavanger_Danmark.csv")
#write_to_file(route_stv_Dk,filename_Stv_Dk)
#route_Aber_Faer,route_Faer_Aales,route_Aales_Dk,route_Dk_Amst,route_DK_Stav,route_New_Aber,route_Amst_New = create_routes()
#route_New_Faer = auto_create_Route(Newcastle,Færøyene,Midpoint_Fa_Ab)
#print(route_stv_Dk)

#createmap(route_stv_Dk)
#createmap(route_Faer_Aales)
#createmap(route_Aales_Dk)
#createmap(route_Dk_Amst)
#createmap(route_DK_Stav)
#createmap(route_New_Aber)
#createmap(route_Amst_New)
#createmap(route_New_Faer)

def create_international_routes():

    #Route start and finish
    A_F = [Aalesund,Færøyene]
    F_M = [Færøyene,Midpoint_Fa_Ab]
    M_A = [Midpoint_Fa_Ab,Aberdeen]
    A_N = [Aberdeen,Newcastle]
    N_A = [Newcastle,Amsterdam]
    A_E = [Amsterdam, Danmark]
    E_M1    = [Danmark, Midpoint_1_Dan_Aal]
    M1_M2   = [Midpoint_1_Dan_Aal, Midpoint_2_Dan_Aal]
    M2_A    = [Midpoint_2_Dan_Aal, Aalesund]

    #Distances of routes, divded by 7, to find amount of 7 km segments route consists of

    A_F_dist = geodesic(Aalesund,Færøyene).kilometers/7
    F_M_dist = geodesic(Færøyene, Midpoint_Fa_Ab).kilometers/7
    M_A_dist = geodesic(Midpoint_Fa_Ab, Aberdeen).kilometers/7
    A_N_dist = geodesic(Aberdeen, Newcastle).kilometers/7
    N_A_dist = geodesic(Newcastle, Amsterdam).kilometers/7
    A_E_dist = geodesic(Amsterdam, Danmark).kilometers / 7
    E_M1_dist   = geodesic(Danmark, Midpoint_1_Dan_Aal).kilometers / 7
    M1_M2_dist  = geodesic(Midpoint_1_Dan_Aal, Midpoint_2_Dan_Aal).kilometers / 7
    M2_A_dist   = geodesic(Midpoint_2_Dan_Aal, Aalesund).kilometers / 7


    #generating routes with this amount of spits:

    A_F_route   = generate_intricate_route(A_F,int(np.floor(A_F_dist)))
    F_M_route   = generate_intricate_route(F_M,int(np.floor(F_M_dist)))
    M_A_route   = generate_intricate_route(M_A,int(np.floor(M_A_dist)))
    A_N_route   = generate_intricate_route(A_N,int(np.floor(A_N_dist)))
    N_A_route   = generate_intricate_route(N_A,int(np.floor(N_A_dist)))
    A_E_route   = generate_intricate_route(A_E,int(np.floor(A_E_dist)))
    E_M1_route  = generate_intricate_route(E_M1,int(np.floor(E_M1_dist)))
    M1_M2_route = generate_intricate_route(M1_M2,int(np.floor(M1_M2_dist)))
    M2_A_route  = generate_intricate_route(M2_A,int(np.floor(M2_A_dist)))

    return A_F_route,F_M_route,M_A_route,A_N_route,N_A_route,A_E_route, E_M1_route, M1_M2_route, M2_A_route

def create_folder_with_files(Route_Name, directory_path):

    folder_name = str(Route_Name)
    folder_path = os.path.join(directory_path, folder_name)
    os.makedirs(folder_path, exist_ok=True)
    for i in range(1):
        savespeed = "savespeed_" + Route_Name + ".csv"
        saveTWD = "saveTWD_" + Route_Name + ".csv"
        saveTWS = "saveTWS_" + Route_Name + ".csv"
        path_speed = os.path.join(folder_path, savespeed)
        path_TWD = os.path.join(folder_path, saveTWD)
        path_TWS = os.path.join(folder_path, saveTWS)
        open(path_speed, 'a').close()
        open(path_TWD, 'a').close()
        open(path_TWS, 'a').close()
def create_route_files(Route_Name, directory_path):

    folder_path = os.path.join(directory_path, Route_Name)
    os.makedirs(folder_path, exist_ok=True)
    savespeed = "Route_" + Route_Name + ".csv"
    path_speed = os.path.join(folder_path, savespeed)
    open(path_speed, 'a').close()

file_path_routes = "C:\\Users\\haako\\Documents\\10. Semester\\MASTAH\\Master-2.0\\Route_data"
route_directory = "C:\\Users\\haako\Documents\\10. Semester\\MASTAH\\Master-2.0\\Output_files"

def invert_csv(csv_File):

    with open(csv_File, 'r') as textfile:
        for row in reversed(list(csv.reader(textfile))):
            print(', '.join(row))


#to create route:
#1. Run create_folder_with_files(Route_name,directory_path)
#   Directory path is usually = "C:\\Users\\haako\\Documents\\10. Semester\\MASTAH\\Master-2.0\\Output_files"
#3. See if intermediate points are needed
#4. run route_start_end = auto_create_route(Startcoordinate, endcoordinate, optional midcoordinate 1, optional midcoordinate 2)
#5. run write_to_file(route_start_end), this writes route coordinates to file
#6. run main(routename, simulation count)

def auto_create_Route_2(start_point_coordinates, end_point_coordinates, midpoint1 = (), midpoint2 = (),midpoint3 = (),midpoint4 = (),midpoint5 = ()):

    route_step_1 = [start_point_coordinates, midpoint1]  # first segment of 1 step route
    route_step_2 = [midpoint1, midpoint2]  # second segment of 1 step route
    route_step_3 = [midpoint2, midpoint3]  # second segment of 1 step route
    route_step_4 = [midpoint3, midpoint4]  # second segment of 1 step route
    route_step_5 = [midpoint4, end_point_coordinates]  # second segment of 1 step route

    distance_step_1 = geodesic(start_point_coordinates, midpoint1).kilometers  # distance of first step
    distance_step_2 = geodesic(midpoint1, midpoint2).kilometers  # distance of second step
    distance_step_3 = geodesic(midpoint2, midpoint3).kilometers  # distance of second step
    distance_step_4 = geodesic(midpoint3, midpoint4).kilometers  # distance of second step
    distance_step_5 = geodesic(midpoint4, end_point_coordinates).kilometers  # distance of second step

    distance_tot = int(np.floor((distance_step_1 + distance_step_2 + distance_step_3 + distance_step_4 +distance_step_5) // 7))  # distance div by 7

    Route1 = generate_intricate_route(route_step_1, distance_tot)  # route part 1
    Route2 = generate_intricate_route(route_step_2, distance_tot)  # route part 2
    Route3 = generate_intricate_route(route_step_3, distance_tot)  # route part 3
    Route4 = generate_intricate_route(route_step_4, distance_tot)  # route part 4
    Route5 = generate_intricate_route(route_step_5, distance_tot)  # route part 5



    Route = Route1 + Route2 + Route3 + Route4 + Route5
    createmap(Route)
    return Route

def auto_create_Route_3(start_point_coordinates, end_point_coordinates, midpoint1 = (), midpoint2 = (),midpoint3 = (),midpoint4 = (),midpoint5 = ()):

    route_step_1 = [start_point_coordinates, midpoint1]  # first segment of 1 step route
    route_step_2 = [midpoint1, midpoint2]  # second segment of 2. step route
    route_step_3 = [midpoint2, midpoint3]  # second segment of 3. step route
    route_step_4 = [midpoint3, midpoint4]  # second segment of 4. step route
    route_step_5 = [midpoint4, midpoint5]  # second segment of 5. step route
    route_step_6 = [midpoint5, end_point_coordinates]  # second segment of final step route

    distance_step_1 = geodesic(start_point_coordinates, midpoint1).kilometers  # distance of 1. step
    distance_step_2 = geodesic(midpoint1, midpoint2).kilometers  # distance of 2. step
    distance_step_3 = geodesic(midpoint2, midpoint3).kilometers  # distance of 3. step
    distance_step_4 = geodesic(midpoint3, midpoint4).kilometers  # distance of 4. step
    distance_step_5 = geodesic(midpoint4, midpoint5).kilometers  # distance of 5. step
    distance_step_6 = geodesic(midpoint5, end_point_coordinates).kilometers  # distance of final step

    distance_tot = int(np.floor((distance_step_1 + distance_step_2 + distance_step_3 + distance_step_4 + distance_step_5 + distance_step_6) // 7))  # distance div by 7

    Route1 = generate_intricate_route(route_step_1, distance_tot)  # route part 1
    Route2 = generate_intricate_route(route_step_2, distance_tot)  # route part 2
    Route3 = generate_intricate_route(route_step_3, distance_tot)  # route part 3
    Route4 = generate_intricate_route(route_step_4, distance_tot)  # route part 4
    Route5 = generate_intricate_route(route_step_5, distance_tot)  # route part 5
    Route6 = generate_intricate_route(route_step_6, distance_tot)  # route part 6

    Route = Route1 + Route2 + Route3 + Route4 + Route5 + Route6

    #Creates a tripvector of points from start_point to end_point with stepdistance of max 7 km.
    # (meaning weather will always be recalculated when entering new 0.125 degree lat/lon which is weather granularity
    #print("no route found")
    #createmap(Route)
    return Route

#route_poor = [Danmark,Mid_Danmark_stvg,Aalesund,Aberdeen,Newcastle,Amsterdam]
#route_good = [Danmark,Mid_Danmark_stvg,Florø,Færøyene,Midpoint_Fa_Ab,Newcastle,Amsterdam]

#Route_good_complete = auto_create_Route_3(Danmark, Amsterdam,Mid_Danmark_stvg, Midpoint_2_Dan_Aal, Aalesund,Aberdeen,Newcastle)
#Route_bad_complete = auto_create_Route_3(Danmark, Amsterdam,Mid_Danmark_stvg,   Florø,              Færøyene,Midpoint_Fa_Ab,Newcastle)
#create_route_files("Good_Route",file_path_routes)
#create_route_files("Poor_Route",file_path_routes)
#create_folder_with_files("Poor_Route",route_directory)

#write_to_file(Route_good_complete,"../Route_data/Good_Route/Route_Good_Route.csv")
#write_to_file(Route_bad_complete,"../Route_data/Poor_Route/Route_Poor_Route.csv")

def read_route2(csv):
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

#Poor_route = read_route2("../Route_data/Poor_Route/Poor_route.csv")
#createmap(Poor_route)

Good_route = read_route2("../Route_data/Good_Route/Good_route.csv")
createmap(Good_route)