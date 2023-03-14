from file_handling import  mac_windows_file_handle, write_to_file, read_array_from_file, read_datestamp_from_file
import pandas as pd

timestamp_file   = mac_windows_file_handle("Output_files/datestamp9.csv")
timestamp_data   = read_datestamp_from_file(timestamp_file)
data_file        = mac_windows_file_handle("Output_files/Newcastle_Aberdeen/SS_Newcastle_Aberdeen.csv")
data_array       = read_array_from_file(data_file)



def add_to_dataframe(datafile,timestamp_file):
    datetime_array = read_datestamp_from_file(timestamp_file)
    df = pd.DataFrame(datafile)
    df = df.drop(df.index[len(datetime_array):])
    df.insert(0, "timestamp", datetime_array)
    print(df)

    return df


a = add_to_dataframe(data_array,timestamp_file)
#testdataframe = VS
#df = add_timestamp_to_dataframe(testdataframe)
#write_to_file(df,testvs)
