import os
import pandas as pd
import functools
import operator
import numpy as np

### Constants

proc_path = os.path.expanduser("~/dados-ic/output/")
default_time_deltas = [60 * 5, 60 * 10, 60 * 15, 60 * 30, 60 * 60, 60 * 90,
                       60 * 120, 60 * 150, 60 * 180, 60*210, 60 * 240,
                       60 * 270, 60 * 300, 60 * 330, 60 * 360]
default_path_dict = {"tsi": os.path.join(proc_path, "tsi.csv"),
                     "rad": os.path.join(proc_path, "rad.csv"),
                     "goes": os.path.join(proc_path, "goes.csv"),
                     "mod": os.path.join(proc_path, "modis_terra.csv"),
                     "myd": os.path.join(proc_path, "modis_aqua.csv")}

default_output_dict = {"tsi": os.path.join(proc_path, "tsi-centered.csv"),
                     "rad": os.path.join(proc_path, "rad-centered.csv"),
                     "goes": os.path.join(proc_path, "goes-centered.csv")}

def center_data(data1, data2, time_deltas=default_time_deltas):
    """
    Averages data2 content near data1.Time points using
    given time_delta in seconds. The output Time column
    are given from data1.Time points and the N
    column are the size of data2 points used for
    the output row.
    
    Keyword Arguments:
    data1 -- Reference dataset (eg. MODIS data)
    data2 -- Dataset to average (eg. TSI data)
    time_delta -- Seconds apart from data1 points to average    
    """

    time_deltas = sorted(time_deltas, reverse=True)
    data3 = pd.DataFrame()
    data3_times = []
    data3_N = []
    data3_time_delta = []
    
    total = len(data1)
    i = 0
    
    for item_iter in data1.iterrows():
        filtered_data2 = data2.copy()
        
        for time_delta in time_deltas:
            item = item_iter[1]
            item_time = item.Time
            min_time = item_time - pd.Timedelta(seconds=time_delta)
            max_time = item_time + pd.Timedelta(seconds=time_delta)

            select_filter = filtered_data2.Time > min_time
            select_filter &= filtered_data2.Time < max_time

            filtered_data2 = filtered_data2[select_filter]
            N = len(filtered_data2)
            avg = filtered_data2.mean()

            data3_times.append(item.Time)
            data3_N.append(N)
            data3_time_delta.append(time_delta)
            data3 = data3.append(avg, ignore_index=True)
            
        i += 1
        print("\r%s / %s points processed" % (i, total), end="")
    data3["Timedelta"] = data3_time_delta
    data3["Time"] = data3_times
    data3["N"] = data3_N
    
    return data3

def goes_center_data(data1, goes_data, time_deltas=default_time_deltas):    
    """
    Similiar to center_data, but this function takes into account
    that GOES have spatial deltas.
    """
    delta_list = np.unique(goes_data.Delta)    
    output = []    
    for delta in delta_list:
        filter_indices = (goes_data.Delta == delta)
        data2 = goes_data[filter_indices]
        data3 = center_data(data1, data2, time_deltas)
        output.append(data3)
    
    dataframe = pd.concat(output)
    
    return dataframe.reset_index()

def process_all(modis_data, path_dict=default_path_dict,
                output_dict=default_output_dict,
                time_deltas=default_time_deltas):
    """
    
    """
    data1 = modis_data
    
    deltas = np.unique(data1.Delta)
    gambi_delta = deltas[np.argmax(deltas)]
    data1 = data1[data1.Delta == gambi_delta]

    for k in output_dict.keys():
        print("Processing %s" % (k))
        data2 = pd.read_csv(path_dict[k], na_values=["--"])
        data2.Time = pd.to_datetime(data2.Time)
        if k == "goes":
            data3 = goes_center_data(data1, data2, time_deltas)
        else:
            data3 = center_data(data1, data2)
        
        data3 = data3[data3.N > 0]
        data3 = data3[pd.notnull(data3.Time)]
        data3.to_csv(output_dict[k], index=False)           
    
    
def process_all_lazy():
    
    
    print("Centering on GOES Data")
    modis_path = default_path_dict["goes"]
    modis_data = pd.read_csv(modis_path, na_values=["--"])
    modis_data.Time = pd.to_datetime(modis_data.Time)
    output_dict = default_output_dict.copy()
    del output_dict["goes"]
    
    for k in output_dict.keys():
        output_dict[k] = output_dict[k].replace("centered", "centered-goes")
    process_all(modis_data, output_dict=output_dict)
    
    
    print("Centering on MODIS Terra data")
    modis_path = default_path_dict["mod"]
    modis_data = pd.read_csv(modis_path, na_values=["--"])
    modis_data.Time = pd.to_datetime(modis_data.Time)
    output_dict = default_output_dict.copy()
    
    for k in output_dict.keys():
        output_dict[k] = output_dict[k].replace("centered", "centered-mod")
    process_all(modis_data, output_dict=output_dict)
    
        
    print("Centering on MODIS Aqua Data")
    modis_path = default_path_dict["myd"]
    modis_data = pd.read_csv(modis_path, na_values=["--"])
    modis_data.Time = pd.to_datetime(modis_data.Time)
    output_dict = default_output_dict.copy()
    
    for k in output_dict.keys():
        output_dict[k] = output_dict[k].replace("centered", "centered-myd")
    process_all(modis_data, output_dict=output_dict)

    
####### Average modis times
def average_modis_path(modis_path, minutes=5):
    #modis_path = default_path_dict["mod"]
    modis_data = pd.read_csv(modis_path, na_values=["--"])
    modis_data.Time = pd.to_datetime(modis_data.Time)
    data = modis_data


    for it in data.iterrows():
        item = it[1]
        item_time = item.Time
        min_time = item_time - pd.Timedelta(minutes=minutes)
        max_time = item_time + pd.Timedelta(minutes=minutes)

        select_filter = data.Time > min_time
        select_filter &= data.Time < max_time

        avg_time = np.mean(data[select_filter].Time.values.astype(np.int64))
        siz = len(data[select_filter])
        if siz > 1:
            data.loc[select_filter, "Time"] = avg_time
            print("\r%s     " % (len(data) - it[0]), end="")


    data = data[pd.notnull(data.Time)]
    data.to_csv(modis_path, index=False)
    return data

def main():
    average_modis_path(default_path_dict["mod"])
    average_modis_path(default_path_dict["myd"])
    process_all_lazy()
