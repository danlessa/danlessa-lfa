"""
Functions for processing MODIS L6 data
Author: Danilo Lessa Bernardineli (danlessa@if.usp.br)
"""
import numpy as np
from pyhdf.SD import *
import time
import os
import pandas as pd
import sys
import common
import datetime as dt
import pyhdf

lat_mao = common.lat_mao
lon_mao = common.lon_mao

# MODIS base time
initial_time = time.mktime(time.strptime("00:00 01/01/1993", "%H:%M %d/%m/%Y"))

# Default settings for data filtering and post-processing output
default_delta = 0.5
default_variables = ["Cloud_Fraction", "Solar_Zenith"]

default_data_folder = os.path.expanduser("~/dados-ic/MODIS_Terra")
default_output_path = os.path.expanduser("~/dados-ic/output/modis_terra.csv")


def get_modis_filepaths(folder_path):
    """
    Returns a list of tuple MODIS filepaths. An element can contain two
    elements if there is more than a file associated with the same instant.
    """

    # Getting list of files
    files = os.listdir(folder_path)

    files = [f[10:-22] for f in files]  # clean useless filename text
    files.sort()
    files
    date_hour = [[f[:-5], f[-4:]] for f in files]  # getting the date and time
    date_hour = pd.DataFrame(date_hour, columns=["date", "time"], dtype=int)

    # Deciding what files are 'duplicates'
    duplicates = {}
    duplicates["date"] = []
    duplicates["time1"] = []
    duplicates["time2"] = []

    unique_dates = np.unique(date_hour.date)

    for unique_date in unique_dates:

        # Get only same-day files
        inds1 = date_hour.date == unique_date
        date_hour_day = date_hour[inds1]

        unique_times = date_hour_day.time

        for time in unique_times:
            delta = np.abs(unique_times - time)

            # We assume that "duplicate" files are within 30
            inds = delta < 30
            same_times = unique_times[inds]

            if len(same_times) > 1:
                a = []
                for t in same_times:
                    a.append(t)
                duplicates["date"].append(unique_date)
                duplicates["time1"].append(a[0])
                duplicates["time2"].append(a[1])

    dup = pd.DataFrame(duplicates)
    dup = dup.drop_duplicates()

    # Creating an nice list of duplicate filenames
    files = os.listdir(folder_path)
    dup_files = []
    duplicates = []

    for duplicate in dup.iterrows():
        duplicate = duplicate[1]
        pseudo_path1 = "%s.%04d" % (duplicate.date, duplicate.time1)
        pseudo_path2 = "%s.%04d" % (duplicate.date, duplicate.time2)

        file1 = None
        file2 = None
        for f in files:
            if f.find(pseudo_path1) > -1:
                file1 = f
            if f.find(pseudo_path2) > -1:
                file2 = f
        dup_files.append((file1, file2))

    dup_filepaths = [(folder_path + "/" + f[0], folder_path + "/" + f[1])
                     for f in dup_files]

    filepaths = [folder_path + "/" + f for f in files]
    un_filepaths = []
    for file in filepaths:
        if file not in np.array(dup_filepaths).flatten():
            un_filepaths.append(file)

    output_filepaths = []

    for filepath in un_filepaths:
        output_filepaths.append([filepath])

    for filepath in dup_filepaths:
        output_filepaths.append(filepath)

    return output_filepaths


def load_dataset(dataset, variables=default_variables):
    """
    Loads MODIS dataset and retrieve pertinent variables.
    Output is an dictionary.

    Keyword arguments:
    dataset -- MODIS dataset (pyHDF4)
    variables -- List containing vars to output from HDF file.
    """

    int_time = time.mktime(time.strptime("00:00 01/01/1993", "%H:%M %d/%m/%Y"))

    latitude = np.squeeze(dataset.select("Latitude")[:])
    longitude = np.squeeze(dataset.select("Longitude")[:])
    t = np.squeeze(dataset.select("Scan_Start_Time")[:]) + int_time

    output = {}
    output["Latitude"] = latitude
    output["Longitude"] = longitude
    output["Time"] = t

    for variable in variables:
        variable_data = np.squeeze(dataset.select(variable)[:])
        output[variable] = variable_data

    return output


def concatenate_datadict(datadict1, datadict2):
    """
    Concatenates two datadicts.
    """

    datadict_out = {}
    for k in datadict1.keys():
        datadict_out[k] = np.concatenate((datadict1[k], datadict2[k]))

    return datadict_out


def filter_modis_datadict(datadict, delta=default_delta,
                          variables=default_variables):
    """
    Filter pertinent variables from a load_dataset(..) datadict
    and returns a filtered version

    Keyword arguments:
    dataset -- MODIS dataset (pyHDF4)
    delta -- Distance to filter from MAO site (degrees, default=0.5)
    variables -- List containing vars to output from HDF file.
    """

    latitude = datadict["Latitude"]
    longitude = datadict["Longitude"]
    t = datadict["Time"]

    latitude_delta = latitude - lat_mao
    longitude_delta = longitude - lon_mao

    # Filter variables according to
    filter_indices = latitude_delta < delta
    filter_indices &= longitude_delta < delta
    filter_indices &= latitude_delta > -delta
    filter_indices &= longitude_delta > -delta

    filtered_datadict = {}

    filtered_datadict["Latitude"] = latitude[filter_indices]
    filtered_datadict["Longitude"] = longitude[filter_indices]
    filtered_datadict["Time"] = t[filter_indices]

    for variable in variables:
        filtered_datadict[variable] = datadict[variable][filter_indices]

    filtered_datadict["Delta"] = delta

    return filtered_datadict


def read_modis_datasets(filepaths, delta_list=[default_delta],
                        variables=default_variables):
    """Reads and filters MODIS HDF4 filepaths and concatenates them together.

    Keyword arguments:
    filepaths -- List containing filepaths for the HDF4 datasets.
    delta -- Distance to filter from MAO site (degrees, default=0.5)
    variables -- List containing vars to output from HDF file.
    """
    delta_list.sort(reverse=True)

    output = []
    datadict = {}

    i = 0
    for filepath in filepaths:
        i += 1
        dataset = SD(filepath)
        loaded_datadict = load_dataset(dataset, variables)
        if (i > 1):
            datadict = concatenate_datadict(datadict, loaded_datadict)
        else:
            datadict = loaded_datadict

        if (i == len(filepaths)):
            for delta in delta_list:
                datadict = filter_modis_datadict(datadict, delta, variables)
                output.append(datadict)
    return output


def post_process(datadicts):
    """Do some post-processing work on the filtered pandas datadicts."""

    output_list = []
    for datadict in datadicts:
        output = {}
        output["Cloud_Fraction"] = np.nanmean(datadict["Cloud_Fraction"]) / 100
        t = np.nanmean(datadict["Time"])
        if np.isfinite(t):
            output["Time"] = dt.datetime.utcfromtimestamp(t)
        else:
            output["Time"] = np.nan
        output["Count"] = np.sum(np.isfinite(datadict["Cloud_Fraction"]))
        output["Solar_Zenith"] = np.nanmean(datadict["Solar_Zenith"])
        output["Delta"] = datadict["Delta"]
        output_list.append(output)
    return output_list


def generate_output(data_folder=default_data_folder,
                    delta_list=[0.25, 0.5, 1.0],
                    variables=default_variables):
    """Process all MODIS files on given folder.

    Keyword arguments:
    data_folder -- location to MODIS HDF4 files
    delta -- Distance to filter from MAO site (degrees, default=0.5)
    variables -- List containing vars to output from HDF file.
    """
    delta_list.sort(reverse=True)

    filepaths = get_modis_filepaths(data_folder)
    outputs = []

    N = len(filepaths)
    i = 1

    for filepath in filepaths:
        try:
            print("\r%d/%d - %s  " % (i, N, filepath[0]), end="")
            i += 1
            dataset = read_modis_datasets(filepath, delta_list, variables)
            processed_dataset = post_process(dataset)
            outputs.append(processed_dataset)
        except KeyboardInterrupt:
            raise
        except pyhdf.error.HDF4Error:
            print("[Error]")

    flat_output = [item for sublist in outputs for item in sublist]
    output_data = pd.DataFrame(flat_output)
    output_data = output_data.sort_values(by=["Time"])
    return output_data


def main():
    """
    Generic do-all function
    """
    delta_list = np.arange(0, 2, 0.1) + 0.1

    data_folder = os.path.expanduser("~/dados-ic/MODIS_Aqua")
    output_path = os.path.expanduser("~/dados-ic/output/modis_aqua.csv")
    output = generate_output(data_folder, delta_list)
    output.to_csv(output_path, index=False)

    data_folder = os.path.expanduser("~/dados-ic/MODIS_Terra")
    output_path = os.path.expanduser("~/dados-ic/output/modis_terra.csv")
    output = generate_output(data_folder, delta_list)
    output.to_csv(output_path, index=False)

if __name__ == "__main__":
    main()
