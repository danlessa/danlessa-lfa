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

# MODIS base time
initial_time = time.mktime(time.strptime("00:00 01/01/1993", "%H:%M %d/%m/%Y"))

# Default settings for data filtering and post-processing output
default_delta = 0.5
default_variables = ["Cloud_Fraction", "Solar_Zenith"]

default_data_folder = os.path.expanduser("~/dados-ic/MODIS_Aqua")
default_output_path = os.path.expanduser("~/dados-ic/output/modis_aqua.csv")


def get_modis_filepaths(folder_path):
    """
    Returns a list of tuple MODIS filepaths. An element can contain two
    elements if there is more than a file associated with the same instant.
    """

    # Getting list of files
    files = os.listdir(data_path)

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
    files = os.listdir(data_path)
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

    dup_filepaths = [(data_path + "/" + f[0], data_path + "/" + f[1])
                     for f in dup_files]

    filepaths = [data_path + "/" + f for f in files]
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


def filter_modis_dataset(dataset, delta=default_delta,
                         variables=default_variables):
    """
    Filter pertinent variables from a MODIS HDF4 dataset
    and returns a pandas dataset.

    Keyword arguments:
    dataset -- MODIS dataset (pyHDF4)
    delta -- Distance to filter from MAO site (degrees, default=0.5)
    variables -- List containing vars to output from HDF file.
    """

    # MODIS06_L2 base time

    _time = time.mktime(time.strptime("00:00 01/01/1993", "%H:%M %d/%m/%Y"))

    latitude = np.array(dataset.select("Latitude")[:])
    longitude = np.array(dataset.select("Longitude")[:])

    t = np.array(dataset.select("Scan_Start_Time")[:]) + initial_time

    latitude_delta = latitude - lat_mao
    longitude_delta = longitude - lon_mao

    # Filter variables according to
    filter_indices = latitude_delta < delta
    filter_indices &= longitude_delta < delta
    filter_indices &= latitude_delta > -delta
    filter_indices &= longitude_delta > -delta

    filtered_dataset = {}

    filtered_dataset["Latitude"] = latitude[filter_indices]
    filtered_dataset["Longitude"] = longitude[filter_indices]
    filtered_dataset["Time"] = t[filter_indices]
    for variable in variables:
        variable_data = np.array(dataset.select(variable)[:])
        variable_data = variable_data[filter_indices]
        filtered_dataset[variable] = variable_data

    return pd.DataFrame(filtered_dataset)


def read_modis_datasets(filepaths, delta=default_delta,
                        variables=default_variables):
    """Reads and filters MODIS HDF4 filepaths and concatenates them together.

    Keyword arguments:
    filepaths -- List containing filepaths for the HDF4 datasets.
    delta -- Distance to filter from MAO site (degrees, default=0.5)
    variables -- List containing vars to output from HDF file.
    """
    datasets = []

    for filepath in filepaths:
        dataset = SD(filepath)
        processed_dataset = filter_modis_dataset(dataset, delta, variables)
        datasets.append(processed_dataset)

    return pd.concat(datasets)


def post_process(dataset):
    """Do some post-processing work on the filtered pandas dataset."""
    output = {}
    output["Cloud_Fraction"] = np.nanmean(dataset.Cloud_Fraction)
    output["Time"] = np.nanmean(dataset.Time)
    output["Count"] = np.sum(np.isfinite(dataset.Cloud_Fraction))
    output["Solar_Zenith"] = np.nanmean(dataset.Solar_Zenith)
    return output


def generate_output(data_folder=default_data_folder, delta=default_delta,
                    variables=default_variables):
    """Process all MODIS files on given folder.

    Keyword arguments:
    data_folder -- location to MODIS HDF4 files
    delta -- Distance to filter from MAO site (degrees, default=0.5)
    variables -- List containing vars to output from HDF file.
    """
    filepaths = get_modis_filepaths(data_folder)
    outputs = []

    N = len(filepaths)
    i = 1

    for filepath in filepaths:
        try:
            dataset = read_modis_datasets(filepath, delta, variables)
            processed_dataset = post_process(dataset)
            outputs.append(processed_dataset)
        except Exception as e:
            print("[Erro]")

        print("\r%d/%d - %s\t" % (i, N, filepath[0]), end="")
        i += 1

    output_data = pd.DataFrame(outputs)
    output_data = output_data.sort_values(by=["Time"])
    return output_data
