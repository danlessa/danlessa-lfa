import os
import pandas as pd
import numpy as np

lat_mao = -(3 + 12/60 + 46.70/3600)
lon_mao = -(60 + 35/60 + 53/3600)


def load_data(path, data_origin):

        if data_origin == "modis-aqua" or data_origin == "modis-terra":
            raw = pd.read_csv(path, na_values="--", parse_dates=["Time"])
        elif data_origin == "tsi":
            raw = pd.read_csv(path, na_values="--", parse_dates=["Time"])
        elif data_origin == "rad":
            raw = pd.read_csv(path, na_values="--", parse_dates=["Time"])
        else:
            raise Exception("Data type not specified")

        return raw


def get_files(folder_path, extension):
    """Get all files with given extension on given folder path."""
    files = os.listdir(folder_path)
    files = [f for f in files if f[-len(extension):] == extension]
    return files


def get_filepaths(folder_path, extension):
    files = get_files(folder_path, extension)
    paths = [folder_path + f for f in files]
    return paths


def higienize_data(dataset):
    """Coerces all cols on dataset to numeric and
        sort them according to time"""
    print("Higienizando dados")
    for column in dataset.columns:
        dataset[column] = pd.to_numeric(dataset[column], errors="coerce")
    dataset = dataset.sort_values("Time",
                                  ascending=True).reset_index(drop=True)
    dataset.Time = pd.to_datetime(dataset.Time, unit='s')
    return dataset


def validate_data(data_dict):
    """Given an dict, return an subset of it containing only finite elements"""
    for k in data_dict.keys():
        data = data_dict[k]
        for c in data.columns:
            if data[c].dtype == np.float64:
                inds = np.isfinite(data[c])
                data = data[inds]
    return data_dict


def convert_time(dataset):
    """Wrapper function"""
    dataset.Time = pd.to_datetime(dataset.Time)
    return dataset


def load():
    """Load data from standard csv files."""
    raw_XL = pd.read_csv("CF-XL.csv", na_values="--")
    raw_TSI = pd.read_csv("CF-TSI.csv", na_values="--")
    raw_terra = pd.read_csv("modis-terra.csv", na_values="--")
    raw_aqua = pd.read_csv("modis-aqua.csv", na_values="--")
    data = {"rad": raw_XL, "tsi": raw_TSI, "terra": raw_terra,
            "aqua": raw_aqua}

    for k in data.keys():
        data[k] = convert_time(data[k])

    return data
