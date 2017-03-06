import os
import pandas as pd
import numpy as np

mao_lat = -(3 + 12/60 + 46.70/3600)
mao_lon = -(60 + 35/60 + 53/3600)

def get_files(folder_path, extension):
    files = os.listdir(folder_path)
    files = [f for f in files if f[-len(extension):] == extension] 
    return files

def process(folder_path, type):
    if data_origin == "MODIS":
        extension = ".hdf"
    elif data_origin == "TSI":
        extension = ".cdf"
    elif data_origin == "RAD":
        extension = ".nc"
    else:
        raise Exception("Data type not specified")
             
def higienize_data(dataset):
    print("Higienizando dados")
    for column in dataset.columns:
        dataset[column] = pd.to_numeric(dataset[column], errors="coerce")
    dataset = dataset.sort_values("Time", ascending=True).reset_index(drop=True)
    dataset.Time = pd.to_datetime(dataset.Time, unit='s')
    return dataset

def validate_data(data_dict):
    for k in data_dict.keys():
        data = data_dict[k]
        for c in data.columns:
            if data[c].dtype == np.float64:
                inds = np.isfinite(data[c])
                data = data[inds]
    return data_dict

def convert_time(dataset):
    dataset.Time = pd.to_datetime(dataset.Time)
    return dataset


def load():
    raw_XL = pd.read_csv("CF-XL.csv", na_values="--")
    raw_TSI = pd.read_csv("CF-TSI.csv", na_values="--")
    raw_terra = pd.read_csv("modis-terra.csv", na_values="--")
    raw_aqua = pd.read_csv("modis-aqua.csv", na_values="--")    
    data = {"rad": raw_XL, "tsi": raw_TSI, "terra": raw_terra, "aqua": raw_aqua}
    
    for k in data.keys():
        data[k] = convert_time(data[k])
        
    return data
