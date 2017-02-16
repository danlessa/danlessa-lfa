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
    print("Salvando")
    return dataset
