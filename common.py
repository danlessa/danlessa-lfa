import os

mao_lat = -(3 + 12/60 + 46.70/3600)
mao_lon = -(60 + 35/60 + 53/3600)


def get_files(folder_path, extension):
    files = os.listdir(folder_path)
    files = [f for f in files if f[-len(extension):] == extension
    
    return files

def process(folder_path, type):
    if data_origin == "MODIS":
        extension = ".hdf"
    else if data_origin == "TSI":
        extension = ".cdf"
    else if data_origin == "RAD":
        extension = ".nc"
    else:
        raise Exception("Data type not specified")
