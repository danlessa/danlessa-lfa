import numpy as np
from pyhdf.SD import *
import time
import os
import pandas as pd

lat_mao = -(3 + 12/60 + 46.70/3600)
lon_mao = -(60 + 35/60 + 53/3600)
output_path = "modis.csv"
data_folder = "/home/danilo/dados-ic/modis-aqua/"

initial_time = time.mktime(time.strptime("00:00 01/01/1993", "%H:%M %d/%m/%Y"))

def do_work(path):
    sd = SD(path)
    lat = np.array(sd.select("Latitude")[:])
    lat_dist = np.abs(lat - lat_mao)

    lon = np.array(sd.select("Longitude")[:])
    lon_dist = np.abs(lon - lon_mao)

    cf = np.array(sd.select("Cloud_Fraction")[:])
    t = np.array(sd.select("Scan_Start_Time")[:]) + initial_time

    inds = (lat_dist < 0.5) & (lon_dist < 0.5)

    cf_avg = np.nanmean(cf[inds]) / 100
    t_avg = np.nanmean(t[inds])
    point_count = len(cf[inds])

    return {"CloudFraction": cf_avg, "Time": t_avg, "Count": point_count}



def main():
    files = os.listdir(data_folder)
    files = [file for file in files if file[-4:] == ".hdf"]

    i = 0
    files_count = len(files)
    cf = np.zeros(files_count)
    t = np.zeros(files_count)
    count = np.zeros(files_count)

    for file in files:

        status = "%s/%s - %s" % (i, files_count, file)
        print(status)
        file_path = data_folder + file
        output = do_work(file_path)
        cf[i] = output["CloudFraction"]
        t[i] = output["Time"]
        count[i] = output["Count"]
        i += 1


    print("CFs obtidas")
    data = pd.DataFrame({"CloudFraction": cf, "Time": t, "Count": count})
    data.to_csv("modis.csv", index=False)
    print("Arquivo salvo em modis.csv")




if __name__ == "__main__":
    main()
