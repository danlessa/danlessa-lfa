import numpy as np
from pyhdf.SD import *
import time
import os

output_path = "modis.csv"
lat_mara = -(3 + 12/60 + 46.70/3600)
lon_mara = -(60 + 35/60 + 53/3600)
initial_time = time.mktime(time.strptime("00:00 01/01/1993", "%H:%M %d/%m/%Y"))


def distance(lat, lon):
    return np.sqrt(np.power(lat-lat_mara, 2) + np.power(lon - lon_mara, 2))


def weighted_average(x, d):
    return np.sum(x * d)


def work(file):
    sd = SD(file)
    lat = np.array(sd.select("Latitude")[:])
    lon = np.array(sd.select("Longitude")[:])
    cf = np.array(sd.select("Cloud_Fraction")[:])
    t = np.array(sd.select("Scan_Start_Time")[:]) + initial_time
    dist = distance(lat, lon)
    ind = np.where(dist == dist.min())
    x = ind[0]
    y = ind[1]

    lat = lat[x-1:x+2, y-1:y+2]
    lon = lon[x-1:x+2, y-1:y+2]
    cf = cf[x-1:x+2, y-1:y+2]
    t = t[x-1:x+2, y-1:y+2]
    dist = dist[x-1:x+2, y-1:y+2]
    d = dist / np.sum(dist)

    cf = weighted_average(cf, d)
    lat = weighted_average(lat, d)
    lon = weighted_average(lon, d)
    t = weighted_average(t, d)

    return (cf, lat, lon, t)


def main():
    files = os.listdir()
    files = [file for file in files if file[-4:] == ".hdf"]
    with open(output_path, "w") as out:
        out.write("CloudFraction,Latitude,Longitude,ScanTime\n")
        i = 0
        j = len(files)
        for file in files:
            row = work(file)
            out.write("%f,%f,%f,%d\n" % row)
            print("%d/%d\t" % (i, j) + file + "\t" + str(row))
            i += 1

if __name__ == "__main__":
    main()
