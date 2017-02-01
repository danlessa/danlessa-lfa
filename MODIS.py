import numpy as np
from pyhdf.SD import *
import time
import os
from geopy.distance import vincenty

lat_mao = -(3 + 12/60 + 46.70/3600)
lon_mao = -(60 + 35/60 + 53/3600)
coord_mao = (lat_mao, lon_mao)
output_path = "modis.csv"

initial_time = time.mktime(time.strptime("00:00 01/01/1993", "%H:%M %d/%m/%Y"))

def do_work(path):
    sd = SD(path)
    lat = np.array(sd.select("Latitude")[:])
    lon = np.array(sd.select("Longitude")[:])
    cf = np.array(sd.select("Cloud_Fraction")[:])
    t = np.array(sd.select("Scan_Start_Time")[:]) + initial_time
    coord_pt = (lat, lon)
    get_dist = np.vectorize(vincenty, excluded=(0, 1))
    dist = get_dist(coord_mao, coord_pt)
    ind = np.where(dist == np.min(dist))
    x = ind[0]
    y = ind[1]

    lat = lat[x-1:x+2, y-1:y+2]
    lon = lon[x-1:x+2, y-1:y+2]
    cf = cf[x-1:x+2, y-1:y+2]
    t = t[x-1:x+2, y-1:y+2]
    dist = dist[x-1:x+2, y-1:y+2]

    d = dist / np.sum(dist)

    def weighted_average(x):
        return np.sum(x * d)

    total_dist = np.sum(dist)
    cf = weighted_average(cf)
    lat = weighted_average(lat)
    lon = weighted_average(lon)
    t = weighted_average(t)

    return (cf, lat, lon, t)


def main():
    files = os.listdir()
    files = [file for file in files if file[-4:] == ".hdf"]
    with open(output_path, "w") as out:
        out.write("CloudFraction,Latitude,Longitude,ScanTime\n")
        i = 0
        j = len(files)
        for file in files:
            row = do_work(file)
            out.write("%f,%f,%f,%d\n" % row)
            print("%d/%d\t" % (i, j) + file + "\t" + str(row))
            i += 1

if __name__ == "__main__":
    main()
