"""
Functions for processing MODIS L6 data
Author: Danilo Lessa Bernardineli
"""
import common
import numpy as np
from pyhdf.SD import *
import time
import os
import pandas as pd

lat_mao = common.lat_mao
lon_mao = common.lon_mao

# MODIS base time
initial_time = time.mktime(time.strptime("00:00 01/01/1993", "%H:%M %d/%m/%Y"))


def do_work(path):
    """Extract variables from file."""
    sd = SD(path)
    lat = np.array(sd.select("Latitude")[:])
    lat_delta = lat - lat_mao

    lon = np.array(sd.select("Longitude")[:])
    lon_delta = lon - lon_mao

    cf = np.array(sd.select("Cloud_Fraction")[:])
    t = np.array(sd.select("Scan_Start_Time")[:]) + initial_time

    d = 0.5
    inds = (lat_delta < d) & (lon_delta < d)
    inds &= (lat_delta > -d) & (lon_delta > -d)

    lat_delta = lat_delta[inds]
    lon_delta = lon_delta[inds]

    lat_mean = np.nanmean(lat_delta)
    lon_mean = np.nanmean(lon_delta)

    cut = 0.03
    cutted = False

    if (np.abs(lat_mean) > cut) or (np.abs(lon_mean) > cut):
        cutted = True

    cf_avg = np.nanmean(cf[inds]) / 100
    t_avg = np.nanmean(t[inds])
    point_count = len(cf[inds])

    return (cf_avg, t_avg, point_count, cutted)


def work(path_list):
    """Work through path list and consolidate in a list."""
    result = {"CloudFraction": [], "Time": [], "Count": [], "Outlier": []}

    total_files = len(path_list)
    i = 0

    for path in path_list:
        i += 1
        print("%s/%s - %s" % (i, total_files, path))
        res = do_work(path)
        result["CloudFraction"].append(res[0])
        result["Time"].append(res[1])
        result["Count"].append(res[2])
        result["Outlier"].append(res[3])

    return result
