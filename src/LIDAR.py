"""
Functions for processing GOAmazon TSI data
Author: Danilo Lessa Bernardineli
"""
import common
import netCDF4 as nc
import os
import numpy as np
import pandas as pd


def do_work(path):
    """Extract variables from file."""
    dataset = nc.Dataset(path)

    cf = dataset["dl_cloud_frequency"]
    time = dataset["time_offset"][:] + dataset["base_time"][0]

    return (cf, time)


def work(path_list):
    """Work through path list and consolidate in a list."""

    result = {"Lidar_CloudFraction": [], "Time": []}

    total_files = len(path_list)
    i = 0

    for path in path_list:
        i += 1
        print("\r%d/%d - %s\t" % (i, total_files, path), end="")
        res = do_work(path)
        result["Lidar_CloudFraction"].extend(res[0])
        result["Time"].extend(res[1])

    return result
