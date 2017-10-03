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

    cf_thick = dataset["percent_opaque"][:] / 100
    cf_thin = dataset["percent_thin"][:] / 100
    time = dataset["time_offset"][:] + dataset["base_time"][0]

    cf = cf_thick + cf_thin

    return (cf, time)


def work(path_list):
    """Work through path list and consolidate in a list."""

    result = {"TSI_CloudFraction": [], "Time": []}

    total_files = len(path_list)
    i = 0

    for path in path_list:
        i += 1
        print("\r%d/%d - %s\t" % (i, total_files, path), end="")
        res = do_work(path)
        result["TSI_CloudFraction"].extend(res[0])
        result["Time"].extend(res[1])

    return result
