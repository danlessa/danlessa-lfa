#!/usr/bin/env python3
"""
This file contains an implementation of the functions and procedures described
in Xie & Liu 2013 paper for obtaining cloud fraction values
"""
import common
import numpy as np
import netCDF4 as nc
import os
import pandas as pd
import scipy.optimize as opt
import time
import datetime
import pysolar

g = 0.87
e = 1e-14


def vector_solve(A, B, f):
    if len(A) != len(B):
        return None
    n = len(A)
    mtr = np.zeros(n)
    for i in range(0, n):
        mtr[i] = opt.brentq(f, 0 + e, 1 - e, args=(A[i], B[i]))

    return mtr


def get_B(Fdn_clr, Fdn_all, Fup_all, Fdn_clr_d, Fdn_all_d, T):
    """Obtains B1 and B2 parameter used in Xie & Liu 2013 paper (eq 13).
    Keyword arguments:
    Fdn_clr -- Downward flux in clear sky
    Fdn_all -- Downward flux
    Fup_all -- Upward flux
    Fdn_clr_d -- Diffuse downward flux in clear sky
    Fdn_all_d -- Diffuse downward flux
    T -- Cloud transmissibility"""
    invalid_inds = np.ma.getmask(Fdn_clr)
    invalid_inds |= np.ma.getmask(Fdn_all)
    invalid_inds |= np.ma.getmask(Fup_all)
    invalid_inds |= np.ma.getmask(Fdn_clr_d)
    invalid_inds |= np.ma.getmask(Fdn_all_d)
    invalid_inds |= np.ma.getmask(T)
    B1 = (Fdn_clr - Fdn_all) / (Fdn_clr - Fup_all * T * T)
    B2 = (Fdn_clr_d - Fdn_all_d) / Fdn_clr_d

    B1[invalid_inds] = np.nan
    B2[invalid_inds] = np.nan
    B1 = np.array(B1)
    B2 = np.array(B2)
    return (B1, B2)


def cloud_albedo_func(alfa_r, B1, B2):
    """Cloud albedo function -- used for obtaining the cloud albedo through some
    root finding method.
    Keyword Arguments:
    alfa_r -- Cloud albedo
    B1 -- ...
    B2 -- ..."""
    termo = 1 - np.exp(-2 * alfa_r / ((1 - alfa_r) * (1 - g)))
    return (alfa_r / termo) - B1 / B2
cloud_albedo_func = np.vectorize(cloud_albedo_func)


def get_cloud_albedo(B1, B2):
    """Calculates an aproximate cloud albedo value using the development in
     eq. 14 from Xie & Liu 2013 article.
    Keyword arguments:
    B1 -- B1 parameter from get_B function
    B2 -- B2 parameter from get_B function"""
    invalid_inds = (B1 < 0) | (B1 > 1) | (B2 < 0) | (
        B2 > 1) | (B1 / B2 > 1.0) | (
        B1 / B2 < 0.07) | ~np.isfinite(B1) | ~np.isfinite(B2)
    zero_inds = ((B2 == 0) | (B1 == 0)) & (~invalid_inds)
    one_inds = (B1 == B2) & (~invalid_inds)
    calc_inds = (~one_inds) & (~zero_inds) & (~invalid_inds)

    cloud_albedo = np.zeros(len(B1))
    cloud_albedo[invalid_inds] = np.nan
    cloud_albedo[zero_inds] = 0.0
    cloud_albedo[one_inds] = 1.0

    B1t = B1[calc_inds]
    B2t = B2[calc_inds]

    cloud_albedo[calc_inds] = vector_solve(B1t, B2t, cloud_albedo_func)

    return cloud_albedo


def cloud_albedo_correction(cloud_albedo, solar_zenith_cosine):
    """Cloud albedo correction for cloud absorption (eq. 17c). """
    return cloud_albedo / (1.0537 + 0.0788 * solar_zenith_cosine)


def get_cloud_fraction(cloud_albedo, B1):
    """Formula for obtaining cloud fraction."""
    return B1 / cloud_albedo


def do_work(path):
    """Do all hard work for a file."""
    # Variable loading
    dataset = nc.Dataset(path, "r")

    base_time = dataset["base_time"][0]
    time_offset = dataset["time_offset"][:]
    timestamp = time_offset + base_time
    Fdn_clr = dataset["clearsky_downwelling_shortwave"][:]
    Fdn_all = dataset["downwelling_shortwave"][:]
    Fdn_clr_d = dataset["clearsky_direct_downwelling_shortwave"][:]
    rad_cf_short = dataset["cloudfraction_shortwave"][:]
    rad_cf_long = dataset["cloudfraction_longwave"][:]
    Fup_all = dataset["upwelling_shortwave"][:]
    Fdn_all_d = dataset["direct_downwelling_shortwave"][:]
    cod = dataset["visible_cloud_optical_depth"][:]
    T = dataset["cloud_transmissivity_shortwave"][:]
    solar_zen_cosine = dataset["cosine_zenith"][:]

    # Obtaining cloud albedo
    (B1, B2) = get_B(Fdn_clr, Fdn_all, Fup_all, Fdn_clr_d, Fdn_all_d, T)
    cloud_albedo = get_cloud_albedo(B1, B2)
    cloud_albedo_cor = cloud_albedo_correction(cloud_albedo, solar_zen_cosine)

    # Obtaining CF
    CF1 = get_cloud_fraction(cloud_albedo, B1)
    CF2 = get_cloud_fraction(cloud_albedo_cor, B1)

    return (cloud_albedo, cloud_albedo_cor, CF1, CF2, timestamp, rad_cf_short,
            rad_cf_long)


def work(path_list):
    """Work through path list and consolidate in a list."""

    result = {"Long_SW_CloudFraction": [], "Long_LW_CloudFraction": [],
              "XL_CloudAlbedo": [],
              "XL_CloudAlbedoCorrected": [], "XL_CloudFraction": [],
              "XL_CloudFractionCorrected": [], "Time": []}

    total_files = len(path_list)
    i = 0

    for path in path_list:
        print("%s/%s - %s" % (i, total_files, path))
        i += 1
        res = do_work(path)
        result["XL_CloudAlbedo"].extend(res[0])
        result["XL_CloudAlbedoCorrected"].extend(res[1])
        result["XL_CloudFraction"].extend(res[2])
        result["XL_CloudFractionCorrected"].extend(res[3])
        result["Time"].extend(res[4])
        result["Long_SW_CloudFraction"].extend(res[5])
        result["Long_LW_CloudFraction"].extend(res[6])

    return result
