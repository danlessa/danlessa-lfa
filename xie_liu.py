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

__author__ = "Danilo Lessa Bernardineli"
__email__ = "danilo.bernardineli@usp.br"

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
    B1 = (Fdn_clr - Fdn_all) / (Fdn_clr - Fup_all * T * T)
    B2 = (Fdn_clr_d - Fdn_all_d) / Fdn_clr_d

    #B1[B1.mask] = np.nan
    #B2[B2.mask] = np.nan
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
        B2 > 1) | (B1 / B2 > 1.0) | (B1 / B2 < 0.068)
    zero_inds = ((B2 == 0) | (B1 == 0)) & (~invalid_inds)
    one_inds = (B1 == B2) & (~invalid_inds)
    calc_inds = (~one_inds) & (~zero_inds) & (~invalid_inds)

    cloud_albedo = np.zeros(len(B1))
    cloud_albedo[zero_inds] = 0.0
    cloud_albedo[one_inds] = 1.0
    #cloud_albedo[invalid_inds] = np.nan

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
    Fdn_clr = dataset["clearsky_downwelling_shortwave"][:]
    Fdn_all = dataset["downwelling_shortwave"][:]
    Fup_all = dataset["upwelling_shortwave"][:]
    Fdn_clr_d = dataset["clearsky_direct_downwelling_shortwave"][:]
    Fdn_all_d = dataset["direct_downwelling_shortwave"][:]
    rad_cf_short = dataset["cloudfraction_shortwave"][:]
    rad_cf_long = dataset["cloudfraction_longwave"][:]
    rad_cf_short_status = dataset["cloudfraction_shortwave_status"][:]
    T = dataset["cloud_transmissivity_shortwave"][:]
    solar_zen_cosine = dataset["cosine_zenith"][:]

    # Obtaining cloud albedo
    (B1, B2) = get_B(Fdn_clr, Fdn_all, Fup_all, Fdn_clr_d, Fdn_all_d, T)
    cloud_albedo = get_cloud_albedo(B1, B2)
    cloud_albedo_cor = cloud_albedo_correction(cloud_albedo, solar_zen_cosine)

    # Obtaining CF
    CF1 = get_cloud_fraction(cloud_albedo, B1)
    CF2 = get_cloud_fraction(cloud_albedo_cor, B1)
    time = time_offset + base_time

    return (cloud_albedo, cloud_albedo_cor, CF1, CF2, time, rad_cf_short, rad_cf_long, rad_cf_short_status)


def main():

    data_folder = "dados/maoradflux1long/"

    print("Adquirindo o albedo da nuvem e CF pelo algoritimo de Xie & Liu")
    print("Data path:" + data_folder)
    files = os.listdir(data_folder)
    files = [f for f in files if f[-3:] == ".nc"]

    result = {"Shortwave_CloudFraction": [], "Longwave_CloudFraction": [],
              "Shortwave_CloudFraction_Status": [],  "XL_CloudAlbedo": [],
              "XL_CloudAlbedoCorrected": [], "XL_CloudFraction": [],
              "XL_CloudFractionCorrected": [], "Time": []}

    i = 0
    n = len(files)
    for filename in files:
        path = data_folder + filename

        res = do_work(path)
        result["XL_CloudAlbedo"].extend(res[0])
        result["XL_CloudAlbedoCorrected"].extend(res[1])
        result["XL_CloudFraction"].extend(res[2])
        result["XL_CloudFractionCorrected"].extend(res[3])
        result["Time"].extend(res[4])
        result["Shortwave_CloudFraction"].extend(res[5])
        result["Longwave_CloudFraction"].extend(res[6])
        result["Shortwave_CloudFraction_Status"].extend(res[7])

        i += 1
        print("%s/%s\t%s" % (i, n, path))

    df = pd.DataFrame(result)
    df = common.higienize_data(df)
    df.to_csv("CF-XL.csv", index=False)


if __name__ == "__main__":
    main()
