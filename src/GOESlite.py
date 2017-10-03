"""
Optimized version of the GOES code used at LFA-USP
Compiled by Danilo Lessa Bernardineli
Author: Alexandre Lima Correia
"""
import common
import numpy as np
import netCDF4 as nc
import os
import pandas as pd
import time
import datetime
import matplotlib.pyplot as plt

######################################################
"""
Constants
Obtained from physconstants.py and goesconstants.py
 """
######################################################

Bfunction = [3.22252E-004, 2.29538E-049, 1.99043E+001, -4.15215E-006,
             1.32596E-008, 8.85164E-001, -3.57508E-002, 5.90154E-004,
             -5.12428E-006, 2.47901E-008, -6.35867E-011, 6.77820E-014,
             3.25625E-007, 5.56472E-052, 4.61646E-004, 2.56534E-009,
             1.26165E-011, 4.35725E-006, 3.09653E-008, 1.29105E-010,
             5.14528E-013, 1.99084E-015, 7.02852E-018, 2.04987E-020]

X0 = 29.0
kVIS = 0.001160
a4 = -0.52227011
C = 1.255
C1 = 1.191066e-5
C2 = 1.438833
b4 = 1.0023802
bb2 = 68.2167
bb4 = 15.6854
g4 = -2.0798856e-6
m2 = 227.3889
m4 = 5.2285
n4 = 937.23449
factor = 232.4 / 354.1531901968
F0 = 9.74384716655927
t0 = 0.75
cloudVISthreshold = 0.125
cloudTempthreshold = 300.0

###########################################
"""
Copy-pasted functions
Obtained from 01-goes-to-pixeldata.py
"""
###########################################


def sunzen(time, lat, lon):
    # function to calculate the sun zenith angle and the Sun-Earth distance
    # input: datetime object at UTC
    # input: lat, lon
    # output: sun zenith angle in degrees
    # output: Sun-Earth distance in astronomical units

    days = time.date() - datetime.date(1900, 1, 1)
    days = days.days + 2
    timehr = (((time.hour) / 24.) +
              (time.minute / 60. / 24.) + (time.second / 60. / 60. / 24.))

    # Calculate Julian Day and Julian Century
    JD = days + 2415018.5 + timehr / 24.
    JC = (JD - 2451545.0) / 36525.0

    # Calculate Geometric Mean Anomaly Sun (deg) and Geometric Mean Lon Sun
    # (deg)
    GMAS = 357.52911 + JC * (35999.05029 - (0.0001537) * JC)
    GMLS = np.mod(280.46646 + JC * (36000.76983 + JC * 0.0003032), 360)

    # Calculate Sun Eq of Ctr, Sun True Long (deg),  Sun True Anom (deg) and
    # Sun App Long (deg)
    SEC = (np.sin(np.deg2rad(GMAS)) *
           (1.914602 - JC * (0.004817 + 0.000014 * JC)) +
           np.sin(np.deg2rad(2 * GMAS)) * (0.019993 - 0.000101 * JC) +
           np.sin(np.deg2rad(3 * GMAS)) * 0.000289)

    STL = GMLS + SEC
    STA = GMAS + SEC
    SAL = STL - 0.00569 - 0.00478 * np.sin(np.deg2rad(125.04 - 1934.136 * JC))

    # Calculate Mean Obliq Ecliptic (deg), Obliq Correction (deg), Eccent
    # Earth Orbit
    MOE = (23 + (26 +
                 ((21.448 - JC * (46.815 + JC *
                  (0.00059 - JC * 0.001813)))) / 60) / 60)

    OC = MOE + 0.00256 * np.cos(np.deg2rad(125.04 - 1934.136 * JC))

    EEO = 0.016708634 - JC * (0.000042037 + 0.0000001267 * JC)

    # Calculate the Equation of Time (min) and the True Solar Time (min)
    vary = np.tan(np.deg2rad(OC / 2)) * np.tan(np.deg2rad(OC / 2))
    EOT = 4 * np.rad2deg(vary *
                         np.sin(2 * np.deg2rad(GMLS)) -
                         2 * EEO * np.sin(np.deg2rad(GMAS)) +
                         4 * EEO * vary * np.sin(np.deg2rad(GMAS)) *
                         np.cos(2 * np.deg2rad(GMLS)) -
                         0.5 * vary * vary * np.sin(4 * np.deg2rad(GMLS)) -
                         1.25 * EEO * EEO * np.sin(2 * np.deg2rad(GMAS)))

    TST = np.mod(timehr * 1440 + EOT + 4 * lon, 1440)

    # Calculate the Hour Angle (deg) and the Sun Declination (deg)
    HA = TST / 4 - 180
    id = TST < 0
    HA[id] = TST[id] / 4 + 180
    SD = np.rad2deg(np.arcsin(np.sin(np.deg2rad(OC)) *
                    np.sin(np.deg2rad(SAL))))

    # Calculate Sun-Earth distance (AUs)
    SED = (1.000001018 * (1 - EEO * EEO)) / (1 + EEO * np.cos(np.deg2rad(STA)))

    # Calculate the Sun Zenith Angle (deg)
    SZA = np.rad2deg(np.arccos(np.sin(np.deg2rad(lat)) *
                               np.sin(np.deg2rad(SD)) +
                               np.cos(np.deg2rad(lat)) *
                               np.cos(np.deg2rad(SD)) *
                               np.cos(np.deg2rad(HA))))

    # Retorna o ângulo zenital solar "SZA" e a distância Terra-Sol "Sun-Earth
    # Distance, ou SED"
    return {"Solar_Zenith": SZA, "Sun_Earth_Distance": SED}


def viewzen(lat, lon, satlat=-0.125907, satlon=-74.574043, ER=6371.0,
            SR=42164.17478):
    """
    Obtained from 01-goes-to-pixeldata.py.
    """
    # function to calculate the satellite view zenith angle
    #
    # input: lat, lon such as -3.0, -60.0 for Manaus
    # output: sattelite zenith angle in degrees

    # Surface polar angle (latitude) and equatorial angle (longitude), in
    # radians
    PA = np.deg2rad(90.0 - lat)
    EA = np.deg2rad(lon)

    # Sattelite polar angle (latitude) and equatorial angle (longitude), in
    # radians
    SPA = np.deg2rad(90.0 - satlat)
    SEA = np.deg2rad(satlon)

    # X,Y,Z projections of R1 vector, from Earth center to surface pixel
    # ER = Earth radius in km, defined in physconstants.py
    R1X = ER * np.sin(PA) * np.cos(EA)
    R1Y = ER * np.sin(PA) * np.sin(EA)
    R1Z = ER * np.cos(PA)

    # X,Y,Z projections of R2 vector, from Earth center to sattelite
    # SR = Sattelite orbit radius in km, defined in physconstants.py
    R2X = SR * np.sin(SPA) * np.cos(SEA)
    R2Y = SR * np.sin(SPA) * np.sin(SEA)
    R2Z = SR * np.cos(SPA)

    # X,Y,Z projections of R3 vector, from surface to sattelite
    R3X = R2X - R1X
    R3Y = R2Y - R1Y
    R3Z = R2Z - R1Z

    # X,Y,Z projections of U, unitary vertical surface vector
    UX = np.sin(PA) * np.cos(EA)
    UY = np.sin(PA) * np.sin(EA)
    UZ = np.cos(PA)

    # Calculate the cosine of VZA, then VZA in degrees
    CVZA = (R3X * UX + R3Y * UY + R3Z * UZ) / \
        np.sqrt(R3X * R3X + R3Y * R3Y + R3Z * R3Z)
    VZA = np.rad2deg(np.arccos(CVZA))

    # Retorna o ângulo zenital de visada do satélite,
    # ou 'sattelite view zenith angle = VZA'
    return VZA


def getemission(Temp, Bfunction):
    """
    Obtained from 01-goes-to-pixeldata.py.
    """
    myB = -1.0 * Temp
    myB[myB < 0] = np.nan  # -9999

    (l0, l1, l2, l3, l4, h0, h1, h2, h3, h4, h5, h6) = Bfunction[:12]

    lotemp = (Temp > 154.0) & (Temp <= 190.0)
    hitemp = (Temp > 190.0)
    myB[lotemp] = l0 + l1 * (np.power(Temp[lotemp], l2)) + \
        l3 * Temp[lotemp] + l4 * (np.power(Temp[lotemp], 2))

    myB[hitemp] = (h0 + h1 * Temp[hitemp] + h2 * (np.power(Temp[hitemp], 2)) +
                   h3 * (np.power(Temp[hitemp], 3)) + h4 *
                   (np.power(Temp[hitemp], 4)) + h5 *
                   (np.power(Temp[hitemp], 5)) + h6 *
                   (np.power(Temp[hitemp], 6)))

    return myB


###########################################
"""
Modified functions
Originals on 01-goes-to-pixeldata.py
"""
###########################################


def reflectance(time, lat, lon, data, temperature, procmode="IR"):
    """
    Calculates reflectance on GOES 1/4 band
    Lat/lon/data/temperature must have same shape.

    Keyword Arguments:

    time -- Time when data was adquired.
    lat -- GOES latitude matrix.
    lon -- GOES longitude matrix.
    data -- GOES band data matrix.
    temperature -- Output from temperature function.
    procmode -- Must be 'IR' (band 4) or 'VIS' (band 1).
    """
    lat = lat
    lon = lon
    Raw = data / 32

    sz_out = sunzen(time, lat, lon)
    solar_zenith = sz_out["Solar_Zenith"]
    solarearth_distance = sz_out["Sun_Earth_Distance"]
    view_zenith = viewzen(lat, lon)
    mu0 = np.cos(np.deg2rad(solar_zenith))

    if procmode == "VIS":
        cs = kVIS * (Raw - X0) * C
    elif procmode == "IR":
        RadIR = (Raw - bb2) / m2 * factor
        Emission = getemission(temperature, Bfunction)
        alpha = (t0 * F0 * mu0 / (np.pi * solarearth_distance ** 2)) - Emission
        cs = (RadIR - Emission) / alpha

    return cs


def temperature(data):
    """
    Calculates temperature using GOES-13 band 4 on Kelvin.

    Keyword Arguments:
    data -- GOES-13 band 4 data matrix.
    """
    Raw = data / 32.  # Convert 16bit to 10bit
    RadIR = (Raw - bb4) / m4
    Teff4 = RadIR
    Teff4 = Teff4 - RadIR + np.nan  # -9999.
    TIR = Teff4
    Teff4 = (C2 * n4) / (np.log(1. + (C1 * np.power(n4, 3.) / RadIR)))
    TIR = a4 + b4 * Teff4 + g4 * (np.power(Teff4, 2.))
    return TIR  # Retorna a temperatura de brilho (Canal 4, far infra-red)


def cloudfraction(refl, temp):
    """
    Obtains cloud fraction value based on reflectance and temperature data.

    Keyword Arguments:
    refl -- output from reflectance(..) function.
    temperature -- output from temperature(..) function.
    """

    cf = np.mean((refl > cloudVISthreshold) & (temp < cloudTempthreshold))
    return cf


"""
Processing functions
"""


def process(datapoint, d=0.5):
    """
    Outputs some variables using an 2-tuple containing IR and VIS
    GOES retrieval data for the same instant.

    Keyword Arguments:
    datapoint -- 2-tuple containing IR/VIS filepaths
    d -- Distance filter (degrees from MAO site)
    """

    IR_path = datapoint[0]
    VIS_path = datapoint[1]

    raw_IR = nc.Dataset(IR_path)
    time = datetime.datetime.fromtimestamp(raw_IR["time"][0])

    lat = np.squeeze(raw_IR["lat"])
    lon = np.squeeze(raw_IR["lon"])
    data = np.squeeze(raw_IR["data"])

    filter_indices = (lat - lat_mao < d)
    filter_indices &= (lat - lat_mao > -d)
    filter_indices &= (lon - lon_mao < d)
    filter_indices &= (lon - lon_mao > -d)

    lat = lat[filter_indices]
    lon = lon[filter_indices]
    data = data[filter_indices]

    temp = temperature(data)
    refl_IR = reflectance(time, lat, lon, data, temp, "IR")
    cf_IR = cloudfraction(refl_IR, temp)

    raw_VIS = nc.Dataset(VIS_path)
    data = np.squeeze(raw_VIS["data"])[filter_indices]
    refl_VIS = reflectance(time, lat, lon, data, temp, "VIS")
    cf_VIS = cloudfraction(refl_VIS, temp)

    return {"CloudFraction_VIS": cf_VIS, "CloudFraction_IR": cf_IR,
            "Refl_VIS": np.mean(refl_VIS), "Refl_IR": np.mean(refl_IR),
            "temp": np.mean(temp), "N_IR": len(refl_IR),
            "N_VIS": len(refl_VIS), "Time": time}


def do_work(folder_paths, d=0.5):
    """
    Apply process function to given folders and return its output

    Keyword Arguments:
    folder_paths -- list containg folders to process
    d -- Distance filter (degrees from MAO site)
    """

    datapoints = []

    for folder_path in folder_paths:
        filelist = os.listdir(folder_path)
        truncated_filelist = [f[:-11] for f in filelist]
        datapoint_list = np.unique(truncated_filelist)
        for datapoint in datapoint_list:
            filepath_VIS = folder_path + "/" + datapoint + ".BAND_01.nc"
            filepath_IR = folder_path + "/" + datapoint + ".BAND_04.nc"
            datapoints.append([filepath_IR, filepath_VIS])

    output = []

    i = 1
    N = len(datapoints)

    for datapoint in datapoints:
        try:
            print("\r%s - %d/%d\t" % (datapoint[0][:-11], i, N), end="")
            output.append(process(datapoint))
            i += 1
        except:
            print("[Error]")

    output = pd.DataFrame(output)
    return output


def main():
    """
    Do the generic thing.
    """

    folder_paths = [os.path.expanduser(
        "~/dados-ic/GOES/2014"), os.path.expanduser("~/dados-ic/GOES/2015")]
    output = do_work(folder_paths, d=0.5)

    output_path = os.expanduser("~/dados-ic/processed/")
    output_path = os.path.join(output_path, "goes.csv")
    output.to_csv(output_path, index=False)


if __name__ == "__init__":
    main()
