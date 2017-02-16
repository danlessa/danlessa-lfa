"""

"""
import common
import netCDF4 as nc
import os
import numpy as np
import pandas as pd

data_folder = "dados/maotsi/"


def do_work(path):
    dataset = nc.Dataset(path)

    cf_thick = dataset["percent_opaque"][:] / 100
    cf_thick_qc = dataset["qc_percent_opaque"][:]

    cf_thin = dataset["percent_thin"][:] / 100
    cf_thin_qc = dataset["qc_percent_thin"][:]

    time = dataset["time_offset"][:] + dataset["base_time"][0]

    cf = cf_thick + cf_thin

    return (cf, time)


def main():

    print("Adquirindo o CF do TSI")
    print("Data path:" + data_folder)
    files = os.listdir(data_folder)
    files = [f for f in files if f[-4:] == ".cdf"]

    result = {"CloudFraction": [], "Time": []}

    i = 0
    n = len(files)
    for filename in files:
        path = data_folder + filename

        res = do_work(path)
        result["CloudFraction"].extend(res[0])
        result["Time"].extend(res[1])

        i += 1
        print("%s/%s\t%s" % (i, n, path))

    df = pd.DataFrame(result)
    df = common.higienize_data(df)
    df.to_csv("CF-TSI.csv", index=False)

if __name__ == "__main__":
    main()
