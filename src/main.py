import common
import os
import TSI
import MODIS
import XL
import LIDAR
import numpy as np
import pandas as pd

def process(folder_path, data_origin):
    """Process all files on given folder"""

    work_function = None
    if data_origin == "modis-aqua" or data_origin == "modis-terra":
        extension = ".hdf"
        work_function = MODIS.work
    elif data_origin == "tsi":
        extension = ".cdf"
        work_function = TSI.work
    elif data_origin == "rad":
        extension = ".nc"
        work_function = XL.work
    elif data_origin == "lidar":
        extension = ".nc"
        work_function = LIDAR.work
    else:
        raise Exception("Data type not specified")

    paths = common.get_filepaths(folder_path, extension)

    return work_function(paths)


def process_all():
    """Process all files using standard folders."""
    data_paths = {
                  "tsi": "~/dados-ic/maotsiskycoverM1.b1/",
                  "rad": "~/dados-ic/maoradflux1longM1.c2/",
                  "lidar": "~/dados-ic/dlprofwstats4news/",
                  "modis-aqua": "~/dados-ic/MODIS_Aqua/",
                  "modis-terra": "~/dados-ic/MODIS_Terra/"}

    for key in data_paths.keys():
        data_origin = key
        folder_path = os.path.expanduser(data_paths[key])

        print("Working on %s dataset" % data_origin)
        raw_output = pd.DataFrame(process(folder_path, data_origin))
        clean_output = common.higienize_data(raw_output)
        clean_output.to_csv(data_origin + ".csv", index=False)


def main():
    process_all()

if __name__ == "__main__":
    main()
