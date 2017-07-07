This repository contain my fiddling and code related to my scientific internship at LFA-IFUSP.

My project consists of comparing and validating several measured cloud fraction between different instruments at the GOAmazon site. Specifically, I will compare the following CF sources:

* CF (Long algorithm) from the GOAmazon radiative flux analyzer
* Calculated Xie & Liu CF algorithm using the GOAmazon radiative flux analyzer data
* GOAmazon TSI
* MODIS Aqua and Terra CF product, channel 6 - level 2
* GOES 12 e 13

## Requeriments

All code in this repository was writen and tested in Python 3.5.3, with the interactive notebooks written in Jupyter (>=1.0.0) and IPython (>=6.0.0).

I've used the following Python packages in this project:

* numpy >= 1.12.1
* pandas >= 0.20.1
* scipy >= 0.19.0
* matplotlib >= 2.0.2
* netCDF4 >= 1.2.7
* python-hdf4 >= 0.9

## How to use




## Listings


### Folders

* /var/ - Configuration
* /notebooks/ - Several interactive notebooks
* /img/ - Some output images
* /data/ - Processed data

### Source files

* /src/main.py - 
* /src/common.py -
* /src/XL.py - GOAmazon radiative flux analyzer CFs
* /src/MODIS.py - MODIS Aqua and Terra CFs
* /src/TSI.py - TSI CF
* /src/LIDAR.py - LIDAR CF