#!/usr/bin/python2.7
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 19 14:15:03 2015

@author: acorreia
"""


#
# Usage:
#  $ python ./07-prepare-reff.py /home/acorreia/datafolder inputfilename outputfilename
#
# Example:
# $ python ./07-prepare-reff.py  /home/acorreia/Desktop/phase/data/test/MN/DRY  /home/acorreia/Desktop/phase/data/test/Thi_Tcw_test.csv  /home/acorreia/Desktop/phase/data/test/MN-DRY-reffready.csv

import os
import sys
import pandas as pd
import glob
import numpy as np

datadir = sys.argv[1]
curdir = os.getcwd()
infile = sys.argv[2]
outfile = sys.argv[3]
inheader = ['asat', 'asite', 'aseason', 'ayear', 'ajulian', 'ahhmmss', 'aboxsize', 'aNtot', 'aNcld', 'acf', 'aClRefAvg', 'aClRefStd', 'aClTmpAvg', 'aClTmpStd', 'aNCRefAvg', 'aNCRefStd', 'aNCTmpAvg', 'aNCTmpStd', 'aN_cw', 'amin_temp_cw', 'amax_temp_cw', 'amedian_temp_cw', 'amean_temp_cw', 'astdev_temp_cw',
            'aN_hi', 'amin_temp_hi', 'amax_temp_hi', 'amedian_temp_hi', 'amean_temp_hi', 'astdev_temp_hi', 'aN_mix', 'amin_temp_mix', 'amax_temp_mix', 'amedian_temp_mix', 'amean_temp_mix', 'astdev_temp_mix', 'aN_hw', 'amin_temp_hw', 'amax_temp_hw', 'amedian_temp_hw', 'amean_temp_hw', 'astdev_temp_hw']
cwhi = pd.read_csv(infile, index_col=False, header=0, names=inheader)
flist = sorted(glob.glob(datadir + '/*_stat.csv'))
frame = pd.DataFrame()
list = []
count = 0
while count < len(flist):
    filename = flist[count]
    filename = filename.split('/')[len(filename.split('/')) - 1]
    sat = filename.split('.')[0]
    year = filename.split('.')[1]
    julian = filename.split('.')[2]
    last = filename.split('.')[3]
    hhmmss = last.split('_')[0]
    boxsize = last.split('_')[1]

    data = pd.read_csv(flist[count], index_col=False, header=0, names=['sat', 'site', 'season', 'year',
                                                                       'julian', 'hhmmss', 'sza', 'vza', 'boxsize', 'ref063', 'ref390', 'temp', 'ratio', 'phase'])

    site = data.site[0]
    num = np.shape(data)[0]
    # Remove last 17 lines of dataframe
    data = data.drop(data.index[list(range(num - 17 - 1, num))])
    num = np.shape(data)[0]
    wl = np.ones(num) * 3.9
    data.insert(1, 'wavelength', wl)

    id = (cwhi.asat == sat) & (cwhi.asite == site) & (cwhi.ayear == int(year)) & (
        cwhi.ajulian == int(julian)) & (cwhi.ahhmmss == int(hhmmss)) & (cwhi.aboxsize == boxsize)
    tcw1 = np.squeeze(cwhi.amean_temp_cw[id])

    if type(tcw1).__name__ != 'float64':
        tcw1 = np.average(tcw1)
    tcw = np.ones(num) * tcw1
    data.insert(np.shape(data)[1], 'tcw', tcw)

    thi = cwhi.amean_temp_hi[id]
    if np.shape(thi)[0] > 1:
        thi = np.average(thi)
    thi = np.ones(num) * np.squeeze(cwhi.amean_temp_hi[id])
    data.insert(np.shape(data)[1], 'thi', thi)

    data.drop(['boxsize', 'ratio'], inplace=True, axis=1)

    list.append(data)
    count = count + 1

if len(list) > 0:
    frame = pd.concat(list)
    headerflag = True
    if os.path.exists(outfile):
        headerflag = False
    csv_out = open(outfile, 'ab+')
    frame.to_csv(csv_out, index=False, header=headerflag)
    csv_out.close()
os.chdir(curdir)
exit()
