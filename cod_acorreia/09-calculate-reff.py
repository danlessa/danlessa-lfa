#!/usr/bin/python2.7
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 19 14:15:03 2015

@author: acorreia
"""
# 
# Usage: 
#  $ python ./09-calculate-reff.py /home/acorreia/outputfolder/filename /home/acorreia/LUTfile   /home/acorreia/Desktop/phase/data/test/reffxThistogram.png
#
# Example:
# $ python ./09-calculate-reff.py /home/acorreia/Desktop/phase/data/test/Thi_Tcw_test-reffready.csv /home/acorreia/Desktop/phase/code/production2/goes12.effective-radius-vs-reflectance-lut.csv  /home/acorreia/Desktop/phase/data/test/reffxThistogram.png

import sys
import pandas as pd
import numpy as np
from calcradius2 import calcreff
from makehist2d import hist2d
from physconstants import *

datafile=sys.argv[1]
lutfile=sys.argv[2]
histsavename=''
if len(sys.argv)>3 :
    histsavename=sys.argv[3]
outfile=datafile[:-4]+'-final.csv'
reff=calcreff(datafile,lutfile)
frame = pd.DataFrame()
frame  = pd.read_csv(datafile,index_col=False,header=0,names=['sat','wavelength','site','season','year','julian','hhmmss','sza','vza','ref063','ref390','temp','phase','tcw','thi'])
frame.insert(12,'reff',reff)
csv_out = open(outfile, 'ab+')
frame.to_csv(csv_out,index=False)
csv_out.close()

if (histsavename !=''):
    temp=np.array(frame.temp)
    valid=temp>0
    if np.any(valid):
        temp=temp[valid]
        temp=temp-zeroc
        reff=np.array(reff)
        reff=reff[valid]
        valid2=reff>0
        if np.any(valid2):
            reff=reff[valid2]
            temp=temp[valid2]
            hist2d(reff,temp,savename=histsavename)

exit()

