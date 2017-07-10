#!/usr/bin/python2.7
# -*- coding: utf-8 -*-
"""
Created on Sat May  6 17:03:14 2017

@author: acorreia
"""
#
# Usage:
#  $ python ./08-compile-reff.py /home/acorreia/datafolder outputfilename
#
# Example:
# $ python ./08-compile-reff.py  /home/acorreia/Desktop/phase/data/test/ /home/acorreia/Desktop/phase/data/test/allsites-reffready.csv

import os
import sys
import pandas as pd
import glob

datadir = sys.argv[1]
curdir = os.getcwd()
outfile = sys.argv[2]

slash = ''
if datadir[len(datadir) - 1] != '/':
    slash = '/'

flist = sorted(glob.glob(datadir + slash + '*reffready.csv'))
frame = pd.DataFrame()
lista = []
for fil in flist:
    str = ''
    df = pd.read_csv(fil, index_col=None, header=0)
    lista.append(df)
frame = pd.concat(lista)
csv_out = open(outfile, 'wb')
frame.to_csv(csv_out, index=False)
csv_out.close()

exit()
