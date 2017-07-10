#!/usr/bin/python2.7
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 18 16:22:06 2014

@author: acorreia
"""

#Usage: 
#  $ python ./compile-data.py /home/acorreia/datafolder /home/acorreia/outputfolder

import os
import sys
import glob

datadir=sys.argv[1]
curdir=os.getcwd()
outdir=curdir
if len(sys.argv)>2 :
    outdir=sys.argv[2]
os.chdir(outdir)

listvis=sorted(glob.glob(datadir+'/*cloudreflectance063.csv'))
listird=sorted(glob.glob(datadir+'/*cloudreflectance390.csv'))
listtmp=sorted(glob.glob(datadir+'/*cloudtemperature.csv'))
listsza=sorted(glob.glob(datadir+'/*sza.csv'))
listvza=sorted(glob.glob(datadir+'/*vza.csv'))
cflist=sorted(glob.glob(datadir+'/*cloudfraction.csv'))

count=0
while count < len(listvis):
    if (os.stat(listvis[count]).st_size == 0):
        os.remove(listvis[count])
        os.remove(listird[count])
        os.remove(listtmp[count])
        os.remove(listsza[count])
        os.remove(listvza[count])
        os.remove(cflist[count])
    else:
        fh1, fh2, fh3 = open(listvis[count], 'r'), open(listird[count], 'r'), open(listtmp[count], 'r')
        fh4, fh5 = open(listsza[count], 'r'), open(listvza[count], 'r')
        outfile = listvis[count]
        name=listvis[count]
        outfile=name.split('/')[len(name.split('/'))-1]
        outfile=outfile[:-23]+'alldata.csv'
        fh6 = open(outfile, 'w')
        for line in fh1:
            d1=line[:-2]
            d2=fh2.readline()[:-1]
            d3=fh3.readline()[:-1]
            d4=fh4.readline()[:-1]
            d5=fh5.readline()[:-1]
            outstr=d1+','+d2+','+d3+','+d4+','+d5+'\n'
            fh6.writelines(outstr)
        fh1.close()
        fh2.close()
        fh3.close()
        fh4.close()
        fh5.close()
        fh6.close()        
    count = count + 1
os.chdir(curdir)
exit()


