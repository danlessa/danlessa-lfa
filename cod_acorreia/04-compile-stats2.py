#!/usr/bin/python2.7
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 19 14:15:03 2015

@author: acorreia
"""

# Feb19/2015 This script reads CSV files in a given folder and calculates statistics about them. The input files are supposed to have 3 numeric columns 
# containing real numbers. The output files will contain all the original data and statistics about each column.  
#
# Usage: 
#  $ python ./compile-stats2.py /home/acorreia/datafolder /home/acorreia/outputfolder DRY AH
#

import os
import sys
import glob
import pandas as pd
import numpy as np

from physconstants import *

datadir=sys.argv[1]
curdir=os.getcwd()
outdir=sys.argv[2]
os.chdir(outdir)
season=sys.argv[3]
site=sys.argv[4]

maxratio,maxphase='---','---'
minratio,minphase='---','---'
medratio,medphase='---','---'
avgratio,avgphase='---','---'

flist = sorted(glob.glob(datadir+'/*alldata.csv'))
cflist = sorted(glob.glob(datadir+'/*cloudfraction.csv'))

count=0
while count < len(flist):
    filename = flist[count]
    filename = filename.split('/')[len(filename.split('/'))-1]
    outfile  = filename[:-4]+'_stat.csv'
    sat      = filename.split('.')[0]
    year     = filename.split('.')[1]
    julian   = filename.split('.')[2]
    last     = filename.split('.')[3]
    hhmmss   = last.split('_')[0]
    boxsize  = last.split('_')[1]
    data     = pd.read_csv(flist[count],header=None,names=['vis','ir','t','sza','vza'])
    fcf = open(cflist[count], 'r')
    line=fcf.readline()
    fcf.close()
    ntot,nc,cf=line.split(',')[0],line.split(',')[1],line.split(',')[2]
    crefavg,crefstd,ctmpavg,ctmpstd=line.split(',')[3],line.split(',')[4],line.split(',')[5],line.split(',')[6]
    ncrefavg,ncrefstd,nctmpavg,nctmpstd=line.split(',')[7],line.split(',')[8],line.split(',')[9],line.split(',')[10]
    data.fillna(mynan)
    rvis,rir,temp=data.vis,data.ir,data.t
    rvis,rir,temp=np.float64(rvis),np.float64(rir),np.float64(temp)
    ntot,nc,cf=np.int(ntot),np.int(nc),np.float64(cf)
    crefavg,crefstd,ctmpavg,ctmpstd=np.float64(crefavg),np.float64(crefstd),np.float64(ctmpavg),np.float64(ctmpstd)
    ncrefavg,ncrefstd,nctmpavg,nctmpstd=np.float64(ncrefavg),np.float64(ncrefstd),np.float64(nctmpavg),np.float64(nctmpstd)
    rvis,rir,temp=np.array(rvis),np.array(rir),np.array(temp)    
    ntot,nc,cf=str(np.array(ntot)),str(np.array(nc)),str(np.array(cf))
    crefavg,crefstd,ctmpavg,ctmpstd=str(np.array(crefavg)),str(np.array(crefstd)),str(np.array(ctmpavg)),str(np.array(ctmpstd))
    ncrefavg,ncrefstd,nctmpavg,nctmpstd=str(np.array(ncrefavg)),str(np.array(ncrefstd)),str(np.array(nctmpavg)),str(np.array(nctmpstd))    
    nonzero=np.array(rir>0.0)
    ratio=rvis[nonzero]/rir[nonzero]
    icelim = icedry
    waterlim = waterdry
    if season == 'WET':
        icelim = icewet
        waterlim = waterwet    
    nv=str(rvis.size)
    vismin,vismax,vismed,visavg,visstd=str(rvis),str(rvis),str(rvis),str(rvis),'0'       
    if rvis.size>1:
        vismin=str(min(rvis))
        vismax=str(max(rvis))
        vismed=str(np.median(rvis))
        visavg=str(np.average(rvis))
        visstd=str(np.std(rvis))
    ni=str(rir.size)
    irmin,irmax,irmed,iravg,irstd=str(rir),str(rir),str(rir),str(rir),'0'
    if rir.size>1:
        irmin=str(min(rir))
        irmax=str(max(rir))
        irmed=str(np.median(rir))
        iravg=str(np.average(rir))
        irstd=str(np.std(rir))        
    nt=str(temp.size)
    tmin,tmax,tmed,tavg,tstd=str(temp),str(temp),str(temp),str(temp),'0'
    if temp.size>1:
        tmin=str(min(temp))
        tmax=str(max(temp))
        tmed=str(np.median(temp))
        tavg=str(np.average(temp))
        tstd=str(np.std(temp))    
    minratio,minphase='---','---'
    if rir.size>1:        
        if (min(rir)>0.0):
            minphase='mix'
            minratio=float(vismin)/float(irmin)
            if ( (minratio>icelim) and (float(vismin)>cloudVISthreshold) and (float(tmin)<cloudTempthreshold) and (float(tmin)<zeroc) ):   # APR07 2016
                minphase='ice'
            if ( (minratio<waterlim) and (float(vismin)>cloudVISthreshold) and (float(tmin)<cloudTempthreshold) and (float(tmin)>=zeroc)):
                minphase='water'
        minratio=str(minratio)
        maxratio,maxphase='---','---'
        if (max(rir)>0.0):
            maxphase='mix'
            maxratio=float(vismax)/float(irmax)
            if ( (maxratio>icelim) and (float(vismax)>cloudVISthreshold) and (float(tmax)<cloudTempthreshold) and (float(tmax)<zeroc) ):
                maxphase='ice'
            if ( (maxratio<waterlim) and (float(vismax)>cloudVISthreshold) and (float(tmax)<cloudTempthreshold) and (float(tmax)>=zeroc) ):
                maxphase='water'
        maxratio=str(maxratio)    
        medratio,medphase='---','---'
        if (np.median(rir)>0.0):
            medphase='mix'
            medratio=float(vismed)/float(irmed)
            if ( (medratio>icelim) and (float(vismed)>cloudVISthreshold) and (float(tmed)<cloudTempthreshold) and (float(tmed)<zeroc) ):
                medphase='ice'
            if ( (medratio<waterlim) and (float(vismed)>cloudVISthreshold) and (float(tmed)<cloudTempthreshold) and (float(tmed)>=zeroc)):
                medphase='water'
        medratio=str(medratio)        
        avgratio,avgphase='---','---'
        if (np.average(rir)>0.0):
            avgphase='mix'
            avgratio=float(visavg)/float(iravg)
            if ( (avgratio>icelim) and (float(visavg)>cloudVISthreshold) and (float(tavg)<cloudTempthreshold) and (float(tavg)<zeroc) ):
                avgphase='ice'
            if ( (avgratio<waterlim) and (float(visavg)>cloudVISthreshold) and (float(tavg)<cloudTempthreshold) and (float(tavg)>=zeroc)):
                avgphase='water'
        avgratio=str(avgratio)    
    if os.path.isfile(outfile):
        outfn = open(outfile, 'a')
    else:
        header='sat,site,season,year,julian,hhmmss,sza,vza,boxsize,r_vis_0.63,r_IR_3.90,TEMPtop(K)_11,VIS/IR_ratio,phase\n'
        outfn = open(outfile, 'w')
        outfn.writelines(header)    
    nv=str(rvis.size)
    ni=str(rir.size)
    nt=str(temp.size)
    with open(flist[count],'r') as fh0:
        for line in fh0:
            rv,ri,tp,sza,vza=line.split(',')[0],line.split(',')[1],line.split(',')[2],line.split(',')[3],line.split(',')[4]
            vza=vza[:-1]
            tp=float(tp)
            ratio,phase='---','---'
            if (float(ri)>0.0):
                phase='mix'
                ratio=float(rv)/float(ri)
                if ((float(rv)>cloudVISthreshold) and (float(tp)<cloudTempthreshold) and (float(tp)<tlim)):
                    phase='ice'
                if ( (ratio>icelim) and (float(rv)>cloudVISthreshold) and (float(tp)<cloudTempthreshold) and (float(tp)<zeroc) ):
                    phase='ice'
                if ((float(rv)>cloudVISthreshold) and (float(tp)<cloudTempthreshold) and (float(tp)>=zeroc)):
                    phase='water'
                if ( (ratio<waterlim) and (float(rv)>cloudVISthreshold) and (float(tp)<cloudTempthreshold) and (float(tp)>=zeroc)):
                    phase='water'
            ratio=str(ratio)
            tp=str(tp)
            outstr=sat+','+site+','+season+','+year+','+julian+','+hhmmss+','+sza+','+vza+','+boxsize+','+rv+','+ri+','+tp+','+ratio+','+phase+'\n'                        
            outfn.writelines(outstr)
    outstr = sat+','+site+','+season+','+year+','+julian+','+hhmmss+',N,'+nv+','+ni+','+nt+'\n'
    outfn.writelines(outstr)
    outstr = sat+','+site+','+season+','+year+','+julian+','+hhmmss+',MIN,'+vismin+','+irmin+','+tmin+','+minratio+','+minphase+'\n'
    outfn.writelines(outstr)
    outstr = sat+','+site+','+season+','+year+','+julian+','+hhmmss+',MAX,'+vismax+','+irmax+','+tmax+','+maxratio+','+maxphase+'\n'
    outfn.writelines(outstr)
    outstr = sat+','+site+','+season+','+year+','+julian+','+hhmmss+',AVG,'+visavg+','+iravg+','+tavg+','+avgratio+','+avgphase+'\n'
    outfn.writelines(outstr)
    outstr = sat+','+site+','+season+','+year+','+julian+','+hhmmss+',MED,'+vismed+','+irmed+','+tmed+','+medratio+','+medphase+'\n'
    outfn.writelines(outstr)
    outstr = sat+','+site+','+season+','+year+','+julian+','+hhmmss+',STD,'+visstd+','+irstd+','+tstd+'\n'
    outfn.writelines(outstr)
    outstr = sat+','+site+','+season+','+year+','+julian+','+hhmmss+',Ntot,'+ntot+'\n'
    outfn.writelines(outstr)
    outstr = sat+','+site+','+season+','+year+','+julian+','+hhmmss+',Ncld,'+nc+'\n'
    outfn.writelines(outstr)
    outstr = sat+','+site+','+season+','+year+','+julian+','+hhmmss+',CF,'+cf+'\n'
    outfn.writelines(outstr)
    outstr = sat+','+site+','+season+','+year+','+julian+','+hhmmss+',ClRefAvg,'+crefavg+'\n'
    outfn.writelines(outstr)
    outstr = sat+','+site+','+season+','+year+','+julian+','+hhmmss+',ClRefStd,'+crefstd+'\n'
    outfn.writelines(outstr)
    outstr = sat+','+site+','+season+','+year+','+julian+','+hhmmss+',ClTmpAvg,'+ctmpavg+'\n'
    outfn.writelines(outstr)
    outstr = sat+','+site+','+season+','+year+','+julian+','+hhmmss+',ClTmpStd,'+ctmpstd+'\n'
    outfn.writelines(outstr)
    outstr = sat+','+site+','+season+','+year+','+julian+','+hhmmss+',NCRefAvg,'+ncrefavg+'\n'
    outfn.writelines(outstr)
    outstr = sat+','+site+','+season+','+year+','+julian+','+hhmmss+',NCRefStd,'+ncrefstd+'\n'
    outfn.writelines(outstr)
    outstr = sat+','+site+','+season+','+year+','+julian+','+hhmmss+',NCTmpAvg,'+nctmpavg+'\n'
    outfn.writelines(outstr)
    outstr = sat+','+site+','+season+','+year+','+julian+','+hhmmss+',NCTmpStd,'+nctmpstd+'\n'
    outfn.writelines(outstr+'\n')
    count = count + 1
if count>0:
    outfn.close()
os.chdir(curdir)
exit()

