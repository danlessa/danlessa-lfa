#!/usr/bin/python2.7
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 19 14:15:03 2015

@author: acorreia
"""

# 
# Usage: 
#  $ python ./coldwater-hotice.py /home/acorreia/datafolder /home/acorreia/outputfolder outputfilename DRY AH
#

import os
import sys
import glob
import pandas as pd
import numpy as np

from physconstants import *


def rdvar(myvar,mylabel):
    idvar=(data.asza==mylabel)
    myvar=data[idvar][['avza']].values
    myvar=myvar.flatten()
    myvar=[float(i) for i in myvar]
    myvar=np.array(myvar)
    myvar=str(myvar[0])
    return myvar


datadir=sys.argv[1]
curdir=os.getcwd()
outdir=sys.argv[2]
outfile=sys.argv[3]
season=sys.argv[4]
site=sys.argv[5]
os.chdir(outdir)
header='sat,site,season,year,julian,hhmmss,boxsize,Ntot,Ncld,cf,ClRefAvg,ClRefStd,ClTmpAvg,ClTmpStd,NCRefAvg,NCRefStd,NCTmpAvg,NCTmpStd'
header=header+',N_cw,min_temp_cw,max_temp_cw,median_temp_cw,mean_temp_cw,stdev_temp_cw,N_hi,min_temp_hi,max_temp_hi,median_temp_hi,mean_temp_hi,stdev_temp_hi'
header=header+',N_mix,min_temp_mix,max_temp_mix,median_temp_mix,mean_temp_mix,stdev_temp_mix'
header=header+',N_hw,min_temp_hw,max_temp_hw,median_temp_hw,mean_temp_hw,stdev_temp_hw'
header=header+'\n'
outfn = open(outfile, 'w')
outfn.writelines(header) 
flist = sorted(glob.glob(datadir+'/*_stat.csv'))
count=0
while count < len(flist):
    filename = flist[count]
    filename = filename.split('/')[len(filename.split('/'))-1]
    sat      = filename.split('.')[0]    
    year     = filename.split('.')[1]    
    julian   = filename.split('.')[2]
    last     = filename.split('.')[3]
    hhmmss   = last.split('_')[0]
    boxsize  = last.split('_')[1]    
    data     = pd.read_csv(flist[count],index_col=None,header=0,names=['asat','asite','aseason','ayear','ajulian','ahhmmss','asza','avza','aboxsize','avis','air','atemp','aratio','aphase'])
    watertemp,icetemp,mixtemp = strnan,strnan,strnan #'-9999','-9999','-9999'
    Ntot,Nc,ClRefAvg,ClRefStd,ClTmpAvg,ClTmpStd,NCRefAvg,NCRefStd,NCTmpAvg,NCTmpStd=strnan,strnan,strnan,strnan,strnan,strnan,strnan,strnan,strnan,strnan #'-9999','-9999','-9999','-9999','-9999','-9999','-9999','-9999','-9999','-9999'
    cloudfraction=strnan #'-9999'
    Ntot=rdvar(Ntot,'Ntot')
    Nc=rdvar(Nc,'Ncld')
    ClRefAvg=rdvar(ClRefAvg,'ClRefAvg')
    ClRefStd=rdvar(ClRefStd,'ClRefStd')
    ClTmpAvg=rdvar(ClTmpAvg,'ClTmpAvg')
    ClTmpStd=rdvar(ClTmpStd,'ClTmpStd')
    NCRefAvg=rdvar(NCRefAvg,'NCRefAvg')
    NCRefStd=rdvar(NCRefStd,'NCRefStd')
    NCTmpAvg=rdvar(NCTmpAvg,'NCTmpAvg')
    NCTmpStd=rdvar(NCTmpStd,'NCTmpStd')

    idcf=(data.asza=='CF')
    cloudfraction=data[idcf][['avza']].values
    cloudfraction=cloudfraction.flatten()
    cloudfraction=[float(i) for i in cloudfraction]
    cloudfraction=np.array(cloudfraction)
    cloudfraction=str(cloudfraction[0])
        
    idx=(data.aboxsize==boxsize) & (data.aphase=='water')
    watertemp=data[idx][['atemp']].values
    watertemp=watertemp.flatten()
    watertemp=[float(i) for i in watertemp]
    watertemp=np.array(watertemp)
    valid=watertemp>=tlim
    watertemp=watertemp[valid]
    
    waterr063=data[idx][['avis']].values
    waterr063=waterr063.flatten()
    waterr063=[float(i) for i in waterr063]
    waterr063=np.array(waterr063)
    waterr063=waterr063[valid]
    
    waterr390=data[idx][['air']].values
    waterr390=waterr390.flatten()
    waterr390=[float(i) for i in waterr390]
    waterr390=np.array(waterr390)
    waterr390=waterr390[valid]    
    
    waterphase=data[idx][['aphase']].values
    waterphase=waterphase.flatten()
    waterphase=np.array(waterphase)
    waterphase=waterphase[valid]
    
    idx=(data.aboxsize==boxsize) & (data.aphase=='ice')
    icetemp=data[idx][['atemp']].values
    icetemp=icetemp.flatten()
    icetemp=[float(i) for i in icetemp]
    icetemp=np.array(icetemp)
       
    icer063=data[idx][['avis']].values
    icer063=icer063.flatten()
    icer063=[float(i) for i in icer063]
    icer063=np.array(icer063)
    
    icer390=data[idx][['air']].values
    icer390=icer390.flatten()
    icer390=[float(i) for i in icer390]
    icer390=np.array(icer390)
        
    icephase=data[idx][['aphase']].values
    icephase=icephase.flatten()
    icephase=np.array(icephase)
    
    idx=(data.aboxsize==boxsize) & (data.aphase=='mix')
    mixtemp=data[idx][['atemp']].values
    mixtemp=mixtemp.flatten()
    mixtemp=[float(i) for i in mixtemp]
    mixtemp=np.array(mixtemp)
    valid=mixtemp>=tlim
    mixtemp=mixtemp[valid]
    
    mixr063=data[idx][['avis']].values
    mixr063=mixr063.flatten()
    mixr063=[float(i) for i in mixr063]
    mixr063=np.array(mixr063)
    mixr063=mixr063[valid]
    
    mixr390=data[idx][['air']].values
    mixr390=mixr390.flatten()
    mixr390=[float(i) for i in mixr390]
    mixr390=np.array(mixr390)
    mixr390=mixr390[valid]
    
    mixphase=data[idx][['aphase']].values
    mixphase=mixphase.flatten()
    mixphase=np.array(mixphase)
    mixphase=mixphase[valid]

    nwater=np.shape(watertemp)[0]
    nice=np.shape(icetemp)[0]
    nmix=np.shape(mixtemp)[0]
    
    coldwater,coldwatermed,coldwatersd,coldwatern,coldwatermin,coldwatermax= strnan,strnan,strnan,strnan,strnan,strnan #'-9999','-9999','-9999','-9999','-9999','-9999'
    hotice,hoticemed,hoticesd,hoticen,hoticemin,hoticemax= strnan,strnan,strnan,strnan,strnan,strnan #'-9999','-9999','-9999','-9999','-9999','-9999'
    mixtmp,mixmed,mixsd,mixn,mixmin,mixmax= strnan,strnan,strnan,strnan,strnan,strnan #'-9999','-9999','-9999','-9999','-9999','-9999'   
    hotwater,hotwatermed,hotwatersd,hotwatern,hotwatermin,hotwatermax= strnan,strnan,strnan,strnan,strnan,strnan #'-9999','-9999','-9999','-9999','-9999','-9999'
    
    if (nwater>0):
        sortedwater = np.sort(watertemp,axis=None)
        waterset = sortedwater[0:int(coldwaterpercent*len(sortedwater))+1]
        if (np.shape(waterset)[0]==0):
            waterset=sortedwater
        coldwater=str(np.average(waterset))
        coldwatermed=str(np.median(waterset))
        coldwatermin=str(min(waterset))
        coldwatermax=str(max(waterset))
        coldwatersd=str(np.std(waterset))
        coldwatern=str(len(waterset))
        waterset2   = sortedwater[int((1-hotwaterpercent)*len(sortedwater))-1:len(sortedwater)-1]
        if (np.shape(waterset2)[0]==0):
            waterset2=sortedwater
        hotwater=str(np.average(waterset2))
        hotwatermed=str(np.median(waterset2))
        hotwatermin=str(min(waterset2))
        hotwatermax=str(max(waterset2))
        hotwatersd=str(np.std(waterset2))
        hotwatern=str(len(waterset2))
        
    if (nice>0):
        sortedice   = np.sort(icetemp,axis=None)
        iceset   = sortedice[int((1-hoticepercent)*len(sortedice))-1:len(sortedice)-1]        
        if (np.shape(iceset)[0]==0):
            iceset=sortedice
        valid=iceset>=tlim
        iceset=iceset[valid]
        nval=np.shape(iceset)[0]
        if (nval>0):
            hotice=str(np.average(iceset))
            hoticemed=str(np.median(iceset))
            hoticemin=str(min(iceset))
            hoticemax=str(max(iceset))
            hoticesd=str(np.std(iceset))
            hoticen=str(len(iceset))

    if (nmix>0):
        mixset   = np.sort(mixtemp,axis=None)
        mixtmp=str(np.average(mixset))
        mixmed=str(np.median(mixset))
        mixmin=str(min(mixset))
        mixmax=str(max(mixset))
        mixsd=str(np.std(mixset))
        mixn=str(len(mixset))
        
    outstring = sat+','+site+','+season+','+str(year)+','+str(julian)+','+str(hhmmss)+','+boxsize+','+Ntot+','+Nc+','+cloudfraction+','
    outstring = outstring + ClRefAvg+','+ClRefStd+','+ClTmpAvg+','+ClTmpStd+','+NCRefAvg+','+NCRefStd+','+NCTmpAvg+','+NCTmpStd+','
    outstring = outstring + str(coldwatern)+','+coldwatermin+','+coldwatermax+','+coldwatermed+','+coldwater+','+coldwatersd+','
    outstring = outstring + hoticen+','+hoticemin+','+hoticemax+','+hoticemed+','+hotice+','+hoticesd+','
    outstring = outstring + mixn+','+mixmin+','+mixmax+','+mixmed+','+mixtmp+','+mixsd+','
    outstring = outstring + hotwatern+','+hotwatermin+','+hotwatermax+','+hotwatermed+','+hotwater+','+hotwatersd
    outstring = outstring + '\n'
    outfn.writelines(outstring)
    count = count + 1
outfn.close()
os.chdir(curdir)
exit()



        
