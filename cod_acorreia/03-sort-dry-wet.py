#!/usr/bin/python2.7
# -*- coding: utf-8 -*-
"""
Created on Sun Mar  1 11:40:41 2015

@author: acorreia
"""

# python ./sort-dry-wet.py datadir

import os
import sys
import glob
import calendar

from physconstants import *

datadir=sys.argv[1]
outdry=datadir+'/DRY'
outwet=datadir+'/WET'
outvis=datadir+'/VIS/'
outir=datadir+'/IR/'
outtmp=datadir+'/TMP/'
outphase=datadir+'/phase/'
#outjjl=datadir+'/JJL'
if ( os.path.isdir(outdry)==False ):
    os.makedirs(outdry)
if ( os.path.isdir(outwet)==False ):
    os.makedirs(outwet)
#if ( os.path.isdir(outjjl)==False ):
#    os.makedirs(outjjl)
outdry,outwet=outdry+'/',outwet+'/'
#outjjl=outjjl+'/'
os.chdir(datadir)
flist = sorted(glob.glob(datadir+'/*alldata.csv'))
cflist=sorted(glob.glob(datadir+'/*cloudfraction.csv'))
vislist   = sorted(glob.glob(datadir+'/*r063.png'))
irlist    = sorted(glob.glob(datadir+'/*r390.png'))
tmplist    = sorted(glob.glob(datadir+'/*temp.png'))
phaselist = sorted(glob.glob(datadir+'/*phase.png'))
if len (vislist)>0:
    if ( os.path.isdir(outvis)==False ):
        os.makedirs(outvis)
        count=0
        while count < len (vislist):
            fn = vislist[count]
            filename = fn.split('/')[len(fn.split('/'))-1]
            os.rename(vislist[count],outvis+filename)
            count=count+1
if len (irlist)>0:
    if ( os.path.isdir(outir)==False ):
        os.makedirs(outir)
        count=0
        while count < len (irlist):
            fn = irlist[count]
            filename = fn.split('/')[len(fn.split('/'))-1]
            os.rename(irlist[count],outir+filename)
            count=count+1            
if len (tmplist)>0:
    if ( os.path.isdir(outtmp)==False ):
        os.makedirs(outtmp)
        count=0
        while count < len (tmplist):
            fn = tmplist[count]
            filename = fn.split('/')[len(fn.split('/'))-1]
            os.rename(tmplist[count],outtmp+filename)
            count=count+1            
if len (phaselist)>0:
    if ( os.path.isdir(outphase)==False ):
        os.makedirs(outphase)
        count=0
        while count < len (phaselist):
            fn = phaselist[count]
            filename = fn.split('/')[len(fn.split('/'))-1]
            os.rename(phaselist[count],outphase+filename)
            count=count+1
count=0
while count < len(flist):
    fn = flist[count]
    fn2=cflist[count]
    filename  = fn.split('/')[len(fn.split('/'))-1]
    filename2 = fn2.split('/')[len(fn2.split('/'))-1]
    year     = int(filename.split('.')[1])
    julian   = int(filename.split('.')[2])
    leap=0
    if (calendar.isleap(year)):
        leap=1
    wetend=dryseasonbegin-1+leap
    wetini=dryseasonend+1+leap
    dryend=dryseasonend+leap
    dryini=dryseasonbegin+leap
    #jjlini=152+leap
    #jjlend=212+leap  
    if ((julian<=wetend) or (julian>=wetini)):
        os.rename(filename,outwet+filename)
        os.rename(filename2,outwet+filename2)        
    if ((julian>=dryini) and (julian<=dryend)):
        os.rename(filename,outdry+filename)
        os.rename(filename2,outdry+filename2)        
    #if ((julian>=jjlini) and (julian<=jjlend)):
    #    os.rename(filename,outjjl+filename)
    count=count+1
exit()
