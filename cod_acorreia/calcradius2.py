#!/usr/bin/python2.7
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 17 14:57:15 2014

@author: acorreia
"""

import numpy as np
import pandas as pd
from physconstants import *
from lutparams import *


def calcreff(myfile,mylut):
    # Read the LUT file
    data = pd.read_csv(mylut,index_col=None,header=0,names=lutvars)
    lutpw=np.array(data.pw).astype(np.float64)
    lutphase=np.array(data.phase).astype(np.str)
    lutcod=np.array(data.cod).astype(np.float64)
    lutvza=np.array(data.vza).astype(np.float64)
    lutsza=np.array(data.sza).astype(np.float64)
    lutreff=np.array(data.reff).astype(np.float64)
    lutr30=np.array(data.r30).astype(np.float64)
    lutr150=np.array(data.r150).astype(np.float64)

    # Read the data file    
    data=pd.read_csv(myfile,header=0,names=['asat','awavelength','asite','aseason','ayear','ajulian','ahhmmss','asza','avza','avis','air','atemp','aphase','atcw','athi'])
    hhmmss=np.array(data.ahhmmss).astype(np.float64)
    vza=np.array(data.avza).astype(np.float64)
    sza=np.array(data.asza).astype(np.float64)
    ir=np.array(data.air).astype(np.float64)
    temp=np.array(data.atemp).astype(np.float64)
    phase=np.array(data.aphase).astype(np.str)
    tcw=np.array(data.atcw).astype(np.float64)
    thi=np.array(data.athi).astype(np.float64)
    
    # Initialize effective radius  
    reff=np.zeros(len(vza))+mynan
    reffwater=np.zeros(len(vza))+mynan
    reffice=np.zeros(len(vza))+mynan

    #Initialize waterfrac and icefrac
    icefrac   = np.zeros(len(vza))
    waterfrac = np.zeros(len(vza))
    waterfrac[phase=='water']=1.
    waterfrac[phase=='ice']=0.
    mixid = (phase=='mix')
    if np.any(mixid):
        waterfrac[mixid]=1.0#was 0.5
    wid = (phase=='mix') & (tcw>0) & (temp>0) & (temp>tcw) 
    if np.any(wid):
        waterfrac[wid]=1.
    iid = (phase=='mix') & (thi>0) & (temp>0) & (temp<thi)
    if np.any(iid):
        waterfrac[iid]=0.
    mid = (phase=='mix') & (thi>0) & (thi<zeroc) & (tcw>0) & (temp>0) & (temp>=thi) #& (temp<=tcw)
    if np.any(mid):
        waterfrac[mid]=(temp[mid]-thi[mid])/(tcw[mid]-thi[mid])
    icefrac=1.-waterfrac
    icefrac[icefrac<0]=mynan

    # Select ICE LUT data for a given COD and PW     
    sel=(lutpw==myPW) & (lutcod==myCOD) & (lutphase=='ice')

    # for loop iterating through the ICE LUT
    for id in range(0,len(lutvza[sel]),totIceReff):
        vtop,vbot=lutvza[sel][id]+VZAstep/2.,lutvza[sel][id]-VZAstep/2.
        sztop,szbot=lutsza[sel][id]+SZAstep/2.,lutsza[sel][id]-SZAstep/2.
        if lutvza[sel][id]<=10:
            vbot=0
        if lutsza[sel][id]==0:
            szbot=0
        if lutsza[sel][id]==lutsza[sel][len(lutsza[sel])-1] :
            sztop=lutsza[sel][len(lutsza[sel])-1]

        selectAM = (vza>=vbot) & (vza<vtop) & (sza>=szbot) & (sza<sztop) & (hhmmss<=LocNoon) & ((phase=='ice') | (phase=='mix'))
        selectPM = (vza>=vbot) & (vza<vtop) & (sza>=szbot) & (sza<sztop) & (hhmmss>LocNoon)  & ((phase=='ice') | (phase=='mix')) 

        if np.any(selectAM):
            yAM=ir[selectAM]
            for k in range(id,id+totIceReff):
                if (k-id)==0:
                    smin=yAM>lutr150[sel][k]
                    if np.any(smin):
                        sminid=(selectAM) & (ir > lutr150[sel][k])
                        reffice[sminid]=lutreff[sel][k]-((ir[sminid]-lutr150[sel][k])/(lutr150[sel][k]-lutr150[sel][k+1])*(lutreff[sel][k+1]-lutreff[sel][k]))
                if k!=id+totIceReff-1:
                    smid=(yAM>lutr150[sel][k+1])&(yAM<=lutr150[sel][k])
                    if np.any(smid):
                        smidid=(selectAM)&(ir>lutr150[sel][k+1])&(ir<=lutr150[sel][k])
                        reffice[smidid]=((ir[smidid]-lutr150[sel][k+1])/(lutr150[sel][k]-lutr150[sel][k+1])*(lutreff[sel][k+1]-lutreff[sel][k])) + lutreff[sel][k]
                else:
                    smax=yAM<lutr150[sel][k]
                    if np.any(smax):
                        smaxid=(selectAM)&(ir<lutr150[sel][k])
                        reffice[smaxid]=((lutr150[sel][k]-ir[smaxid])/(lutr150[sel][k-1]-lutr150[sel][k])*(lutreff[sel][k]-lutreff[sel][k-1])) + lutreff[sel][k]

        if np.any(selectPM):
            yPM=ir[selectPM]
            for k in range(id,id+totIceReff):
                if (k-id)==0:
                    smin=yPM>lutr30[sel][k]
                    if np.any(smin):
                        sminid=(selectPM) & (ir > lutr30[sel][k])
                        reffice[sminid]=lutreff[sel][k]-((ir[sminid]-lutr30[sel][k])/(lutr30[sel][k]-lutr30[sel][k+1])*(lutreff[sel][k+1]-lutreff[sel][k]))
                if k!=id+totIceReff-1:
                    smid=(yPM>lutr30[sel][k+1])&(yPM<=lutr30[sel][k])
                    if np.any(smid):
                        smidid=(selectPM)&(ir>lutr30[sel][k+1])&(ir<=lutr30[sel][k])
                        reffice[smidid]=((ir[smidid]-lutr30[sel][k+1])/(lutr30[sel][k]-lutr30[sel][k+1])*(lutreff[sel][k+1]-lutreff[sel][k])) + lutreff[sel][k]
                else:
                    smax=yPM<lutr30[sel][k]
                    if np.any(smax):
                        smaxid=(selectPM)&(ir<lutr30[sel][k])
                        reffice[smaxid]=((lutr30[sel][k]-ir[smaxid])/(lutr30[sel][k-1]-lutr30[sel][k])*(lutreff[sel][k]-lutreff[sel][k-1])) + lutreff[sel][k]
                                                                           
    # Select WATER LUT data for a given COD and PW     
    sel=(lutpw==myPW) & (lutcod==myCOD) & (lutphase=='water')

    # for loop iterating through the WATER LUT
    for id in range(0,len(lutvza[sel]),totWaterReff):
        vtop,vbot=lutvza[sel][id]+VZAstep/2.,lutvza[sel][id]-VZAstep/2.
        sztop,szbot=lutsza[sel][id]+SZAstep/2.,lutsza[sel][id]-SZAstep/2.
        if lutvza[sel][id]<=10:
            vbot=0
        if lutsza[sel][id]==0:
            szbot=0
        if lutsza[sel][id]==lutsza[sel][len(lutsza[sel])-1] :
            sztop=lutsza[sel][len(lutsza[sel])-1]

        selectAM = (vza>=vbot) & (vza<vtop) & (sza>=szbot) & (sza<sztop) & (hhmmss<=LocNoon) & ((phase=='water') | (phase=='mix'))
        selectPM = (vza>=vbot) & (vza<vtop) & (sza>=szbot) & (sza<sztop) & (hhmmss>LocNoon)  & ((phase=='water') | (phase=='mix')) 

        if np.any(selectAM):
            yAM=ir[selectAM]
            for k in range(id,id+totWaterReff):
                if (k-id)==0:
                    smin=yAM>lutr150[sel][k]
                    if np.any(smin):
                        sminid=(selectAM) & (ir > lutr150[sel][k])
                        reffwater[sminid]=lutreff[sel][k]-((ir[sminid]-lutr150[sel][k])/(lutr150[sel][k]-lutr150[sel][k+1])*(lutreff[sel][k+1]-lutreff[sel][k]))
                                                
                if k!=id+totWaterReff-1:
                    smid=(yAM>lutr150[sel][k+1])&(yAM<=lutr150[sel][k])
                    if np.any(smid):
                        smidid=(selectAM)&(ir>lutr150[sel][k+1])&(ir<=lutr150[sel][k])
                        reffwater[smidid]=((ir[smidid]-lutr150[sel][k+1])/(lutr150[sel][k]-lutr150[sel][k+1])*(lutreff[sel][k+1]-lutreff[sel][k])) + lutreff[sel][k]
                else:
                    smax=yAM<lutr150[sel][k]
                    if np.any(smax):
                        smaxid=(selectAM)&(ir<lutr150[sel][k])
                        reffwater[smaxid]=((lutr150[sel][k]-ir[smaxid])/(lutr150[sel][k-1]-lutr150[sel][k])*(lutreff[sel][k]-lutreff[sel][k-1])) + lutreff[sel][k]

        if np.any(selectPM):
            yPM=ir[selectPM]
            for k in range(id,id+totWaterReff):
                if (k-id)==0:
                    smin=yPM>lutr30[sel][k]
                    if np.any(smin):
                        sminid=(selectPM) & (ir > lutr30[sel][k])
                        reffwater[sminid]=lutreff[sel][k]-((ir[sminid]-lutr30[sel][k])/(lutr30[sel][k]-lutr30[sel][k+1])*(lutreff[sel][k+1]-lutreff[sel][k]))
                if k!=id+totWaterReff-1:
                    smid=(yPM>lutr30[sel][k+1])&(yPM<=lutr30[sel][k])
                    if np.any(smid):
                        smidid=(selectPM)&(ir>lutr30[sel][k+1])&(ir<=lutr30[sel][k])
                        reffwater[smidid]=((ir[smidid]-lutr30[sel][k+1])/(lutr30[sel][k]-lutr30[sel][k+1])*(lutreff[sel][k+1]-lutreff[sel][k])) + lutreff[sel][k]
                else:
                    smax=yPM<lutr30[sel][k]
                    if np.any(smax):
                        smaxid=(selectPM)&(ir<lutr30[sel][k])
                        reffwater[smaxid]=((lutr30[sel][k]-ir[smaxid])/(lutr30[sel][k-1]-lutr30[sel][k])*(lutreff[sel][k]-lutreff[sel][k-1])) + lutreff[sel][k]
                                      
    reff = (waterfrac*reffwater) + (icefrac*reffice)
    reff[reff<=0]=mynan
       
    return reff

