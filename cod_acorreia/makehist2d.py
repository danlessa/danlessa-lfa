#!/usr/bin/python2.7
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 17 14:57:15 2014

@author: acorreia
"""
import numpy as np
import matplotlib.pyplot as plt

def hist2d(reff,temp,savename=None):
    """Makes a 2d histogram of temperature in Celsius vs effective radius in um.
    'reff' and 'temp' must be numpy arrays. If 'savename' is specified, a png
    file with the histogram will be saved. 

    @type reff: numpy array
    @param reff: effective radius data array
    @type temp: numpy array
    @param temp: temperature in C data array
    @type savename: string
    @param savename: name of a optional png output file
        
    """
    # Define graph limits
    nzlims = [10.0, -80.0] # inverted temperature limits
    rlims = [0.0,30.0]     # reff (micrometers) limits
    
    # Find the min/max of the data to define the aspect ratio for the plot
    nzmax = max(nzlims)
    nzmin = min(nzlims)
    rmax  = max (rlims)
    rmin  = min (rlims)
    aspectratioYZ  = (1.0*rmax-1.0*rmin)/(1.0*nzmax-1.0*nzmin)
         
    # Set up x, y and colorbar labels
    # rlabel = '$\mathrm{Reflectance\\ at\\ 0.63\\ \mu m}$'
    # nzlabel = '$\mathrm{Reflectance\\ at\\ 3.90\\ \mu m}$'
    rlabel = '$\mathrm{Effective\\ radius\\ (\mu m)}$'  # X-Axis
    nzlabel = '$\mathrm{Temperature\\ (\degree C)}$'    # Y-Axis
    hlabel = '$\mathrm{Relative\\ frequency\\ ( \% )}$' # Colorbar 
    
    # Define the locations for the axes and the geometry of the plot
    left, width = 0.1, 0.85
    bottom, height = 0.1, 0.85
    rect_temperature = [left, bottom, width, height] 
     
    # Define the number of bins
    #nybins =  50 #  50 is fine for 'small' statistics
    #nzbins =  50 #  50 is fine for 'small' statistics
    nybins = 100 # 100 is fine for 'large' statistics
    nzbins = 100 # 100 is fine for 'large' statistics
     
    # Build bin arrays
    zbins = np.linspace(start = nzmin, stop = nzmax, num = nzbins)
    rbins = np.linspace(start = rmin, stop = rmax, num = nybins)
 
    # Calculate the histogram and normalize it    
    nHYZ, nxyzedges,nyyzedges = np.histogram2d(temp,reff,bins=(zbins,rbins),normed=True)
    nHYZ=100.0*nHYZ/(1.0*np.sum(np.abs(nHYZ)))
    
    # Set up the colorbar
    plt.jet()
    colorbarmin=0.0
    #colorbarmax=0.7 # 0.7 is fine for 'small' statistics with nybins,nzbins=50
    colorbarmax=0.2 # 0.2 is fine for 'large' statistics with nybins,nzbins=100
    #colorbarticknum=7 # 7 is fine for 'small' statistics with nybins,nzbins=50
    colorbarticknum=10 # 10 is fine for 'large' statistics with nybins,nzbins=100
    nHYZ = nHYZ.clip(min=colorbarmin,max=colorbarmax)
    nHYZ[0,1] = colorbarmin
    nHYZ[0,0] = colorbarmax
    
    # Set up the size of the figure and add axes
    fig = plt.figure(3, figsize=(10.0,8.0))
    axTemperature = plt.axes(rect_temperature) 
    
    # Plot the temperature data
    cax = (axTemperature.imshow(nHYZ, extent=[rmin,rmax,nzmax,nzmin], interpolation='None', origin='upper',aspect=aspectratioYZ))
    
    # Configure colorbar
    cbaxes = fig.add_axes([.9, 0.1, 0.03, 0.85])
    cbar=plt.colorbar(cax,cax=cbaxes)
    cbar.ax.tick_params(labelsize=10)
    cbar.ax.set_title(hlabel,rotation='vertical',position=(-0.6,0.7),fontsize=14)
    cbar.set_clim(colorbarmin,colorbarmax)
    cbar.set_ticks(np.linspace(colorbarmin,colorbarmax,num=colorbarticknum+1))
    
    # Plot the axes labels
    axTemperature.set_xlabel(rlabel,fontsize=16)
    axTemperature.set_ylabel(nzlabel,fontsize=16)
    #start, end = axTemperature.get_xlim()
    #axTemperature.xaxis.set_ticks(np.arange(start, end, 1))
     
    # Make the ticklabels pretty
    ticklabels = axTemperature.get_xticklabels()
    for label in ticklabels:
        label.set_fontsize(14)
        label.set_family('serif')
    ticklabels = axTemperature.get_yticklabels()
    for label in ticklabels:
        label.set_fontsize(14)
        label.set_family('serif')
     
    #Set up the plot limits
    axTemperature.set_xlim(rlims)
    axTemperature.set_ylim(nzlims)
     
    #Show the plot
    plt.draw()
     
    # Save to a File
    if savename :
        plt.savefig(savename,format = 'png', transparent=False)
    
    return
    