#!/usr/bin/python2.7
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 18 16:22:06 2014

@author: acorreia
"""

# These are physical constants and definitions needed to calculate GOES reflectances
# Usage in a separate code:
#
# from physconstants import *
#

# Nominal Earth radius in km
ER = 6371.0

# Nominal Geosynchronous orbit (sattelite radius) in km
SR = 42164.17478

# Planck's function:
h  = 6.62606896e-34 #   Planck's constant in  J.s
ls = 299792458     #   light speed in  m/s
k  = 1.3806504e-23  #   Boltzmann constant in J/K

t0 = 0.75  #0.75      # Bidirectional transmission function (at 3.9um) between GOES - cloud tops - GOES, t0=0.75 by Kaufman and Nakajima 1993 
F0 = 9.68465   # Solar descending irradiance at TOA for 3.9um, by Platnick and Fontenla 2008, in W/[m2 um]. This is only indicative, the F0 to be used is convoluted with sensor response function (see goesconstants.py)

# Define VIS/IR thresholds to detect water, ice or mixed phase pixels
cloudVISthreshold  = 0.125  # The limit of 0.125 (for the Amazon) works well when comparing to ISCCP data. Previously it was 0.15   # VIS reflectance above this limit indicates cloud pixel 
cloudTempthreshold = 300.0  # The limit of 300 K (for the Amazon) works well when comparing to ISCCP data. Previously it was 283K   # Temp in Kelvin below this limit indicates cloud pixel
#icedry   = 17.2399  # VIS/IR above this limit indicates ice pixel
#icewet   = 20.2301  # VIS/IR above this limit indicates ice píxel
#waterdry =  4.6801  # VIS/IR below this limit indicates water pixel
#waterwet =  7.3172  # VIS/IR below this limit indicates water pixel
icedry   =  7.6176
icewet   =  9.6324
waterdry =  3.9559
waterwet =  4.2353

# Define lower limit for the temperature of liquid supercooled water or hottest ice in K
# -38.1C = 235.05K
tlim = 235.05

# Define percentage of pixels to be considered as 'hot ice', e.g. 0.10 means 10% hottest ice pixels
hoticepercent = 0.050

# Define percentage of pixels to be considered as 'cold water', e.g. 0.10 means 10% coldest water pixels
coldwaterpercent = 0.050

# Define percentage of pixels to be considered as 'hot water', e.g. 0.10 means 10% hottest water pixels
hotwaterpercent = 0.050

# Definition of DRY Season limits in Julian day, for a non-leap year
#dryseasonbegin = 213   # Aug  1st
dryseasonbegin = 152   # Jun  1st
dryseasonend   = 304   # Oct 31st

# Define local noon in UTC time. Example: in Manaus, time is UTC-4, so LocNoon=160000
LocNoon=160000

# Define 0 deg C in Kelvin
zeroc=273.15

# Definition of NaN (not-a-number) or missing data
mynan=-9999
strnan='-9999'
