#!/usr/bin/python2.7
# -*- coding: utf-8 -*-
"""
Created on Sat May  6 16:22:38 2017

@author: acorreia
"""

# List of parameters used when the ice/water look up table was built
# Usage in a separate code:
#
# from lutparams import *
#

# Define steps for VZA and SZA in the LUT file. For example: if the LUT file has VZA=10,20,30, use VZAstep=10; if the LUT file has SZA=0,5,10,15,etc. use SZAstep=5
VZAstep=10
SZAstep=5

#Define total number of Reff values in the ice and water LUTs. For example: if Reff varies from 5 to 20 um, in 5 um steps, total number is 4.
totIceReff=63
totWaterReff=86

# Define general parameters to be selected in the LUT
myPW=6.0   # precipitable water in cm
myCOD=50   # COD=50 : only thick clouds are considered 

# LUT file variables in the header (1 header line exactly)
lutvars=['sat','wl','pw','phase','cod','vza','sza','reff','topdown','l30','l150','r30','r150']