#!/usr/bin/python2.7
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 25 16:38:33 2015

@author: acorreia
"""

import os
#######################################
#
# CONFIGURE HERE WHAT TO RUN
#
######################################

# Define the sites to be used in the calculations. Other sites can be defined in the 01-goes-to-pixeldata.py code 
# sitelist=['AF','AH','BT','CB','CG','JP','MN','RB']
sitelist=['AF','AH','CB','JP','RB']

# Define the seasons to be run. Season definition dates in physconstants.py
seasonlist=['DRY','WET']

# Define the half size of a box around each site location. Nominal pixel size is 4km, so a box size of 12 pixels means ~96 x 96 km^2 box.
boxsizelist=['12']

# Define if thermodynamic phase plots, reflectance and brightness temperature plots are to be saved.
phaseopt='No'

# Define the folder path to where the data is located. This folder should contain subfolders named as '2010', '2011', '2012', etc.
# Each subfolder should contain all the NETCDF GOES data to be analysed.  
#dataroot='/media/acorreia/6d49ef17-ce96-4c84-9aea-e0b790bfafd1/'
dataroot='/home/acorreia/Desktop/phase/data/test-may06-2017/'

# Define the start and end years to be used in the calculations. These are subfolders in the dataroot folder.
startyear = 2010
endyear   = 2012

# Define where the results will be stored
#outroot='/home/acorreia/Desktop/phase/output/18z-isccp-hotwater/'
outroot=dataroot

# Define the name of the output file. This should be only the "rootname", like 'analysis-'. The final filenames will be like 'analysis-2010.csv', 'analysis-2011.csv', etc.
outfileroot=outroot+'allsites-isscp-hotwater-'

# Define the name for a 2D histogram to be plotted at the end
png2dhistogram=outroot+'Temp-Reff.png'

# Define the location of the scripts that are part of the processing chain
# Don't change this unless you know what you are doing.
pexec='/usr/bin/python2.7'
codedir='/home/acorreia/Desktop/phase/code/production5-reff/'
p1=codedir+'01-goes-to-pixeldata.py'
p2=codedir+'02-compile-data.py'
p3=codedir+'03-sort-dry-wet.py'
p4=codedir+'04-compile-stats2.py'
p5=codedir+'05-coldwater-hotice3.py'
p6=codedir+'06-joincsv.py'
p7=codedir+'07-prepare-reff3.py'
p8=codedir+'08-compile-reff.py'
p9=codedir+'09-calculate-reff.py'
mylut=codedir+'lutfile'


##########################################
# DON'T MESS ANYTHING BELOW THIS POINT
##########################################
for fdir in range(startyear, endyear+1):
    datadir=dataroot+str(fdir)
    outdir=outroot+str(fdir)
    outfilename=outfileroot+str(fdir)+'.csv'

    # Check if output dirs exist, create them otherwise
    i=0
    while (i<len(sitelist)):
        sitedir=outdir+'/'+sitelist[i]
        if ( os.path.isdir(sitedir)==False ):
            os.makedirs(sitedir)
        j=0
        while (j<len(seasonlist)):
            seasondir=sitedir+'/'+seasonlist[j]
            if ( os.path.isdir(seasondir)==False ):
                os.makedirs(seasondir)
            j=j+1
        i=i+1
    
    # Loop sitelist, seasonlist, boxsizelist
    i=0
    while (i<len(sitelist)):
        j=0
        while (j<len(boxsizelist)):
            cmd1=pexec+' '+p1+' '+datadir+' '+outdir+'/'+sitelist[i]+' '+sitelist[i]+' '+boxsizelist[j]+' phase='+phaseopt
            cmd2=pexec+' '+p2+' '+outdir+'/'+sitelist[i]+' '+outdir+'/'+sitelist[i]
            cmd3=pexec+' '+p3+' '+outdir+'/'+sitelist[i]
            print()
            print(cmd1)
            os.system(cmd1)
            print()
            print(cmd2)
            os.system(cmd2)
            print()
            print(cmd3)
            os.system(cmd3)        
            k=0
            while (k<len(seasonlist)):            
                cmd4=pexec+' '+p4+' '+outdir+'/'+sitelist[i]+'/'+seasonlist[k]+' '+outdir+'/'+sitelist[i]+'/'+seasonlist[k]+' '+seasonlist[k]+' '+sitelist[i]
                cmd5=pexec+' '+p5+' '+outdir+'/'+sitelist[i]+'/'+seasonlist[k]+' '+outdir+' '+seasonlist[k]+'-'+sitelist[i]+'.csv '+seasonlist[k]+' '+sitelist[i]                     
                print()            
                print(cmd4)     
                os.system(cmd4)
                print()
                print(cmd5)
                os.system(cmd5)                       
                k=k+1     
            j=j+1    
        i=i+1
    
    cmd6=pexec+' '+p6+' '+outdir+' '+outfilename
    print()
    print(cmd6)
    os.system(cmd6)
    
    # Loop sitelist, seasonlist
    i=0
    while (i<len(sitelist)):
        k=0
        while (k<len(seasonlist)):
            cmd7=pexec+' '+p7+' '+outdir+'/'+sitelist[i]+'/'+seasonlist[k]+' '+outfilename+' '+outroot+sitelist[i]+'-'+seasonlist[k]+'-reffready.csv'
            print()
            print(cmd7)
            os.system(cmd7)                         
            k=k+1
        i=i+1

cmd8=pexec+' '+p8+' '+outroot+' '+outfileroot+'reffready.csv'
cmd9=pexec+' '+p9+' '+outfileroot+'reffready.csv '+mylut+' '+png2dhistogram
print()
print(cmd8)
os.system(cmd8)
print()
print(cmd9)
os.system(cmd9)       
    
exit() 
    
