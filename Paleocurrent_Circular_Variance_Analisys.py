# -*- coding: utf-8 -*-
# Version 08.02.2021
# Calculate circular variance from file .KML

import sys
import numpy as np
import matplotlib.pyplot as plt
from Circular_Variance_functions import *
import os
from astropy.stats import circvar
from astropy import units as u
from geopy.distance import geodesic
from scipy.interpolate import griddata
import random
from scipy.stats import vonmises
import math
from astropy.stats import median_absolute_deviation


# PART 1
#########################################################
#Files oppening

Dir = os.getcwd()		# Actual Directory

Arquivo='Brahmaputra_Pandu_scrolls_var'

Extension='.kml'

f=Dir+"/"+Arquivo+Extension

saida=Dir+"/"+Arquivo+'_dadospalc.txt'

# Read in the file
with open(f, 'r') as file :
  filedata = file.read()

# Find the word "coordinates" in the .kml file and thus define where the data is to be cut
linecut=[]
lookup = 'coordinates'

with open(f, 'r') as myFile:
    for num, line in enumerate(myFile, 1):
        
	    linecut.extend([num]) # Adds a value to the "linecut" array which is used to cut only the data and not the header
		

# Cut .kml data
with open(f, 'r') as file :
	filedata = file.readlines()[linecut[0]:linecut[1]]
	file.close()
filedata=str(filedata[0]) 

filedata = filedata.replace(',0 ', '\n')
filedata = filedata.replace('\t\t\t\t-', '-')
filedata = filedata.replace(',', ' ')

with open(saida, 'w') as file:
 file.write(filedata)

# PART 2
#########################################################

# Creates latitude (lat),longitude (lon) and azimute arrays with the two columns of the .txt file

Lon,Lat = np.loadtxt(saida, unpack='True') # Creates two arrays with latitude and longitude values (unpack divides the file into columns in arrays)

N_points = len(Lon)

if N_points%5 != 0:
	
	print("ERROR: Number of measures is not a multiple of 5")

count1 = 0

Coord_Azim = []

while count1 < N_points:
	
	count2 = 0
	
	while count2 < 4:
		
		Lat1 = Lat[count1 + count2]
		
		Lon1 = Lon[count1 + count2]
		
		Lat2 = Lat[count1 + count2 + 1]
		
		Lon2 = Lon[count1 + count2 + 1]
		
		Azim = Azimuth_2points(Lat1, Lon1, Lat2, Lon2)
		
		Coord_Azim.append([Lat1, Lon1, np.rad2deg(Azim)])
		
		count2 += 1
	
	count1 += 5
	
np.savetxt(Arquivo+"_Streams"+".txt", Coord_Azim)

#########################################################

Dir = os.getcwd()		# Actual Directory

Arquivo= Arquivo + '_Streams'

Extension='.txt'

f=Dir+"/"+Arquivo+Extension

Lat,Lon,PLS = np.loadtxt(f, unpack='True')

# PART 3
#########################################################
# Interpolation

[Interp_Lat, Interp_Lon, Interp_PLS, Perce] = Interp_PLS(Lat, Lon, PLS)

# PART 4
#########################################################
# Grid sampling

N_Vars = 1000				# Subsamples Numbers

Grid_step_len_1 = Perce*5.0	# Meters

#[Circ_Var_Ativo_Grid, X_Grid, Y_Grid, PLS_Grid, Circ_Var_out] = Subsample_Grid(Lat, Lon, PLS, Grid_step_len, N_Vars)
[Circ_Var_Ativo_Grid_1, X_Grid_1, Y_Grid_1, PLS_Grid_1, Circ_Var_out_1] = Subsample_Grid(Interp_Lat, Interp_Lon, Interp_PLS, Grid_step_len_1, N_Vars)

Circ_Var_Ativo = []

Circ_Var_Ativo.append(circvar(PLS*u.deg))

# PART 5
#########################################################

Circ_Var_Ativo_Interp = []

Circ_Var_Ativo_Interp.append(circvar(Interp_PLS*u.deg))

print("Circular Variance with Interpolation (2th Percetile)")
print(Circ_Var_Ativo_Interp)
print("Total number of flow directions : %d" %len(Interp_PLS))

np.savetxt("Circular_Variance_Interp.txt", Circ_Var_Ativo_Interp)

#########################################################
# Rosette angular distribution

PalcRose=np.copy(PLS_Grid_1)

PalcRose[(PalcRose<0)]=360+PalcRose[(PalcRose<0)]

n_numbers = 100
bins_number = 50  # the [0, 360) interval will be subdivided into this number of equal bins
bins = np.linspace(0.0, 2 * np.pi, bins_number + 1)
angles = PalcRose*np.pi/180 #2 * np.pi * np.random.rand(n_numbers)
n, _, _ = plt.hist(angles, bins)

plt.clf()
width = 2 * np.pi / bins_number
ax = plt.subplot(1, 1, 1, projection='polar')
ax.set_theta_offset(np.pi/2) 
ax.set_theta_direction(-1) 
bars = ax.bar(bins[:bins_number], n, width=width, bottom=0.0)


for bar in bars:
    bar.set_alpha(0.5)

plt.savefig("Rosette of the Variancia Circular Median Grid Subsample %d meters.png" %Grid_step_len_1)

plt.clf()

MAD_1 = median_absolute_deviation(Circ_Var_Ativo_Grid_1)
plt.hist(Circ_Var_Ativo_Grid_1)
plt.title("Bootstrap %d samples; Median Circular Variance : %f \n MAD : %f" %(N_Vars,Circ_Var_out_1, MAD_1))
plt.savefig("Histogram Variance Circular Grid %d meters.png" %Grid_step_len_1)

np.savetxt("Bootstrap_Circular_Variances_%d.txt" %Grid_step_len_1, Circ_Var_Ativo_Grid_1)

####################################################################

