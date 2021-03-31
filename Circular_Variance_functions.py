# -*- coding: utf-8 -*-
# Version 10.06.2019 
# FUNCTIONS: Calculate paleocurrent directions from paleo scrolls data from active rivers

import sys
import numpy as np
import matplotlib.pyplot as plt
from geopy.distance import geodesic
import random
from astropy.stats import circvar
from astropy import units as u


#########################################################
#Function to calculate azimuth between two points with coordinates

def Azimuth_2points(Lat1, Lon1, Lat2, Lon2):
	
	A = (Lat1, Lon1)
		
	B = (Lat2 , Lon2)
		
	dist1 = geodesic(A, B).meters
	
	B2 = (Lat1, Lon2)
	
	dist2 = geodesic(A, B2).meters
	
	if dist2 > dist1:
		
		dist2 = dist1
	
	Theta = np.arccos(dist2/dist1)
	
	if Lat1 <= Lat2 and Lon1 <= Lon2:
		
		Azimuth = np.pi/2 - Theta
		
	elif Lat1 > Lat2 and Lon1 <= Lon2:
		
		Azimuth = np.pi/2 + Theta
		
	elif Lat1 >= Lat2 and Lon1 > Lon2:
		
		Azimuth = 3*np.pi/2 - Theta
		
	elif Lat1 < Lat2 and Lon1 > Lon2:
		
		Azimuth = 3*np.pi/2 + Theta
		
	return (Azimuth)


#########################################################
# Interpolation

def Interp_PLS(Lat, Lon, PLS):	
	
	Dist_Azim = []

	N_Azim = len(PLS)

	Count_Azim = 0

	Count_Salto = 0

	while Count_Azim < N_Azim - 1:

		A = (Lat[Count_Azim], Lon[Count_Azim])

		B = (Lat[Count_Azim + 1] , Lon[Count_Azim + 1])

		Dist = geodesic(A, B).meters

		Dist_Azim = np.append(Dist_Azim, Dist)

		if Count_Salto == 2:

			Count_Salto = 0

			Count_Azim += 2

		else:

			Count_Salto += 1

			Count_Azim += 1


	Perce = np.percentile(Dist_Azim, 2)

	Interp_PLS = []
	Interp_Lat = []
	Interp_Lon = []

	Count_Azim = 0

	Count_Salto = 0

	while Count_Azim < N_Azim - 1:

		A = (Lat[Count_Azim], Lon[Count_Azim])

		B = (Lat[Count_Azim + 1] , Lon[Count_Azim + 1])

		Dist = geodesic(A, B).meters

		Interp_PLS = np.append(Interp_PLS, PLS[Count_Azim])
		Interp_Lat = np.append(Interp_Lat, Lat[Count_Azim])
		Interp_Lon = np.append(Interp_Lon, Lon[Count_Azim])

		N_Int = (Dist // Perce)*1 + 1

		if N_Int == 1:

			N_Int += 1

		New_Pls = np.interp(np.linspace(0,Dist, N_Int), [0, Dist], [PLS[Count_Azim], PLS[Count_Azim + 1]])
		New_Lat = np.interp(np.linspace(0,Dist, N_Int), [0, Dist], [Lat[Count_Azim], Lat[Count_Azim + 1]])
		New_Lon = np.interp(np.linspace(0,Dist, N_Int), [0, Dist], [Lon[Count_Azim], Lon[Count_Azim + 1]])

		New_Pls = np.delete(New_Pls, len(New_Pls) - 1)
		New_Lat = np.delete(New_Lat, len(New_Lat) - 1)
		New_Lon = np.delete(New_Lon, len(New_Lon) - 1)

		New_Pls = np.delete(New_Pls, 0)
		New_Lat = np.delete(New_Lat, 0)
		New_Lon = np.delete(New_Lon, 0)


		Interp_PLS = np.append(Interp_PLS, New_Pls)
		Interp_Lat = np.append(Interp_Lat, New_Lat)
		Interp_Lon = np.append(Interp_Lon, New_Lon)


		if Count_Salto == 2:

			Interp_PLS = np.append(Interp_PLS, PLS[Count_Azim + 1])
			Interp_Lat = np.append(Interp_Lat, Lat[Count_Azim + 1])
			Interp_Lon = np.append(Interp_Lon, Lon[Count_Azim + 1])

			Count_Salto = 0

			Count_Azim += 2

		else:

			Count_Salto += 1

			Count_Azim += 1
	
	return(Interp_Lat, Interp_Lon, Interp_PLS, Perce)

#########################################################
# Subsampling Grid Data

def Subsample_Grid(Lat, Lon, PLS, Grid_step_len, N_Vars):

	N_elem = np.linspace(0, len(PLS) - 1, len(PLS))

	Zero_Lat = min(Lat)
	Zero_Lon = min(Lon)

	Coord_X = []
	Coord_Y = []

	Count_Coord = 0

	while Count_Coord < len(N_elem):

		P_Lat = Lat[Count_Coord]
		P_Lon = Lon[Count_Coord]


		X = geodesic([Zero_Lat, Zero_Lon], [Zero_Lat, P_Lon]).meters

		Y = geodesic([Zero_Lat, Zero_Lon], [P_Lat, Zero_Lon]).meters

		Coord_X = np.append(Coord_X, X)
		Coord_Y = np.append(Coord_Y, Y)

		Count_Coord += 1

	Circ_Var_Ativo_Grid = []
	
	Circ_Var_out = 0

	Count_vars = 0

	while Count_vars < N_Vars:


		X_Grid = []
		Y_Grid = []
		PLS_Grid = []

		Copy_Coord_X = np.copy(Coord_X)
		Copy_Coord_Y = np.copy(Coord_Y)
		Copy_Interp_PLS = np.copy(PLS)

		#print(len(Copy_Coord_X))

		while len(Copy_Coord_X) > 0:

			List_elem_temp = np.linspace(0, len(Copy_Coord_X) - 1, len(Copy_Coord_X))

			Seed = random.randint(0, len(Copy_Coord_X) - 1)

			X_Grid = np.append(X_Grid, Copy_Coord_X[Seed])
			Y_Grid = np.append(Y_Grid, Copy_Coord_Y[Seed])
			PLS_Grid = np.append(PLS_Grid, Copy_Interp_PLS[Seed])

			Temp = [i for i,v in enumerate(Copy_Coord_X) if v > Copy_Coord_X[Seed] - Grid_step_len and v < Copy_Coord_X[Seed] + Grid_step_len]

			Temp2 = [i for i,v in enumerate(Copy_Coord_Y) if v > Copy_Coord_Y[Seed] - Grid_step_len and v < Copy_Coord_Y[Seed] + Grid_step_len]

			Delete_Elem = [value for value in Temp if value in Temp2] 

			Copy_Coord_X = np.delete(Copy_Coord_X, Delete_Elem)
			Copy_Coord_Y = np.delete(Copy_Coord_Y, Delete_Elem)
			Copy_Interp_PLS = np.delete(Copy_Interp_PLS, Delete_Elem)


		Circ_Var_Ativo_Grid.append(circvar(PLS_Grid*u.deg))
		
		Circ_Var_mean = np.median(Circ_Var_Ativo_Grid)
		
		Circ_Var_actual = circvar(PLS_Grid*u.deg)
		
		if abs(Circ_Var_out - Circ_Var_mean) > abs(Circ_Var_actual - Circ_Var_mean):
			
			Circ_Var_out = Circ_Var_actual
			
			X_Grid_out = X_Grid
			Y_Grid_out = Y_Grid
			PLS_Grid_out = PLS_Grid
		
		
		print(Count_vars)
		Count_vars += 1

	return(Circ_Var_Ativo_Grid, X_Grid_out, Y_Grid_out, PLS_Grid_out, Circ_Var_out)