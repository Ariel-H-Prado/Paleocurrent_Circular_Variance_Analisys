Progam to calculate the circular variance from a .kml file

There are two .py codes. Paleocurrent_Circular_Variance_Analisys.py is the code to be used, Circular_Variance_functions.py is auxiliary.

 1. Place the two codes and the .kml file in the same folder.
 2. Input the name of the .kml file on the line 26 of the code.
 3. Run the code.

Outputs:

"Original_name_dadospalc.txt": contains the latitude and longitude extracted from the .kml file.
"Original_name_Streams.txt": contains the latitude, longitude and azimuth direction of the paleostreams (in degrees).
"Circular_Variance_Interp.txt": contains the total circular variance of the interpolated points, without random selection.
"Bootstrap_Circular_Variances_%d.txt": contains the circular variance after random selection.