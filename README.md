# 2019-exoplanets

python3 is used
From terminal , run the following command. Assumes the script is in the same directory as Swift.
Otherwise, give the absolute path to Swift. The script runs for 7 hours, so make sure not to let the computer
to sleep in "Energy Saver" setting
python read_fits.py ./Swift

The script will produce the following files:
coordinates1.txt, uvotsource1.txt
coordinates2.txt, uvotsource2.txt
coordinates3.txt, uvotsource3.txt
coordinates4.txt, uvotsource4.txt
coordinates5.txt, uvotsource5.txt
coordinates#.txt only contains RA and DEC coordinates, and is used to look up Kepler Target Search webpage.
uvotsource#.txt contains more columns of information for the source.
Line X in coordinates#.txt and line X in uvotsource#.txt are for same source.
Each file has 9000 or less lines, as there is a 10000 line per file restriction in the webpage.

# In a web broswer, go to URL http://archive.stsci.edu/kepler/kepler_fov/search.php?form=fuf
# In "Local File Name (required)", choose coordinates1.txt. Repeat it for coordinates[2-5].txt
# In "Delimiter", choose ","
# In "RA,Target or Data ID Column #", choose "1"
# In "Dec Column #(if used)", choose "2"
# In "File Contents", choose "Coordinates"
# In "Radius (arcmin), put 0.1
# In "Equinox", choose J2000
# In "Sort By:", let "ang_sep(')" be the top choice
# In "Search Output Format", choose "File: comma-separated values"
# In "Search Output Columns", choose the columns in the following order:
# Kepler_ID,RA (J2000),Dec (J2000),GALEX FUV,GALEX NUV,g,r,i,z,J,H,K,Teff,Data Availability,Ang Sep (')
Click "Search" button with the choices above

# When the file is downloaded, run the following command to move the downloaded file to the same directory as the script and rename it.
# Do this when searching coordinates[2-5].txt, and rename them to kepler_fov_search_[2-5].txt  
mv ~/Downloads/kepler_fov_search.txt kepler_fov_search_1.txt

# Run following commands
python find_kics.py kepler_fov_search_1.txt 
python find_kics.py kepler_fov_search_2.txt 
python find_kics.py kepler_fov_search_3.txt 
python find_kics.py kepler_fov_search_4.txt 
python find_kics.py kepler_fov_search_5.txt 
# There will be 5 output files: uvotsource_result_1.csv, uvotsource_result_2.csv,
# uvotsource_result_3.csv, uvotsource_result_4.csv, uvotsource_result_5.csv
# 1st column is Kepler_ID.
