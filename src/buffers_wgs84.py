# Cut environmental variables for each individual within the cgs-WGS84 location

# This script runs within GRASS GIS environment, in a location with datum WGS84, RCS lat-lon,
# with all environmental/background data imported into it.

# bernardo niebuhr
# 2019-07-02

# 1st: create UTM zone locations, datum wgs84
# mannualy
#+proj=utm +zone=12 +datum=WGS84 +units=m +no_defs
#+proj=utm +zone=15 +datum=WGS84 +units=m +no_defs
#+proj=utm +zone=16 +datum=WGS84 +units=m +no_defs
#+proj=utm +zone=20 +south +datum=WGS84 +units=m +no_defs
#+proj=utm +zone=21 +south +datum=WGS84 +units=m +no_defs
#+proj=utm +zone=22 +south +datum=WGS84 +units=m +no_defs
#+proj=utm +zone=23 +south +datum=WGS84 +units=m +no_defs

# open Python
python

# import modules
import os, fnmatch, imp
import subprocess
import grass.script as grass
from grass.pygrass.modules.shortcuts import general as g
from grass.pygrass.modules.shortcuts import vector as v
from grass.pygrass.modules.shortcuts import raster as r

#---------------------------------------
# Create a single mapset for each environmental variables to be cut, still in lat-lon, before being reprojected to UTM
# This means that all the data for all correspondent UTm zones will be in a single mapset, in lat-lon
# The next step will be to reproject these
g.mapset(mapset = "variables_cut", flags = "c")

# resturn to PERMANENT
g.mapset(mapset = "PERMANENT")

#---------------------------------------
# Import buffers

folder_path = r'E:\_neojaguardatabase\Buffer70_zones\buffers'
os.chdir(folder_path) # Change to this folder

files = os.listdir(folder_path) # List files in the folder
for i in files:
    if i.endswith('gpkg'): # Select tif files
        print i
        name = i.replace('.gpkg', '')
        v.in_ogr(input = i, output = name, overwrite = True) # Import maps


#------
# cut for buffers

maps = ['drainage_15s_binary_tif_exp', ]

# List of buffers
list_buffers = grass.list_grouped('vect', pattern = '*buffer*')['PERMANENT']

# Map used as the base to align the other variables
map_for_define_region = 'Neotropic_Hansen_percenttreecoverd_2000_wgs84'

# For each buffer/individual
for i in list_buffers:
    
    print i
    
    # define region
    grass.run_command('g.region', vect = i, res = 30,
                      align = map_for_define_region, flags = 'p')
    
    # creates polygon around the region we are working in
    v.in_region(output = 'region_'+i)
    
    # mask
    grass.run_command('r.mask', vector = i) # Mask for the buffer for individual i
    
    # cut variables using the mask
    for mm in maps:
        expr = i+'_'+mm+' = '+mm+'@PERMANENT'
        print expr
        r.mapcalc(expr, overwrite = True)
    
    # Remove mask
    grass.run_command('r.mask', flags = 'r')


# remove all these, they have error 0,1,128
Neotropic_Hansen_forest1_0_95percenttreecover_2000_wgs84_30m_tif
