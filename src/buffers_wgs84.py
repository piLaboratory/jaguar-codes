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
import numpy as np
import grass.script as grass
from grass.pygrass.modules.shortcuts import general as g
from grass.pygrass.modules.shortcuts import vector as v
from grass.pygrass.modules.shortcuts import raster as r

#---------------------------------------
# Create a single mapset for each environmental variables to be cut, still in lat-lon, before being reprojected to UTM
# This means that all the data for all correspondent UTm zones will be in a single mapset, in lat-lon
# The next step will be to reproject these
g.mapset(mapset = "variables_cut", flags = "c")

#---------------------------------------
# Import buffers
# folder_path = r'E:\_neojaguardatabase\Buffer70_zones\buffers'
# os.chdir(folder_path) # Change to this folder
# 
# files = os.listdir(folder_path) # List files in the folder
# for i in files:
#     if i.endswith('gpkg'): # Select tif files
#         print i
#         name = i.replace('.gpkg', '')
#         v.in_ogr(input = i, output = name, overwrite = True) # Import maps

#------
# cut for buffers

# variables to be exported
pa = r'E:\_neojaguardatabase'
os.chdir(pa)
arq_maps = np.genfromtxt('Environmental spatial data_UTM_2019_teste.csv', delimiter = ';', dtype = None, skip_header = 1)
maps_to_export = [arq_maps[index][0] for index in range(len(arq_maps))]
names_maps = [arq_maps[index][1] for index in range(len(arq_maps))]

# years for forest prop
arq = np.genfromtxt('YEARS_JAGUAR_FINAL.csv', delimiter = ';', dtype = None, skip_header = 1)
years = [arq[index][4] for index in range(len(arq))]
ind_cod = [arq[index][0] for index in range(len(arq))]

# List of buffers
list_buffers = grass.list_grouped('vect', pattern = '*buffer*')['variables_cut']

# Map used as the base to align the other variables
map_for_define_region = 'Neotropic_Hansen_percenttreecoverd_2000_wgs84@PERMANENT'

# For each buffer/individual
for i in list_buffers:
    
    print i
    
    ind_codmatch = i.split('buffer_')[1][:-8]
    which_ind = np.where(np.char.find(ind_cod, ind_codmatch) == 0)[0][0]
    yy = years[which_ind]
    
    # define region
    grass.run_command('g.region', vect = i, res = '0:00:00.9',
                      align = map_for_define_region, flags = 'p')
    
    # creates polygon around the region we are working in
    v.in_region(output = 'region_'+i, overwrite = True)
    
    # mask
    grass.run_command('r.mask', vector = i) # Mask for the buffer for individual i
    
    # cut variables using the mask
    counter = 0
    for mm in maps_to_export:
        expr = i+'_'+names_maps[counter]+' = '+mm+'@PERMANENT'
        print expr, counter
        r.mapcalc(expr, overwrite = True)
        counter = counter + 1
    
    # Loop for Hansen
    
    # thresholds for binary values of natural vegetation
    thresholds = [0, 30, 50, 80]
    
    # loop to cut for each one and account for deforestation
    for tr in thresholds:
        # Hansen bin
        r.mapcalc(i+'_Treecover_Hansen_2000_30m_gt'+str(tr)+'_binary = if('+i+'_Treecover_Hansen_2000_30m > '+str(tr)+', 1, 0)', overwrite = True)
        # Hansen correspondent year
        expr = i+'_Treecover_Hansen_2000_30m_gt'+str(tr)+'_binary_year'+str(yy)+' = if('+i+'_Treecoverloss_Hansen_2017_30m > 0 && '+ \
        i+'_Treecoverloss_Hansen_2017_30m < '+str(yy)+', 0, '+i+'_Treecover_Hansen_2000_30m_gt'+str(tr)+'_binary)'
        try:
            r.mapcalc(expr, overwrite = True)
        except:
            pass
    
    # Remove mask
    grass.run_command('r.mask', flags = 'r')




# remove all these, they have error 0,1,128
Neotropic_Hansen_forest1_0_95percenttreecover_2000_wgs84_30m_tif

g.remove raster pattern=* -f
