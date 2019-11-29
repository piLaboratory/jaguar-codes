# reproject layers to UTM, calculate areas and distances, export

# bernardo niebuhr
# 2019-07-02

# This script runs within GRASS GIS environment.

# open Python
python

# import modules
import os, fnmatch, imp
import subprocess
import grass.script as grass
from grass.pygrass.modules.shortcuts import general as g
from grass.pygrass.modules.shortcuts import vector as v
from grass.pygrass.modules.shortcuts import raster as r

#------------------------------------------
# Location utm 12N

# list maps in cgs-location
all_maps = grass.read_command('r.proj', location = 'newLocation_wgs84', mapset = 'PERMANENT', flags = 'l').replace('\r', '').split('\n')
all_regions = grass.read_command('v.proj', location = 'newLocation_wgs84', mapset = 'PERMANENT', flags = 'l').replace('\r', '').split('\n')

# get only those for zone 12
maps = [i for i in all_maps if i.endswith('zone_12')]
regions = [i for i in all_regions if i.endswith('zone_12') and i.startswith('region')]

# loop to reproject regions
for i in regions:
    v.proj(location = 'newLocation_wgs84', mapset = 'PERMANENT', input = i, output = i)

g.region(vector = (regions))

# loop to reproject
for i in maps:
    
    if 'forest1_0_95' in i:
        method = 'nearest'
    else:
        method = 'bilinear'
    
    r.proj(location = 'newLocation_wgs84', mapset = 'PERMANENT', input = i, output = i, method = method, resolution = 30)

#-------
# distance to edges
maps_forest = grass.list_grouped(type = 'raster', pattern = 'Neotropic_Hansen_percenttreecover*')['PERMANENT']

for i in maps_forest:
    
    # region
    g.region(raster = i, res = 30, flags = 'ap')
    # binary 0.95
    forest_bin = i.replace('Neotropic_Hansen_percenttreecover', 'Neotropic_Hansen_forest1_0_95percenttreecover')
    expr = forest_bin+' = if('+i+' >= 95, 1, 0)'
    r.mapcalc(expr, overwrite = True)
    
    # calculate patch size and distance to forest edges
    # distance to forest edges: Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp
    lsmetrics_dir = r'F:\_neojaguardatabase\LS_METRICS\_LSMetrics_v1_0_0'
    os.chdir(lsmetrics_dir)
    
    # Run LSMetrics
    from LSMetrics_v1_0_0 import dist_edge, rulesreclass, patch_size
    
    name_forest = 'dist2forestedges_'
    dist_edge(input_maps = [forest_bin],
        classify_edge_as_zero = False,
        prepare_biodim = False, remove_trash = True,
        prefix = name_forest, add_counter_name = False, export = False, dirout = '')
    
    # patch size
    name_patchsize = 'local_'
    patch_size(input_maps = [forest_bin],
        prefix = name_patchsize)
    
    #------
    # export
    
    output_folder = r'F:\_neojaguardatabase\Buffer70_zones\variables_30m'
    
    # Outpt folder
    os.chdir(output_folder)
    
    ind = i.split('buffer')[1]
    buffer_dir = 'ind' + ind
    
    # Create output folder
    if not os.path.exists(buffer_dir):
        os.mkdir(buffer_dir)
        print "Directory " + buffer_dir + " Created "
    else:    
        print "Directory " + buffer_dir + " already exists"
    
    os.chdir(buffer_dir)
    
    # List of raster to be exported
    list_rast_export = grass.list_grouped(type = 'rast', pattern = '*'+ind+'*')['PERMANENT']
    
    # Loop to export all raster
    for j in list_rast_export:
        print j
    
        r.out_gdal(input = j, output = j+'.tif', createopt = "TFW=YES,COMPRESS=DEFLATE",
            overwrite = True)
