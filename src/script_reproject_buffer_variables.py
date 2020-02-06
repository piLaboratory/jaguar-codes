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
all_maps = grass.read_command('r.proj', location = 'newLocation_wgs84', mapset = 'variables_cut', flags = 'l').replace('\r', '').split('\n')
all_regions = grass.read_command('v.proj', location = 'newLocation_wgs84', mapset = 'variables_cut', flags = 'l').replace('\r', '').split('\n')

# change location
g.mapset(mapset = "PERMANENT", location = "location_N_zone_12")

# get only those for zone 12
maps = [i for i in all_maps if 'zone_21' in i]
regions = [i for i in all_regions if i.endswith('zone_21') and i.startswith('region')]

# loop to reproject regions
for i in regions:
    v.proj(location = 'newLocation_wgs84', mapset = 'variables_cut', input = i, output = i, overwrite = True)

g.region(vector = (regions))

categorical = ['Drainage', 'Ecoregions', 'Landcover', 'Protected_areas', 
'roads', 'Treecover', 'Tree_plantations', 'Water_presence']

#len(grass.list_grouped(type = "raster", pattern = "*"+"101"+"*")["PERMANENT"])

# loop to reproject
for reg in regions[70:]:
    
    # define region
    g.region(vector = reg, res = 30, flags = 'ap')
    # get prefix for maps
    reg_pref = reg.split("buffer")[1]
    # select maps in this region
    maps_reg = [m for m in maps if reg_pref in m]
    
    # loop for reprojecting maps
    for i in maps_reg:
        
        print i
        
        if any(word in i for word in categorical):
            method = 'nearest'
        else:
            method = 'bilinear'
        
        r.proj(location = 'newLocation_wgs84', mapset = 'variables_cut', input = i, output = i, method = method, resolution = 30, overwrite = True)
        
#-------
# distance to edges, patch size, distance to roads, distance to water bodies

# forest related metrics
maps_forest = grass.list_grouped(type = 'raster', pattern = '*year*')['PERMANENT']
maps_forest = [i for i in maps_forest if ('dist' not in i and 'AreaHA' not in i and 'pid' not in i)]
#maps_forest = [i for i in maps_forest if ('dist' not in i and '99' in i and 'gt80' in i)]

# calculate patch size and distance to forest edges
# distance to forest edges: Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp
lsmetrics_dir = r'E:\_neojaguardatabase\scripts\LS_METRICS\_LSMetrics_v1_0_0'
os.chdir(lsmetrics_dir)
    
# Run LSMetrics
from LSMetrics_v1_0_0 import dist_edge, rulesreclass, patch_size

# distance to forest edges
name_forest = ''
dist_edge(input_maps = maps_forest[73:],
    classify_edge_as_zero = False,
    prepare_biodim = False, remove_trash = True,
    prefix = name_forest, add_counter_name = False, export = False, dirout = '')
    
maps_forest[47]

# patch size
patch_size(input_maps = maps_forest[47:], remove_trash = True, export = False)

# distance to roads
road_maps = grass.list_grouped(type = 'raster', pattern = '*roads*')['PERMANENT']
road_maps = [i for i in road_maps if 'dist' not in i]

for i in road_maps:
    
    g.region(raster = i, flags = 'p')
    r.grow_distance(input = i, distance = i + '_road_dist_m', overwrite = True)

# distance to water
water_maps1 = grass.list_grouped(type = 'raster', pattern = '*Drainage*')['PERMANENT']
water_maps2 = grass.list_grouped(type = 'raster', pattern = '*Water_presence*')['PERMANENT']

water_maps = water_maps1 + water_maps2

for i in water_maps:
    
    print i
    
    g.region(raster = i, flags = 'p')
    if 'Drainage' in i:
        r.mapcalc(i + '_1null = if('+i+' == 1, 1, null())', overwrite = True)
    else:
        r.mapcalc(i + '_1null = if('+i+' == 2, 1, null())', overwrite = True)
    
    r.grow_distance(input = i + '_1null', distance = i + '_water_dist_m', overwrite = True)
    g.remove(type = 'raster', name = i + '_1null', flags = 'f')

#------
# export

output_folder = r'E:\_neojaguardatabase\Buffer70_zones\variables_30m'

# Outpt folder
os.chdir(output_folder)

for i in regions:
    
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


#--------------------------------------------
# Loop for all the other zones 

# list maps in cgs-location
all_maps = grass.read_command('r.proj', location = 'newLocation_wgs84', mapset = 'variables_cut', flags = 'l').replace('\r', '').split('\n')
all_regions = grass.read_command('v.proj', location = 'newLocation_wgs84', mapset = 'variables_cut', flags = 'l').replace('\r', '').split('\n')

# zones
zones = ['zone_15', 'zone_16', 'zone_20', 'zone_21', 'zone_22', 'zone_23']
zone_locations = ['N_zone_15', 'N_zone_16', 'S_zone_20', 'S_zone_21', 'S_zone_22', 'S_zone_23']

zones = ['zone_21', 'zone_22', 'zone_23']
zone_locations = ['S_zone_21', 'S_zone_22', 'S_zone_23']

for zz in range(len(zones)):
    
    # change location
    g.mapset(mapset = "PERMANENT", location = "location_" + zone_locations[zz])
    
    # get only those for zone 12
    maps = [i for i in all_maps if zones[zz] in i]
    regions = [i for i in all_regions if i.endswith(zones[zz]) and i.startswith('region')]
    
    # loop to reproject regions
    for i in regions:
        v.proj(location = 'newLocation_wgs84', mapset = 'variables_cut', input = i, output = i, overwrite = True)
    
    g.region(vector = (regions))
    
    categorical = ['Drainage', 'Ecoregions', 'Landcover', 'Protected_areas', 
    'roads', 'Treecover', 'Tree_plantations', 'Water_presence']
    
    # loop to reproject
    for reg in regions:
        
        # define region
        g.region(vector = reg, res = 30, flags = 'ap')
        # get prefix for maps
        reg_pref = reg.split("buffer_")[1]
        # select maps in this region
        maps_reg = [m for m in maps if reg_pref in m]
        
        # loop for reprojecting maps
        for i in maps_reg:
            
            if any(word in i for word in categorical):
                method = 'nearest'
            else:
                method = 'bilinear'
            
            r.proj(location = 'newLocation_wgs84', mapset = 'variables_cut', input = i, output = i, method = method, resolution = 30, overwrite = True)
            
    #-------
    # distance to edges, patch size, distance to roads, distance to water bodies
    
    # forest related metrics
    maps_forest = grass.list_grouped(type = 'raster', pattern = '*'+zones[zz]+'*year*')['PERMANENT']
    maps_forest = [i for i in maps_forest if 'dist' not in i]
    
    # calculate patch size and distance to forest edges
    # distance to forest edges: Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp
    lsmetrics_dir = r'E:\_neojaguardatabase\scripts\LS_METRICS\_LSMetrics_v1_0_0'
    os.chdir(lsmetrics_dir)
        
    # Run LSMetrics
    from LSMetrics_v1_0_0 import dist_edge, rulesreclass, patch_size
        
    # distance to forest edges
    name_forest = ''
    dist_edge(input_maps = maps_forest,
        classify_edge_as_zero = False,
        prepare_biodim = False, remove_trash = True,
        prefix = name_forest, add_counter_name = False, export = False, dirout = '')
        
    # patch size
    patch_size(input_maps = maps_forest, remove_trash = True, export = False)
    
    # distance to roads
    road_maps = grass.list_grouped(type = 'raster', pattern = '*'+zones[zz]+'*roads*')['PERMANENT']
    road_maps = [i for i in road_maps if 'dist' not in i]
    
    for i in road_maps:
        
        g.region(raster = i, flags = 'p')
        r.grow_distance(input = i, distance = i + '_road_dist_m', overwrite = True)
    
    # distance to water
    water_maps1 = grass.list_grouped(type = 'raster', pattern = '*'+zones[zz]+'*Drainage*')['PERMANENT']
    water_maps2 = grass.list_grouped(type = 'raster', pattern = '*'+zones[zz]+'*Water_presence*')['PERMANENT']
    
    water_maps = water_maps1 + water_maps2
    
    for i in water_maps:
        
        print i
        
        g.region(raster = i, flags = 'p')
        if 'Drainage' in i:
            r.mapcalc(i + '_1null = if('+i+' == 1, 1, null())', overwrite = True)
        else:
            r.mapcalc(i + '_1null = if('+i+' == 2, 1, null())', overwrite = True)
        
        r.grow_distance(input = i + '_1null', distance = i + '_water_dist_m', overwrite = True)
        g.remove(type = 'raster', name = i + '_1null', flags = 'f')
    
    #------
    # export
    
    output_folder = r'E:\_neojaguardatabase\Buffer70_zones\variables_30m'
    
    # Outpt folder
    os.chdir(output_folder)
    
    for i in regions:
        
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

