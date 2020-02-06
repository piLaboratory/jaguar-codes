# buffers in WGS84 - akde 50 and akde 95

# bernardo niebuhr
# 2019-07-02

# open Python
python

# import modules
import os, fnmatch, imp
import numpy as np
import subprocess
import grass.script as grass
from grass.pygrass.modules.shortcuts import general as g
from grass.pygrass.modules.shortcuts import vector as v
from grass.pygrass.modules.shortcuts import raster as r


#---------------------------------------
# Change to a new mapset

g.mapset(mapset = 'akde', flags = 'c') #-c to create

#---------------------------------------
# Import buffers

# akde 95
folder_path = r'E:\_neojaguardatabase\other_limits_cut_variables\akde\95_per'
os.chdir(folder_path) # Change to this folder

files = os.listdir(folder_path) # List files in the folder
for i in files:
    if i.endswith('shp'): # Select tif files
        print i
        name = 'akde_ind'+i.replace('.shp', '').replace('.', '_').replace('%', '_')
        grass.run_command('v.import', input = i, output = name, overwrite = True) # Import maps

# akde 50
folder_path = r'E:\_neojaguardatabase\other_limits_cut_variables\akde\50_per'
os.chdir(folder_path) # Change to this folder

files = os.listdir(folder_path) # List files in the folder
for i in files:
    if i.endswith('shp'): # Select tif files
        print i
        name = 'akde_ind'+i.replace('.shp', '').replace('.', '_').replace('%', '_')
        grass.run_command('v.import', input = i, output = name, overwrite = True) # Import maps

#------
# cut for buffers

# Variables
# HFP2009_wgs84_1km_neotropic_tif_exp - ok
# population density - ok
# distance from roads - r.grow.distance -m input=sea distance=dist_from_sea_geodetic metric=geodesic - calc
# road density buffer 10km - calc
# PAs - prop - ok
# Hansen - bin 0/1 > 0, calc

# Hansen patch size - calc
# Hansen edge area, core area, edge/core areas
# productivity

maps = ['HFP2009_wgs84_1km_neotropic_tif_exp', 'Livestock_Cattle_CC2006_AD_1km_neotropic_tif_exp',
#'Neotropic_Earthenv_dem90m_wgs84_new', 'Neotropic_Earthenv_slope_from_dem90m_wgs84_new',
'Neotropic_Hansen_percenttreecoverd_2000_wgs84', 'Neotropical_Hansen_treecoverlossperyear_wgs84_2017',
'wc2.0_bio_30s_12_neotropic_tif_exp', 'wc2.0_bio_30s_15_neotropic_tif_exp', 'MOD17A3_Science_NPP_mean_00_15',
'Population_density_gpw_v4_rev10_2015_1km_neotropic_tif_exp',
'protected_areas_2018_bin_tif_exp', 'gROADS_v1_americas_rast']

# years for forest prop
pa = r'E:\_neojaguardatabase\other_limits_cut_variables\akde'
os.chdir(pa)
#arq = np.genfromtxt('HR_year2dcod.csv', skip_header = 1, delimiter = ';')
arq = np.genfromtxt('HR_year_2d.csv', delimiter = ';', dtype = None, skip_header = 1)#, encoding = "utf-8")
years = [arq[index][8] for index in range(len(arq))]
ind_cod = [str(arq[index][0]) for index in range(len(arq))]

# List of buffers
list_buffers = grass.list_grouped('vect', pattern = 'akde*')['akde']

map_for_define_region = 'Neotropic_Hansen_percenttreecoverd_2000_wgs84@PERMANENT'

# output file
names = ['HFP2009', 'livestock_cattle_density', 'net_primary_productivity', 'forest_year', 
        'population_density', 'road_dist', 'road_density_km_per100km2', 'protected_areas',
        'bio12_annual_precipitation', 'bio15_precipitation_seasonality']

pa = r'E:\_neojaguardatabase\other_limits_cut_variables\akde\maps'
os.chdir(pa)

f = open('table_akde.csv', 'a+')
f.write('id;kernel;'+';'.join(names)+'\n')
f.close()

# individuals with bad akde fit
not_run = ['103', '59_2', '60_2', '25_50_']

#### FALTA O INDIVIDUO 122 44_50!!!

# For each buffer
for i in list_buffers[-2:]:
    
    print i
    
    if 'J' in i:
        ind_codmatch = i.split('ind')[1][:-7]
        kern = i.split('ind')[1][-6:-4]
    else:
        ind_codmatch = i.split('ind')[1][:-3]
        kern = i.split('ind')[1][-2:]
    
    if ind_codmatch not in not_run:
        
        which_ind = np.where(np.char.find(ind_cod, ind_codmatch) == 0)[0][0]
        yy = years[which_ind]
        
        # Define region
        grass.run_command('g.region', vect = i, res = '00:00:00.9',
                          align = map_for_define_region, flags = 'p')
        # Use vector as a mask
        grass.run_command('r.mask', vector = i, cat = 2, overwrite = True) # Mask for the buffer for individual i
        
        # Cut maps
        for mm in maps:
            expr = i+'_'+mm+' = '+mm+'@PERMANENT'
            print expr
            r.mapcalc(expr, overwrite = True)
         
        # distance from roads
        r.grow_distance(input = i+'_'+'gROADS_v1_americas_rast', 
            distance = i+'_'+'gROADS_v1_dist_from_roads_meters_exp', metric = 'geodesic', 
            flags = 'm', overwrite = True)
        
        # road density
        size = 9990/30 # 10 km = 333 pixels
        r.neighbors(input = i+'_'+'gROADS_v1_americas_rast', 
            output = i+'_'+'gROADS_v1_road_density_pixels_per_100km2', method = 'sum', 
            size = size, overwrite = True)
        r.mapcalc(i+'_'+'gROADS_v1_road_density_km_per_100km2_exp = '+i+'_'+'gROADS_v1_road_density_pixels_per_100km2 * 30 / 1000', overwrite = True)
        
        # Hansen bin
        r.mapcalc(i+'_Neotropic_Hansen_percenttreecoverd_2000_wgs84_gt0_binary = if(Neotropic_Hansen_percenttreecoverd_2000_wgs84 > 0, 1, 0)', overwrite = True)
        
        # Hansen correspondent year
        expr = i+'_Neotropic_Hansen_percenttreecoverd_2000_wgs84_gt0_binary_year_exp = if('+i+'_Neotropical_Hansen_treecoverlossperyear_wgs84_2017 > 0 && '+ \
        i+'_Neotropical_Hansen_treecoverlossperyear_wgs84_2017 < '+str(yy)+', 0, '+i+'_Neotropic_Hansen_percenttreecoverd_2000_wgs84_gt0_binary)'
        r.mapcalc(expr, overwrite = True)
         
        # NPP ext in the end
        r.mapcalc(i+'_MOD17A3_Science_NPP_mean_00_15_exp = '+i+'_MOD17A3_Science_NPP_mean_00_15', overwrite = True)
        
        # export maps
        maps_exp = grass.list_grouped(type = 'raster', pattern = i+'*exp')['akde']
        names = ['HFP2009', 'livestock_cattle_density', 'net_primary_productivity', 'forest_year', 
        'population_density', 'road_dist', 'road_density_km_per100km2', 'protected_areas',
        'bio12_annual_precipitation', 'bio15_precipitation_seasonality']
        names_exp = [i+'_'+nm for nm in names]
        
        pa = r'E:\_neojaguardatabase\other_limits_cut_variables\akde\maps'
        os.chdir(pa)
        
        for mm in range(len(maps_exp)):
            if 'NPP' in maps_exp[mm]:
                r.out_gdal(input = maps_exp[mm], output = names_exp[mm], 
                    createopt = 'COMPRESS=DEFLATE', overwrite = True, flags = 'f')
            else:
                r.out_gdal(input = maps_exp[mm], output = names_exp[mm], 
                    createopt = 'COMPRESS=DEFLATE', overwrite = True)
        
        # save statistics
        f = open('table_akde.csv', 'a+')
        f.write(ind_codmatch+';'+kern)
        for mm in range(len(maps_exp)):
            var = grass.parse_command('r.univar', map = maps_exp[mm], flags = 'g')
            if len(var) == 0:
                val = 'NA'
            else:
                val = var['mean'].encode("utf-8")
            f.write(';'+val)
        
        f.write('\n')
        f.close()
        
        # Remove mask
        grass.run_command('r.mask', flags = 'r')

# same thing, but once things are calculated, just to write averages in a table file
for i in list_buffers[121:124]:
    
    print i
    
    if 'J' in i:
        ind_codmatch = i.split('ind')[1][:-7]
        kern = i.split('ind')[1][-6:-4]
    else:
        ind_codmatch = i.split('ind')[1][:-3]
        kern = i.split('ind')[1][-2:]
    
    if ind_codmatch not in not_run:
        
        which_ind = np.where(np.char.find(ind_cod, ind_codmatch) == 0)[0][0]
        yy = years[which_ind]
        
        # Define region
        grass.run_command('g.region', vect = i, res = '00:00:00.9',
                          align = map_for_define_region, flags = 'p')
        # Use vector as a mask
        grass.run_command('r.mask', vector = i, cat = 2, overwrite = True) # Mask for the buffer for individual i
        
        # export maps
        maps_exp = grass.list_grouped(type = 'raster', pattern = i+'*exp')['akde']
        names = ['HFP2009', 'livestock_cattle_density', 'net_primary_productivity', 'forest_year', 
        'population_density', 'road_dist', 'road_density_km_per100km2', 'protected_areas',
        'bio12_annual_precipitation', 'bio15_precipitation_seasonality']
        names_exp = [i+'_'+nm for nm in names]
        
        pa = r'E:\_neojaguardatabase\other_limits_cut_variables\akde\maps'
        os.chdir(pa)
        
        for mm in range(len(maps_exp)):
            if 'NPP' in maps_exp[mm]:
                r.out_gdal(input = maps_exp[mm], output = names_exp[mm], 
                    createopt = 'COMPRESS=DEFLATE', overwrite = True, flags = 'f')
            else:
                r.out_gdal(input = maps_exp[mm], output = names_exp[mm], 
                    createopt = 'COMPRESS=DEFLATE', overwrite = True)
        
        # save statistics
        f = open('table_akde.csv', 'a+')
        f.write(ind_codmatch+';'+kern)
        for mm in range(len(maps_exp)):
            var = grass.parse_command('r.univar', map = maps_exp[mm], flags = 'g')
            if len(var) == 0:
                val = 'NA'
            else:
                val = var['mean'].encode("utf-8")
            f.write(';'+val)
        
        f.write('\n')
        f.close()
        
        # Remove mask
        grass.run_command('r.mask', flags = 'r')


# function to define region 2 degrees
def reg_2deg(input = '', size = 2):

    # get region center
    reg = grass.parse_command('g.region', input = input, flags = 'c')
    lat = float()
    lon = float()

    g.region(n = lat + size/2, s = lat - size/2, w = lon - size/2, e = lon + size/2,
        align = input, flags = 'ap')

reg_2deg(input = 'MASK', size = 2)
