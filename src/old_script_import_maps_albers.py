#---------------------------------------------------------------------------------------
"""
 Neo Jaguar Environmental Database - Albers

 Script to import, re-project, organize and cut the spatial variables for analysis of the
 Neotropical Jaguar movement database
 
 Bernardo B. S. Niebuhr - bernardo_brandaum@yahoo.com.br
 Julia E. F. Oshima - juliaoshima@yahoo.com.br
 Vanesa Bejarano - vanesa.bejarano@gmail.com
 
 Laboratorio de Ecologia Espacial e Conservacao
 Universidade Estadual Paulista - UNESP
 Rio Claro - SP - Brasil
 
 License GPLv2 - feel free to use, modify, and share, but cite us and make your code free.
"""
#---------------------------------------------------------------------------------------

#---------------------------------------
# Prepare the environment

# Call python environment
python

# Load modules
import os
import grass.script as grass
import subprocess

# Define region - the Neotropical region
# Here we use the tree cover map from Hansen, already cut to the Neotropics, as a layer to define the region
map_for_define_region = 'Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp'

#---------------------------------------
# Import variables

# 1. Water Frequency 2000 - 30m

# Import map
folder_path = r'H:\_neojaguardatabase\Envdatabase\30m\Neotropic\Water frequency\2010'
os.chdir(folder_path) # Change to this folder
files = os.listdir(folder_path) # List files in the folder
for i in files:
    if i[-3:] == 'tif': # Select tif files
        print i
        name = i.replace('.tif', '_rast')  
        grass.run_command('r.import', input = i, output = name, overwrite = True) # Import maps

# Mosaic of water frequency maps

# List of maps
maps_water = grass.list_grouped('rast', pattern = 'p*2010_rast')['PERMANENT']

# Region of study
grass.run_command('g.region', rast = map_for_define_region, res = 30, flags = 'ap')

# Combine maps
water_map_mosaic = 'water_frequency_2010_30m_tif_exp'
grass.run_command('r.patch', input = maps_water, output = water_map_mosaic, overwrite = True)

# Delete input maps
grass.run_command('g.remove', type = 'raster', pattern = 'p*2010_rast', flags = 'f')

# 2. Ecoregions 2017 - vector

# Import vector
folder_path = r'H:\_neojaguardatabase\Envdatabase\Vetores\World\Ecoregions'
os.chdir(folder_path) # Change to this folder
grass.run_command('v.import', input = 'Ecoregions2017.shp', output = 'Ecoregions2017_shp', overwrite = True)

# Region of study
grass.run_command('g.region', rast = map_for_define_region, res = 30, flags = 'ap')

# Rasterize using biome number, that indicates biome type (column BIOME_NUM) - 30m
grass.run_command('v.to.rast', input = 'Ecoregions2017_shp',
                  output = 'Ecoregions2017_biome_type_rast_exp', use = 'attr', attribute_column = 'BIOME_NUM', overwrite = True)

# Rasterize using the Ecoregion code (column ECO_ID) - 30m
grass.run_command('v.to.rast', input = 'Ecoregions2017_shp',
                  output = 'Ecoregions2017_ecoregion_code_rast_exp', use = 'attr', attribute_column = 'ECO_ID', overwrite = True)


# 3. Roads 1980-2010 - vector

# Import
folder_path = r'H:\_neojaguardatabase\Envdatabase\Vetores\World\Roads\groads-v1-americas-shp'
os.chdir(folder_path) # Change to this folder
grass.run_command('v.import', input = 'gROADS-v1-americas.shp', output = 'gROADS_v1_americas_shp', overwrite = True)

# Region of study
grass.run_command('g.region', rast = map_for_define_region, res = 30, flags = 'ap')

# Rasterize (1 = roads, null = non-roads) - 30m
grass.run_command('v.to.rast', input = 'gROADS_v1_americas_shp',
                  output = 'gROADS_v1_americas_rast', use = 'val', value = 1, overwrite = True)
###Falta acrescentar o calculo da distancia euclidiana - ver se vamos ter pra todo neotropico ou só para cada área individual

#---------------------------------------
# Import raster maps that were already cut for the Neotropical region outside GRASS

# 4. Percent tree cover Hansen 2000 - 30m
folder_path = r'H:\_neojaguardatabase\Envdatabase\30m\Neotropic\Hansen_perctreecover_2000'
os.chdir(folder_path) # Change to this folder
grass.run_command('r.import', input = 'Neotropic_Hansen_percenttreecoverd_2000_albers.tif', output = 'Neotropic_Hansen_percenttreecover_2000_30m_tif_exp', overwrite = True)

# 4.5 Hansen Forest patches 2000 (tree cover >= 95%) - 30m

# Region of study
grass.run_command('g.region', rast = map_for_define_region, res = 30, flags = 'ap')

# Transform Tree cover map into a binary forest non-forest map (forest = 1, non forest = 0)
output = 'Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp'
grass.mapcalc(output + ' = if(Neotropic_Hansen_percenttreecover_2000_30m_tif_exp >= 95, 1, 0)', overwrite = True)

# 5. Treecover loss Hansen 2001-2017 - 30m
folder_path = r'H:\_neojaguardatabase\Envdatabase\30m\Neotropic\Hansen_treecoverlossperyear_2017'
os.chdir(folder_path) # Change to this folder
grass.run_command('r.import', input = 'Neotropical_Hansen_treecoverlossperyear_albers_2017.tif', output = 'Neotropical_Hansen_treecoverlossperyear_2017_30m_tif_exp', overwrite = True)

# 5.5 Hansen patches of forest loss 2001-2017 - 30m

# Region of study
grass.run_command('g.region', rast = map_for_define_region, res = 30, flags = 'ap')

# Transform Tree cover map into a binary forest non-forest map (forest = 1, non forest = 0)
output = 'Neotropical_Hansen_treecoverlossperyear_binary_2017_30m_tif_exp'
grass.mapcalc(output + ' = if(Neotropical_Hansen_treecoverlossperyear_2017_30m_tif_exp >= 1, 1, 0)', overwrite = True)

# 6. Elevation - DEM 90m
folder_path = r'H:\_neojaguardatabase\Envdatabase\90m\Neotropic\Earthenv_elevation\Earthenv_dem90_Neotr'
os.chdir(folder_path) # Change to this folder
grass.run_command('r.import', input = 'Neotropic_Earthenv_dem90m_wgs84.tif', output = 'Neotropic_Earthenv_dem90m_tif_exp', overwrite = True)

#---------------------------------------
# Now we are going to reproject and import maps with global coverage.
# First, we created another GRASS GIS location, with Datum WGS84 and Projection lat-lon, and imported/reprojected these
# maps to this location. Finally we cut these maps to the Neotropical region, using the same Hansen tree cover map to define
# the working region.

# Below we reproject this maps from this other location.

# 7. Slope - 1km

# Import (re-project from the WGS84 location)
# I am almost sure this line below is irrelavant for the process, it does not affect the reprojection.
# Still, just to make it sure, we'll keep it
grass.run_command('g.region', rast = map_for_define_region, res = 1000, flags = 'ap')

grass.run_command('r.proj', location = 'newLocation_wgs84', mapset = 'PERMANENT',
                  input = 'Slope_mn_SRTM_1km_neotropic_tif_exp', output = 'Slope_mn_SRTM_1km_neotropic_albers_tif_exp',
                  overwrite = True)
grass.run_command('r.proj', location = 'newLocation_wgs84', mapset = 'PERMANENT',
                  input = 'Slope_md_SRTM_1km_neotropic_tif_exp', output = 'Slope_md_SRTM_1km_neotropic_albers_tif_exp',
                  overwrite = True)

# 8. Human footprint - 1km
grass.run_command('g.region', rast = map_for_define_region, res = 1000, flags = 'ap')

grass.run_command('r.proj', location = 'newLocation_wgs84', mapset = 'PERMANENT',
                  input = 'HFP2009_wgs84_1km_neotropic_tif_exp', output = 'HFP2009_1km_neotropic_albers_tif_exp',
                  overwrite = True)
grass.run_command('r.proj', location = 'newLocation_wgs84', mapset = 'PERMANENT',
                  input = 'HFP1993_wgs84_1km_neotropic_tif_exp', output = 'HFP1993_1km_neotropic_albers_tif_exp',
                  overwrite = True)

# 9. Elevation - 1km
grass.run_command('g.region', rast = map_for_define_region, res = 1000, flags = 'ap')

grass.run_command('r.proj', location = 'newLocation_wgs84', mapset = 'PERMANENT',
                  input = 'Elevation_mn_SRTM_1km_neotropic_tif_exp', output = 'Elevation_mn_SRTM_1km_neotropic_albers_tif_exp',
                  overwrite = True)
grass.run_command('r.proj', location = 'newLocation_wgs84', mapset = 'PERMANENT',
                  input = 'Elevation_md_SRTM_1km_neotropic_tif_exp', output = 'Elevation_md_SRTM_1km_neotropic_albers_tif_exp',
                  overwrite = True)

# 10. Livestock - 1km
grass.run_command('g.region', rast = map_for_define_region, res = 1000, flags = 'ap')

grass.run_command('r.proj', location = 'newLocation_wgs84', mapset = 'PERMANENT',
                  input = 'Livestock_Cattle_CC2006_AD_1km_neotropic_tif_exp', output = 'Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp',
                  overwrite = True)

# 11. Landcover - 300m
grass.run_command('g.region', rast = map_for_define_region, res = 300, flags = 'ap')

grass.run_command('r.proj', location = 'newLocation_wgs84', mapset = 'PERMANENT',
                  input = 'Landcover_ESACCI_2015_300m_neotropic_tif_exp', output = 'Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp',
                  overwrite = True)

# 12. Climate - 1km
grass.run_command('g.region', rast = map_for_define_region, res = 1000, flags = 'ap')

grass.run_command('r.proj', location = 'newLocation_wgs84', mapset = 'PERMANENT',
                  input = 'wc2.0_bio_30s_01_neotropic_tif_exp', output = 'wc2.0_bio_30s_01_neotropic_albers_tif_exp', overwrite = True)
grass.run_command('r.proj', location = 'newLocation_wgs84', mapset = 'PERMANENT',
                  input = 'wc2.0_bio_30s_02_neotropic_tif_exp', output = 'wc2.0_bio_30s_02_neotropic_albers_tif_exp', overwrite = True)
grass.run_command('r.proj', location = 'newLocation_wgs84', mapset = 'PERMANENT',
                  input = 'wc2.0_bio_30s_03_neotropic_tif_exp', output = 'wc2.0_bio_30s_03_neotropic_albers_tif_exp', overwrite = True)
grass.run_command('r.proj', location = 'newLocation_wgs84', mapset = 'PERMANENT',
                  input = 'wc2.0_bio_30s_04_neotropic_tif_exp', output = 'wc2.0_bio_30s_04_neotropic_albers_tif_exp', overwrite = True)
grass.run_command('r.proj', location = 'newLocation_wgs84', mapset = 'PERMANENT',
                  input = 'wc2.0_bio_30s_05_neotropic_tif_exp', output = 'wc2.0_bio_30s_05_neotropic_albers_tif_exp', overwrite = True)
grass.run_command('r.proj', location = 'newLocation_wgs84', mapset = 'PERMANENT',
                  input = 'wc2.0_bio_30s_06_neotropic_tif_exp', output = 'wc2.0_bio_30s_06_neotropic_albers_tif_exp', overwrite = True)
grass.run_command('r.proj', location = 'newLocation_wgs84', mapset = 'PERMANENT',
                  input = 'wc2.0_bio_30s_07_neotropic_tif_exp', output = 'wc2.0_bio_30s_07_neotropic_albers_tif_exp', overwrite = True)
grass.run_command('r.proj', location = 'newLocation_wgs84', mapset = 'PERMANENT',
                  input = 'wc2.0_bio_30s_08_neotropic_tif_exp', output = 'wc2.0_bio_30s_08_neotropic_albers_tif_exp', overwrite = True)
grass.run_command('r.proj', location = 'newLocation_wgs84', mapset = 'PERMANENT',
                  input = 'wc2.0_bio_30s_09_neotropic_tif_exp', output = 'wc2.0_bio_30s_09_neotropic_albers_tif_exp', overwrite = True)
grass.run_command('r.proj', location = 'newLocation_wgs84', mapset = 'PERMANENT',
                  input = 'wc2.0_bio_30s_10_neotropic_tif_exp', output = 'wc2.0_bio_30s_10_neotropic_albers_tif_exp', overwrite = True)
grass.run_command('r.proj', location = 'newLocation_wgs84', mapset = 'PERMANENT',
                  input = 'wc2.0_bio_30s_11_neotropic_tif_exp', output = 'wc2.0_bio_30s_11_neotropic_albers_tif_exp', overwrite = True)
grass.run_command('r.proj', location = 'newLocation_wgs84', mapset = 'PERMANENT',
                  input = 'wc2.0_bio_30s_12_neotropic_tif_exp', output = 'wc2.0_bio_30s_12_neotropic_albers_tif_exp', overwrite = True)
grass.run_command('r.proj', location = 'newLocation_wgs84', mapset = 'PERMANENT',
                  input = 'wc2.0_bio_30s_13_neotropic_tif_exp', output = 'wc2.0_bio_30s_13_neotropic_albers_tif_exp', overwrite = True)
grass.run_command('r.proj', location = 'newLocation_wgs84', mapset = 'PERMANENT',
                  input = 'wc2.0_bio_30s_14_neotropic_tif_exp', output = 'wc2.0_bio_30s_14_neotropic_albers_tif_exp', overwrite = True)
grass.run_command('r.proj', location = 'newLocation_wgs84', mapset = 'PERMANENT',
                  input = 'wc2.0_bio_30s_15_neotropic_tif_exp', output = 'wc2.0_bio_30s_15_neotropic_albers_tif_exp', overwrite = True)
grass.run_command('r.proj', location = 'newLocation_wgs84', mapset = 'PERMANENT',
                  input = 'wc2.0_bio_30s_16_neotropic_tif_exp', output = 'wc2.0_bio_30s_16_neotropic_albers_tif_exp', overwrite = True)
grass.run_command('r.proj', location = 'newLocation_wgs84', mapset = 'PERMANENT',
                  input = 'wc2.0_bio_30s_17_neotropic_tif_exp', output = 'wc2.0_bio_30s_17_neotropic_albers_tif_exp', overwrite = True)
grass.run_command('r.proj', location = 'newLocation_wgs84', mapset = 'PERMANENT',
                  input = 'wc2.0_bio_30s_18_neotropic_tif_exp', output = 'wc2.0_bio_30s_18_neotropic_albers_tif_exp', overwrite = True)
grass.run_command('r.proj', location = 'newLocation_wgs84', mapset = 'PERMANENT',
                  input = 'wc2.0_bio_30s_19_neotropic_tif_exp', output = 'wc2.0_bio_30s_19_neotropic_albers_tif_exp', overwrite = True)


# 13. Population density 2015 - 1km
grass.run_command('g.region', rast = map_for_define_region, res = 1000, flags = 'ap')

grass.run_command('r.proj', location = 'newLocation_wgs84', mapset = 'PERMANENT',
                  input = 'Population_density_gpw_v4_rev10_2015_1km_neotropic_tif_exp', output = 'Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp',
                  overwrite = True)

# 14. Heterogeneity - 1km
grass.run_command('g.region', rast = map_for_define_region, res = 1000, flags = 'ap')

grass.run_command('r.proj', location = 'newLocation_wgs84', mapset = 'PERMANENT',
                  input = 'Homogeneity_01_05_uint16_1km_neotropic_tif_exp', output = 'Homogeneity_01_05_uint16_1km_neotropic_albers_tif_exp',
                  overwrite = True)
grass.run_command('r.proj', location = 'newLocation_wgs84', mapset = 'PERMANENT',
                  input = 'Std_01_05_1km_uint16_neotropic_tif_exp', output = 'Std_01_05_1km_uint16_neotropic_albers_tif_exp',
                  overwrite = True)

# 15. Water frequency 2000 - 1km

grass.run_command('g.region', rast = map_for_define_region, res = 1000, flags = 'ap')

grass.run_command('r.proj', location = 'newLocation_wgs84', mapset = 'PERMANENT',
                  input = 'GIW_sin_1km_bilinear_wgs84_1km_neotropic_tif_exp', output = 'GIW_sin_1km_bilinear_1km_neotropic_albers_tif_exp',
                  overwrite = True)


#---------------------------------------
# Import vectors and rasterize them 
# Import (re-project from the WGS84 location)
# I am almost sure this line below is irrelavant for the process, it does not affect the reprojection.
# Still, just to make it sure, we'll keep it

##Vetores que tivemos que importar no Grass wgs84 mesmo
# 16. Protected Areas 2018 - vector
# We tried to do it in the WGS84 location and then just import the final raster versions
# But the shapefile WDPA_Aug2018-shapefile-polygons.shp have topology issues so we could not import it...
# Protected areas polygons overlap

# Ju: Ber corrigi a topologia e estou tentando importar pra fazer a reprojeção on the fly mesmo, vamos ver se funciona com esse shape novo, criei uma coluna numerica também
# Nome do shape alterado aqui

# Import vector
folder_path = r'H:\_neojaguardatabase\Envdatabase\Vetores\World\Protected areas'
os.chdir(folder_path) # Change to this folder
grass.run_command('v.import', input = 'Neot_WDPA_Aug2018_wgs84_overlapeliminated_numadd.shp', output = 'protected_areas_2018_shp', overwrite = True)

# Region of study
grass.run_command('g.region', rast = map_for_define_region, flags = 'ap')

# Cut - não rolou essa linha
grass.mapcalc('protected_areas_2018_neot_shp = protected_areas_2018_shp', overwrite = True)

# Rasterize separating protected (value=1) and not protected area(value=null) - 30m 
grass.run_command('v.to.rast', input = 'protected_areas_2018_shp',
                  output = 'protected_areas_2018_bin_tif_exp', use = 'attr', attribute_column = 'PA_num', overwrite = True)

# Transform in binary (protected area = 1, non protected = 0) - 30m 
grass.run_command('r.null', map = 'protected_areas_2018_bin_tif_exp', null = 0)
 
# Rasterize using the IUCN category (column IUCN_CAT) - 30m

# Ju: Ber I created a new numeric columm with values based on the IUCN_CAT (Ia = 10, Ib = 11, II = 2, III = 3, IV = 4, V = 5, VI = 6, Not Applicable = 7, Not Assigned = 8, Not Reported = 9)
grass.run_command('v.to.rast', input = 'protected_areas_2018_shp',
                  output = 'protected_areas_2018_iucn_category_tif_exp', use = 'attr', attribute_column = 'IUCN_numer', overwrite = True)


##17. Tree plantations 2013-2014 - vector
# Import vector
folder_path = r'H:\_neojaguardatabase\Envdatabase\Vetores\World\Tree plantations'
os.chdir(folder_path) # Change to this folder
grass.run_command('v.import', input = 'Tree_plantations.shp', output = 'tree_plantations_2013_shp', overwrite = True)

# Region of study
grass.run_command('g.region', rast = map_for_define_region, res = 30, flags = 'ap')

# Rasterize using a binary classification (tree plantation = 1, non tree plantation = null) - 30m
grass.run_command('v.to.rast', input = 'tree_plantations_2013_shp',
                  output = 'tree_plantations_binary_tif_exp', use = 'attr', attribute_column = 'TREEPL_num', overwrite = True)

# Transform in binary (tree plantation = 1, non tree plantation = 0) - 30m
grass.run_command('r.null', map = 'tree_plantations_binary_tif_exp', null = 0)


# Rasterize using the type of tree plantation (column spec_1) - 30m 

# We created a numerice columm for each category of the columm spec_simp: Fruit = 1, Fruit mix = 2, Oil palm = 3, Oil palm mix = 4,
# Other = 5, Other mix = 6, Recently cleared = 7, Rubber = 8,Rubber mix = 9, Unknown = 10, Wood fiber / timber = 11,
# Wood fiber / timber mix = 12
grass.run_command('v.to.rast', input = 'tree_plantations_2013_shp',
                  output = 'tree_plantation_type_tif_exp', use = 'attr', attribute_column = 'SPEC_num', overwrite = True)


# 18. Drainage

# Import vector
folder_path = r'H:\_neojaguardatabase\Envdatabase\Vetores\World\Hydrosheds\River networks'
os.chdir(folder_path) # Change to this folder
grass.run_command('v.import', input = 'ca_sa_riv_merged_15s_albers.shp', output = 'drainage_15s_shp', overwrite = True)

# Region of study
grass.run_command('g.region', rast = map_for_define_region, res = 30, flags = 'ap')


#Oi Ber terminou de rodar v.to.rast de drenagem mas também não rolou pq esse valor de 1 no script não fazia sentido mesmo
#pro shape de linha da drenagem - vou tentar gerar no ArcGIS também e depois importamos

#Ju: Ber criei a coluna no shape DREN_num com valor 1, e também o raster binário no ArcGIS. Vou tentar importar o shape se não der importamos o raster direto
# corrigi o nome aqui no script

# Rasterize using a binary classification (drainage = 1, non protected = null) - 30m
grass.run_command('v.to.rast', input = 'drainage_15s_shp',
                  output = 'drainage_15s_binary_tif_exp', use = 'attr', attribute_column = 'DREN_num', overwrite = True)

# Transform in binary (drainage = 1, non drainage = 0) - 30m 
grass.run_command('r.null', map = 'drainage_15s_binary_tif_exp', null = 0)

####Rodei at'e aqui :)


###Falta acrescentar o calculo da distancia euclidiana - ver se vamos ter pra todo neotropico ou só para cada área individual


# 13. Major Dams
# Are we keeping this one? In my opinion, it is not worth it!! We have points, not polygons, so it does not make sense
# for fine scale analyses such as step selection functions


# 44. Patch size

grass.run_command('r.null', map = map_for_define_region, null = '0')

lsmetrics_dir = r'H:\_neojaguardatabase\LS_METRICS\_LSMetrics_v1_0_0'
os.chdir(lsmetrics_dir)

# Run LSMetrics
subprocess.call('python LSMetrics_v1_0_0.py', shell=True) # runs and wait

#---------------------------------------
# Import buffers

folder_path = r'H:\_neojaguardatabase\Buffer 70'
os.chdir(folder_path) # Change to this folder
files = os.listdir(folder_path) # List files in the folder
for i in files:
    if i[-3:] == 'shp': # Select tif files
        print i
        name = 'b' + i.replace('.shp', '_shp')
        grass.run_command('v.import', input = i, output = name, overwrite = True) # Import maps




####Ju: Ber, rodei até aqui, aparentemente importou os buffers certinho

#---------------------------------------
# Cortes

# Region of study
grass.run_command('g.region', rast = map_for_define_region, res = 30, flags = 'ap')

# List of buffers
list_buffers = grass.list_grouped('vect', pattern = '*Buffer*')['PERMANENT']

# Root of output folder
output_folder = r'H:\_neojaguardatabase\Envdatabase\Basesfinais_res30m'

# Mask
for i in list_buffers:
    print i
    grass.run_command('g.region', vect = i, res = 30,
                      align = map_for_define_region, flags = 'p')
    grass.run_command('r.mask', vector = i) # Mask for the buffer for individual i
    #grass.run_command('g.region', vect = i, res = 30,
    #                  align = map_for_define_region, flags = 'p')
    
    # Aqui devemos calcular as distancias ou outros procedimentos - ROADS e RIVERS
    # !!!! overwrite = True
    # roads: gROADS_v1_americas_30m_tif@PERMANENT
    name_roads = 'dist2roads_exp'
    grass.run_command('r.grow.distance', input = 'gROADS_v1_americas_30m_tif', distance = name_roads, overwrite = True)
    
    # drainage: drainage_15s_2018_tif
    name_drainage = 'dist2drainage_exp'
    grass.run_command('r.grow.distance', input = 'drainage_15s_2018_tif', distance = name_drainage, overwrite = True)
    
    # distance to forest edges: Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp
    lsmetrics_dir = r'H:\_neojaguardatabase\LS_METRICS\_LSMetrics_v1_0_0'
    os.chdir(lsmetrics_dir)
    
    # Run LSMetrics
    from LSMetrics_v1_0_0 import dist_edge
    
    name_forest = 'dist2forestedges_'
    
    dist_edge(input_maps = ['Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp'],
              classify_edge_as_zero = False,
              prepare_biodim = False, remove_trash = True,
              prefix = name_forest, add_counter_name = False, export = False, dirout = '')
    
    # distance to water: water_frequency_2010_30m_tif_exp
    water_1_null = 'waterfreq_1_null'
    grass.mapcalc(water_1_null + ' = if(water_frequency_2010_30m_tif_exp == 1, null(), 2)', overwrite = True)
    name_water = 'dist2waterbodies_exp'
    grass.run_command('r.grow.distance', input = water_1_null, distance = name_water, overwrite = True)
    
    # patch size
    lsmetrics_dir = r'H:\_neojaguardatabase\LS_METRICS\_LSMetrics_v1_0_0'
    os.chdir(lsmetrics_dir)
    
    # Run LSMetrics
    from LSMetrics_v1_0_0 import rulesreclass, patch_size
        
    name_patchsize = 'local_'
    
    patch_size(input_maps = ['Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp'],
              prefix = name_patchsize)
        
    # List of raster to be exported
    list_rast_export = grass.list_grouped(type = 'rast', pattern = '*exp*')['PERMANENT']
    
    # Outpt folder
    os.chdir(output_folder)
    
    buffer_dir = i.replace('_shp', '')
    
    # Create output folder
    if not os.path.exists(buffer_dir):
        os.mkdir(buffer_dir)
        print "Directory " + buffer_dir + " Created "
    else:    
        print "Directory " + buffer_dir + " already exists"
    
    os.chdir(buffer_dir)
    
    grass.run_command('g.region', vect = i, res = 30,
                      align = map_for_define_region, flags = 'p')
    
    # Loop to export all raster
    for j in list_rast_export:
        print j
        
        prefix = i[0:4] + '_'
        
        # Define no data value
        try:
            
            grass.run_command('r.out.gdal', input = j, output = prefix+j+'.tif', createopt = "TFW=YES,COMPRESS=DEFLATE",
                              overwrite = True)
            
        except:
            max_val = float(grass.parse_command('r.info', map = j, flags = 'r')['max'])
            maxvals = [255, 32767, 2e6, 4e6, 3e38]
            
            nodata_val = [x for x in maxvals if max_val <= x][0]
            
            if j == 'Neotropic_Earthenv_dem90m_tif_exp':
                grass.run_command('r.out.gdal', input = j, output = prefix+j+'.tif', createopt = "TFW=YES,COMPRESS=DEFLATE",
                                  nodata = 32767, type = 'Int16', overwrite = True)
            else:
                grass.run_command('r.out.gdal', input = j, output = prefix+j+'.tif', createopt = "TFW=YES,COMPRESS=DEFLATE",
                              nodata = nodata_val, overwrite = True)
            ## This nodata value above may create problems. Let's keep an eye on that just to make sure.
    
    # Remove mask
    grass.run_command('r.mask', flags = 'r')
    grass.run_command('g.region', rast = map_for_define_region, res = 30, flags = 'ap')
    
    # Remove temp distance maps
    grass.run_command('g.remove', type = 'raster', pattern = 'dist2*', flags = 'f')
    grass.run_command('g.remove', type = 'raster', pattern = name_patchsize+'*Hansen*', flags = 'f')
    grass.run_command('g.remove', type = 'raster', name = 'waterfreq_1_null', flags = 'f')
    

