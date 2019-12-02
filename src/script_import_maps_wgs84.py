# Neo Jaguar Database - WGS84

# os.chdir(r'H:\_neojaguardatabase\Envdatabase\30m\Neotropic\Water frequency\2010')
# grass.run_command('r.import', input = 'p001r050_WF_2010.tif', output = 'p001r050_WF_2010')
python

# Load modules
import os
import subprocess
import grass.script as grass
from grass.pygrass.modules.shortcuts import general as g
from grass.pygrass.modules.shortcuts import vector as v
from grass.pygrass.modules.shortcuts import raster as r

# 0. Hansen tree cover
folder_path = r'H:\_neojaguardatabase\Envdatabase\30m\Neotropic\Hansen_perctreecover_2000'
os.chdir(folder_path) # Change to this folder
grass.run_command('r.import', input = 'Neotropic_Hansen_forest1_0_95percenttreecover_2000_wgs84.tif',
                  output = 'Neotropic_Hansen_forest1_0_95percenttreecover_2000_wgs84_30m_tif', overwrite = True) # Import maps

# Define region
map_for_define_region = 'Neotropic_Hansen_percenttreecoverd_2000_wgs84'
# Region of study

# Ju: Ber essa linha abaixo acho que ta errada pq ta em metros né? Tava res = 30, seria res = '0:00:01'  ?

# grass.run_command('g.region', rast = map_for_define_region,  res = 30, flags = 'ap')
grass.run_command('g.region', rast = map_for_define_region, flags = 'ap') #Não coloquei nada ele pegou a resolução original

# 1. Outros Hansens
# 1.2 Treecover loss Hansen 2001-2017 - 30m
folder_path = r'F:\_neojaguardatabase\Envdatabase\30m\Neotropic\Hansen_treecoverlossperyear_2017'
os.chdir(folder_path) # Change to this folder

ff = os.listdir('.')

for i in ff:
    if i.startswith('Hansen'):
        print i
        name = i.replace('.tif', '')
        try:
            grass.run_command('r.in.gdal', input = i, output = i, overwrite = True)

grass.run_command('r.import', input = 'Neotropical_Hansen_treecoverlossperyear_wgs84_2017.tif', 
	output = 'Neotropical_Hansen_treecoverlossperyear_wgs84_2017', overwrite = True)


# 2. Water Frequency 2000 - 30m

# Import map
folder_path = r'F:\_neojaguardatabase\Envdatabase\30m\Neotropic\Water frequency\2010'
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
grass.run_command('g.region', rast = map_for_define_region, flags = 'ap')

# Combine maps
water_map_mosaic = 'water_frequency_2010_30m_tif_exp'
grass.run_command('r.patch', input = maps_water, output = water_map_mosaic, overwrite = True)

# Delete input maps
grass.run_command('g.remove', type = 'raster', pattern = 'p*2010_rast', flags = 'f')

# transformar em 1/0 e 1/null

# Cortar todos os rasters mundiais só para região neotropical

# 2. Slope
# Import
folder_path = r'H:\_neojaguardatabase\Envdatabase\1km\World\Earthenv_topography\Slope'
os.chdir(folder_path) # Change to this folder
grass.run_command('r.import', input = 'slope_1KMmn_SRTM.tif', output = 'Slope_mn_SRTM_1km_tif_exp', overwrite = True)
grass.run_command('r.import', input = 'slope_1KMmd_SRTM.tif', output = 'Slope_md_SRTM_1km_tif_exp', overwrite = True)

# Region of study
grass.run_command('g.region', rast = map_for_define_region, res = '0:00:30', flags = 'ap') # res = 1km

# Cut
grass.mapcalc('Slope_mn_SRTM_1km_neotropic_tif_exp = Slope_mn_SRTM_1km_tif_exp', overwrite = True)
grass.mapcalc('Slope_md_SRTM_1km_neotropic_tif_exp = Slope_md_SRTM_1km_tif_exp', overwrite = True)

#3. Human footprint
# Import
folder_path = r'H:\_neojaguardatabase\Envdatabase\1km\World\WCS_Humanfootprint\HumanFootprintv2\Dryadv3\Maps'
os.chdir(folder_path) # Change to this folder
grass.run_command('r.import', input = 'HFP2009_wgs84.tif', output = 'HFP2009_wgs84_1km_tif_exp', overwrite = True)
grass.run_command('r.import', input = 'HFP1993_wgs84.tif', output = 'HFP1993_wgs84_1km_tif_exp', overwrite = True)

# Region of study
grass.run_command('g.region', rast = map_for_define_region, res = '0:00:30', flags = 'ap')

# Cut
grass.mapcalc('HFP2009_wgs84_1km_neotropic_tif_exp = HFP2009_wgs84_1km_tif_exp', overwrite = True)
grass.mapcalc('HFP1993_wgs84_1km_neotropic_tif_exp = HFP1993_wgs84_1km_tif_exp', overwrite = True)

# 5. Elevation
# Import
folder_path = r'H:\_neojaguardatabase\Envdatabase\1km\World\Earthenv_topography\Elevation'
os.chdir(folder_path) # Change to this folder
grass.run_command('r.import', input = 'elevation_1KMmn_SRTM.tif', output = 'Elevation_mn_SRTM_1km_tif_exp', overwrite = True)
grass.run_command('r.import', input = 'elevation_1KMmd_SRTM.tif', output = 'Elevation_md_SRTM_1km_tif_exp', overwrite = True)

# Region of study
grass.run_command('g.region', rast = map_for_define_region, res = '0:00:30', flags = 'ap')

# Cut
grass.mapcalc('Elevation_mn_SRTM_1km_neotropic_tif_exp = Elevation_mn_SRTM_1km_tif_exp', overwrite = True)
grass.mapcalc('Elevation_md_SRTM_1km_neotropic_tif_exp = Elevation_md_SRTM_1km_tif_exp', overwrite = True)

# 6. Livestock
# Import
folder_path = r'H:\_neojaguardatabase\Envdatabase\1km\World\Livestock\livestock_CATTLE\CATTLE'
os.chdir(folder_path) # Change to this folder
grass.run_command('r.import', input = 'Glb_Cattle_CC2006_AD.tif', output = 'Livestock_Cattle_CC2006_AD_1km_tif_exp', overwrite = True)

# Region of study
grass.run_command('g.region', rast = map_for_define_region, res = '0:00:30', flags = 'ap')

# Cut
grass.mapcalc('Livestock_Cattle_CC2006_AD_1km_neotropic_tif_exp = Livestock_Cattle_CC2006_AD_1km_tif_exp', overwrite = True)


# 7. Landcover 300m  -
# Import
folder_path = r'H:\_neojaguardatabase\Envdatabase\300m\Landcover ESA Climate Change Initiative - led by UCLouvain\ESACCI-LC-L4-LCCS-Map-300m-P1Y-2015-v2.0.7\product'
os.chdir(folder_path) # Change to this folder
grass.run_command('r.import', input = 'ESACCI-LC-L4-LCCS-Map-300m-P1Y-2015-v2.0.7.tif', output = 'Landcover_ESACCI_2015_300m_tif_exp', overwrite = True)

# Region of study
grass.run_command('g.region', rast = map_for_define_region, res = '0:00:10', flags = 'ap')

# Cut - NÃO ROLOU....BER pode ver se entende o erro?
grass.mapcalc('Landcover_ESACCI-2015_1km_neotropic_tif_exp = Landcover_ESACCI_2015_300m_tif_exp', overwrite = True)


# 8. Climate - ESSA PARTE JA FOI FEITA
# Import 
folder_path = r'H:\_neojaguardatabase\Envdatabase\1km\World\Worldclim_v2\Bioclim_30s\wc2.0_30s_bio'
os.chdir(folder_path) # Change to this folder
grass.run_command('r.import', input = 'wc2.0_bio_30s_01.tif', output = 'wc2.0_bio_30s_01_tif_exp', overwrite = True)
grass.run_command('r.import', input = 'wc2.0_bio_30s_02.tif', output = 'wc2.0_bio_30s_02_tif_exp', overwrite = True)
grass.run_command('r.import', input = 'wc2.0_bio_30s_03.tif', output = 'wc2.0_bio_30s_03_tif_exp', overwrite = True)
grass.run_command('r.import', input = 'wc2.0_bio_30s_04.tif', output = 'wc2.0_bio_30s_04_tif_exp', overwrite = True)
grass.run_command('r.import', input = 'wc2.0_bio_30s_05.tif', output = 'wc2.0_bio_30s_05_tif_exp', overwrite = True)
grass.run_command('r.import', input = 'wc2.0_bio_30s_06.tif', output = 'wc2.0_bio_30s_06_tif_exp', overwrite = True)
grass.run_command('r.import', input = 'wc2.0_bio_30s_07.tif', output = 'wc2.0_bio_30s_07_tif_exp', overwrite = True)
grass.run_command('r.import', input = 'wc2.0_bio_30s_08.tif', output = 'wc2.0_bio_30s_08_tif_exp', overwrite = True)
grass.run_command('r.import', input = 'wc2.0_bio_30s_09.tif', output = 'wc2.0_bio_30s_09_tif_exp', overwrite = True)
grass.run_command('r.import', input = 'wc2.0_bio_30s_10.tif', output = 'wc2.0_bio_30s_10_tif_exp', overwrite = True)
grass.run_command('r.import', input = 'wc2.0_bio_30s_11.tif', output = 'wc2.0_bio_30s_11_tif_exp', overwrite = True)
grass.run_command('r.import', input = 'wc2.0_bio_30s_12.tif', output = 'wc2.0_bio_30s_12_tif_exp', overwrite = True)
grass.run_command('r.import', input = 'wc2.0_bio_30s_13.tif', output = 'wc2.0_bio_30s_13_tif_exp', overwrite = True)
grass.run_command('r.import', input = 'wc2.0_bio_30s_14.tif', output = 'wc2.0_bio_30s_14_tif_exp', overwrite = True)
grass.run_command('r.import', input = 'wc2.0_bio_30s_15.tif', output = 'wc2.0_bio_30s_15_tif_exp', overwrite = True)
grass.run_command('r.import', input = 'wc2.0_bio_30s_16.tif', output = 'wc2.0_bio_30s_16_tif_exp', overwrite = True)
grass.run_command('r.import', input = 'wc2.0_bio_30s_17.tif', output = 'wc2.0_bio_30s_17_tif_exp', overwrite = True)
grass.run_command('r.import', input = 'wc2.0_bio_30s_18.tif', output = 'wc2.0_bio_30s_18_tif_exp', overwrite = True)
grass.run_command('r.import', input = 'wc2.0_bio_30s_19.tif', output = 'wc2.0_bio_30s_19_tif_exp', overwrite = True)

# Region of study 
grass.run_command('g.region', rast = map_for_define_region, res = '0:00:30', flags = 'ap')

# Cut
grass.mapcalc('wc2.0_bio_30s_01_neotropic_tif_exp = wc2.0_bio_30s_01_tif_exp', overwrite = True)
grass.mapcalc('wc2.0_bio_30s_02_neotropic_tif_exp = wc2.0_bio_30s_02_tif_exp', overwrite = True)
grass.mapcalc('wc2.0_bio_30s_03_neotropic_tif_exp = wc2.0_bio_30s_03_tif_exp', overwrite = True)
grass.mapcalc('wc2.0_bio_30s_04_neotropic_tif_exp = wc2.0_bio_30s_04_tif_exp', overwrite = True)
grass.mapcalc('wc2.0_bio_30s_05_neotropic_tif_exp = wc2.0_bio_30s_05_tif_exp', overwrite = True)
grass.mapcalc('wc2.0_bio_30s_06_neotropic_tif_exp = wc2.0_bio_30s_06_tif_exp', overwrite = True)
grass.mapcalc('wc2.0_bio_30s_07_neotropic_tif_exp = wc2.0_bio_30s_07_tif_exp', overwrite = True)
grass.mapcalc('wc2.0_bio_30s_08_neotropic_tif_exp = wc2.0_bio_30s_08_tif_exp', overwrite = True)
grass.mapcalc('wc2.0_bio_30s_09_neotropic_tif_exp = wc2.0_bio_30s_09_tif_exp', overwrite = True)
grass.mapcalc('wc2.0_bio_30s_10_neotropic_tif_exp = wc2.0_bio_30s_10_tif_exp', overwrite = True)
grass.mapcalc('wc2.0_bio_30s_11_neotropic_tif_exp = wc2.0_bio_30s_11_tif_exp', overwrite = True)
grass.mapcalc('wc2.0_bio_30s_12_neotropic_tif_exp = wc2.0_bio_30s_12_tif_exp', overwrite = True)
grass.mapcalc('wc2.0_bio_30s_13_neotropic_tif_exp = wc2.0_bio_30s_13_tif_exp', overwrite = True)
grass.mapcalc('wc2.0_bio_30s_14_neotropic_tif_exp = wc2.0_bio_30s_14_tif_exp', overwrite = True)
grass.mapcalc('wc2.0_bio_30s_15_neotropic_tif_exp = wc2.0_bio_30s_15_tif_exp', overwrite = True)
grass.mapcalc('wc2.0_bio_30s_16_neotropic_tif_exp = wc2.0_bio_30s_16_tif_exp', overwrite = True)
grass.mapcalc('wc2.0_bio_30s_17_neotropic_tif_exp = wc2.0_bio_30s_17_tif_exp', overwrite = True)
grass.mapcalc('wc2.0_bio_30s_18_neotropic_tif_exp = wc2.0_bio_30s_18_tif_exp', overwrite = True)
grass.mapcalc('wc2.0_bio_30s_19_neotropic_tif_exp = wc2.0_bio_30s_19_tif_exp', overwrite = True)


# 9. Population density 1km - Raster esta com varias áreas de NoData na Amazonia,
# ver se precisa transformar em 0 na hora das análises de SSF e RSF
# Import
folder_path = r'H:\_neojaguardatabase\Envdatabase\1km\World\Population Density\gpw v4 2015'
os.chdir(folder_path) # Change to this folder
grass.run_command('r.import', input = 'gpw_v4_population_density_rev10_2015_30_sec.tif', output = 'Population_density_gpw_v4_rev10_2015_1km_tif_exp', overwrite = True)

# Region of study
grass.run_command('g.region', rast = map_for_define_region, res = '0:00:30', flags = 'ap')

# Cut
grass.mapcalc('Population_density_gpw_v4_rev10_2015_1km_neotropic_tif_exp = Population_density_gpw_v4_rev10_2015_1km_tif_exp', overwrite = True)


# 10. Heterogeneity
# Import
folder_path = r'H:\_neojaguardatabase\Envdatabase\1km\World\Earthenv_heterogeneity'
os.chdir(folder_path) # Change to this folder
grass.run_command('r.import', input = 'Homogeneity_01_05_1km_uint16.tif', output = 'Homogeneity_01_05_uint16_1km_tif_exp', overwrite = True)
grass.run_command('r.import', input = 'std_01_05_1km_uint16.tif', output = 'Std_01_05_1km_uint16_tif_exp', overwrite = True)

# Region of study
grass.run_command('g.region', rast = map_for_define_region, res = '0:00:30', flags = 'ap')

# Cut
grass.mapcalc('Homogeneity_01_05_uint16_1km_neotropic_tif_exp = Homogeneity_01_05_uint16_1km_tif_exp', overwrite = True)
grass.mapcalc('Std_01_05_1km_uint16_neotropic_tif_exp = Std_01_05_1km_uint16_tif_exp', overwrite = True)


# 11. Water frequency 1km
# Import
folder_path = r'H:\_neojaguardatabase\Envdatabase\1km\World\Water frequency'
os.chdir(folder_path) # Change to this folder
grass.run_command('r.import', input = 'GIW_sin_1km_bilinear_wgs84.tif', output = 'GIW_sin_1km_bilinear_wgs84_1km_tif_exp', overwrite = True)

# Region of study
grass.run_command('g.region', rast = map_for_define_region, res = '0:00:30', flags = 'ap')

# Cut
grass.mapcalc('GIW_sin_1km_bilinear_wgs84_1km_neotropic_tif_exp = GIW_sin_1km_bilinear_wgs84_1km_tif_exp', overwrite = True)


# 2. Ecoregions 2017 - vector

# Import vector
folder_path = r'E:\_neojaguardatabase\Envdatabase\Vetores\World\Ecoregions'
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
folder_path = r'F:\_neojaguardatabase\Envdatabase\Vetores\World\Roads\groads-v1-americas-shp'
os.chdir(folder_path) # Change to this folder
grass.run_command('v.import', input = 'gROADS-v1-americas.shp', output = 'gROADS_v1_americas_shp', overwrite = True)

# Region of study
grass.run_command('g.region', rast = map_for_define_region, #res = '00:00:30', 
	flags = 'ap')

# Rasterize (1 = roads, null = non-roads) - 30m
grass.run_command('v.to.rast', input = 'gROADS_v1_americas_shp',
                  output = 'gROADS_v1_americas_rast', use = 'val', value = 1, overwrite = True)
###Falta acrescentar o calculo da distancia euclidiana - ver se vamos ter pra todo neotropico ou só para cada área individual


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
folder_path = r'F:\_neojaguardatabase\Envdatabase\Vetores\World\Protected areas'
os.chdir(folder_path) # Change to this folder
grass.run_command('v.import', input = 'Neot_WDPA_Aug2018_wgs84_overlapeliminated_numadd.shp', output = 'protected_areas_2018_shp', overwrite = True)

# Region of study
grass.run_command('g.region', rast = map_for_define_region, flags = 'ap')

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
folder_path = r'E:\_neojaguardatabase\Envdatabase\Vetores\World\Tree plantations'
os.chdir(folder_path) # Change to this folder
grass.run_command('v.import', input = 'Tree_plantations.shp', output = 'tree_plantations_2013_shp', overwrite = True)

# Region of study
grass.run_command('g.region', rast = map_for_define_region, res = '00:00:00.9', flags = 'ap')

# Rasterize using a binary classification (tree plantation = 1, non tree plantation = null) - 30m
grass.run_command('v.to.rast', input = 'tree_plantations_2013_shp',
                  output = 'tree_plantations_binary_exp', use = 'val', value = '1', overwrite = True)

# Transform in binary (tree plantation = 1, non tree plantation = 0) - 30m
grass.run_command('r.null', map = 'tree_plantations_binary_exp', null = 0)


# Rasterize using the type of tree plantation (column spec_1) - 30m 

# We created a numerice columm for each category of the columm spec_simp: Fruit = 1, Fruit mix = 2, Oil palm = 3, Oil palm mix = 4,
# Other = 5, Other mix = 6, Recently cleared = 7, Rubber = 8,Rubber mix = 9, Unknown = 10, Wood fiber / timber = 11,
# Wood fiber / timber mix = 12
grass.run_command('v.to.rast', input = 'tree_plantations_2013_shp',
                  output = 'tree_plantation_type_exp', use = 'attr', attribute_column = 'SPEC_num', overwrite = True)


# 18. Drainage

# Import vector
folder_path = r'E:\_neojaguardatabase\Envdatabase\Vetores\World\Hydrosheds\River networks'
os.chdir(folder_path) # Change to this folder
grass.run_command('v.import', input = 'ca_sa_riv_merged_15s_albers.shp', output = 'drainage_15s_shp', overwrite = True)

# Region of study
grass.run_command('g.region', rast = map_for_define_region, res = '00:00:00.9', flags = 'ap')


# Rasterize using a binary classification (drainage = 1, non protected = null) - 30m
grass.run_command('v.to.rast', input = 'drainage_15s_shp',
                  output = 'drainage_15s_binary_tif_exp', use = 'attr', attribute_column = 'DREN_num', overwrite = True)

# Transform in binary (drainage = 1, non drainage = 0) - 30m 
grass.run_command('r.null', map = 'drainage_15s_binary_tif_exp', null = 0)
