#---------------------------------------------------------------------------------------
"""
 Neo Jaguar Environmental Database - WGS84

 Script to import spatial layers with global coverage in WGS84 Datum and clip them from
 the Neotropics to reproject them into Albers-SAD69.
 This is part of the Neotropical Jaguar movement research agenda
 
 Bernardo B. S. Niebuhr - bernardo_brandaum@yahoo.com.br
 Julia E. F. Oshima - juliaoshima@yahoo.com.br
 Vanesa Bejarano - vanesa.bejarano@gmail.com
 
 Laboratorio de Ecologia Espacial e Conservacao
 Universidade Estadual Paulista - UNESP
 Rio Claro - SP - Brasil
 
 License GPLv2 - feel free to use, modify, and share, but cite us and make your code free.
"""
#---------------------------------------------------------------------------------------

# Load modules
import os
import grass.script as grass
import subprocess

# 0. Hansen tree cover
folder_path = r'H:\_neojaguardatabase\Envdatabase\30m\Neotropic\Hansen_perctreecover_2000'
os.chdir(folder_path) # Change to this folder
grass.run_command('r.import', input = 'Neotropic_Hansen_forest1_0_95percenttreecover_2000_wgs84.tif',
                  output = 'Neotropic_Hansen_forest1_0_95percenttreecover_2000_wgs84_30m_tif', overwrite = True) # Import maps

# Define region
map_for_define_region = 'Neotropic_Hansen_forest1_0_95percenttreecover_2000_wgs84_30m_tif'
# Region of study

# Ju: Ber essa linha abaixo acho que ta errada pq ta em metros né? Tava res = 30, seria res = '0:00:01'  ?

# grass.run_command('g.region', rast = map_for_define_region,  res = 30, flags = 'ap')
grass.run_command('g.region', rast = map_for_define_region, flags = 'ap') #Não coloquei nada ele pegou a resolução original


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



