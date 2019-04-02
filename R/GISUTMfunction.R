#' #  **GIS bases in UTM**
#' 
#' #### *Alan E. de Barros, Bernardo Niebuhr, Vanesa Bejarano, Julia Oshima,Claudia Kanda, Milton Ribeiro, Ronaldo Morato,Paulo Prado*
#' date: "March, 21 2019"
#' Run JaguarDataPrep first !!! 
#' 
#' Load packages a few more packages                      
if(!require(install.load)) install.packages('install.load'); library(install.load)
install.load::install_load("raster","RCurl","rts","rvest","rgl","lubridate","lattice","rgdal","sp","stringr","methods",
                           "maptools", "vegan","sp","spatialEco","installr")

lcover <-"D:/GISbases/b001JaguarAlbersBuffer70/b001_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif"
land_cover=raster(lcover)
crs(land_cover)   

newproj <- "+proj=longlat +datum=WGS84"              
land_cover<- projectRaster(land_cover, crs=newproj)


aea_to_longlat =  function(tiff.path, output.path, newproj = "+proj=longlat +datum=WGS84", Jaguar.id){
  land_cover=raster(tiff.path)
  land_cover <- projectRaster(land_cover, crs=newproj1)
  writeRaster(land_cover, filename=output.path, format="GTiff", overwrite=TRUE)
  #comand para exportar, provavelmente precisa de um nome para o arquivo. Ja reservei um argumento output.path para isto
}

## Um loop para transformar em todos os arquivos de uma pasta
## Path do doretorio onde estao todos os arquivos de cada bicho
path1 = "D:/GISbases/"
## nomes dos diretorios de cada bicho
dir2 = dir(path1, full.names=TRUE, pattern ="^b[0-9]")
## Indice inicial do vetor de nomes de projecoes (um nome para cada bicho)
j = 1
for(path2 in dir2){
  ## nomes de todos os arquivos de uma pasta de um bicho
  dir3 = dir(path2, full.names=TRUE, pattern = "tif$")
  for(path3 in dir3){
    out.name = paste(gsub(".tif", "_WGS.tif", path3, fixed=TRUE), sep="")
    
    cat(dir1[i], "\n")
    cat(out.name, "\n \n")
    #aea_to_longlat(tiff.path = path3, output.path = out.name , new.proj2=vetor.nomes.proj[j])
  }
  j=j+1
}












#tiff_to_utm =  function(tiff.path, newproj1 = "+proj=longlat +datum=WGS84", newproj2, Jaguar.id){
lcover <-"D:/GISbases/b001JaguarAlbersBuffer70/b001_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif"
land_cover=raster(lcover)
crs(land_cover)   
plot(land_cover)
newproj <- "+proj=longlat +datum=WGS84"
pr1 <- projectRaster(land_cover, crs=newproj)
newproj2 <- "+proj=utm +zone=21K +south +datum=WGS84 +units=m +no_defs"
land_cover <- projectRaster(land_cover, crs=newproj2)
coord.UTM  = SpatialPoints(cbind(J1$utm_x ,J1$utm_y),CRS("+proj=utm +zone=21K +south +datum=WGS84 +units=m +no_defs"))     
coord.UTM 
compareCRS(coord.UTM,land_cover)
plot(land_cover)
#test <- writeRaster(land_cover, filename="test.tif", format="GTiff", overwrite=TRUE)
#test
plot(coord.UTM,add=T)


t1 <-"D:/Documents/GitHub/jaguar-codes/R/test.tif"
t1=raster(t1)
crs(t1)   
plot(t1)

newproj <- "+proj=longlat +datum=WGS84"
pr1 <- projectRaster(land_cover, crs=newproj)
newproj2 <- "+proj=utm +zone=21K +south +datum=WGS84 +units=m +no_defs"
land_cover <- projectRaster(land_cover, crs=newproj2)

tiff_to_utm =  function(tiff.path, output.path, newproj1 = "+proj=longlat +datum=WGS84", newproj2, Jaguar.id){
  land_cover=raster(tiff.path)
  pr1 <- projectRaster(land_cover, crs=newproj1)
  land_cover <- projectRaster(land_cover, crs=newproj2)
  #coord.UTM  = SpatialPoints(cbind(J1$utm_x ,J1$utm_y),CRS(newproj2))     
  writeRaster(land_cover, filename="test.tif", format="GTiff", overwrite=TRUE)
  #comand para exportar, provavelmente precisa de um nome para o arquivo. Ja reservei um argumento output.path para isto
}

## Um loop para transformar em todos os arquivos de uma pasta
## Path do doretorio onde estao todos os arquivos de cada bicho
path1 = "D:/GISbases/"
## nomes dos diretorios de cada bicho
dir2 = dir(path1, full.names=TRUE, pattern ="^b[0-9]")
## Indice inicial do vetor de nomes de projecoes (um nome para cada bicho)
j = 1
for(path2 in dir2){
## nomes de todos os arquivos de uma pasta de um bicho
  dir3 = dir(path2, full.names=TRUE, pattern = "tif$")
  for(path3 in dir3){
    out.name = paste(gsub(".tif", "_UTM.tif", path3, fixed=TRUE), sep="")
    
    #cat(dir1[i], "\n")
    #cat(out.name, "\n \n")
    #tiff_to_utm(tiff.path = path3, output.path = out.name , new.proj2=vetor.nomes.proj[j])
  }
  j=j+1
}


# "D:/GISbases/J1/Landcover.tif" ### individual folders?
# export LD_LIBRARY_PATH=/lib:/usr/lib:/usr/local/lib





#############################################################################################################################



#' Anthropic
#b001_dist2roads_exp
#' J1
#' b001_dist2roads_exp
dr <-"D:/GISbases/b001JaguarAlbersBuffer70/b001_dist2roads_exp.tif"
dr=raster(dr)
#crs(dr)   
#plot(dr)
dr <- projectRaster(dr, crs="+proj=longlat +datum=WGS84")
#crs(dr)   
dr <- projectRaster(dr, crs=get_crs(J001trk))
#crs(dr)   
#compareCRS(dr,get_crs(J001trk))
writeRaster(dr,filename=file.path("D:/GISUTM/J001/UTM_b0001_dist2roads_exp"), format="GTiff", overwrite=TRUE)

test <-"D:/GISUTM/J93/UTM_b093_dist2roads_exp.tif";(test=raster(test))
x11()
plot(test)

#' Anthropic
#b001_dist2roads_exp        <--------------------------------------------------------------
#b001_human_footprint_1993_1km_tif_exp
#b001_human_footprint_2009_1km_tif_exp    <--------------------------------------------------------------

b001_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp    <--------------------------------------------------------------
b001_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp    <--------------------------------------------------------------
#' Classification
b001_Ecoregions2017_biome_type_rast_exp
b001_Ecoregions2017_ecoregion_code_rast_exp
b001_protected_areas_2018_bin_tif_exp
b001_protected_areas_2018_iucn_category_tif_exp
#' Climatic
b001_wc2.0_bio_30s_01_neotropic_albers_tif_exp
b001_wc2.0_bio_30s_02_neotropic_albers_tif_exp
b001_wc2.0_bio_30s_03_neotropic_albers_tif_exp
b001_wc2.0_bio_30s_04_neotropic_albers_tif_exp
b001_wc2.0_bio_30s_05_neotropic_albers_tif_exp
b001_wc2.0_bio_30s_06_neotropic_albers_tif_exp
b001_wc2.0_bio_30s_07_neotropic_albers_tif_exp
b001_wc2.0_bio_30s_08_neotropic_albers_tif_exp
b001_wc2.0_bio_30s_09_neotropic_albers_tif_exp
b001_wc2.0_bio_30s_10_neotropic_albers_tif_exp
b001_wc2.0_bio_30s_11_neotropic_albers_tif_exp
b001_wc2.0_bio_30s_12_neotropic_albers_tif_exp
bb001_wc2.0_bio_30s_13_neotropic_albers_tif_exp
bXXX_wc2.0_bio_30s_14_neotropic_albers_tif_exp
bXXX_wc2.0_bio_30s_15_neotropic_albers_tif_exp
bXXX_wc2.0_bio_30s_16_neotropic_albers_tif_exp
bXXX_wc2.0_bio_30s_17_neotropic_albers_tif_exp
bXXX_wc2.0_bio_30s_18_neotropic_albers_tif_exp
bXXX_wc2.0_bio_30s_19_neotropic_albers_tif_exp
# Hydrography
bXXX_water_frequency_2010_1km_tif_exp
#bXXX_water_frequency_2010_30m_tif_exp      <--------------------------------------------------------------
# Landscape
bXXX_dist2drainage_exp        <--------------------------------------------------------------
bXXX_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist   <--------------------------------------------------------------
bXXX_dist2waterbodies_exp     <--------------------------------------------------------------
bXXX_Homogeneity_01_05_uint16_1km_neotropic_albers_tif_exp
#bXXX_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp   <--------------------------------------------------------------
bXXX_local_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_patch_AreaHA
bXXX_local_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_pid
bXXX_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp
#bXXX_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp     <--------------------------------------------------------------
bXXX_Neotropical_Hansen_treecoverlossperyear_2017_30m_tif_exp
bXXX_Neotropical_Hansen_treecoverlossperyear_binary_2017_30m_tif_exp
bXXX_Std_01_05_1km_uint16_neotropic_albers_tif_exp
bXXX_tree_plantation_type_tif_exp
bXXX_tree_plantations_binary_tif_exp
# Topography
bXXX_Elevation_md_SRTM_1km_neotropic_albers_tif_exp
bXXX_Elevation_mn_SRTM_1km_neotropic_albers_tif_exp
bXXX_Neotropic_Earthenv_dem90m_tif_exp
bXXX_Slope_md_SRTM_1km_neotropic_albers_tif_exp
bXXX_Slope_mn_SRTM_1km_neotropic_albers_tif_exp
bXXX_drainage_15s_binary_tif_exp







#' Anthropic
#b001_dist2roads_exp
#' J001
#' b0001_dist2roads_exp
dr <-"D:/GISbases/b0001JaguarAlbersBuffer70/b0001_dist2roads_exp.tif"
dr=raster(dr)
#crs(dr)   
#plot(dr)
dr <- projectRaster(dr, crs="+proj=longlat +datum=WGS84")
#crs(dr)   
dr <- projectRaster(dr, crs=get_crs(J001trk))
#crs(dr)   
#compareCRS(dr,get_crs(J001trk))
writeRaster(dr,filename=file.path("D:/GISUTM/J001/UTM_b0001_dist2roads_exp"), format="GTiff", overwrite=TRUE)
#test <-"D:/GISUTM/J001/UTM_b001_dist2roads_exp.tif";(test=raster(test))
#x11()
#plot(test)


lcover <-"c:/SIG/b115_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif"
land_cover=raster(lcover)
crs(land_cover)   
plot(land_cover)
newproj <- "+proj=longlat +datum=WGS84"
pr1 <- projectRaster(land_cover, crs=newproj)
newproj <- "+proj=utm +zone=21K +south +datum=WGS84 +units=m +no_defs"
land_cover <- projectRaster(land_cover, crs=newproj)
coord.UTM  = SpatialPoints(cbind(SaoBento$utm_x ,SaoBento$utm_y),CRS("+proj=utm +zone=21K +south +datum=WGS84 +units=m +no_defs"))
coord.UTM 
compareCRS(coord.UTM,land_cover)
plot(land_cover)
plot(coord.UTM,add=T)

tc <-"c:/SIG/b115_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif"
tree_cover=raster(tc)
crs(tree_cover)   
plot(tree_cover)
newproj <- "+proj=longlat +datum=WGS84"
pr1 <- projectRaster(tree_cover, crs=newproj)
newproj <- "+proj=utm +zone=21K +south +datum=WGS84 +units=m +no_defs"
tree_cover <- projectRaster(tree_cover, crs=newproj)
coord.UTM = SpatialPoints(cbind(SaoBento$utm_x ,SaoBento$utm_y),CRS("+proj=utm +zone=21K +south +datum=WGS84 +units=m +no_defs"))
coord.UTM 
compareCRS(coord.UTM,tree_cover)
plot(tree_cover)
plot(coord.UTM,add=T)


liv <-"c:/SIG/b115_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif"
livestock=raster(liv)
crs(livestock)   
plot(livestock)
newproj <- "+proj=longlat +datum=WGS84"
pr1 <- projectRaster(livestock, crs=newproj)
newproj <- "+proj=utm +zone=21K +south +datum=WGS84 +units=m +no_defs"
livestock <- projectRaster(livestock, crs=newproj)
coord.UTM = SpatialPoints(cbind(SaoBento$utm_x ,SaoBento$utm_y),CRS("+proj=utm +zone=21K +south +datum=WGS84 +units=m +no_defs"))
coord.UTM
compareCRS(coord.UTM,livestock)
plot(livestock)
plot(coord.UTM,add=T)


hf<-"c:/SIG/b115_human_footprint_2009_1km_tif_exp.tif"
human_foot=raster(hf)
crs(human_foot)   
plot(human_foot)
newproj <- "+proj=longlat +datum=WGS84"
pr1 <- projectRaster(human_foot, crs=newproj)
newproj <- "+proj=utm +zone=21K +south +datum=WGS84 +units=m +no_defs"
human_foot <- projectRaster(human_foot, crs=newproj)
coord.UTM = SpatialPoints(cbind(SaoBento$utm_x ,SaoBento$utm_y),CRS("+proj=utm +zone=21K +south +datum=WGS84 +units=m +no_defs"))
coord.UTM
compareCRS(coord.UTM,human_foot)
plot(human_foot)
plot(coord.UTM,add=T)



el<-"c:/SIG/b115_Elevation_md_SRTM_1km_neotropic_albers_tif_exp.tif"
elevat=raster(el)
crs(elevat)   
plot(elevat)
newproj <- "+proj=longlat +datum=WGS84"
pr1 <- projectRaster(elevat, crs=newproj)
newproj <- "+proj=utm +zone=21K +south +datum=WGS84 +units=m +no_defs"
elevat <- projectRaster(elevat, crs=newproj)
coord.UTM = SpatialPoints(cbind(SaoBento$utm_x ,SaoBento$utm_y),CRS("+proj=utm +zone=21K +south +datum=WGS84 +units=m +no_defs"))
coord.UTM
compareCRS(coord.UTM,elevat)
plot(elevat)
plot(coord.UTM,add=T)
elevat->elevation


w30 <-"c:/SIG/b115_water_frequency_2010_30m_tif_exp.tif"
water30m=raster(w30)
crs(water30m)   
plot(water30m)
newproj <- "+proj=longlat +datum=WGS84"
pr1 <- projectRaster(water30m, crs=newproj)
newproj <- "+proj=utm +zone=21K +south +datum=WGS84 +units=m +no_defs"
water30m <- projectRaster(water30m, crs=newproj)
coord.UTM = SpatialPoints(cbind(SaoBento$utm_x ,SaoBento$utm_y),CRS("+proj=utm +zone=21K +south +datum=WGS84 +units=m +no_defs"))
coord.UTM
compareCRS(coord.UTM,water30m)
plot(water30m)
plot(coord.UTM,add=T)


popD <-"c:/SIG/b115_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif"
pop=raster(popD)
crs(pop)   
plot(pop)
newproj <- "+proj=longlat +datum=WGS84"
pr1 <- projectRaster(pop, crs=newproj)
newproj <- "+proj=utm +zone=21K +south +datum=WGS84 +units=m +no_defs"
pop <- projectRaster(pop, crs=newproj)
coord.UTM = SpatialPoints(cbind(SaoBento$utm_x ,SaoBento$utm_y),CRS("+proj=utm +zone=21K +south +datum=WGS84 +units=m +no_defs"))
coord.UTM
compareCRS(coord.UTM,pop)
plot(pop)
plot(coord.UTM,add=T)


dn<-"c:/SIG/b115_dist2drainage_exp.tif"
distdrain=raster(dn)
crs(distdrain)   
plot(distdrain)
newproj <- "+proj=longlat +datum=WGS84"
pr1 <- projectRaster(distdrain, crs=newproj)
newproj <- "+proj=utm +zone=21K +south +datum=WGS84 +units=m +no_defs"
distdrain <- projectRaster(distdrain, crs=newproj)
coord.UTM = SpatialPoints(cbind(SaoBento$utm_x ,SaoBento$utm_y),CRS("+proj=utm +zone=21K +south +datum=WGS84 +units=m +no_defs"))
coord.UTM
compareCRS(coord.UTM,distdrain)
plot(distdrain)
plot(coord.UTM,add=T)


#dw<-"c:/SIG/b115_dist2waterbodies_exp.tif"
#distwater=raster(dw)
#crs(distwater)   
#plot(distwater)


#dr<-"c:/SIG/b115_dist2roads_exp.tif"
#distroad=raster(dr)
#crs(distroad)   
#plot(distroad)

#dt<-"c:/SIG/b115_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif"
#disttreecover=raster(dt)
#crs(disttreecover)   
#plot(disttreecover)





################################################################################################################################
# dist2roads_exp

#' J81    b081JaguarAlbersBuffer70
#' b081_dist2roads_exp
dr <-"D:/GISbases/b081JaguarAlbersBuffer70/b081_dist2roads_exp.tif"
dr=raster(dr)
#crs(dr)   
#plot(dr)
dr <- projectRaster(dr, crs="+proj=longlat +datum=WGS84")
#crs(dr)   
dr <- projectRaster(dr, crs=get_crs(J81trk))
#crs(dr)   
#compareCRS(dr,get_crs(J81trk))
writeRaster(dr,filename=file.path("D:/GISUTM/J81/UTM_b081_dist2roads_exp"), format="GTiff", overwrite=TRUE)


#' J82
#' b082_dist2roads_exp
dr <-"D:/GISbases/b082JaguarAlbersBuffer70/b082_dist2roads_exp.tif"
dr=raster(dr)
#crs(dr)   
#plot(dr)
dr <- projectRaster(dr, crs="+proj=longlat +datum=WGS84")
#crs(dr)   
dr <- projectRaster(dr, crs=get_crs(J82trk))
#crs(dr)   
#compareCRS(dr,get_crs(J82trk))
writeRaster(dr,filename=file.path("D:/GISUTM/J82/UTM_b082_dist2roads_exp"), format="GTiff", overwrite=TRUE)

#b083_dist2roads_exp
#' J83
#' b083_dist2roads_exp
dr <-"D:/GISbases/b083JaguarAlbersBuffer70/b083_dist2roads_exp.tif"
dr=raster(dr)
#crs(dr)   
#plot(dr)
dr <- projectRaster(dr, crs="+proj=longlat +datum=WGS84")
#crs(dr)   
dr <- projectRaster(dr, crs=get_crs(J83trk))
#crs(dr)   
#compareCRS(dr,get_crs(J83trk))
writeRaster(dr,filename=file.path("D:/GISUTM/J83/UTM_b083_dist2roads_exp"), format="GTiff", overwrite=TRUE)

#b084_dist2roads_exp
#' J84
#' b084_dist2roads_exp
dr <-"D:/GISbases/b084JaguarAlbersBuffer70/b084_dist2roads_exp.tif"
dr=raster(dr)
#crs(dr)   
#plot(dr)
dr <- projectRaster(dr, crs="+proj=longlat +datum=WGS84")
#crs(dr)   
dr <- projectRaster(dr, crs=get_crs(J84trk))
#crs(dr)   
#compareCRS(dr,get_crs(J84trk))
writeRaster(dr,filename=file.path("D:/GISUTM/J84/UTM_b084_dist2roads_exp"), format="GTiff", overwrite=TRUE)

#b085_dist2roads_exp
#' J85
#' b085_dist2roads_exp
dr <-"D:/GISbases/b085JaguarAlbersBuffer70/b085_dist2roads_exp.tif"
dr=raster(dr)
#crs(dr)   
#plot(dr)
dr <- projectRaster(dr, crs="+proj=longlat +datum=WGS84")
#crs(dr)   
dr <- projectRaster(dr, crs=get_crs(J85trk))
#crs(dr)   
#compareCRS(dr,get_crs(J85trk))
writeRaster(dr,filename=file.path("D:/GISUTM/J85/UTM_b085_dist2roads_exp"), format="GTiff", overwrite=TRUE)

#b086_dist2roads_exp
#' J86
#' b086_dist2roads_exp
dr <-"D:/GISbases/b086JaguarAlbersBuffer70/b086_dist2roads_exp.tif"
dr=raster(dr)
#crs(dr)   
#plot(dr)
dr <- projectRaster(dr, crs="+proj=longlat +datum=WGS84")
#crs(dr)   
dr <- projectRaster(dr, crs=get_crs(J86trk))
#crs(dr)   
#compareCRS(dr,get_crs(J86trk))
writeRaster(dr,filename=file.path("D:/GISUTM/J86/UTM_b086_dist2roads_exp"), format="GTiff", overwrite=TRUE)

#b087_dist2roads_exp
#' J87
#' b087_dist2roads_exp
dr <-"D:/GISbases/b087JaguarAlbersBuffer70/b087_dist2roads_exp.tif"
dr=raster(dr)
#crs(dr)   
#plot(dr)
dr <- projectRaster(dr, crs="+proj=longlat +datum=WGS84")
#crs(dr)   
dr <- projectRaster(dr, crs=get_crs(J87trk))
#crs(dr)   
#compareCRS(dr,get_crs(J87trk))
writeRaster(dr,filename=file.path("D:/GISUTM/J87/UTM_b087_dist2roads_exp"), format="GTiff", overwrite=TRUE)

#b088_dist2roads_exp
#' J88
#' b088_dist2roads_exp
dr <-"D:/GISbases/b088JaguarAlbersBuffer70/b088_dist2roads_exp.tif"
dr=raster(dr)
#crs(dr)   
#plot(dr)
dr <- projectRaster(dr, crs="+proj=longlat +datum=WGS84")
#crs(dr)   
dr <- projectRaster(dr, crs=get_crs(J88trk))
#crs(dr)   
#compareCRS(dr,get_crs(J88trk))
writeRaster(dr,filename=file.path("D:/GISUTM/J88/UTM_b088_dist2roads_exp"), format="GTiff", overwrite=TRUE)


#b089_dist2roads_exp
#' J89
#' b089_dist2roads_exp
dr <-"D:/GISbases/b089JaguarAlbersBuffer70/b089_dist2roads_exp.tif"
dr=raster(dr)
#crs(dr)   
#plot(dr)
dr <- projectRaster(dr, crs="+proj=longlat +datum=WGS84")
#crs(dr)   
dr <- projectRaster(dr, crs=get_crs(J89trk))
#crs(dr)   
#compareCRS(dr,get_crs(J89trk))
writeRaster(dr,filename=file.path("D:/GISUTM/J89/UTM_b089_dist2roads_exp"), format="GTiff", overwrite=TRUE)


#b091_dist2roads_exp
#' J91
#' b091_dist2roads_exp
dr <-"D:/GISbases/b091JaguarAlbersBuffer70/b091_dist2roads_exp.tif"
dr=raster(dr)
#crs(dr)   
#plot(dr)
dr <- projectRaster(dr, crs="+proj=longlat +datum=WGS84")
#crs(dr)   
dr <- projectRaster(dr, crs=get_crs(J91trk))
#crs(dr)   
#compareCRS(dr,get_crs(J91trk))
writeRaster(dr,filename=file.path("D:/GISUTM/J91/UTM_b091_dist2roads_exp"), format="GTiff", overwrite=TRUE)

#' J92
#' b092_dist2roads_exp
dr <-"D:/GISbases/b092JaguarAlbersBuffer70/b092_dist2roads_exp.tif"
dr=raster(dr)
#crs(dr)   
#plot(dr)
dr <- projectRaster(dr, crs="+proj=longlat +datum=WGS84")
#crs(dr)   
dr <- projectRaster(dr, crs=get_crs(J92trk))
#crs(dr)   
#compareCRS(dr,get_crs(J92trk))
writeRaster(dr,filename=file.path("D:/GISUTM/J92/UTM_b092_dist2roads_exp"), format="GTiff", overwrite=TRUE)

#' J93
#' b093_dist2roads_exp
dr <-"D:/GISbases/b093JaguarAlbersBuffer70/b093_dist2roads_exp.tif"
dr=raster(dr)
crs(dr)   
x11()
plot(dr)
dr <- projectRaster(dr, crs="+proj=longlat +datum=WGS84")
#crs(dr)   
dr <- projectRaster(dr, crs=get_crs(J93trk))
#crs(dr)   
#compareCRS(dr,get_crs(J93trk))
writeRaster(dr,filename=file.path("D:/GISUTM/J93/UTM_b093_dist2roads_exp"), format="GTiff", overwrite=TRUE)

#' J94
#' b094_dist2roads_exp
dr <-"D:/GISbases/b094JaguarAlbersBuffer70/b094_dist2roads_exp.tif"
dr=raster(dr)
crs(dr)   
plot(dr)
dr <- projectRaster(dr, crs="+proj=longlat +datum=WGS84")
#crs(dr)   
dr <- projectRaster(dr, crs=get_crs(J94trk))
#crs(dr)   
#compareCRS(dr,get_crs(J94trk))
writeRaster(dr,filename=file.path("D:/GISUTM/J94/UTM_b094_dist2roads_exp"), format="GTiff", overwrite=TRUE)

#' J95
#' b095_dist2roads_exp
dr <-"D:/GISbases/b095JaguarAlbersBuffer70/b095_dist2roads_exp.tif"
dr=raster(dr)
#crs(dr)   
#plot(dr)
dr <- projectRaster(dr, crs="+proj=longlat +datum=WGS84")
#crs(dr)   
dr <- projectRaster(dr, crs=get_crs(J95trk))
#crs(dr)   
#compareCRS(dr,get_crs(J95trk))
writeRaster(dr,filename=file.path("D:/GISUTM/J95/UTM_b095_dist2roads_exp"), format="GTiff", overwrite=TRUE)

#' J96
#' b096_dist2roads_exp
dr <-"D:/GISbases/b096JaguarAlbersBuffer70/b096_dist2roads_exp.tif"
dr=raster(dr)
#crs(dr)   
#plot(dr)
dr <- projectRaster(dr, crs="+proj=longlat +datum=WGS84")
#crs(dr)   
dr <- projectRaster(dr, crs=get_crs(J96trk))
#crs(dr)   
#compareCRS(dr,get_crs(J96trk))
writeRaster(dr,filename=file.path("D:/GISUTM/J96/UTM_b096_dist2roads_exp"), format="GTiff", overwrite=TRUE)


#' J97
#' b097_dist2roads_exp
dr <-"D:/GISbases/b097JaguarAlbersBuffer70/b097_dist2roads_exp.tif"
dr=raster(dr)
#crs(dr)   
#plot(dr)
dr <- projectRaster(dr, crs="+proj=longlat +datum=WGS84")
#crs(dr)   
dr <- projectRaster(dr, crs=get_crs(J97trk))
#crs(dr)   
#compareCRS(dr,get_crs(J97trk))
writeRaster(dr,filename=file.path("D:/GISUTM/J97/UTM_b097_dist2roads_exp"), format="GTiff", overwrite=TRUE)

#' J98
#' b098_dist2roads_exp
dr <-"D:/GISbases/b098JaguarAlbersBuffer70/b098_dist2roads_exp.tif"
dr=raster(dr)
#crs(dr)   
#plot(dr)
dr <- projectRaster(dr, crs="+proj=longlat +datum=WGS84")
#crs(dr)   
dr <- projectRaster(dr, crs=get_crs(J98trk))
#crs(dr)   
#compareCRS(dr,get_crs(J98trk))
writeRaster(dr,filename=file.path("D:/GISUTM/J98/UTM_b098_dist2roads_exp"), format="GTiff", overwrite=TRUE)


#' J99
#' b099_dist2roads_exp
dr <-"D:/GISbases/b099JaguarAlbersBuffer70/b099_dist2roads_exp.tif"
dr=raster(dr)
#crs(dr)   
#plot(dr)
dr <- projectRaster(dr, crs="+proj=longlat +datum=WGS84")
#crs(dr)   
dr <- projectRaster(dr, crs=get_crs(J99trk))
#crs(dr)   
#compareCRS(dr,get_crs(J99trk))
writeRaster(dr,filename=file.path("D:/GISUTM/J99/UTM_b099_dist2roads_exp"), format="GTiff", overwrite=TRUE)

#' J100
#' b100_dist2roads_exp
dr <-"D:/GISbases/b100JaguarAlbersBuffer70/b100_dist2roads_exp.tif"
dr=raster(dr)
#crs(dr)   
#plot(dr)
dr <- projectRaster(dr, crs="+proj=longlat +datum=WGS84")
#crs(dr)   
dr <- projectRaster(dr, crs=get_crs(J100trk))
#crs(dr)   
#compareCRS(dr,get_crs(J100trk))
writeRaster(dr,filename=file.path("D:/GISUTM/J100/UTM_b100_dist2roads_exp"), format="GTiff", overwrite=TRUE)


#' J101
#' b101_dist2roads_exp
dr <-"D:/GISbases/b101JaguarAlbersBuffer70/b101_dist2roads_exp.tif"
dr=raster(dr)
#crs(dr)   
#plot(dr)
dr <- projectRaster(dr, crs="+proj=longlat +datum=WGS84")
#crs(dr)   
dr <- projectRaster(dr, crs=get_crs(J101trk))
#crs(dr)   
#compareCRS(dr,get_crs(J101trk))
writeRaster(dr,filename=file.path("D:/GISUTM/J101/UTM_b101_dist2roads_exp"), format="GTiff", overwrite=TRUE)


#' J102
#' b102_dist2roads_exp
dr <-"D:/GISbases/b102JaguarAlbersBuffer70/b102_dist2roads_exp.tif"
dr=raster(dr)
#crs(dr)   
#plot(dr)
dr <- projectRaster(dr, crs="+proj=longlat +datum=WGS84")
#crs(dr)   
dr <- projectRaster(dr, crs=get_crs(J102trk))
#crs(dr)   
#compareCRS(dr,get_crs(J102trk))
writeRaster(dr,filename=file.path("D:/GISUTM/J102/UTM_b102_dist2roads_exp"), format="GTiff", overwrite=TRUE)

#' J103
#' b103_dist2roads_exp
dr <-"D:/GISbases/b103JaguarAlbersBuffer70/b103_dist2roads_exp.tif"
dr=raster(dr)
#crs(dr)   
#plot(dr)
dr <- projectRaster(dr, crs="+proj=longlat +datum=WGS84")
#crs(dr)   
dr <- projectRaster(dr, crs=get_crs(J103trk))
#crs(dr)   
#compareCRS(dr,get_crs(J103trk))
writeRaster(dr,filename=file.path("D:/GISUTM/J103/UTM_b103_dist2roads_exp"), format="GTiff", overwrite=TRUE)

#' J104
#' b104_dist2roads_exp
dr <-"D:/GISbases/b104JaguarAlbersBuffer70/b104_dist2roads_exp.tif"
dr=raster(dr)
#crs(dr)   
#plot(dr)
dr <- projectRaster(dr, crs="+proj=longlat +datum=WGS84")
#crs(dr)   
dr <- projectRaster(dr, crs=get_crs(J104trk))
#crs(dr)   
#compareCRS(dr,get_crs(J104trk))
writeRaster(dr,filename=file.path("D:/GISUTM/J104/UTM_b104_dist2roads_exp"), format="GTiff", overwrite=TRUE)

#' J105
#' b105_dist2roads_exp
dr <-"D:/GISbases/b105JaguarAlbersBuffer70/b105_dist2roads_exp.tif"
dr=raster(dr)
#crs(dr)   
#plot(dr)
dr <- projectRaster(dr, crs="+proj=longlat +datum=WGS84")
#crs(dr)   
dr <- projectRaster(dr, crs=get_crs(J105trk))
#crs(dr)   
#compareCRS(dr,get_crs(J105trk))
writeRaster(dr,filename=file.path("D:/GISUTM/J105/UTM_b105_dist2roads_exp"), format="GTiff", overwrite=TRUE)

#' J106
#' b106_dist2roads_exp
dr <-"D:/GISbases/b106JaguarAlbersBuffer70/b106_dist2roads_exp.tif"
dr=raster(dr)
#crs(dr)   
#plot(dr)
dr <- projectRaster(dr, crs="+proj=longlat +datum=WGS84")
#crs(dr)   
dr <- projectRaster(dr, crs=get_crs(J106trk))
#crs(dr)   
#compareCRS(dr,get_crs(J106trk))
writeRaster(dr,filename=file.path("D:/GISUTM/J106/UTM_b106_dist2roads_exp"), format="GTiff", overwrite=TRUE)

#' J107
#' b107_dist2roads_exp
dr <-"D:/GISbases/b107JaguarAlbersBuffer70/b107_dist2roads_exp.tif"
dr=raster(dr)
#crs(dr)   
#plot(dr)
dr <- projectRaster(dr, crs="+proj=longlat +datum=WGS84")
#crs(dr)   
dr <- projectRaster(dr, crs=get_crs(J107trk))
#crs(dr)   
#compareCRS(dr,get_crs(J107trk))
writeRaster(dr,filename=file.path("D:/GISUTM/J107/UTM_b107_dist2roads_exp"), format="GTiff", overwrite=TRUE)


#' J108
#' b108_dist2roads_exp
dr <-"D:/GISbases/b108JaguarAlbersBuffer70/b108_dist2roads_exp.tif"
dr=raster(dr)
#crs(dr)   
#plot(dr)
dr <- projectRaster(dr, crs="+proj=longlat +datum=WGS84")
#crs(dr)   
dr <- projectRaster(dr, crs=get_crs(J108trk))
#crs(dr)   
#compareCRS(dr,get_crs(J108trk))
writeRaster(dr,filename=file.path("D:/GISUTM/J108/UTM_b108_dist2roads_exp"), format="GTiff", overwrite=TRUE)


#' J109
#' b109_dist2roads_exp
dr <-"D:/GISbases/b109JaguarAlbersBuffer70/b109_dist2roads_exp.tif"
dr=raster(dr)
#crs(dr)   
#plot(dr)
dr <- projectRaster(dr, crs="+proj=longlat +datum=WGS84")
#crs(dr)   
dr <- projectRaster(dr, crs=get_crs(J109trk))
#crs(dr)   
#compareCRS(dr,get_crs(J109trk))
writeRaster(dr,filename=file.path("D:/GISUTM/J109/UTM_b109_dist2roads_exp"), format="GTiff", overwrite=TRUE)

#' J110
#' b110_dist2roads_exp
dr <-"D:/GISbases/b110JaguarAlbersBuffer70/b110_dist2roads_exp.tif"
dr=raster(dr)
#crs(dr)   
#plot(dr)
dr <- projectRaster(dr, crs="+proj=longlat +datum=WGS84")
#crs(dr)   
dr <- projectRaster(dr, crs=get_crs(J110trk))
#crs(dr)   
#compareCRS(dr,get_crs(J110trk))
writeRaster(dr,filename=file.path("D:/GISUTM/J110/UTM_b110_dist2roads_exp"), format="GTiff", overwrite=TRUE)

#' J111
#' b111_dist2roads_exp
dr <-"D:/GISbases/b111JaguarAlbersBuffer70/b111_dist2roads_exp.tif"
dr=raster(dr)
#crs(dr)   
#plot(dr)
dr <- projectRaster(dr, crs="+proj=longlat +datum=WGS84")
#crs(dr)   
dr <- projectRaster(dr, crs=get_crs(J111trk))
#crs(dr)   
#compareCRS(dr,get_crs(J111trk))
writeRaster(dr,filename=file.path("D:/GISUTM/J111/UTM_b111_dist2roads_exp"), format="GTiff", overwrite=TRUE)


#' J112
#' b112_dist2roads_exp
dr <-"D:/GISbases/b112JaguarAlbersBuffer70/b112_dist2roads_exp.tif"
dr=raster(dr)
#crs(dr)   
#plot(dr)
dr <- projectRaster(dr, crs="+proj=longlat +datum=WGS84")
#crs(dr)   
dr <- projectRaster(dr, crs=get_crs(J112trk))
#crs(dr)   
#compareCRS(dr,get_crs(J112trk))
writeRaster(dr,filename=file.path("D:/GISUTM/J112/UTM_b112_dist2roads_exp"), format="GTiff", overwrite=TRUE)


#' J113
#' b113_dist2roads_exp
dr <-"D:/GISbases/b113JaguarAlbersBuffer70/b113_dist2roads_exp.tif"
dr=raster(dr)
#crs(dr)   
#plot(dr)
dr <- projectRaster(dr, crs="+proj=longlat +datum=WGS84")
#crs(dr)   
dr <- projectRaster(dr, crs=get_crs(J113trk))
#crs(dr)   
#compareCRS(dr,get_crs(J113trk))
writeRaster(dr,filename=file.path("D:/GISUTM/J113/UTM_b113_dist2roads_exp"), format="GTiff", overwrite=TRUE)

#' J114
#' b114_dist2roads_exp
dr <-"D:/GISbases/b114JaguarAlbersBuffer70/b114_dist2roads_exp.tif"
dr=raster(dr)
#crs(dr)   
#plot(dr)
dr <- projectRaster(dr, crs="+proj=longlat +datum=WGS84")
#crs(dr)   
dr <- projectRaster(dr, crs=get_crs(J114trk))
#crs(dr)   
#compareCRS(dr,get_crs(J114trk))
writeRaster(dr,filename=file.path("D:/GISUTM/J114/UTM_b114_dist2roads_exp"), format="GTiff", overwrite=TRUE)

#' J115
#' b115_dist2roads_exp
dr <-"D:/GISbases/b115JaguarAlbersBuffer70/b115_dist2roads_exp.tif"
dr=raster(dr)
#crs(dr)   
#plot(dr)
dr <- projectRaster(dr, crs="+proj=longlat +datum=WGS84")
#crs(dr)   
dr <- projectRaster(dr, crs=get_crs(J115trk))
#crs(dr)   
#compareCRS(dr,get_crs(J115trk))
writeRaster(dr,filename=file.path("D:/GISUTM/J115/UTM_b115_dist2roads_exp"), format="GTiff", overwrite=TRUE)

#' J116
#' b116_dist2roads_exp
dr <-"D:/GISbases/b116JaguarAlbersBuffer70/b116_dist2roads_exp.tif"
dr=raster(dr)
#crs(dr)   
#plot(dr)
dr <- projectRaster(dr, crs="+proj=longlat +datum=WGS84")
#crs(dr)   
dr <- projectRaster(dr, crs=get_crs(J116trk))
#crs(dr)   
#compareCRS(dr,get_crs(J116trk))
writeRaster(dr,filename=file.path("D:/GISUTM/J116/UTM_b116_dist2roads_exp"), format="GTiff", overwrite=TRUE)

#' J117
#' b117_dist2roads_exp
dr <-"D:/GISbases/b117JaguarAlbersBuffer70/b117_dist2roads_exp.tif"
dr=raster(dr)
#crs(dr)   
#plot(dr)
dr <- projectRaster(dr, crs="+proj=longlat +datum=WGS84")
#crs(dr)   
dr <- projectRaster(dr, crs=get_crs(J117trk))
#crs(dr)   
#compareCRS(dr,get_crs(J117trk))
writeRaster(dr,filename=file.path("D:/GISUTM/J117/UTM_b117_dist2roads_exp"), format="GTiff", overwrite=TRUE)


###########################################################################################################################
# human_footprint_2009_1km_tif_exp

#' J1
#' b001_human_footprint_2009_1km_tif_exp
hf <-"D:/GISbases/b001JaguarAlbersBuffer70/b001_human_footprint_2009_1km_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J1trk))
#crs(hf)   
#compareCRS(hf,get_crs(J1trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J1/UTM_b001_human_footprint_2009_1km_tif_exp"), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J1/UTM_b001_human_footprint_2009_1km_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J2
#' b002_human_footprint_2009_1km_tif_exp
hf <-"D:/GISbases/b002JaguarAlbersBuffer70/b002_human_footprint_2009_1km_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J2trk))
#crs(hf)   
#compareCRS(hf,get_crs(J2trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J2/UTM_b002_human_footprint_2009_1km_tif_exp"), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J2/UTM_b002_human_footprint_2009_1km_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J3
#' b003_human_footprint_2009_1km_tif_exp
hf <-"D:/GISbases/b003JaguarAlbersBuffer70/b003_human_footprint_2009_1km_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J3trk))
#crs(hf)   
#compareCRS(hf,get_crs(J3trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J3/UTM_b003_human_footprint_2009_1km_tif_exp"), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J3/UTM_b003_human_footprint_2009_1km_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J4
#' b004_human_footprint_2009_1km_tif_exp
hf <-"D:/GISbases/b004JaguarAlbersBuffer70/b004_human_footprint_2009_1km_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J4trk))
#crs(hf)   
#compareCRS(hf,get_crs(J4trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J4/UTM_b004_human_footprint_2009_1km_tif_exp"), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J4/UTM_b004_human_footprint_2009_1km_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J5
#' b005_human_footprint_2009_1km_tif_exp
hf <-"D:/GISbases/b005JaguarAlbersBuffer70/b005_human_footprint_2009_1km_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J5trk))
#crs(hf)   
#compareCRS(hf,get_crs(J5trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J5/UTM_b005_human_footprint_2009_1km_tif_exp"), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J5/UTM_b005_human_footprint_2009_1km_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J6
#' b006_human_footprint_2009_1km_tif_exp
hf <-"D:/GISbases/b006JaguarAlbersBuffer70/b006_human_footprint_2009_1km_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J6trk))
#crs(hf)   
#compareCRS(hf,get_crs(J6trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J6/UTM_b006_human_footprint_2009_1km_tif_exp"), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J6/UTM_b006_human_footprint_2009_1km_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J7
#' b007_human_footprint_2009_1km_tif_exp
hf <-"D:/GISbases/b007JaguarAlbersBuffer70/b007_human_footprint_2009_1km_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J7trk))
#crs(hf)   
#compareCRS(hf,get_crs(J7trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J7/UTM_b007_human_footprint_2009_1km_tif_exp"), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J7/UTM_b007_human_footprint_2009_1km_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J8
#' b008_human_footprint_2009_1km_tif_exp
hf <-"D:/GISbases/b008JaguarAlbersBuffer70/b008_human_footprint_2009_1km_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J8trk))
#crs(hf)   
#compareCRS(hf,get_crs(J8trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J8/UTM_b008_human_footprint_2009_1km_tif_exp"), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J8/UTM_b008_human_footprint_2009_1km_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J9
#' b009_human_footprint_2009_1km_tif_exp
hf <-"D:/GISbases/b009JaguarAlbersBuffer70/b009_human_footprint_2009_1km_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J9trk))
#crs(hf)   
#compareCRS(hf,get_crs(J9trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J9/UTM_b009_human_footprint_2009_1km_tif_exp"), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J9/UTM_b009_human_footprint_2009_1km_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J10
#' b010_human_footprint_2009_1km_tif_exp
hf <-"D:/GISbases/b010JaguarAlbersBuffer70/b010_human_footprint_2009_1km_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J10trk))
#crs(hf)   
#compareCRS(hf,get_crs(J10trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J10/UTM_b010_human_footprint_2009_1km_tif_exp"), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J10/UTM_b010_human_footprint_2009_1km_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J11
#' b011_human_footprint_2009_1km_tif_exp
hf <-"D:/GISbases/b011JaguarAlbersBuffer70/b011_human_footprint_2009_1km_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J11trk))
#crs(hf)   
#compareCRS(hf,get_crs(J11trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J11/UTM_b011_human_footprint_2009_1km_tif_exp"), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J11/UTM_b011_human_footprint_2009_1km_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J12
#' b012_human_footprint_2009_1km_tif_exp
hf <-"D:/GISbases/b012JaguarAlbersBuffer70/b012_human_footprint_2009_1km_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J12trk))
#crs(hf)   
#compareCRS(hf,get_crs(J12trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J12/UTM_b012_human_footprint_2009_1km_tif_exp"), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J12/UTM_b012_human_footprint_2009_1km_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J13
#' b013_human_footprint_2009_1km_tif_exp
hf <-"D:/GISbases/b013JaguarAlbersBuffer70/b013_human_footprint_2009_1km_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J13trk))
#crs(hf)   
#compareCRS(hf,get_crs(J13trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J13/UTM_b013_human_footprint_2009_1km_tif_exp"), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J13/UTM_b013_human_footprint_2009_1km_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J14
#' b014_human_footprint_2009_1km_tif_exp
hf <-"D:/GISbases/b014JaguarAlbersBuffer70/b014_human_footprint_2009_1km_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J14trk))
#crs(hf)   
#compareCRS(hf,get_crs(J14trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J14/UTM_b014_human_footprint_2009_1km_tif_exp"), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J14/UTM_b014_human_footprint_2009_1km_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J15
#' b015_human_footprint_2009_1km_tif_exp
hf <-"D:/GISbases/b015JaguarAlbersBuffer70/b015_human_footprint_2009_1km_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J15trk))
#crs(hf)   
#compareCRS(hf,get_crs(J15trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J15/UTM_b015_human_footprint_2009_1km_tif_exp"), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J15/UTM_b015_human_footprint_2009_1km_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J16
#' b016_human_footprint_2009_1km_tif_exp
hf <-"D:/GISbases/b016JaguarAlbersBuffer70/b016_human_footprint_2009_1km_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J16trk))
#crs(hf)   
#compareCRS(hf,get_crs(J16trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J16/UTM_b016_human_footprint_2009_1km_tif_exp"), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J16/UTM_b016_human_footprint_2009_1km_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J17
#' b017_human_footprint_2009_1km_tif_exp
hf <-"D:/GISbases/b017JaguarAlbersBuffer70/b017_human_footprint_2009_1km_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J17trk))
#crs(hf)   
#compareCRS(hf,get_crs(J17trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J17/UTM_b017_human_footprint_2009_1km_tif_exp"), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J17/UTM_b017_human_footprint_2009_1km_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J18
#' b018_human_footprint_2009_1km_tif_exp
hf <-"D:/GISbases/b018JaguarAlbersBuffer70/b018_human_footprint_2009_1km_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J18trk))
#crs(hf)   
#compareCRS(hf,get_crs(J18trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J18/UTM_b018_human_footprint_2009_1km_tif_exp"), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J18/UTM_b018_human_footprint_2009_1km_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J19
#' b019_human_footprint_2009_1km_tif_exp
hf <-"D:/GISbases/b019JaguarAlbersBuffer70/b019_human_footprint_2009_1km_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J19trk))
#crs(hf)   
#compareCRS(hf,get_crs(J19trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J19/UTM_b019_human_footprint_2009_1km_tif_exp"), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J19/UTM_b019_human_footprint_2009_1km_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J20
#' b020_human_footprint_2009_1km_tif_exp
hf <-"D:/GISbases/b020JaguarAlbersBuffer70/b020_human_footprint_2009_1km_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J20trk))
#crs(hf)   
#compareCRS(hf,get_crs(J20trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J20/UTM_b020_human_footprint_2009_1km_tif_exp"), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J20/UTM_b020_human_footprint_2009_1km_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J21
#' b021_human_footprint_2009_1km_tif_exp
hf <-"D:/GISbases/b021JaguarAlbersBuffer70/b021_human_footprint_2009_1km_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J21trk))
#crs(hf)   
#compareCRS(hf,get_crs(J21trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J21/UTM_b021_human_footprint_2009_1km_tif_exp"), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J21/UTM_b021_human_footprint_2009_1km_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J22
#' b022_human_footprint_2009_1km_tif_exp
hf <-"D:/GISbases/b022JaguarAlbersBuffer70/b022_human_footprint_2009_1km_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J22trk))
#crs(hf)   
#compareCRS(hf,get_crs(J22trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J22/UTM_b022_human_footprint_2009_1km_tif_exp"), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J22/UTM_b022_human_footprint_2009_1km_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J23
#' b023_human_footprint_2009_1km_tif_exp
hf <-"D:/GISbases/b023JaguarAlbersBuffer70/b023_human_footprint_2009_1km_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J23trk))
#crs(hf)   
#compareCRS(hf,get_crs(J23trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J23/UTM_b023_human_footprint_2009_1km_tif_exp"), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J23/UTM_b023_human_footprint_2009_1km_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J24
#' b024_human_footprint_2009_1km_tif_exp
hf <-"D:/GISbases/b024JaguarAlbersBuffer70/b024_human_footprint_2009_1km_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J24trk))
#crs(hf)   
#compareCRS(hf,get_crs(J24trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J24/UTM_b024_human_footprint_2009_1km_tif_exp"), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J24/UTM_b024_human_footprint_2009_1km_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J25
#' b025_human_footprint_2009_1km_tif_exp
hf <-"D:/GISbases/b025JaguarAlbersBuffer70/b025_human_footprint_2009_1km_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J25trk))
#crs(hf)   
#compareCRS(hf,get_crs(J25trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J25/UTM_b025_human_footprint_2009_1km_tif_exp"), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J25/UTM_b025_human_footprint_2009_1km_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J26
#' b026_human_footprint_2009_1km_tif_exp
hf <-"D:/GISbases/b026JaguarAlbersBuffer70/b026_human_footprint_2009_1km_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J26trk))
#crs(hf)   
#compareCRS(hf,get_crs(J26trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J26/UTM_b026_human_footprint_2009_1km_tif_exp"), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J26/UTM_b026_human_footprint_2009_1km_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J27
#' b027_human_footprint_2009_1km_tif_exp
hf <-"D:/GISbases/b027JaguarAlbersBuffer70/b027_human_footprint_2009_1km_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J27trk))
#crs(hf)   
#compareCRS(hf,get_crs(J27trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J27/UTM_b027_human_footprint_2009_1km_tif_exp"), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J27/UTM_b027_human_footprint_2009_1km_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J28
#' b028_human_footprint_2009_1km_tif_exp
hf <-"D:/GISbases/b028JaguarAlbersBuffer70/b028_human_footprint_2009_1km_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J28trk))
#crs(hf)   
#compareCRS(hf,get_crs(J28trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J28/UTM_b028_human_footprint_2009_1km_tif_exp"), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J28/UTM_b028_human_footprint_2009_1km_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J29
#' b029_human_footprint_2009_1km_tif_exp
hf <-"D:/GISbases/b029JaguarAlbersBuffer70/b029_human_footprint_2009_1km_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J29trk))
#crs(hf)   
#compareCRS(hf,get_crs(J29trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J29/UTM_b029_human_footprint_2009_1km_tif_exp"), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J29/UTM_b029_human_footprint_2009_1km_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J30
#' b030_human_footprint_2009_1km_tif_exp
hf <-"D:/GISbases/b030JaguarAlbersBuffer70/b030_human_footprint_2009_1km_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J30trk))
#crs(hf)   
#compareCRS(hf,get_crs(J30trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J30/UTM_b030_human_footprint_2009_1km_tif_exp"), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J30/UTM_b030_human_footprint_2009_1km_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J30
#' b030_human_footprint_2009_1km_tif_exp
hf <-"D:/GISbases/b030JaguarAlbersBuffer70/b030_human_footprint_2009_1km_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J30trk))
#crs(hf)   
#compareCRS(hf,get_crs(J30trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J30/UTM_b030_human_footprint_2009_1km_tif_exp"), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J30/UTM_b030_human_footprint_2009_1km_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J31
#' b031_human_footprint_2009_1km_tif_exp
hf <-"D:/GISbases/b031JaguarAlbersBuffer70/b031_human_footprint_2009_1km_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J31trk))
#crs(hf)   
#compareCRS(hf,get_crs(J31trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J31/UTM_b031_human_footprint_2009_1km_tif_exp"), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J31/UTM_b031_human_footprint_2009_1km_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J32
#' b032_human_footprint_2009_1km_tif_exp
hf <-"D:/GISbases/b032JaguarAlbersBuffer70/b032_human_footprint_2009_1km_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J32trk))
#crs(hf)   
#compareCRS(hf,get_crs(J32trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J32/UTM_b032_human_footprint_2009_1km_tif_exp"), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J32/UTM_b032_human_footprint_2009_1km_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J33
#' b033_human_footprint_2009_1km_tif_exp
hf <-"D:/GISbases/b033JaguarAlbersBuffer70/b033_human_footprint_2009_1km_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J33trk))
#crs(hf)   
#compareCRS(hf,get_crs(J33trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J33/UTM_b033_human_footprint_2009_1km_tif_exp"), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J33/UTM_b033_human_footprint_2009_1km_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J34
#' b034_human_footprint_2009_1km_tif_exp
hf <-"D:/GISbases/b034JaguarAlbersBuffer70/b034_human_footprint_2009_1km_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J34trk))
#crs(hf)   
#compareCRS(hf,get_crs(J34trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J34/UTM_b034_human_footprint_2009_1km_tif_exp"), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J34/UTM_b034_human_footprint_2009_1km_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J35
#' b035_human_footprint_2009_1km_tif_exp
hf <-"D:/GISbases/b035JaguarAlbersBuffer70/b035_human_footprint_2009_1km_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J35trk))
#crs(hf)   
#compareCRS(hf,get_crs(J35trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J35/UTM_b035_human_footprint_2009_1km_tif_exp"), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J35/UTM_b035_human_footprint_2009_1km_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J36
#' b036_human_footprint_2009_1km_tif_exp
hf <-"D:/GISbases/b036JaguarAlbersBuffer70/b036_human_footprint_2009_1km_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J36trk))
#crs(hf)   
#compareCRS(hf,get_crs(J36trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J36/UTM_b036_human_footprint_2009_1km_tif_exp"), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J36/UTM_b036_human_footprint_2009_1km_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J37
#' b037_human_footprint_2009_1km_tif_exp
hf <-"D:/GISbases/b037JaguarAlbersBuffer70/b037_human_footprint_2009_1km_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J37trk))
#crs(hf)   
#compareCRS(hf,get_crs(J37trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J37/UTM_b037_human_footprint_2009_1km_tif_exp"), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J37/UTM_b037_human_footprint_2009_1km_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J38
#' b038_human_footprint_2009_1km_tif_exp
hf <-"D:/GISbases/b038JaguarAlbersBuffer70/b038_human_footprint_2009_1km_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J38trk))
#crs(hf)   
#compareCRS(hf,get_crs(J38trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J38/UTM_b038_human_footprint_2009_1km_tif_exp"), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J38/UTM_b038_human_footprint_2009_1km_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J39
#' b039_human_footprint_2009_1km_tif_exp
hf <-"D:/GISbases/b039JaguarAlbersBuffer70/b039_human_footprint_2009_1km_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J39trk))
#crs(hf)   
#compareCRS(hf,get_crs(J39trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J39/UTM_b039_human_footprint_2009_1km_tif_exp"), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J39/UTM_b039_human_footprint_2009_1km_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J40
#' b040_human_footprint_2009_1km_tif_exp
hf <-"D:/GISbases/b040JaguarAlbersBuffer70/b040_human_footprint_2009_1km_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J40trk))
#crs(hf)   
#compareCRS(hf,get_crs(J40trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J40/UTM_b040_human_footprint_2009_1km_tif_exp"), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J40/UTM_b040_human_footprint_2009_1km_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J41
#' b041_human_footprint_2009_1km_tif_exp
hf <-"D:/GISbases/b041JaguarAlbersBuffer70/b041_human_footprint_2009_1km_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J41trk))
#crs(hf)   
#compareCRS(hf,get_crs(J41trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J41/UTM_b041_human_footprint_2009_1km_tif_exp"), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J41/UTM_b041_human_footprint_2009_1km_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J42
#' b042_human_footprint_2009_1km_tif_exp
hf <-"D:/GISbases/b042JaguarAlbersBuffer70/b042_human_footprint_2009_1km_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J42trk))
#crs(hf)   
#compareCRS(hf,get_crs(J42trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J42/UTM_b042_human_footprint_2009_1km_tif_exp"), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J42/UTM_b042_human_footprint_2009_1km_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J43
#' b043_human_footprint_2009_1km_tif_exp
hf <-"D:/GISbases/b043JaguarAlbersBuffer70/b043_human_footprint_2009_1km_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J43trk))
#crs(hf)   
#compareCRS(hf,get_crs(J43trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J43/UTM_b043_human_footprint_2009_1km_tif_exp"), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J43/UTM_b043_human_footprint_2009_1km_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J44
#' b044_human_footprint_2009_1km_tif_exp
hf <-"D:/GISbases/b044JaguarAlbersBuffer70/b044_human_footprint_2009_1km_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J44trk))
#crs(hf)   
#compareCRS(hf,get_crs(J44trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J44/UTM_b044_human_footprint_2009_1km_tif_exp"), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J44/UTM_b044_human_footprint_2009_1km_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J45
#' b045_human_footprint_2009_1km_tif_exp
hf <-"D:/GISbases/b045JaguarAlbersBuffer70/b045_human_footprint_2009_1km_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J45trk))
#crs(hf)   
#compareCRS(hf,get_crs(J45trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J45/UTM_b045_human_footprint_2009_1km_tif_exp"), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J45/UTM_b045_human_footprint_2009_1km_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J46
#' b046_human_footprint_2009_1km_tif_exp
hf <-"D:/GISbases/b046JaguarAlbersBuffer70/b046_human_footprint_2009_1km_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J46trk))
#crs(hf)   
#compareCRS(hf,get_crs(J46trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J46/UTM_b046_human_footprint_2009_1km_tif_exp"), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J46/UTM_b046_human_footprint_2009_1km_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J47
#' b047_human_footprint_2009_1km_tif_exp
hf <-"D:/GISbases/b047JaguarAlbersBuffer70/b047_human_footprint_2009_1km_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J47trk))
#crs(hf)   
#compareCRS(hf,get_crs(J47trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J47/UTM_b047_human_footprint_2009_1km_tif_exp"), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J47/UTM_b047_human_footprint_2009_1km_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J48
#' b048_human_footprint_2009_1km_tif_exp
hf <-"D:/GISbases/b048JaguarAlbersBuffer70/b048_human_footprint_2009_1km_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J48trk))
#crs(hf)   
#compareCRS(hf,get_crs(J48trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J48/UTM_b048_human_footprint_2009_1km_tif_exp"), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J48/UTM_b048_human_footprint_2009_1km_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J49
#' b049_human_footprint_2009_1km_tif_exp
hf <-"D:/GISbases/b049JaguarAlbersBuffer70/b049_human_footprint_2009_1km_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J49trk))
#crs(hf)   
#compareCRS(hf,get_crs(J49trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J49/UTM_b049_human_footprint_2009_1km_tif_exp"), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J49/UTM_b049_human_footprint_2009_1km_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J50
#' b050_human_footprint_2009_1km_tif_exp
hf <-"D:/GISbases/b050JaguarAlbersBuffer70/b050_human_footprint_2009_1km_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J50trk))
#crs(hf)   
#compareCRS(hf,get_crs(J50trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J50/UTM_b050_human_footprint_2009_1km_tif_exp"), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J50/UTM_b050_human_footprint_2009_1km_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J51
#' b051_human_footprint_2009_1km_tif_exp
hf <-"D:/GISbases/b051JaguarAlbersBuffer70/b051_human_footprint_2009_1km_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J51trk))
#crs(hf)   
#compareCRS(hf,get_crs(J51trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J51/UTM_b051_human_footprint_2009_1km_tif_exp"), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J51/UTM_b051_human_footprint_2009_1km_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J52
#' b052_human_footprint_2009_1km_tif_exp
hf <-"D:/GISbases/b052JaguarAlbersBuffer70/b052_human_footprint_2009_1km_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J52trk))
#crs(hf)   
#compareCRS(hf,get_crs(J52trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J52/UTM_b052_human_footprint_2009_1km_tif_exp"), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J52/UTM_b052_human_footprint_2009_1km_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J53
#' b053_human_footprint_2009_1km_tif_exp
hf <-"D:/GISbases/b053JaguarAlbersBuffer70/b053_human_footprint_2009_1km_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J53trk))
#crs(hf)   
#compareCRS(hf,get_crs(J53trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J53/UTM_b053_human_footprint_2009_1km_tif_exp"), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J53/UTM_b053_human_footprint_2009_1km_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J54
#' b054_human_footprint_2009_1km_tif_exp
hf <-"D:/GISbases/b054JaguarAlbersBuffer70/b054_human_footprint_2009_1km_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J54trk))
#crs(hf)   
#compareCRS(hf,get_crs(J54trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J54/UTM_b054_human_footprint_2009_1km_tif_exp"), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J54/UTM_b054_human_footprint_2009_1km_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J55
#' b055_human_footprint_2009_1km_tif_exp
hf <-"D:/GISbases/b055JaguarAlbersBuffer70/b055_human_footprint_2009_1km_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J55trk))
#crs(hf)   
#compareCRS(hf,get_crs(J55trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J55/UTM_b055_human_footprint_2009_1km_tif_exp"), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J55/UTM_b055_human_footprint_2009_1km_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J56
#' b056_human_footprint_2009_1km_tif_exp
hf <-"D:/GISbases/b056JaguarAlbersBuffer70/b056_human_footprint_2009_1km_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J56trk))
#crs(hf)   
#compareCRS(hf,get_crs(J56trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J56/UTM_b056_human_footprint_2009_1km_tif_exp"), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J56/UTM_b056_human_footprint_2009_1km_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J57
#' b057_human_footprint_2009_1km_tif_exp
hf <-"D:/GISbases/b057JaguarAlbersBuffer70/b057_human_footprint_2009_1km_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J57trk))
#crs(hf)   
#compareCRS(hf,get_crs(J57trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J57/UTM_b057_human_footprint_2009_1km_tif_exp"), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J57/UTM_b057_human_footprint_2009_1km_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J58
#' b058_human_footprint_2009_1km_tif_exp
hf <-"D:/GISbases/b058JaguarAlbersBuffer70/b058_human_footprint_2009_1km_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J58trk))
#crs(hf)   
#compareCRS(hf,get_crs(J58trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J58/UTM_b058_human_footprint_2009_1km_tif_exp"), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J58/UTM_b058_human_footprint_2009_1km_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J59
#' b059_human_footprint_2009_1km_tif_exp
hf <-"D:/GISbases/b059JaguarAlbersBuffer70/b059_human_footprint_2009_1km_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J59trk))
#crs(hf)   
#compareCRS(hf,get_crs(J59trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J59/UTM_b059_human_footprint_2009_1km_tif_exp"), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J59/UTM_b059_human_footprint_2009_1km_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J60
#' b060_human_footprint_2009_1km_tif_exp
hf <-"D:/GISbases/b060JaguarAlbersBuffer70/b060_human_footprint_2009_1km_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J60trk))
#crs(hf)   
#compareCRS(hf,get_crs(J60trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J60/UTM_b060_human_footprint_2009_1km_tif_exp"), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J60/UTM_b060_human_footprint_2009_1km_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J61
#' b061_human_footprint_2009_1km_tif_exp
hf <-"D:/GISbases/b061JaguarAlbersBuffer70/b061_human_footprint_2009_1km_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J61trk))
#crs(hf)   
#compareCRS(hf,get_crs(J61trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J61/UTM_b061_human_footprint_2009_1km_tif_exp"), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J61/UTM_b061_human_footprint_2009_1km_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J62
#' b062_human_footprint_2009_1km_tif_exp
hf <-"D:/GISbases/b062JaguarAlbersBuffer70/b062_human_footprint_2009_1km_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J62trk))
#crs(hf)   
#compareCRS(hf,get_crs(J62trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J62/UTM_b062_human_footprint_2009_1km_tif_exp"), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J62/UTM_b062_human_footprint_2009_1km_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J63
#' b063_human_footprint_2009_1km_tif_exp
hf <-"D:/GISbases/b063JaguarAlbersBuffer70/b063_human_footprint_2009_1km_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J63trk))
#crs(hf)   
#compareCRS(hf,get_crs(J63trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J63/UTM_b063_human_footprint_2009_1km_tif_exp"), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J63/UTM_b063_human_footprint_2009_1km_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J64
#' b064_human_footprint_2009_1km_tif_exp
hf <-"D:/GISbases/b064JaguarAlbersBuffer70/b064_human_footprint_2009_1km_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#X11();plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J64trk))
#crs(hf)   
#compareCRS(hf,get_crs(J64trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J64/UTM_b064_human_footprint_2009_1km_tif_exp"), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J64/UTM_b064_human_footprint_2009_1km_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J65
#' b065_human_footprint_2009_1km_tif_exp
hf <-"D:/GISbases/b065JaguarAlbersBuffer70/b065_human_footprint_2009_1km_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J65trk))
#crs(hf)   
#compareCRS(hf,get_crs(J65trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J65/UTM_b065_human_footprint_2009_1km_tif_exp"), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J65/UTM_b065_human_footprint_2009_1km_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J66
#' b066_human_footprint_2009_1km_tif_exp
hf <-"D:/GISbases/b066JaguarAlbersBuffer70/b066_human_footprint_2009_1km_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J66trk))
#crs(hf)   
#compareCRS(hf,get_crs(J66trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J66/UTM_b066_human_footprint_2009_1km_tif_exp"), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J66/UTM_b066_human_footprint_2009_1km_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J67
#' b067_human_footprint_2009_1km_tif_exp
hf <-"D:/GISbases/b067JaguarAlbersBuffer70/b067_human_footprint_2009_1km_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J67trk))
#crs(hf)   
#compareCRS(hf,get_crs(J67trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J67/UTM_b067_human_footprint_2009_1km_tif_exp"), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J67/UTM_b067_human_footprint_2009_1km_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J68
#' b068_human_footprint_2009_1km_tif_exp
hf <-"D:/GISbases/b068JaguarAlbersBuffer70/b068_human_footprint_2009_1km_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J68trk))
#crs(hf)   
#compareCRS(hf,get_crs(J68trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J68/UTM_b068_human_footprint_2009_1km_tif_exp"), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J68/UTM_b068_human_footprint_2009_1km_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J69
#' b069_human_footprint_2009_1km_tif_exp
hf <-"D:/GISbases/b069JaguarAlbersBuffer70/b069_human_footprint_2009_1km_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J69trk))
#crs(hf)   
#compareCRS(hf,get_crs(J69trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J69/UTM_b069_human_footprint_2009_1km_tif_exp"), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J69/UTM_b069_human_footprint_2009_1km_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J70
#' b070_human_footprint_2009_1km_tif_exp
hf <-"D:/GISbases/b070JaguarAlbersBuffer70/b070_human_footprint_2009_1km_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J70trk))
#crs(hf)   
#compareCRS(hf,get_crs(J70trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J70/UTM_b070_human_footprint_2009_1km_tif_exp"), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J70/UTM_b070_human_footprint_2009_1km_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J71
#' b071_human_footprint_2009_1km_tif_exp
hf <-"D:/GISbases/b071JaguarAlbersBuffer70/b071_human_footprint_2009_1km_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J71trk))
#crs(hf)   
#compareCRS(hf,get_crs(J71trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J71/UTM_b071_human_footprint_2009_1km_tif_exp"), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J71/UTM_b071_human_footprint_2009_1km_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J72
#' b072_human_footprint_2009_1km_tif_exp
hf <-"D:/GISbases/b072JaguarAlbersBuffer70/b072_human_footprint_2009_1km_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J72trk))
#crs(hf)   
#compareCRS(hf,get_crs(J72trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J72/UTM_b072_human_footprint_2009_1km_tif_exp"), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J72/UTM_b072_human_footprint_2009_1km_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J73
#' b073_human_footprint_2009_1km_tif_exp
hf <-"D:/GISbases/b073JaguarAlbersBuffer70/b073_human_footprint_2009_1km_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J73trk))
#crs(hf)   
#compareCRS(hf,get_crs(J73trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J73/UTM_b073_human_footprint_2009_1km_tif_exp"), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J73/UTM_b073_human_footprint_2009_1km_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J74
#' b074_human_footprint_2009_1km_tif_exp
hf <-"D:/GISbases/b074JaguarAlbersBuffer70/b074_human_footprint_2009_1km_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J74trk))
#crs(hf)   
#compareCRS(hf,get_crs(J74trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J74/UTM_b074_human_footprint_2009_1km_tif_exp"), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J74/UTM_b074_human_footprint_2009_1km_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J75
#' b075_human_footprint_2009_1km_tif_exp
hf <-"D:/GISbases/b075JaguarAlbersBuffer70/b075_human_footprint_2009_1km_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J75trk))
#crs(hf)   
#compareCRS(hf,get_crs(J75trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J75/UTM_b075_human_footprint_2009_1km_tif_exp"), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J75/UTM_b075_human_footprint_2009_1km_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J76
#' b076_human_footprint_2009_1km_tif_exp
hf <-"D:/GISbases/b076JaguarAlbersBuffer70/b076_human_footprint_2009_1km_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J76trk))
#crs(hf)   
#compareCRS(hf,get_crs(J76trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J76/UTM_b076_human_footprint_2009_1km_tif_exp"), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J76/UTM_b076_human_footprint_2009_1km_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J77
#' b077_human_footprint_2009_1km_tif_exp
hf <-"D:/GISbases/b077JaguarAlbersBuffer70/b077_human_footprint_2009_1km_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J77trk))
#crs(hf)   
#compareCRS(hf,get_crs(J77trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J77/UTM_b077_human_footprint_2009_1km_tif_exp"), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J77/UTM_b077_human_footprint_2009_1km_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J78
#' b078_human_footprint_2009_1km_tif_exp
hf <-"D:/GISbases/b078JaguarAlbersBuffer70/b078_human_footprint_2009_1km_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J78trk))
#crs(hf)   
#compareCRS(hf,get_crs(J78trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J78/UTM_b078_human_footprint_2009_1km_tif_exp"), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J78/UTM_b078_human_footprint_2009_1km_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J79
#' b079_human_footprint_2009_1km_tif_exp
hf <-"D:/GISbases/b079JaguarAlbersBuffer70/b079_human_footprint_2009_1km_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J79trk))
#crs(hf)   
#compareCRS(hf,get_crs(J79trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J79/UTM_b079_human_footprint_2009_1km_tif_exp"), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J79/UTM_b079_human_footprint_2009_1km_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J80
#' b080_human_footprint_2009_1km_tif_exp
hf <-"D:/GISbases/b080JaguarAlbersBuffer70/b080_human_footprint_2009_1km_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J80trk))
#crs(hf)   
#compareCRS(hf,get_crs(J80trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J80/UTM_b080_human_footprint_2009_1km_tif_exp"), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J80/UTM_b080_human_footprint_2009_1km_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J81
#' b081_human_footprint_2009_1km_tif_exp
hf <-"D:/GISbases/b081JaguarAlbersBuffer70/b081_human_footprint_2009_1km_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J81trk))
#crs(hf)   
#compareCRS(hf,get_crs(J81trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J81/UTM_b081_human_footprint_2009_1km_tif_exp"), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J81/UTM_b081_human_footprint_2009_1km_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J82
#' b082_human_footprint_2009_1km_tif_exp
hf <-"D:/GISbases/b082JaguarAlbersBuffer70/b082_human_footprint_2009_1km_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J82trk))
#crs(hf)   
#compareCRS(hf,get_crs(J82trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J82/UTM_b082_human_footprint_2009_1km_tif_exp"), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J82/UTM_b082_human_footprint_2009_1km_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J83
#' b083_human_footprint_2009_1km_tif_exp
hf <-"D:/GISbases/b083JaguarAlbersBuffer70/b083_human_footprint_2009_1km_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J83trk))
#crs(hf)   
#compareCRS(hf,get_crs(J83trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J83/UTM_b083_human_footprint_2009_1km_tif_exp"), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J83/UTM_b083_human_footprint_2009_1km_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J84
#' b084_human_footprint_2009_1km_tif_exp
hf <-"D:/GISbases/b084JaguarAlbersBuffer70/b084_human_footprint_2009_1km_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J84trk))
#crs(hf)   
#compareCRS(hf,get_crs(J84trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J84/UTM_b084_human_footprint_2009_1km_tif_exp"), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J84/UTM_b084_human_footprint_2009_1km_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J85
#' b085_human_footprint_2009_1km_tif_exp
hf <-"D:/GISbases/b085JaguarAlbersBuffer70/b085_human_footprint_2009_1km_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J85trk))
#crs(hf)   
#compareCRS(hf,get_crs(J85trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J85/UTM_b085_human_footprint_2009_1km_tif_exp"), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J85/UTM_b085_human_footprint_2009_1km_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J86
#' b086_human_footprint_2009_1km_tif_exp
hf <-"D:/GISbases/b086JaguarAlbersBuffer70/b086_human_footprint_2009_1km_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J86trk))
#crs(hf)   
#compareCRS(hf,get_crs(J86trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J86/UTM_b086_human_footprint_2009_1km_tif_exp"), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J86/UTM_b086_human_footprint_2009_1km_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J87
#' b087_human_footprint_2009_1km_tif_exp
hf <-"D:/GISbases/b087JaguarAlbersBuffer70/b087_human_footprint_2009_1km_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J87trk))
#crs(hf)   
#compareCRS(hf,get_crs(J87trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J87/UTM_b087_human_footprint_2009_1km_tif_exp"), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J87/UTM_b087_human_footprint_2009_1km_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J88
#' b088_human_footprint_2009_1km_tif_exp
hf <-"D:/GISbases/b088JaguarAlbersBuffer70/b088_human_footprint_2009_1km_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J88trk))
#crs(hf)   
#compareCRS(hf,get_crs(J88trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J88/UTM_b088_human_footprint_2009_1km_tif_exp"), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J88/UTM_b088_human_footprint_2009_1km_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J89
#' b089_human_footprint_2009_1km_tif_exp
hf <-"D:/GISbases/b089JaguarAlbersBuffer70/b089_human_footprint_2009_1km_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J89trk))
#crs(hf)   
#compareCRS(hf,get_crs(J89trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J89/UTM_b089_human_footprint_2009_1km_tif_exp"), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J89/UTM_b089_human_footprint_2009_1km_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J90
#' b090_human_footprint_2009_1km_tif_exp
hf <-"D:/GISbases/b090JaguarAlbersBuffer70/b090_human_footprint_2009_1km_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J90trk))
#crs(hf)   
#compareCRS(hf,get_crs(J90trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J90/UTM_b090_human_footprint_2009_1km_tif_exp"), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J90/UTM_b090_human_footprint_2009_1km_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J91
#' b091_human_footprint_2009_1km_tif_exp
hf <-"D:/GISbases/b091JaguarAlbersBuffer70/b091_human_footprint_2009_1km_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J91trk))
#crs(hf)   
#compareCRS(hf,get_crs(J91trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J91/UTM_b091_human_footprint_2009_1km_tif_exp"), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J91/UTM_b091_human_footprint_2009_1km_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J92
#' b092_human_footprint_2009_1km_tif_exp
hf <-"D:/GISbases/b092JaguarAlbersBuffer70/b092_human_footprint_2009_1km_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J92trk))
#crs(hf)   
#compareCRS(hf,get_crs(J92trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J92/UTM_b092_human_footprint_2009_1km_tif_exp"), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J92/UTM_b092_human_footprint_2009_1km_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J93
#' b093_human_footprint_2009_1km_tif_exp
hf <-"D:/GISbases/b093JaguarAlbersBuffer70/b093_human_footprint_2009_1km_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J93trk))
#crs(hf)   
#compareCRS(hf,get_crs(J93trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J93/UTM_b093_human_footprint_2009_1km_tif_exp"), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J93/UTM_b093_human_footprint_2009_1km_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J94
#' b094_human_footprint_2009_1km_tif_exp
hf <-"D:/GISbases/b094JaguarAlbersBuffer70/b094_human_footprint_2009_1km_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J94trk))
#crs(hf)   
#compareCRS(hf,get_crs(J94trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J94/UTM_b094_human_footprint_2009_1km_tif_exp"), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J94/UTM_b094_human_footprint_2009_1km_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J95
#' b095_human_footprint_2009_1km_tif_exp
hf <-"D:/GISbases/b095JaguarAlbersBuffer70/b095_human_footprint_2009_1km_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J95trk))
#crs(hf)   
#compareCRS(hf,get_crs(J95trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J95/UTM_b095_human_footprint_2009_1km_tif_exp"), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J95/UTM_b095_human_footprint_2009_1km_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J96
#' b096_human_footprint_2009_1km_tif_exp
hf <-"D:/GISbases/b096JaguarAlbersBuffer70/b096_human_footprint_2009_1km_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J96trk))
#crs(hf)   
#compareCRS(hf,get_crs(J96trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J96/UTM_b096_human_footprint_2009_1km_tif_exp"), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J96/UTM_b096_human_footprint_2009_1km_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J97
#' b097_human_footprint_2009_1km_tif_exp
hf <-"D:/GISbases/b097JaguarAlbersBuffer70/b097_human_footprint_2009_1km_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J97trk))
#crs(hf)   
#compareCRS(hf,get_crs(J97trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J97/UTM_b097_human_footprint_2009_1km_tif_exp"), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J97/UTM_b097_human_footprint_2009_1km_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J98
#' b098_human_footprint_2009_1km_tif_exp
hf <-"D:/GISbases/b098JaguarAlbersBuffer70/b098_human_footprint_2009_1km_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J98trk))
#crs(hf)   
#compareCRS(hf,get_crs(J98trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J98/UTM_b098_human_footprint_2009_1km_tif_exp"), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J98/UTM_b098_human_footprint_2009_1km_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J99
#' b099_human_footprint_2009_1km_tif_exp
hf <-"D:/GISbases/b099JaguarAlbersBuffer70/b099_human_footprint_2009_1km_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J99trk))
#crs(hf)   
#compareCRS(hf,get_crs(J99trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J99/UTM_b099_human_footprint_2009_1km_tif_exp"), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J99/UTM_b099_human_footprint_2009_1km_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J100
#' b100_human_footprint_2009_1km_tif_exp
hf <-"D:/GISbases/b100JaguarAlbersBuffer70/b100_human_footprint_2009_1km_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J100trk))
#crs(hf)   
#compareCRS(hf,get_crs(J100trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J100/UTM_b100_human_footprint_2009_1km_tif_exp"), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J100/UTM_b010_human_footprint_2009_1km_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J101
#' b101_human_footprint_2009_1km_tif_exp
hf <-"D:/GISbases/b101JaguarAlbersBuffer70/b101_human_footprint_2009_1km_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J101trk))
#crs(hf)   
#compareCRS(hf,get_crs(J101trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J101/UTM_b101_human_footprint_2009_1km_tif_exp"), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J101/UTM_b101_human_footprint_2009_1km_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J102
#' b102_human_footprint_2009_1km_tif_exp
hf <-"D:/GISbases/b102JaguarAlbersBuffer70/b102_human_footprint_2009_1km_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J102trk))
#crs(hf)   
#compareCRS(hf,get_crs(J102trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J102/UTM_b102_human_footprint_2009_1km_tif_exp"), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J102/UTM_b102_human_footprint_2009_1km_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J103
#' b103_human_footprint_2009_1km_tif_exp
hf <-"D:/GISbases/b103JaguarAlbersBuffer70/b103_human_footprint_2009_1km_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J103trk))
#crs(hf)   
#compareCRS(hf,get_crs(J103trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J103/UTM_b103_human_footprint_2009_1km_tif_exp"), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J103/UTM_b103_human_footprint_2009_1km_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J104
#' b104_human_footprint_2009_1km_tif_exp
hf <-"D:/GISbases/b104JaguarAlbersBuffer70/b104_human_footprint_2009_1km_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J104trk))
#crs(hf)   
#compareCRS(hf,get_crs(J104trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J104/UTM_b104_human_footprint_2009_1km_tif_exp"), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J104/UTM_b104_human_footprint_2009_1km_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J105
#' b105_human_footprint_2009_1km_tif_exp
hf <-"D:/GISbases/b105JaguarAlbersBuffer70/b105_human_footprint_2009_1km_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J105trk))
#crs(hf)   
#compareCRS(hf,get_crs(J105trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J105/UTM_b105_human_footprint_2009_1km_tif_exp"), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J105/UTM_b105_human_footprint_2009_1km_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J106
#' b106_human_footprint_2009_1km_tif_exp
hf <-"D:/GISbases/b106JaguarAlbersBuffer70/b106_human_footprint_2009_1km_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J106trk))
#crs(hf)   
#compareCRS(hf,get_crs(J106trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J106/UTM_b106_human_footprint_2009_1km_tif_exp"), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J106/UTM_b106_human_footprint_2009_1km_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J107
#' b107_human_footprint_2009_1km_tif_exp
hf <-"D:/GISbases/b107JaguarAlbersBuffer70/b107_human_footprint_2009_1km_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J107trk))
#crs(hf)   
#compareCRS(hf,get_crs(J107trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J107/UTM_b107_human_footprint_2009_1km_tif_exp"), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J107/UTM_b107_human_footprint_2009_1km_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J108
#' b108_human_footprint_2009_1km_tif_exp
hf <-"D:/GISbases/b108JaguarAlbersBuffer70/b108_human_footprint_2009_1km_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J108trk))
#crs(hf)   
#compareCRS(hf,get_crs(J108trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J108/UTM_b108_human_footprint_2009_1km_tif_exp"), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J108/UTM_b108_human_footprint_2009_1km_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J109
#' b109_human_footprint_2009_1km_tif_exp
hf <-"D:/GISbases/b109JaguarAlbersBuffer70/b109_human_footprint_2009_1km_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J109trk))
#crs(hf)   
#compareCRS(hf,get_crs(J109trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J109/UTM_b109_human_footprint_2009_1km_tif_exp"), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J109/UTM_b109_human_footprint_2009_1km_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J110
#' b110_human_footprint_2009_1km_tif_exp
hf <-"D:/GISbases/b110JaguarAlbersBuffer70/b110_human_footprint_2009_1km_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J110trk))
#crs(hf)   
#compareCRS(hf,get_crs(J110trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J110/UTM_b110_human_footprint_2009_1km_tif_exp"), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J110/UTM_b110_human_footprint_2009_1km_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J111
#' b111_human_footprint_2009_1km_tif_exp
hf <-"D:/GISbases/b111JaguarAlbersBuffer70/b111_human_footprint_2009_1km_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J111trk))
#crs(hf)   
#compareCRS(hf,get_crs(J111trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J111/UTM_b111_human_footprint_2009_1km_tif_exp"), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J111/UTM_b111_human_footprint_2009_1km_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J112
#' b112_human_footprint_2009_1km_tif_exp
hf <-"D:/GISbases/b112JaguarAlbersBuffer70/b112_human_footprint_2009_1km_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J112trk))
#crs(hf)   
#compareCRS(hf,get_crs(J112trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J112/UTM_b112_human_footprint_2009_1km_tif_exp"), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J112/UTM_b112_human_footprint_2009_1km_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J113
#' b113_human_footprint_2009_1km_tif_exp
hf <-"D:/GISbases/b113JaguarAlbersBuffer70/b113_human_footprint_2009_1km_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J113trk))
#crs(hf)   
#compareCRS(hf,get_crs(J113trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J113/UTM_b113_human_footprint_2009_1km_tif_exp"), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J113/UTM_b113_human_footprint_2009_1km_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J114
#' b114_human_footprint_2009_1km_tif_exp
hf <-"D:/GISbases/b114JaguarAlbersBuffer70/b114_human_footprint_2009_1km_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J114trk))
#crs(hf)   
#compareCRS(hf,get_crs(J114trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J114/UTM_b114_human_footprint_2009_1km_tif_exp"), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J114/UTM_b114_human_footprint_2009_1km_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J115
#' b115_human_footprint_2009_1km_tif_exp
hf <-"D:/GISbases/b115JaguarAlbersBuffer70/b115_human_footprint_2009_1km_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J115trk))
#crs(hf)   
#compareCRS(hf,get_crs(J115trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J115/UTM_b115_human_footprint_2009_1km_tif_exp"), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J115/UTM_b115_human_footprint_2009_1km_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J116
#' b116_human_footprint_2009_1km_tif_exp
hf <-"D:/GISbases/b116JaguarAlbersBuffer70/b116_human_footprint_2009_1km_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J116trk))
#crs(hf)   
#compareCRS(hf,get_crs(J116trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J116/UTM_b116_human_footprint_2009_1km_tif_exp"), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J116/UTM_b116_human_footprint_2009_1km_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J117
#' b117_human_footprint_2009_1km_tif_exp
hf <-"D:/GISbases/b117JaguarAlbersBuffer70/b117_human_footprint_2009_1km_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J117trk))
#crs(hf)   
#compareCRS(hf,get_crs(J117trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J117/UTM_b117_human_footprint_2009_1km_tif_exp"), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J117/UTM_b117_human_footprint_2009_1km_tif_exp.tif";(test=raster(test)); x11(); plot(test)

############################################################################################################################
#b115_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif

#' J1
#' Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif
hf <-"D:/GISbases/b001JaguarAlbersBuffer70/b001_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J1trk))
#crs(hf)   
#compareCRS(hf,get_crs(J1trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J1/UTM_b001_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J1/UTM_b001_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J2
#' b002_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.
hf <-"D:/GISbases/b002JaguarAlbersBuffer70/b002_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J2trk))
#crs(hf)   
#compareCRS(hf,get_crs(J2trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J2/UTM_b002_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp"), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J2/UTM_b002_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J3
#' b003_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.
hf <-"D:/GISbases/b003JaguarAlbersBuffer70/b003_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J3trk))
#crs(hf)   
#compareCRS(hf,get_crs(J3trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J3/UTM_b003_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J3/UTM_b003_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J4
#' b004_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.
hf <-"D:/GISbases/b004JaguarAlbersBuffer70/b004_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J4trk))
#crs(hf)   
#compareCRS(hf,get_crs(J4trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J4/UTM_b004_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J4/UTM_b004_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J5
#' b005_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.
hf <-"D:/GISbases/b005JaguarAlbersBuffer70/b005_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J5trk))
#crs(hf)   
#compareCRS(hf,get_crs(J5trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J5/UTM_b005_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J5/UTM_b005_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J6
#' b006_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.
hf <-"D:/GISbases/b006JaguarAlbersBuffer70/b006_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J6trk))
#crs(hf)   
#compareCRS(hf,get_crs(J6trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J6/UTM_b006_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J6/UTM_b006_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J7
#' b007_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.
hf <-"D:/GISbases/b007JaguarAlbersBuffer70/b007_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J7trk))
#crs(hf)   
#compareCRS(hf,get_crs(J7trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J7/UTM_b007_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J7/UTM_b007_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J8
#' b008_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.
hf <-"D:/GISbases/b008JaguarAlbersBuffer70/b008_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J8trk))
#crs(hf)   
#compareCRS(hf,get_crs(J8trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J8/UTM_b008_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J8/UTM_b008_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J9
#' b009_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.
hf <-"D:/GISbases/b009JaguarAlbersBuffer70/b009_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J9trk))
#crs(hf)   
#compareCRS(hf,get_crs(J9trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J9/UTM_b009_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J9/UTM_b009_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J10
#' b010_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.
hf <-"D:/GISbases/b010JaguarAlbersBuffer70/b010_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J10trk))
#crs(hf)   
#compareCRS(hf,get_crs(J10trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J10/UTM_b010_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J10/UTM_b010_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J11
#' b011_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.
hf <-"D:/GISbases/b011JaguarAlbersBuffer70/b011_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J11trk))
#crs(hf)   
#compareCRS(hf,get_crs(J11trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J11/UTM_b011_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J11/UTM_b011_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J12
#' b012_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.
hf <-"D:/GISbases/b012JaguarAlbersBuffer70/b012_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J12trk))
#crs(hf)   
#compareCRS(hf,get_crs(J12trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J12/UTM_b012_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J12/UTM_b012_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J13
#' b013_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.
hf <-"D:/GISbases/b013JaguarAlbersBuffer70/b013_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J13trk))
#crs(hf)   
#compareCRS(hf,get_crs(J13trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J13/UTM_b013_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J13/UTM_b013_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J14
#' b014_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.
hf <-"D:/GISbases/b014JaguarAlbersBuffer70/b014_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J14trk))
#crs(hf)   
#compareCRS(hf,get_crs(J14trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J14/UTM_b014_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J14/UTM_b014_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J15
#' b015_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.
hf <-"D:/GISbases/b015JaguarAlbersBuffer70/b015_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J15trk))
#crs(hf)   
#compareCRS(hf,get_crs(J15trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J15/UTM_b015_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J15/UTM_b015_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J16
#' b016_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.
hf <-"D:/GISbases/b016JaguarAlbersBuffer70/b016_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J16trk))
#crs(hf)   
#compareCRS(hf,get_crs(J16trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J16/UTM_b016_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J16/UTM_b016_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J17
#' b017_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.
hf <-"D:/GISbases/b017JaguarAlbersBuffer70/b017_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J17trk))
#crs(hf)   
#compareCRS(hf,get_crs(J17trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J17/UTM_b017_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J17/UTM_b017_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J18
#' b018_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.
hf <-"D:/GISbases/b018JaguarAlbersBuffer70/b018_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J18trk))
#crs(hf)   
#compareCRS(hf,get_crs(J18trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J18/UTM_b018_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J18/UTM_b018_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J19
#' b019_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.
hf <-"D:/GISbases/b019JaguarAlbersBuffer70/b019_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J19trk))
#crs(hf)   
#compareCRS(hf,get_crs(J19trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J19/UTM_b019_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J19/UTM_b019_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J20
#' b020_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.
hf <-"D:/GISbases/b020JaguarAlbersBuffer70/b020_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J20trk))
#crs(hf)   
#compareCRS(hf,get_crs(J20trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J20/UTM_b020_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J20/UTM_b020_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J21
#' b021_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.
hf <-"D:/GISbases/b021JaguarAlbersBuffer70/b021_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J21trk))
#crs(hf)   
#compareCRS(hf,get_crs(J21trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J21/UTM_b021_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J21/UTM_b021_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J22
#' b022_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.
hf <-"D:/GISbases/b022JaguarAlbersBuffer70/b022_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J22trk))
#crs(hf)   
#compareCRS(hf,get_crs(J22trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J22/UTM_b022_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J22/UTM_b022_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J23
#' b023_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.
hf <-"D:/GISbases/b023JaguarAlbersBuffer70/b023_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J23trk))
#crs(hf)   
#compareCRS(hf,get_crs(J23trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J23/UTM_b023_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J23/UTM_b023_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J24
#' b024_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.
hf <-"D:/GISbases/b024JaguarAlbersBuffer70/b024_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J24trk))
#crs(hf)   
#compareCRS(hf,get_crs(J24trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J24/UTM_b024_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J24/UTM_b024_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J25
#' b025_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.
hf <-"D:/GISbases/b025JaguarAlbersBuffer70/b025_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J25trk))
#crs(hf)   
#compareCRS(hf,get_crs(J25trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J25/UTM_b025_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J25/UTM_b025_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J26
#' b026_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.
hf <-"D:/GISbases/b026JaguarAlbersBuffer70/b026_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J26trk))
#crs(hf)   
#compareCRS(hf,get_crs(J26trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J26/UTM_b026_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J26/UTM_b026_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J27
#' b027_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.
hf <-"D:/GISbases/b027JaguarAlbersBuffer70/b027_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J27trk))
#crs(hf)   
#compareCRS(hf,get_crs(J27trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J27/UTM_b027_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J27/UTM_b027_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J28
#' b028_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.
hf <-"D:/GISbases/b028JaguarAlbersBuffer70/b028_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J28trk))
#crs(hf)   
#compareCRS(hf,get_crs(J28trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J28/UTM_b028_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J28/UTM_b028_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J29
#' b029_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.
hf <-"D:/GISbases/b029JaguarAlbersBuffer70/b029_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J29trk))
#crs(hf)   
#compareCRS(hf,get_crs(J29trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J29/UTM_b029_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J29/UTM_b029_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J30
#' b030_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.
hf <-"D:/GISbases/b030JaguarAlbersBuffer70/b030_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J30trk))
#crs(hf)   
#compareCRS(hf,get_crs(J30trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J30/UTM_b030_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J30/UTM_b030_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J30
#' b030_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.
hf <-"D:/GISbases/b030JaguarAlbersBuffer70/b030_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J30trk))
#crs(hf)   
#compareCRS(hf,get_crs(J30trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J30/UTM_b030_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J30/UTM_b030_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J31
#' b031_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.
hf <-"D:/GISbases/b031JaguarAlbersBuffer70/b031_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J31trk))
#crs(hf)   
#compareCRS(hf,get_crs(J31trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J31/UTM_b031_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J31/UTM_b031_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J32
#' b032_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.
hf <-"D:/GISbases/b032JaguarAlbersBuffer70/b032_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J32trk))
#crs(hf)   
#compareCRS(hf,get_crs(J32trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J32/UTM_b032_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J32/UTM_b032_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J33
#' b033_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.
hf <-"D:/GISbases/b033JaguarAlbersBuffer70/b033_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J33trk))
#crs(hf)   
#compareCRS(hf,get_crs(J33trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J33/UTM_b033_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J33/UTM_b033_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J34
#' b034_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.
hf <-"D:/GISbases/b034JaguarAlbersBuffer70/b034_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J34trk))
#crs(hf)   
#compareCRS(hf,get_crs(J34trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J34/UTM_b034_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J34/UTM_b034_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J35
#' b035_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.
hf <-"D:/GISbases/b035JaguarAlbersBuffer70/b035_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J35trk))
#crs(hf)   
#compareCRS(hf,get_crs(J35trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J35/UTM_b035_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J35/UTM_b035_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J36
#' b036_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.
hf <-"D:/GISbases/b036JaguarAlbersBuffer70/b036_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J36trk))
#crs(hf)   
#compareCRS(hf,get_crs(J36trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J36/UTM_b036_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J36/UTM_b036_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J37
#' b037_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.
hf <-"D:/GISbases/b037JaguarAlbersBuffer70/b037_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J37trk))
#crs(hf)   
#compareCRS(hf,get_crs(J37trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J37/UTM_b037_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J37/UTM_b037_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J38
#' b038_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.
hf <-"D:/GISbases/b038JaguarAlbersBuffer70/b038_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J38trk))
#crs(hf)   
#compareCRS(hf,get_crs(J38trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J38/UTM_b038_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J38/UTM_b038_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J39
#' b039_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.
hf <-"D:/GISbases/b039JaguarAlbersBuffer70/b039_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J39trk))
#crs(hf)   
#compareCRS(hf,get_crs(J39trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J39/UTM_b039_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J39/UTM_b039_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J40
#' b040_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.
hf <-"D:/GISbases/b040JaguarAlbersBuffer70/b040_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J40trk))
#crs(hf)   
#compareCRS(hf,get_crs(J40trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J40/UTM_b040_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J40/UTM_b040_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J41
#' b041_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.
hf <-"D:/GISbases/b041JaguarAlbersBuffer70/b041_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J41trk))
#crs(hf)   
#compareCRS(hf,get_crs(J41trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J41/UTM_b041_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J41/UTM_b041_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J42
#' b042_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.
hf <-"D:/GISbases/b042JaguarAlbersBuffer70/b042_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J42trk))
#crs(hf)   
#compareCRS(hf,get_crs(J42trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J42/UTM_b042_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J42/UTM_b042_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J43
#' b043_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.
hf <-"D:/GISbases/b043JaguarAlbersBuffer70/b043_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J43trk))
#crs(hf)   
#compareCRS(hf,get_crs(J43trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J43/UTM_b043_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J43/UTM_b043_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J44
#' b044_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.
hf <-"D:/GISbases/b044JaguarAlbersBuffer70/b044_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J44trk))
#crs(hf)   
#compareCRS(hf,get_crs(J44trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J44/UTM_b044_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J44/UTM_b044_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J45
#' b045_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.
hf <-"D:/GISbases/b045JaguarAlbersBuffer70/b045_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J45trk))
#crs(hf)   
#compareCRS(hf,get_crs(J45trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J45/UTM_b045_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J45/UTM_b045_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J46
#' b046_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.
hf <-"D:/GISbases/b046JaguarAlbersBuffer70/b046_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J46trk))
#crs(hf)   
#compareCRS(hf,get_crs(J46trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J46/UTM_b046_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J46/UTM_b046_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J47
#' b047_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.
hf <-"D:/GISbases/b047JaguarAlbersBuffer70/b047_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J47trk))
#crs(hf)   
#compareCRS(hf,get_crs(J47trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J47/UTM_b047_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J47/UTM_b047_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J48
#' b048_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.
hf <-"D:/GISbases/b048JaguarAlbersBuffer70/b048_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J48trk))
#crs(hf)   
#compareCRS(hf,get_crs(J48trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J48/UTM_b048_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J48/UTM_b048_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J49
#' b049_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.
hf <-"D:/GISbases/b049JaguarAlbersBuffer70/b049_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J49trk))
#crs(hf)   
#compareCRS(hf,get_crs(J49trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J49/UTM_b049_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J49/UTM_b049_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J50
#' b050_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.
hf <-"D:/GISbases/b050JaguarAlbersBuffer70/b050_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J50trk))
#crs(hf)   
#compareCRS(hf,get_crs(J50trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J50/UTM_b050_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J50/UTM_b050_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J51
#' b051_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.
hf <-"D:/GISbases/b051JaguarAlbersBuffer70/b051_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J51trk))
#crs(hf)   
#compareCRS(hf,get_crs(J51trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J51/UTM_b051_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J51/UTM_b051_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J52
#' b052_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.
hf <-"D:/GISbases/b052JaguarAlbersBuffer70/b052_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J52trk))
#crs(hf)   
#compareCRS(hf,get_crs(J52trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J52/UTM_b052_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J52/UTM_b052_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J53
#' b053_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.
hf <-"D:/GISbases/b053JaguarAlbersBuffer70/b053_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J53trk))
#crs(hf)   
#compareCRS(hf,get_crs(J53trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J53/UTM_b053_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J53/UTM_b053_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J54
#' b054_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.
hf <-"D:/GISbases/b054JaguarAlbersBuffer70/b054_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J54trk))
#crs(hf)   
#compareCRS(hf,get_crs(J54trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J54/UTM_b054_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J54/UTM_b054_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J55
#' b055_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.
hf <-"D:/GISbases/b055JaguarAlbersBuffer70/b055_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J55trk))
#crs(hf)   
#compareCRS(hf,get_crs(J55trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J55/UTM_b055_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J55/UTM_b055_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J56
#' b056_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.
hf <-"D:/GISbases/b056JaguarAlbersBuffer70/b056_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J56trk))
#crs(hf)   
#compareCRS(hf,get_crs(J56trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J56/UTM_b056_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J56/UTM_b056_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J57
#' b057_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.
hf <-"D:/GISbases/b057JaguarAlbersBuffer70/b057_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J57trk))
#crs(hf)   
#compareCRS(hf,get_crs(J57trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J57/UTM_b057_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J57/UTM_b057_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J58
#' b058_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.
hf <-"D:/GISbases/b058JaguarAlbersBuffer70/b058_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J58trk))
#crs(hf)   
#compareCRS(hf,get_crs(J58trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J58/UTM_b058_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J58/UTM_b058_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J59
#' b059_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.
hf <-"D:/GISbases/b059JaguarAlbersBuffer70/b059_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J59trk))
#crs(hf)   
#compareCRS(hf,get_crs(J59trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J59/UTM_b059_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J59/UTM_b059_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J60
#' b060_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.
hf <-"D:/GISbases/b060JaguarAlbersBuffer70/b060_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J60trk))
#crs(hf)   
#compareCRS(hf,get_crs(J60trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J60/UTM_b060_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J60/UTM_b060_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J61
#' b061_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.
hf <-"D:/GISbases/b061JaguarAlbersBuffer70/b061_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J61trk))
#crs(hf)   
#compareCRS(hf,get_crs(J61trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J61/UTM_b061_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J61/UTM_b061_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J62
#' b062_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.
hf <-"D:/GISbases/b062JaguarAlbersBuffer70/b062_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J62trk))
#crs(hf)   
#compareCRS(hf,get_crs(J62trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J62/UTM_b062_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J62/UTM_b062_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J63
#' b063_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.
hf <-"D:/GISbases/b063JaguarAlbersBuffer70/b063_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J63trk))
#crs(hf)   
#compareCRS(hf,get_crs(J63trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J63/UTM_b063_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J63/UTM_b063_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J64
#' b064_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.
hf <-"D:/GISbases/b064JaguarAlbersBuffer70/b064_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J64trk))
#crs(hf)   
#compareCRS(hf,get_crs(J64trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J64/UTM_b064_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J64/UTM_b064_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J65
#' b065_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.
hf <-"D:/GISbases/b065JaguarAlbersBuffer70/b065_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J65trk))
#crs(hf)   
#compareCRS(hf,get_crs(J65trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J65/UTM_b065_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J65/UTM_b065_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J66
#' b066_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.
hf <-"D:/GISbases/b066JaguarAlbersBuffer70/b066_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J66trk))
#crs(hf)   
#compareCRS(hf,get_crs(J66trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J66/UTM_b066_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J66/UTM_b066_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J67
#' b067_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.
hf <-"D:/GISbases/b067JaguarAlbersBuffer70/b067_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J67trk))
#crs(hf)   
#compareCRS(hf,get_crs(J67trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J67/UTM_b067_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J67/UTM_b067_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J68
#' b068_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.
hf <-"D:/GISbases/b068JaguarAlbersBuffer70/b068_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J68trk))
#crs(hf)   
#compareCRS(hf,get_crs(J68trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J68/UTM_b068_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J68/UTM_b068_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J69
#' b069_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.
hf <-"D:/GISbases/b069JaguarAlbersBuffer70/b069_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J69trk))
#crs(hf)   
#compareCRS(hf,get_crs(J69trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J69/UTM_b069_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J69/UTM_b069_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J70
#' b070_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.
hf <-"D:/GISbases/b070JaguarAlbersBuffer70/b070_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J70trk))
#crs(hf)   
#compareCRS(hf,get_crs(J70trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J70/UTM_b070_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J70/UTM_b070_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J71
#' b071_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.
hf <-"D:/GISbases/b071JaguarAlbersBuffer70/b071_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J71trk))
#crs(hf)   
#compareCRS(hf,get_crs(J71trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J71/UTM_b071_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J71/UTM_b071_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J72
#' b072_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.
hf <-"D:/GISbases/b072JaguarAlbersBuffer70/b072_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J72trk))
#crs(hf)   
#compareCRS(hf,get_crs(J72trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J72/UTM_b072_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J72/UTM_b072_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J73
#' b073_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.
hf <-"D:/GISbases/b073JaguarAlbersBuffer70/b073_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J73trk))
#crs(hf)   
#compareCRS(hf,get_crs(J73trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J73/UTM_b073_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J73/UTM_b073_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J74
#' b074_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.
hf <-"D:/GISbases/b074JaguarAlbersBuffer70/b074_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J74trk))
#crs(hf)   
#compareCRS(hf,get_crs(J74trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J74/UTM_b074_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J74/UTM_b074_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J75
#' b075_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.
hf <-"D:/GISbases/b075JaguarAlbersBuffer70/b075_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J75trk))
#crs(hf)   
#compareCRS(hf,get_crs(J75trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J75/UTM_b075_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J75/UTM_b075_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J76
#' b076_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.
hf <-"D:/GISbases/b076JaguarAlbersBuffer70/b076_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J76trk))
#crs(hf)   
#compareCRS(hf,get_crs(J76trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J76/UTM_b076_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J76/UTM_b076_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J77
#' b077_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.
hf <-"D:/GISbases/b077JaguarAlbersBuffer70/b077_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J77trk))
#crs(hf)   
#compareCRS(hf,get_crs(J77trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J77/UTM_b077_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J77/UTM_b077_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J78
#' b078_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.
hf <-"D:/GISbases/b078JaguarAlbersBuffer70/b078_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J78trk))
#crs(hf)   
#compareCRS(hf,get_crs(J78trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J78/UTM_b078_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J78/UTM_b078_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J79
#' b079_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.
hf <-"D:/GISbases/b079JaguarAlbersBuffer70/b079_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J79trk))
#crs(hf)   
#compareCRS(hf,get_crs(J79trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J79/UTM_b079_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J79/UTM_b079_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J80
#' b080_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.
hf <-"D:/GISbases/b080JaguarAlbersBuffer70/b080_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J80trk))
#crs(hf)   
#compareCRS(hf,get_crs(J80trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J80/UTM_b080_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J80/UTM_b080_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J81
#' b081_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.
hf <-"D:/GISbases/b081JaguarAlbersBuffer70/b081_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J81trk))
#crs(hf)   
#compareCRS(hf,get_crs(J81trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J81/UTM_b081_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J81/UTM_b081_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J82
#' b082_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.
hf <-"D:/GISbases/b082JaguarAlbersBuffer70/b082_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J82trk))
#crs(hf)   
#compareCRS(hf,get_crs(J82trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J82/UTM_b082_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J82/UTM_b082_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J83
#' b083_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.
hf <-"D:/GISbases/b083JaguarAlbersBuffer70/b083_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J83trk))
#crs(hf)   
#compareCRS(hf,get_crs(J83trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J83/UTM_b083_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J83/UTM_b083_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J84
#' b084_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.
hf <-"D:/GISbases/b084JaguarAlbersBuffer70/b084_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J84trk))
#crs(hf)   
#compareCRS(hf,get_crs(J84trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J84/UTM_b084_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J84/UTM_b084_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J85
#' b085_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.
hf <-"D:/GISbases/b085JaguarAlbersBuffer70/b085_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J85trk))
#crs(hf)   
#compareCRS(hf,get_crs(J85trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J85/UTM_b085_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J85/UTM_b085_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J86
#' b086_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.
hf <-"D:/GISbases/b086JaguarAlbersBuffer70/b086_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J86trk))
#crs(hf)   
#compareCRS(hf,get_crs(J86trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J86/UTM_b086_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J86/UTM_b086_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J87
#' b087_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.
hf <-"D:/GISbases/b087JaguarAlbersBuffer70/b087_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J87trk))
#crs(hf)   
#compareCRS(hf,get_crs(J87trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J87/UTM_b087_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J87/UTM_b087_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J88
#' b088_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.
hf <-"D:/GISbases/b088JaguarAlbersBuffer70/b088_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J88trk))
#crs(hf)   
#compareCRS(hf,get_crs(J88trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J88/UTM_b088_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J88/UTM_b088_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J89
#' b089_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.
hf <-"D:/GISbases/b089JaguarAlbersBuffer70/b089_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J89trk))
#crs(hf)   
#compareCRS(hf,get_crs(J89trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J89/UTM_b089_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J89/UTM_b089_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J90
#' b090_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.
hf <-"D:/GISbases/b090JaguarAlbersBuffer70/b090_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J90trk))
#crs(hf)   
#compareCRS(hf,get_crs(J90trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J90/UTM_b090_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J90/UTM_b090_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J91
#' b091_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.
hf <-"D:/GISbases/b091JaguarAlbersBuffer70/b091_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J91trk))
#crs(hf)   
#compareCRS(hf,get_crs(J91trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J91/UTM_b091_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J91/UTM_b091_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J92
#' b092_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.
hf <-"D:/GISbases/b092JaguarAlbersBuffer70/b092_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J92trk))
#crs(hf)   
#compareCRS(hf,get_crs(J92trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J92/UTM_b092_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J92/UTM_b092_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J93
#' b093_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.
hf <-"D:/GISbases/b093JaguarAlbersBuffer70/b093_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J93trk))
#crs(hf)   
#compareCRS(hf,get_crs(J93trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J93/UTM_b093_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J93/UTM_b093_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J94
#' b094_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.
hf <-"D:/GISbases/b094JaguarAlbersBuffer70/b094_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J94trk))
#crs(hf)   
#compareCRS(hf,get_crs(J94trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J94/UTM_b094_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J94/UTM_b094_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J95
#' b095_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.
hf <-"D:/GISbases/b095JaguarAlbersBuffer70/b095_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J95trk))
#crs(hf)   
#compareCRS(hf,get_crs(J95trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J95/UTM_b095_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J95/UTM_b095_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J96
#' b096_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.
hf <-"D:/GISbases/b096JaguarAlbersBuffer70/b096_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J96trk))
#crs(hf)   
#compareCRS(hf,get_crs(J96trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J96/UTM_b096_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J96/UTM_b096_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J97
#' b097_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.
hf <-"D:/GISbases/b097JaguarAlbersBuffer70/b097_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J97trk))
#crs(hf)   
#compareCRS(hf,get_crs(J97trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J97/UTM_b097_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J97/UTM_b097_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J98
#' b098_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.
hf <-"D:/GISbases/b098JaguarAlbersBuffer70/b098_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J98trk))
#crs(hf)   
#compareCRS(hf,get_crs(J98trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J98/UTM_b098_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J98/UTM_b098_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J99
#' b099_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.
hf <-"D:/GISbases/b099JaguarAlbersBuffer70/b099_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J99trk))
#crs(hf)   
#compareCRS(hf,get_crs(J99trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J99/UTM_b099_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J99/UTM_b099_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J100
#' b100_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.
hf <-"D:/GISbases/b100JaguarAlbersBuffer70/b100_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J100trk))
#crs(hf)   
#compareCRS(hf,get_crs(J100trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J100/UTM_b100_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J100/UTM_b010_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J101
#' b101_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.
hf <-"D:/GISbases/b101JaguarAlbersBuffer70/b101_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J101trk))
#crs(hf)   
#compareCRS(hf,get_crs(J101trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J101/UTM_b101_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J101/UTM_b101_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J102
#' b102_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.
hf <-"D:/GISbases/b102JaguarAlbersBuffer70/b102_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J102trk))
#crs(hf)   
#compareCRS(hf,get_crs(J102trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J102/UTM_b102_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J102/UTM_b102_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J103
#' b103_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.
hf <-"D:/GISbases/b103JaguarAlbersBuffer70/b103_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J103trk))
#crs(hf)   
#compareCRS(hf,get_crs(J103trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J103/UTM_b103_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J103/UTM_b103_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J104
#' b104_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.
hf <-"D:/GISbases/b104JaguarAlbersBuffer70/b104_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J104trk))
#crs(hf)   
#compareCRS(hf,get_crs(J104trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J104/UTM_b104_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J104/UTM_b104_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J105
#' b105_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.
hf <-"D:/GISbases/b105JaguarAlbersBuffer70/b105_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J105trk))
#crs(hf)   
#compareCRS(hf,get_crs(J105trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J105/UTM_b105_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J105/UTM_b105_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J106
#' b106_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.
hf <-"D:/GISbases/b106JaguarAlbersBuffer70/b106_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J106trk))
#crs(hf)   
#compareCRS(hf,get_crs(J106trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J106/UTM_b106_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J106/UTM_b106_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J107
#' b107_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.
hf <-"D:/GISbases/b107JaguarAlbersBuffer70/b107_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J107trk))
#crs(hf)   
#compareCRS(hf,get_crs(J107trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J107/UTM_b107_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J107/UTM_b107_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J108
#' b108_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.
hf <-"D:/GISbases/b108JaguarAlbersBuffer70/b108_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J108trk))
#crs(hf)   
#compareCRS(hf,get_crs(J108trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J108/UTM_b108_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J108/UTM_b108_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J109
#' b109_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.
hf <-"D:/GISbases/b109JaguarAlbersBuffer70/b109_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J109trk))
#crs(hf)   
#compareCRS(hf,get_crs(J109trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J109/UTM_b109_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J109/UTM_b109_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J110
#' b110_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.
hf <-"D:/GISbases/b110JaguarAlbersBuffer70/b110_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J110trk))
#crs(hf)   
#compareCRS(hf,get_crs(J110trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J110/UTM_b110_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J110/UTM_b110_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J111
#' b111_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.
hf <-"D:/GISbases/b111JaguarAlbersBuffer70/b111_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J111trk))
#crs(hf)   
#compareCRS(hf,get_crs(J111trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J111/UTM_b111_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J111/UTM_b111_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J112
#' b112_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.
hf <-"D:/GISbases/b112JaguarAlbersBuffer70/b112_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J112trk))
#crs(hf)   
#compareCRS(hf,get_crs(J112trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J112/UTM_b112_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J112/UTM_b112_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J113
#' b113_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.
hf <-"D:/GISbases/b113JaguarAlbersBuffer70/b113_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J113trk))
#crs(hf)   
#compareCRS(hf,get_crs(J113trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J113/UTM_b113_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J113/UTM_b113_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J114
#' b114_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.
hf <-"D:/GISbases/b114JaguarAlbersBuffer70/b114_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J114trk))
#crs(hf)   
#compareCRS(hf,get_crs(J114trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J114/UTM_b114_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J114/UTM_b114_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J115
#' b115_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.
hf <-"D:/GISbases/b115JaguarAlbersBuffer70/b115_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J115trk))
#crs(hf)   
#compareCRS(hf,get_crs(J115trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J115/UTM_b115_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J115/UTM_b115_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J116
#' b116_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.
hf <-"D:/GISbases/b116JaguarAlbersBuffer70/b116_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J116trk))
#crs(hf)   
#compareCRS(hf,get_crs(J116trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J116/UTM_b116_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J116/UTM_b116_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J117
#' b117_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.
hf <-"D:/GISbases/b117JaguarAlbersBuffer70/b117_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J117trk))
#crs(hf)   
#compareCRS(hf,get_crs(J117trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J117/UTM_b117_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J117/UTM_b117_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

###################################################################################################################################
# Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif

#' J1
#' Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif
hf <-"D:/GISbases/b001JaguarAlbersBuffer70/b001_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J1trk))
#crs(hf)   
#compareCRS(hf,get_crs(J1trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J1/UTM_b001_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J1/UTM_b001_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J2
#' b002_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b002JaguarAlbersBuffer70/b002_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J2trk))
#crs(hf)   
#compareCRS(hf,get_crs(J2trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J2/UTM_b002_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J2/UTM_b002_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J3
#' b003_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b003JaguarAlbersBuffer70/b003_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J3trk))
#crs(hf)   
#compareCRS(hf,get_crs(J3trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J3/UTM_b003_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J3/UTM_b003_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J4
#' b004_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b004JaguarAlbersBuffer70/b004_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J4trk))
#crs(hf)   
#compareCRS(hf,get_crs(J4trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J4/UTM_b004_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J4/UTM_b004_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J5
#' b005_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b005JaguarAlbersBuffer70/b005_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J5trk))
#crs(hf)   
#compareCRS(hf,get_crs(J5trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J5/UTM_b005_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J5/UTM_b005_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J6
#' b006_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b006JaguarAlbersBuffer70/b006_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J6trk))
#crs(hf)   
#compareCRS(hf,get_crs(J6trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J6/UTM_b006_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J6/UTM_b006_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J7
#' b007_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b007JaguarAlbersBuffer70/b007_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J7trk))
#crs(hf)   
#compareCRS(hf,get_crs(J7trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J7/UTM_b007_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J7/UTM_b007_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J8
#' b008_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b008JaguarAlbersBuffer70/b008_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J8trk))
#crs(hf)   
#compareCRS(hf,get_crs(J8trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J8/UTM_b008_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J8/UTM_b008_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J9
#' b009_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b009JaguarAlbersBuffer70/b009_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J9trk))
#crs(hf)   
#compareCRS(hf,get_crs(J9trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J9/UTM_b009_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J9/UTM_b009_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J10
#' b010_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b010JaguarAlbersBuffer70/b010_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J10trk))
#crs(hf)   
#compareCRS(hf,get_crs(J10trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J10/UTM_b010_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J10/UTM_b010_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J11
#' b011_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b011JaguarAlbersBuffer70/b011_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J11trk))
#crs(hf)   
#compareCRS(hf,get_crs(J11trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J11/UTM_b011_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J11/UTM_b011_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J12
#' b012_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b012JaguarAlbersBuffer70/b012_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J12trk))
#crs(hf)   
#compareCRS(hf,get_crs(J12trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J12/UTM_b012_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J12/UTM_b012_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J13
#' b013_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b013JaguarAlbersBuffer70/b013_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J13trk))
#crs(hf)   
#compareCRS(hf,get_crs(J13trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J13/UTM_b013_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J13/UTM_b013_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J14
#' b014_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b014JaguarAlbersBuffer70/b014_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J14trk))
#crs(hf)   
#compareCRS(hf,get_crs(J14trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J14/UTM_b014_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J14/UTM_b014_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J15
#' b015_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b015JaguarAlbersBuffer70/b015_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J15trk))
#crs(hf)   
#compareCRS(hf,get_crs(J15trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J15/UTM_b015_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J15/UTM_b015_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J16
#' b016_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b016JaguarAlbersBuffer70/b016_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J16trk))
#crs(hf)   
#compareCRS(hf,get_crs(J16trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J16/UTM_b016_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J16/UTM_b016_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J17
#' b017_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b017JaguarAlbersBuffer70/b017_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J17trk))
#crs(hf)   
#compareCRS(hf,get_crs(J17trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J17/UTM_b017_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J17/UTM_b017_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J18
#' b018_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b018JaguarAlbersBuffer70/b018_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J18trk))
#crs(hf)   
#compareCRS(hf,get_crs(J18trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J18/UTM_b018_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J18/UTM_b018_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J19
#' b019_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b019JaguarAlbersBuffer70/b019_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J19trk))
#crs(hf)   
#compareCRS(hf,get_crs(J19trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J19/UTM_b019_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J19/UTM_b019_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J20
#' b020_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b020JaguarAlbersBuffer70/b020_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J20trk))
#crs(hf)   
#compareCRS(hf,get_crs(J20trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J20/UTM_b020_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J20/UTM_b020_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J21
#' b021_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b021JaguarAlbersBuffer70/b021_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J21trk))
#crs(hf)   
#compareCRS(hf,get_crs(J21trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J21/UTM_b021_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J21/UTM_b021_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J22
#' b022_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b022JaguarAlbersBuffer70/b022_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J22trk))
#crs(hf)   
#compareCRS(hf,get_crs(J22trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J22/UTM_b022_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J22/UTM_b022_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J23
#' b023_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b023JaguarAlbersBuffer70/b023_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J23trk))
#crs(hf)   
#compareCRS(hf,get_crs(J23trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J23/UTM_b023_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J23/UTM_b023_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J24
#' b024_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b024JaguarAlbersBuffer70/b024_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J24trk))
#crs(hf)   
#compareCRS(hf,get_crs(J24trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J24/UTM_b024_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J24/UTM_b024_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J25
#' b025_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b025JaguarAlbersBuffer70/b025_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J25trk))
#crs(hf)   
#compareCRS(hf,get_crs(J25trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J25/UTM_b025_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J25/UTM_b025_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J26
#' b026_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b026JaguarAlbersBuffer70/b026_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J26trk))
#crs(hf)   
#compareCRS(hf,get_crs(J26trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J26/UTM_b026_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J26/UTM_b026_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J27
#' b027_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b027JaguarAlbersBuffer70/b027_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J27trk))
#crs(hf)   
#compareCRS(hf,get_crs(J27trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J27/UTM_b027_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J27/UTM_b027_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J28
#' b028_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b028JaguarAlbersBuffer70/b028_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J28trk))
#crs(hf)   
#compareCRS(hf,get_crs(J28trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J28/UTM_b028_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J28/UTM_b028_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J29
#' b029_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b029JaguarAlbersBuffer70/b029_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J29trk))
#crs(hf)   
#compareCRS(hf,get_crs(J29trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J29/UTM_b029_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J29/UTM_b029_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J30
#' b030_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b030JaguarAlbersBuffer70/b030_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J30trk))
#crs(hf)   
#compareCRS(hf,get_crs(J30trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J30/UTM_b030_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J30/UTM_b030_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J30
#' b030_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b030JaguarAlbersBuffer70/b030_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J30trk))
#crs(hf)   
#compareCRS(hf,get_crs(J30trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J30/UTM_b030_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J30/UTM_b030_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J31
#' b031_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b031JaguarAlbersBuffer70/b031_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J31trk))
#crs(hf)   
#compareCRS(hf,get_crs(J31trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J31/UTM_b031_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J31/UTM_b031_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J32
#' b032_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b032JaguarAlbersBuffer70/b032_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J32trk))
#crs(hf)   
#compareCRS(hf,get_crs(J32trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J32/UTM_b032_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J32/UTM_b032_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J33
#' b033_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b033JaguarAlbersBuffer70/b033_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J33trk))
#crs(hf)   
#compareCRS(hf,get_crs(J33trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J33/UTM_b033_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J33/UTM_b033_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J34
#' b034_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b034JaguarAlbersBuffer70/b034_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J34trk))
#crs(hf)   
#compareCRS(hf,get_crs(J34trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J34/UTM_b034_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J34/UTM_b034_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J35
#' b035_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b035JaguarAlbersBuffer70/b035_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J35trk))
#crs(hf)   
#compareCRS(hf,get_crs(J35trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J35/UTM_b035_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J35/UTM_b035_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J36
#' b036_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b036JaguarAlbersBuffer70/b036_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J36trk))
#crs(hf)   
#compareCRS(hf,get_crs(J36trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J36/UTM_b036_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J36/UTM_b036_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J37
#' b037_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b037JaguarAlbersBuffer70/b037_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J37trk))
#crs(hf)   
#compareCRS(hf,get_crs(J37trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J37/UTM_b037_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J37/UTM_b037_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J38
#' b038_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b038JaguarAlbersBuffer70/b038_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J38trk))
#crs(hf)   
#compareCRS(hf,get_crs(J38trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J38/UTM_b038_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J38/UTM_b038_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J39
#' b039_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b039JaguarAlbersBuffer70/b039_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J39trk))
#crs(hf)   
#compareCRS(hf,get_crs(J39trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J39/UTM_b039_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J39/UTM_b039_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J40
#' b040_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b040JaguarAlbersBuffer70/b040_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J40trk))
#crs(hf)   
#compareCRS(hf,get_crs(J40trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J40/UTM_b040_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J40/UTM_b040_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J41
#' b041_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b041JaguarAlbersBuffer70/b041_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J41trk))
#crs(hf)   
#compareCRS(hf,get_crs(J41trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J41/UTM_b041_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J41/UTM_b041_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J42
#' b042_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b042JaguarAlbersBuffer70/b042_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J42trk))
#crs(hf)   
#compareCRS(hf,get_crs(J42trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J42/UTM_b042_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J42/UTM_b042_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J43
#' b043_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b043JaguarAlbersBuffer70/b043_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J43trk))
#crs(hf)   
#compareCRS(hf,get_crs(J43trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J43/UTM_b043_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J43/UTM_b043_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J44
#' b044_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b044JaguarAlbersBuffer70/b044_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J44trk))
#crs(hf)   
#compareCRS(hf,get_crs(J44trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J44/UTM_b044_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J44/UTM_b044_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J45
#' b045_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b045JaguarAlbersBuffer70/b045_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J45trk))
#crs(hf)   
#compareCRS(hf,get_crs(J45trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J45/UTM_b045_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J45/UTM_b045_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J46
#' b046_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b046JaguarAlbersBuffer70/b046_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J46trk))
#crs(hf)   
#compareCRS(hf,get_crs(J46trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J46/UTM_b046_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J46/UTM_b046_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J47
#' b047_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b047JaguarAlbersBuffer70/b047_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J47trk))
#crs(hf)   
#compareCRS(hf,get_crs(J47trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J47/UTM_b047_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J47/UTM_b047_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J48
#' b048_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b048JaguarAlbersBuffer70/b048_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J48trk))
#crs(hf)   
#compareCRS(hf,get_crs(J48trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J48/UTM_b048_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J48/UTM_b048_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J49
#' b049_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b049JaguarAlbersBuffer70/b049_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J49trk))
#crs(hf)   
#compareCRS(hf,get_crs(J49trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J49/UTM_b049_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J49/UTM_b049_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J50
#' b050_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b050JaguarAlbersBuffer70/b050_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J50trk))
#crs(hf)   
#compareCRS(hf,get_crs(J50trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J50/UTM_b050_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J50/UTM_b050_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J51
#' b051_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b051JaguarAlbersBuffer70/b051_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J51trk))
#crs(hf)   
#compareCRS(hf,get_crs(J51trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J51/UTM_b051_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J51/UTM_b051_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J52
#' b052_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b052JaguarAlbersBuffer70/b052_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J52trk))
#crs(hf)   
#compareCRS(hf,get_crs(J52trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J52/UTM_b052_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J52/UTM_b052_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J53
#' b053_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b053JaguarAlbersBuffer70/b053_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J53trk))
#crs(hf)   
#compareCRS(hf,get_crs(J53trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J53/UTM_b053_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J53/UTM_b053_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J54
#' b054_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b054JaguarAlbersBuffer70/b054_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J54trk))
#crs(hf)   
#compareCRS(hf,get_crs(J54trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J54/UTM_b054_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J54/UTM_b054_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J55
#' b055_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b055JaguarAlbersBuffer70/b055_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J55trk))
#crs(hf)   
#compareCRS(hf,get_crs(J55trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J55/UTM_b055_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J55/UTM_b055_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J56
#' b056_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b056JaguarAlbersBuffer70/b056_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J56trk))
#crs(hf)   
#compareCRS(hf,get_crs(J56trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J56/UTM_b056_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J56/UTM_b056_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J57
#' b057_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b057JaguarAlbersBuffer70/b057_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J57trk))
#crs(hf)   
#compareCRS(hf,get_crs(J57trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J57/UTM_b057_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J57/UTM_b057_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J58
#' b058_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b058JaguarAlbersBuffer70/b058_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J58trk))
#crs(hf)   
#compareCRS(hf,get_crs(J58trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J58/UTM_b058_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J58/UTM_b058_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J59
#' b059_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b059JaguarAlbersBuffer70/b059_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J59trk))
#crs(hf)   
#compareCRS(hf,get_crs(J59trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J59/UTM_b059_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J59/UTM_b059_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J60
#' b060_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b060JaguarAlbersBuffer70/b060_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J60trk))
#crs(hf)   
#compareCRS(hf,get_crs(J60trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J60/UTM_b060_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J60/UTM_b060_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J61
#' b061_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b061JaguarAlbersBuffer70/b061_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J61trk))
#crs(hf)   
#compareCRS(hf,get_crs(J61trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J61/UTM_b061_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J61/UTM_b061_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J62
#' b062_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b062JaguarAlbersBuffer70/b062_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J62trk))
#crs(hf)   
#compareCRS(hf,get_crs(J62trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J62/UTM_b062_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J62/UTM_b062_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J63
#' b063_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b063JaguarAlbersBuffer70/b063_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J63trk))
#crs(hf)   
#compareCRS(hf,get_crs(J63trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J63/UTM_b063_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J63/UTM_b063_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J64
#' b064_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b064JaguarAlbersBuffer70/b064_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J64trk))
#crs(hf)   
#compareCRS(hf,get_crs(J64trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J64/UTM_b064_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J64/UTM_b064_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J65
#' b065_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b065JaguarAlbersBuffer70/b065_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J65trk))
#crs(hf)   
#compareCRS(hf,get_crs(J65trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J65/UTM_b065_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J65/UTM_b065_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J66
#' b066_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b066JaguarAlbersBuffer70/b066_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J66trk))
#crs(hf)   
#compareCRS(hf,get_crs(J66trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J66/UTM_b066_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J66/UTM_b066_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J67
#' b067_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b067JaguarAlbersBuffer70/b067_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J67trk))
#crs(hf)   
#compareCRS(hf,get_crs(J67trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J67/UTM_b067_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J67/UTM_b067_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J68
#' b068_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b068JaguarAlbersBuffer70/b068_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J68trk))
#crs(hf)   
#compareCRS(hf,get_crs(J68trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J68/UTM_b068_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J68/UTM_b068_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J69
#' b069_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b069JaguarAlbersBuffer70/b069_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J69trk))
#crs(hf)   
#compareCRS(hf,get_crs(J69trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J69/UTM_b069_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J69/UTM_b069_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J70
#' b070_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b070JaguarAlbersBuffer70/b070_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J70trk))
#crs(hf)   
#compareCRS(hf,get_crs(J70trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J70/UTM_b070_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J70/UTM_b070_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J71
#' b071_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b071JaguarAlbersBuffer70/b071_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J71trk))
#crs(hf)   
#compareCRS(hf,get_crs(J71trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J71/UTM_b071_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J71/UTM_b071_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J72
#' b072_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b072JaguarAlbersBuffer70/b072_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J72trk))
#crs(hf)   
#compareCRS(hf,get_crs(J72trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J72/UTM_b072_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J72/UTM_b072_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J73
#' b073_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b073JaguarAlbersBuffer70/b073_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J73trk))
#crs(hf)   
#compareCRS(hf,get_crs(J73trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J73/UTM_b073_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J73/UTM_b073_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J74
#' b074_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b074JaguarAlbersBuffer70/b074_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J74trk))
#crs(hf)   
#compareCRS(hf,get_crs(J74trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J74/UTM_b074_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J74/UTM_b074_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J75
#' b075_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b075JaguarAlbersBuffer70/b075_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J75trk))
#crs(hf)   
#compareCRS(hf,get_crs(J75trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J75/UTM_b075_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J75/UTM_b075_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J76
#' b076_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b076JaguarAlbersBuffer70/b076_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J76trk))
#crs(hf)   
#compareCRS(hf,get_crs(J76trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J76/UTM_b076_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J76/UTM_b076_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J77
#' b077_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b077JaguarAlbersBuffer70/b077_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J77trk))
#crs(hf)   
#compareCRS(hf,get_crs(J77trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J77/UTM_b077_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J77/UTM_b077_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J78
#' b078_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b078JaguarAlbersBuffer70/b078_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J78trk))
#crs(hf)   
#compareCRS(hf,get_crs(J78trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J78/UTM_b078_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J78/UTM_b078_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J79
#' b079_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b079JaguarAlbersBuffer70/b079_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J79trk))
#crs(hf)   
#compareCRS(hf,get_crs(J79trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J79/UTM_b079_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J79/UTM_b079_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J80
#' b080_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b080JaguarAlbersBuffer70/b080_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J80trk))
#crs(hf)   
#compareCRS(hf,get_crs(J80trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J80/UTM_b080_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J80/UTM_b080_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J81
#' b081_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b081JaguarAlbersBuffer70/b081_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J81trk))
#crs(hf)   
#compareCRS(hf,get_crs(J81trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J81/UTM_b081_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J81/UTM_b081_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J82
#' b082_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b082JaguarAlbersBuffer70/b082_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J82trk))
#crs(hf)   
#compareCRS(hf,get_crs(J82trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J82/UTM_b082_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J82/UTM_b082_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J83
#' b083_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b083JaguarAlbersBuffer70/b083_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J83trk))
#crs(hf)   
#compareCRS(hf,get_crs(J83trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J83/UTM_b083_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J83/UTM_b083_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J84
#' b084_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b084JaguarAlbersBuffer70/b084_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J84trk))
#crs(hf)   
#compareCRS(hf,get_crs(J84trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J84/UTM_b084_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J84/UTM_b084_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J85
#' b085_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b085JaguarAlbersBuffer70/b085_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J85trk))
#crs(hf)   
#compareCRS(hf,get_crs(J85trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J85/UTM_b085_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J85/UTM_b085_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J86
#' b086_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b086JaguarAlbersBuffer70/b086_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J86trk))
#crs(hf)   
#compareCRS(hf,get_crs(J86trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J86/UTM_b086_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J86/UTM_b086_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J87
#' b087_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b087JaguarAlbersBuffer70/b087_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J87trk))
#crs(hf)   
#compareCRS(hf,get_crs(J87trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J87/UTM_b087_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J87/UTM_b087_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J88
#' b088_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b088JaguarAlbersBuffer70/b088_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J88trk))
#crs(hf)   
#compareCRS(hf,get_crs(J88trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J88/UTM_b088_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J88/UTM_b088_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J89
#' b089_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b089JaguarAlbersBuffer70/b089_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J89trk))
#crs(hf)   
#compareCRS(hf,get_crs(J89trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J89/UTM_b089_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J89/UTM_b089_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J90
#' b090_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b090JaguarAlbersBuffer70/b090_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J90trk))
#crs(hf)   
#compareCRS(hf,get_crs(J90trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J90/UTM_b090_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J90/UTM_b090_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J91
#' b091_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b091JaguarAlbersBuffer70/b091_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J91trk))
#crs(hf)   
#compareCRS(hf,get_crs(J91trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J91/UTM_b091_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J91/UTM_b091_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J92
#' b092_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b092JaguarAlbersBuffer70/b092_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J92trk))
#crs(hf)   
#compareCRS(hf,get_crs(J92trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J92/UTM_b092_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J92/UTM_b092_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J93
#' b093_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b093JaguarAlbersBuffer70/b093_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J93trk))
#crs(hf)   
#compareCRS(hf,get_crs(J93trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J93/UTM_b093_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J93/UTM_b093_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J94
#' b094_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b094JaguarAlbersBuffer70/b094_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J94trk))
#crs(hf)   
#compareCRS(hf,get_crs(J94trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J94/UTM_b094_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J94/UTM_b094_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J95
#' b095_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b095JaguarAlbersBuffer70/b095_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J95trk))
#crs(hf)   
#compareCRS(hf,get_crs(J95trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J95/UTM_b095_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J95/UTM_b095_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J96
#' b096_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b096JaguarAlbersBuffer70/b096_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J96trk))
#crs(hf)   
#compareCRS(hf,get_crs(J96trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J96/UTM_b096_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J96/UTM_b096_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J97
#' b097_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b097JaguarAlbersBuffer70/b097_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J97trk))
#crs(hf)   
#compareCRS(hf,get_crs(J97trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J97/UTM_b097_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J97/UTM_b097_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J98
#' b098_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b098JaguarAlbersBuffer70/b098_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J98trk))
#crs(hf)   
#compareCRS(hf,get_crs(J98trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J98/UTM_b098_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J98/UTM_b098_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J99
#' b099_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b099JaguarAlbersBuffer70/b099_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J99trk))
#crs(hf)   
#compareCRS(hf,get_crs(J99trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J99/UTM_b099_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J99/UTM_b099_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J100
#' b100_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b100JaguarAlbersBuffer70/b100_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J100trk))
#crs(hf)   
#compareCRS(hf,get_crs(J100trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J100/UTM_b100_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J100/UTM_b010_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J101
#' b101_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b101JaguarAlbersBuffer70/b101_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J101trk))
#crs(hf)   
#compareCRS(hf,get_crs(J101trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J101/UTM_b101_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J101/UTM_b101_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J102
#' b102_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b102JaguarAlbersBuffer70/b102_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J102trk))
#crs(hf)   
#compareCRS(hf,get_crs(J102trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J102/UTM_b102_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J102/UTM_b102_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J103
#' b103_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b103JaguarAlbersBuffer70/b103_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J103trk))
#crs(hf)   
#compareCRS(hf,get_crs(J103trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J103/UTM_b103_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J103/UTM_b103_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J104
#' b104_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b104JaguarAlbersBuffer70/b104_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J104trk))
#crs(hf)   
#compareCRS(hf,get_crs(J104trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J104/UTM_b104_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J104/UTM_b104_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J105
#' b105_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b105JaguarAlbersBuffer70/b105_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J105trk))
#crs(hf)   
#compareCRS(hf,get_crs(J105trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J105/UTM_b105_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J105/UTM_b105_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J106
#' b106_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b106JaguarAlbersBuffer70/b106_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J106trk))
#crs(hf)   
#compareCRS(hf,get_crs(J106trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J106/UTM_b106_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J106/UTM_b106_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J107
#' b107_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b107JaguarAlbersBuffer70/b107_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J107trk))
#crs(hf)   
#compareCRS(hf,get_crs(J107trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J107/UTM_b107_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J107/UTM_b107_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J108
#' b108_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b108JaguarAlbersBuffer70/b108_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J108trk))
#crs(hf)   
#compareCRS(hf,get_crs(J108trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J108/UTM_b108_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J108/UTM_b108_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J109
#' b109_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b109JaguarAlbersBuffer70/b109_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J109trk))
#crs(hf)   
#compareCRS(hf,get_crs(J109trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J109/UTM_b109_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J109/UTM_b109_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J110
#' b110_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b110JaguarAlbersBuffer70/b110_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J110trk))
#crs(hf)   
#compareCRS(hf,get_crs(J110trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J110/UTM_b110_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J110/UTM_b110_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J111
#' b111_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b111JaguarAlbersBuffer70/b111_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J111trk))
#crs(hf)   
#compareCRS(hf,get_crs(J111trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J111/UTM_b111_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J111/UTM_b111_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J112
#' b112_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b112JaguarAlbersBuffer70/b112_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J112trk))
#crs(hf)   
#compareCRS(hf,get_crs(J112trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J112/UTM_b112_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J112/UTM_b112_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J113
#' b113_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b113JaguarAlbersBuffer70/b113_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J113trk))
#crs(hf)   
#compareCRS(hf,get_crs(J113trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J113/UTM_b113_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J113/UTM_b113_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J114
#' b114_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b114JaguarAlbersBuffer70/b114_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J114trk))
#crs(hf)   
#compareCRS(hf,get_crs(J114trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J114/UTM_b114_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J114/UTM_b114_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J115
#' b115_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b115JaguarAlbersBuffer70/b115_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J115trk))
#crs(hf)   
#compareCRS(hf,get_crs(J115trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J115/UTM_b115_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J115/UTM_b115_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J116
#' b116_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b116JaguarAlbersBuffer70/b116_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J116trk))
#crs(hf)   
#compareCRS(hf,get_crs(J116trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J116/UTM_b116_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J116/UTM_b116_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J117
#' b117_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b117JaguarAlbersBuffer70/b117_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J117trk))
#crs(hf)   
#compareCRS(hf,get_crs(J117trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J117/UTM_b117_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J117/UTM_b117_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

######################################################################################################################################

# Neotropic_Earthenv_dem90m_tif_exp.
# Neotropic_Earthenv_dem90m_tif_exp.

#' J1
#' Neotropic_Earthenv_dem90m_tif_exp.tif
hf <-"D:/GISbases/b001JaguarAlbersBuffer70/b001_Neotropic_Earthenv_dem90m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J1trk))
#crs(hf)   
#compareCRS(hf,get_crs(J1trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J1/UTM_b001_Neotropic_Earthenv_dem90m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J1/UTM_b001_Neotropic_Earthenv_dem90m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J2
#' b002_Neotropic_Earthenv_dem90m_tif_exp.
hf <-"D:/GISbases/b002JaguarAlbersBuffer70/b002_Neotropic_Earthenv_dem90m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J2trk))
#crs(hf)   
#compareCRS(hf,get_crs(J2trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J2/UTM_b002_Neotropic_Earthenv_dem90m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J2/UTM_b002_Neotropic_Earthenv_dem90m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J3
#' b003_Neotropic_Earthenv_dem90m_tif_exp.
hf <-"D:/GISbases/b003JaguarAlbersBuffer70/b003_Neotropic_Earthenv_dem90m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J3trk))
#crs(hf)   
#compareCRS(hf,get_crs(J3trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J3/UTM_b003_Neotropic_Earthenv_dem90m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J3/UTM_b003_Neotropic_Earthenv_dem90m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J4
#' b004_Neotropic_Earthenv_dem90m_tif_exp.
hf <-"D:/GISbases/b004JaguarAlbersBuffer70/b004_Neotropic_Earthenv_dem90m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J4trk))
#crs(hf)   
#compareCRS(hf,get_crs(J4trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J4/UTM_b004_Neotropic_Earthenv_dem90m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J4/UTM_b004_Neotropic_Earthenv_dem90m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J5
#' b005_Neotropic_Earthenv_dem90m_tif_exp.
hf <-"D:/GISbases/b005JaguarAlbersBuffer70/b005_Neotropic_Earthenv_dem90m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J5trk))
#crs(hf)   
#compareCRS(hf,get_crs(J5trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J5/UTM_b005_Neotropic_Earthenv_dem90m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J5/UTM_b005_Neotropic_Earthenv_dem90m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J6
#' b006_Neotropic_Earthenv_dem90m_tif_exp.
hf <-"D:/GISbases/b006JaguarAlbersBuffer70/b006_Neotropic_Earthenv_dem90m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J6trk))
#crs(hf)   
#compareCRS(hf,get_crs(J6trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J6/UTM_b006_Neotropic_Earthenv_dem90m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J6/UTM_b006_Neotropic_Earthenv_dem90m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J7
#' b007_Neotropic_Earthenv_dem90m_tif_exp.
hf <-"D:/GISbases/b007JaguarAlbersBuffer70/b007_Neotropic_Earthenv_dem90m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J7trk))
#crs(hf)   
#compareCRS(hf,get_crs(J7trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J7/UTM_b007_Neotropic_Earthenv_dem90m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J7/UTM_b007_Neotropic_Earthenv_dem90m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J8
#' b008_Neotropic_Earthenv_dem90m_tif_exp.
hf <-"D:/GISbases/b008JaguarAlbersBuffer70/b008_Neotropic_Earthenv_dem90m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J8trk))
#crs(hf)   
#compareCRS(hf,get_crs(J8trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J8/UTM_b008_Neotropic_Earthenv_dem90m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J8/UTM_b008_Neotropic_Earthenv_dem90m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J9
#' b009_Neotropic_Earthenv_dem90m_tif_exp.
hf <-"D:/GISbases/b009JaguarAlbersBuffer70/b009_Neotropic_Earthenv_dem90m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J9trk))
#crs(hf)   
#compareCRS(hf,get_crs(J9trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J9/UTM_b009_Neotropic_Earthenv_dem90m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J9/UTM_b009_Neotropic_Earthenv_dem90m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J10
#' b010_Neotropic_Earthenv_dem90m_tif_exp.
hf <-"D:/GISbases/b010JaguarAlbersBuffer70/b010_Neotropic_Earthenv_dem90m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J10trk))
#crs(hf)   
#compareCRS(hf,get_crs(J10trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J10/UTM_b010_Neotropic_Earthenv_dem90m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J10/UTM_b010_Neotropic_Earthenv_dem90m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J11
#' b011_Neotropic_Earthenv_dem90m_tif_exp.
hf <-"D:/GISbases/b011JaguarAlbersBuffer70/b011_Neotropic_Earthenv_dem90m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J11trk))
#crs(hf)   
#compareCRS(hf,get_crs(J11trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J11/UTM_b011_Neotropic_Earthenv_dem90m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J11/UTM_b011_Neotropic_Earthenv_dem90m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J12
#' b012_Neotropic_Earthenv_dem90m_tif_exp.
hf <-"D:/GISbases/b012JaguarAlbersBuffer70/b012_Neotropic_Earthenv_dem90m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J12trk))
#crs(hf)   
#compareCRS(hf,get_crs(J12trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J12/UTM_b012_Neotropic_Earthenv_dem90m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J12/UTM_b012_Neotropic_Earthenv_dem90m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J13
#' b013_Neotropic_Earthenv_dem90m_tif_exp.
hf <-"D:/GISbases/b013JaguarAlbersBuffer70/b013_Neotropic_Earthenv_dem90m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J13trk))
#crs(hf)   
#compareCRS(hf,get_crs(J13trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J13/UTM_b013_Neotropic_Earthenv_dem90m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J13/UTM_b013_Neotropic_Earthenv_dem90m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J14
#' b014_Neotropic_Earthenv_dem90m_tif_exp.
hf <-"D:/GISbases/b014JaguarAlbersBuffer70/b014_Neotropic_Earthenv_dem90m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J14trk))
#crs(hf)   
#compareCRS(hf,get_crs(J14trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J14/UTM_b014_Neotropic_Earthenv_dem90m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J14/UTM_b014_Neotropic_Earthenv_dem90m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J15
#' b015_Neotropic_Earthenv_dem90m_tif_exp.
hf <-"D:/GISbases/b015JaguarAlbersBuffer70/b015_Neotropic_Earthenv_dem90m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J15trk))
#crs(hf)   
#compareCRS(hf,get_crs(J15trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J15/UTM_b015_Neotropic_Earthenv_dem90m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J15/UTM_b015_Neotropic_Earthenv_dem90m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J16
#' b016_Neotropic_Earthenv_dem90m_tif_exp.
hf <-"D:/GISbases/b016JaguarAlbersBuffer70/b016_Neotropic_Earthenv_dem90m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J16trk))
#crs(hf)   
#compareCRS(hf,get_crs(J16trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J16/UTM_b016_Neotropic_Earthenv_dem90m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J16/UTM_b016_Neotropic_Earthenv_dem90m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J17
#' b017_Neotropic_Earthenv_dem90m_tif_exp.
hf <-"D:/GISbases/b017JaguarAlbersBuffer70/b017_Neotropic_Earthenv_dem90m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J17trk))
#crs(hf)   
#compareCRS(hf,get_crs(J17trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J17/UTM_b017_Neotropic_Earthenv_dem90m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J17/UTM_b017_Neotropic_Earthenv_dem90m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J18
#' b018_Neotropic_Earthenv_dem90m_tif_exp.
hf <-"D:/GISbases/b018JaguarAlbersBuffer70/b018_Neotropic_Earthenv_dem90m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J18trk))
#crs(hf)   
#compareCRS(hf,get_crs(J18trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J18/UTM_b018_Neotropic_Earthenv_dem90m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J18/UTM_b018_Neotropic_Earthenv_dem90m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J19
#' b019_Neotropic_Earthenv_dem90m_tif_exp.
hf <-"D:/GISbases/b019JaguarAlbersBuffer70/b019_Neotropic_Earthenv_dem90m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J19trk))
#crs(hf)   
#compareCRS(hf,get_crs(J19trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J19/UTM_b019_Neotropic_Earthenv_dem90m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J19/UTM_b019_Neotropic_Earthenv_dem90m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J20
#' b020_Neotropic_Earthenv_dem90m_tif_exp.
hf <-"D:/GISbases/b020JaguarAlbersBuffer70/b020_Neotropic_Earthenv_dem90m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J20trk))
#crs(hf)   
#compareCRS(hf,get_crs(J20trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J20/UTM_b020_Neotropic_Earthenv_dem90m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J20/UTM_b020_Neotropic_Earthenv_dem90m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J21
#' b021_Neotropic_Earthenv_dem90m_tif_exp.
hf <-"D:/GISbases/b021JaguarAlbersBuffer70/b021_Neotropic_Earthenv_dem90m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J21trk))
#crs(hf)   
#compareCRS(hf,get_crs(J21trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J21/UTM_b021_Neotropic_Earthenv_dem90m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J21/UTM_b021_Neotropic_Earthenv_dem90m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J22
#' b022_Neotropic_Earthenv_dem90m_tif_exp.
hf <-"D:/GISbases/b022JaguarAlbersBuffer70/b022_Neotropic_Earthenv_dem90m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J22trk))
#crs(hf)   
#compareCRS(hf,get_crs(J22trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J22/UTM_b022_Neotropic_Earthenv_dem90m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J22/UTM_b022_Neotropic_Earthenv_dem90m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J23
#' b023_Neotropic_Earthenv_dem90m_tif_exp.
hf <-"D:/GISbases/b023JaguarAlbersBuffer70/b023_Neotropic_Earthenv_dem90m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J23trk))
#crs(hf)   
#compareCRS(hf,get_crs(J23trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J23/UTM_b023_Neotropic_Earthenv_dem90m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J23/UTM_b023_Neotropic_Earthenv_dem90m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J24
#' b024_Neotropic_Earthenv_dem90m_tif_exp.
hf <-"D:/GISbases/b024JaguarAlbersBuffer70/b024_Neotropic_Earthenv_dem90m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J24trk))
#crs(hf)   
#compareCRS(hf,get_crs(J24trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J24/UTM_b024_Neotropic_Earthenv_dem90m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J24/UTM_b024_Neotropic_Earthenv_dem90m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J25
#' b025_Neotropic_Earthenv_dem90m_tif_exp.
hf <-"D:/GISbases/b025JaguarAlbersBuffer70/b025_Neotropic_Earthenv_dem90m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J25trk))
#crs(hf)   
#compareCRS(hf,get_crs(J25trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J25/UTM_b025_Neotropic_Earthenv_dem90m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J25/UTM_b025_Neotropic_Earthenv_dem90m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J26
#' b026_Neotropic_Earthenv_dem90m_tif_exp.
hf <-"D:/GISbases/b026JaguarAlbersBuffer70/b026_Neotropic_Earthenv_dem90m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J26trk))
#crs(hf)   
#compareCRS(hf,get_crs(J26trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J26/UTM_b026_Neotropic_Earthenv_dem90m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J26/UTM_b026_Neotropic_Earthenv_dem90m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J27
#' b027_Neotropic_Earthenv_dem90m_tif_exp.
hf <-"D:/GISbases/b027JaguarAlbersBuffer70/b027_Neotropic_Earthenv_dem90m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J27trk))
#crs(hf)   
#compareCRS(hf,get_crs(J27trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J27/UTM_b027_Neotropic_Earthenv_dem90m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J27/UTM_b027_Neotropic_Earthenv_dem90m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J28
#' b028_Neotropic_Earthenv_dem90m_tif_exp.
hf <-"D:/GISbases/b028JaguarAlbersBuffer70/b028_Neotropic_Earthenv_dem90m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J28trk))
#crs(hf)   
#compareCRS(hf,get_crs(J28trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J28/UTM_b028_Neotropic_Earthenv_dem90m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J28/UTM_b028_Neotropic_Earthenv_dem90m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J29
#' b029_Neotropic_Earthenv_dem90m_tif_exp.
hf <-"D:/GISbases/b029JaguarAlbersBuffer70/b029_Neotropic_Earthenv_dem90m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J29trk))
#crs(hf)   
#compareCRS(hf,get_crs(J29trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J29/UTM_b029_Neotropic_Earthenv_dem90m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J29/UTM_b029_Neotropic_Earthenv_dem90m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J30
#' b030_Neotropic_Earthenv_dem90m_tif_exp.
hf <-"D:/GISbases/b030JaguarAlbersBuffer70/b030_Neotropic_Earthenv_dem90m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J30trk))
#crs(hf)   
#compareCRS(hf,get_crs(J30trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J30/UTM_b030_Neotropic_Earthenv_dem90m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J30/UTM_b030_Neotropic_Earthenv_dem90m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J30
#' b030_Neotropic_Earthenv_dem90m_tif_exp.
hf <-"D:/GISbases/b030JaguarAlbersBuffer70/b030_Neotropic_Earthenv_dem90m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J30trk))
#crs(hf)   
#compareCRS(hf,get_crs(J30trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J30/UTM_b030_Neotropic_Earthenv_dem90m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J30/UTM_b030_Neotropic_Earthenv_dem90m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J31
#' b031_Neotropic_Earthenv_dem90m_tif_exp.
hf <-"D:/GISbases/b031JaguarAlbersBuffer70/b031_Neotropic_Earthenv_dem90m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J31trk))
#crs(hf)   
#compareCRS(hf,get_crs(J31trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J31/UTM_b031_Neotropic_Earthenv_dem90m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J31/UTM_b031_Neotropic_Earthenv_dem90m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J32
#' b032_Neotropic_Earthenv_dem90m_tif_exp.
hf <-"D:/GISbases/b032JaguarAlbersBuffer70/b032_Neotropic_Earthenv_dem90m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J32trk))
#crs(hf)   
#compareCRS(hf,get_crs(J32trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J32/UTM_b032_Neotropic_Earthenv_dem90m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J32/UTM_b032_Neotropic_Earthenv_dem90m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J33
#' b033_Neotropic_Earthenv_dem90m_tif_exp.
hf <-"D:/GISbases/b033JaguarAlbersBuffer70/b033_Neotropic_Earthenv_dem90m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J33trk))
#crs(hf)   
#compareCRS(hf,get_crs(J33trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J33/UTM_b033_Neotropic_Earthenv_dem90m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J33/UTM_b033_Neotropic_Earthenv_dem90m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J34
#' b034_Neotropic_Earthenv_dem90m_tif_exp.
hf <-"D:/GISbases/b034JaguarAlbersBuffer70/b034_Neotropic_Earthenv_dem90m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J34trk))
#crs(hf)   
#compareCRS(hf,get_crs(J34trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J34/UTM_b034_Neotropic_Earthenv_dem90m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J34/UTM_b034_Neotropic_Earthenv_dem90m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J35
#' b035_Neotropic_Earthenv_dem90m_tif_exp.
hf <-"D:/GISbases/b035JaguarAlbersBuffer70/b035_Neotropic_Earthenv_dem90m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J35trk))
#crs(hf)   
#compareCRS(hf,get_crs(J35trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J35/UTM_b035_Neotropic_Earthenv_dem90m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J35/UTM_b035_Neotropic_Earthenv_dem90m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J36
#' b036_Neotropic_Earthenv_dem90m_tif_exp.
hf <-"D:/GISbases/b036JaguarAlbersBuffer70/b036_Neotropic_Earthenv_dem90m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J36trk))
#crs(hf)   
#compareCRS(hf,get_crs(J36trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J36/UTM_b036_Neotropic_Earthenv_dem90m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J36/UTM_b036_Neotropic_Earthenv_dem90m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J37
#' b037_Neotropic_Earthenv_dem90m_tif_exp.
hf <-"D:/GISbases/b037JaguarAlbersBuffer70/b037_Neotropic_Earthenv_dem90m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J37trk))
#crs(hf)   
#compareCRS(hf,get_crs(J37trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J37/UTM_b037_Neotropic_Earthenv_dem90m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J37/UTM_b037_Neotropic_Earthenv_dem90m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J38
#' b038_Neotropic_Earthenv_dem90m_tif_exp.
hf <-"D:/GISbases/b038JaguarAlbersBuffer70/b038_Neotropic_Earthenv_dem90m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J38trk))
#crs(hf)   
#compareCRS(hf,get_crs(J38trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J38/UTM_b038_Neotropic_Earthenv_dem90m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J38/UTM_b038_Neotropic_Earthenv_dem90m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J39
#' b039_Neotropic_Earthenv_dem90m_tif_exp.
hf <-"D:/GISbases/b039JaguarAlbersBuffer70/b039_Neotropic_Earthenv_dem90m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J39trk))
#crs(hf)   
#compareCRS(hf,get_crs(J39trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J39/UTM_b039_Neotropic_Earthenv_dem90m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J39/UTM_b039_Neotropic_Earthenv_dem90m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J40
#' b040_Neotropic_Earthenv_dem90m_tif_exp.
hf <-"D:/GISbases/b040JaguarAlbersBuffer70/b040_Neotropic_Earthenv_dem90m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J40trk))
#crs(hf)   
#compareCRS(hf,get_crs(J40trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J40/UTM_b040_Neotropic_Earthenv_dem90m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J40/UTM_b040_Neotropic_Earthenv_dem90m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J41
#' b041_Neotropic_Earthenv_dem90m_tif_exp.
hf <-"D:/GISbases/b041JaguarAlbersBuffer70/b041_Neotropic_Earthenv_dem90m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J41trk))
#crs(hf)   
#compareCRS(hf,get_crs(J41trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J41/UTM_b041_Neotropic_Earthenv_dem90m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J41/UTM_b041_Neotropic_Earthenv_dem90m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J42
#' b042_Neotropic_Earthenv_dem90m_tif_exp.
hf <-"D:/GISbases/b042JaguarAlbersBuffer70/b042_Neotropic_Earthenv_dem90m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J42trk))
#crs(hf)   
#compareCRS(hf,get_crs(J42trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J42/UTM_b042_Neotropic_Earthenv_dem90m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J42/UTM_b042_Neotropic_Earthenv_dem90m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J43
#' b043_Neotropic_Earthenv_dem90m_tif_exp.
hf <-"D:/GISbases/b043JaguarAlbersBuffer70/b043_Neotropic_Earthenv_dem90m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J43trk))
#crs(hf)   
#compareCRS(hf,get_crs(J43trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J43/UTM_b043_Neotropic_Earthenv_dem90m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J43/UTM_b043_Neotropic_Earthenv_dem90m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J44
#' b044_Neotropic_Earthenv_dem90m_tif_exp.
hf <-"D:/GISbases/b044JaguarAlbersBuffer70/b044_Neotropic_Earthenv_dem90m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J44trk))
#crs(hf)   
#compareCRS(hf,get_crs(J44trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J44/UTM_b044_Neotropic_Earthenv_dem90m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J44/UTM_b044_Neotropic_Earthenv_dem90m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J45
#' b045_Neotropic_Earthenv_dem90m_tif_exp.
hf <-"D:/GISbases/b045JaguarAlbersBuffer70/b045_Neotropic_Earthenv_dem90m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J45trk))
#crs(hf)   
#compareCRS(hf,get_crs(J45trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J45/UTM_b045_Neotropic_Earthenv_dem90m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J45/UTM_b045_Neotropic_Earthenv_dem90m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J46
#' b046_Neotropic_Earthenv_dem90m_tif_exp.
hf <-"D:/GISbases/b046JaguarAlbersBuffer70/b046_Neotropic_Earthenv_dem90m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J46trk))
#crs(hf)   
#compareCRS(hf,get_crs(J46trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J46/UTM_b046_Neotropic_Earthenv_dem90m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J46/UTM_b046_Neotropic_Earthenv_dem90m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J47
#' b047_Neotropic_Earthenv_dem90m_tif_exp.
hf <-"D:/GISbases/b047JaguarAlbersBuffer70/b047_Neotropic_Earthenv_dem90m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J47trk))
#crs(hf)   
#compareCRS(hf,get_crs(J47trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J47/UTM_b047_Neotropic_Earthenv_dem90m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J47/UTM_b047_Neotropic_Earthenv_dem90m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J48
#' b048_Neotropic_Earthenv_dem90m_tif_exp.
hf <-"D:/GISbases/b048JaguarAlbersBuffer70/b048_Neotropic_Earthenv_dem90m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J48trk))
#crs(hf)   
#compareCRS(hf,get_crs(J48trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J48/UTM_b048_Neotropic_Earthenv_dem90m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J48/UTM_b048_Neotropic_Earthenv_dem90m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J49
#' b049_Neotropic_Earthenv_dem90m_tif_exp.
hf <-"D:/GISbases/b049JaguarAlbersBuffer70/b049_Neotropic_Earthenv_dem90m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J49trk))
#crs(hf)   
#compareCRS(hf,get_crs(J49trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J49/UTM_b049_Neotropic_Earthenv_dem90m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J49/UTM_b049_Neotropic_Earthenv_dem90m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J50
#' b050_Neotropic_Earthenv_dem90m_tif_exp.
hf <-"D:/GISbases/b050JaguarAlbersBuffer70/b050_Neotropic_Earthenv_dem90m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J50trk))
#crs(hf)   
#compareCRS(hf,get_crs(J50trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J50/UTM_b050_Neotropic_Earthenv_dem90m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J50/UTM_b050_Neotropic_Earthenv_dem90m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J51
#' b051_Neotropic_Earthenv_dem90m_tif_exp.
hf <-"D:/GISbases/b051JaguarAlbersBuffer70/b051_Neotropic_Earthenv_dem90m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J51trk))
#crs(hf)   
#compareCRS(hf,get_crs(J51trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J51/UTM_b051_Neotropic_Earthenv_dem90m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J51/UTM_b051_Neotropic_Earthenv_dem90m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J52
#' b052_Neotropic_Earthenv_dem90m_tif_exp.
hf <-"D:/GISbases/b052JaguarAlbersBuffer70/b052_Neotropic_Earthenv_dem90m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J52trk))
#crs(hf)   
#compareCRS(hf,get_crs(J52trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J52/UTM_b052_Neotropic_Earthenv_dem90m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J52/UTM_b052_Neotropic_Earthenv_dem90m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J53
#' b053_Neotropic_Earthenv_dem90m_tif_exp.
hf <-"D:/GISbases/b053JaguarAlbersBuffer70/b053_Neotropic_Earthenv_dem90m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J53trk))
#crs(hf)   
#compareCRS(hf,get_crs(J53trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J53/UTM_b053_Neotropic_Earthenv_dem90m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J53/UTM_b053_Neotropic_Earthenv_dem90m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J54
#' b054_Neotropic_Earthenv_dem90m_tif_exp.
hf <-"D:/GISbases/b054JaguarAlbersBuffer70/b054_Neotropic_Earthenv_dem90m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J54trk))
#crs(hf)   
#compareCRS(hf,get_crs(J54trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J54/UTM_b054_Neotropic_Earthenv_dem90m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J54/UTM_b054_Neotropic_Earthenv_dem90m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J55
#' b055_Neotropic_Earthenv_dem90m_tif_exp.
hf <-"D:/GISbases/b055JaguarAlbersBuffer70/b055_Neotropic_Earthenv_dem90m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J55trk))
#crs(hf)   
#compareCRS(hf,get_crs(J55trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J55/UTM_b055_Neotropic_Earthenv_dem90m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J55/UTM_b055_Neotropic_Earthenv_dem90m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J56
#' b056_Neotropic_Earthenv_dem90m_tif_exp.
hf <-"D:/GISbases/b056JaguarAlbersBuffer70/b056_Neotropic_Earthenv_dem90m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J56trk))
#crs(hf)   
#compareCRS(hf,get_crs(J56trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J56/UTM_b056_Neotropic_Earthenv_dem90m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J56/UTM_b056_Neotropic_Earthenv_dem90m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J57
#' b057_Neotropic_Earthenv_dem90m_tif_exp.
hf <-"D:/GISbases/b057JaguarAlbersBuffer70/b057_Neotropic_Earthenv_dem90m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J57trk))
#crs(hf)   
#compareCRS(hf,get_crs(J57trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J57/UTM_b057_Neotropic_Earthenv_dem90m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J57/UTM_b057_Neotropic_Earthenv_dem90m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J58
#' b058_Neotropic_Earthenv_dem90m_tif_exp.
hf <-"D:/GISbases/b058JaguarAlbersBuffer70/b058_Neotropic_Earthenv_dem90m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J58trk))
#crs(hf)   
#compareCRS(hf,get_crs(J58trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J58/UTM_b058_Neotropic_Earthenv_dem90m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J58/UTM_b058_Neotropic_Earthenv_dem90m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J59
#' b059_Neotropic_Earthenv_dem90m_tif_exp.
hf <-"D:/GISbases/b059JaguarAlbersBuffer70/b059_Neotropic_Earthenv_dem90m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J59trk))
#crs(hf)   
#compareCRS(hf,get_crs(J59trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J59/UTM_b059_Neotropic_Earthenv_dem90m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J59/UTM_b059_Neotropic_Earthenv_dem90m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J60
#' b060_Neotropic_Earthenv_dem90m_tif_exp.
hf <-"D:/GISbases/b060JaguarAlbersBuffer70/b060_Neotropic_Earthenv_dem90m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J60trk))
#crs(hf)   
#compareCRS(hf,get_crs(J60trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J60/UTM_b060_Neotropic_Earthenv_dem90m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J60/UTM_b060_Neotropic_Earthenv_dem90m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J61
#' b061_Neotropic_Earthenv_dem90m_tif_exp.
hf <-"D:/GISbases/b061JaguarAlbersBuffer70/b061_Neotropic_Earthenv_dem90m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J61trk))
#crs(hf)   
#compareCRS(hf,get_crs(J61trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J61/UTM_b061_Neotropic_Earthenv_dem90m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J61/UTM_b061_Neotropic_Earthenv_dem90m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J62
#' b062_Neotropic_Earthenv_dem90m_tif_exp.
hf <-"D:/GISbases/b062JaguarAlbersBuffer70/b062_Neotropic_Earthenv_dem90m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J62trk))
#crs(hf)   
#compareCRS(hf,get_crs(J62trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J62/UTM_b062_Neotropic_Earthenv_dem90m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J62/UTM_b062_Neotropic_Earthenv_dem90m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J63
#' b063_Neotropic_Earthenv_dem90m_tif_exp.
hf <-"D:/GISbases/b063JaguarAlbersBuffer70/b063_Neotropic_Earthenv_dem90m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J63trk))
#crs(hf)   
#compareCRS(hf,get_crs(J63trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J63/UTM_b063_Neotropic_Earthenv_dem90m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J63/UTM_b063_Neotropic_Earthenv_dem90m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J64
#' b064_Neotropic_Earthenv_dem90m_tif_exp.
hf <-"D:/GISbases/b064JaguarAlbersBuffer70/b064_Neotropic_Earthenv_dem90m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J64trk))
#crs(hf)   
#compareCRS(hf,get_crs(J64trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J64/UTM_b064_Neotropic_Earthenv_dem90m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J64/UTM_b064_Neotropic_Earthenv_dem90m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J65
#' b065_Neotropic_Earthenv_dem90m_tif_exp.
hf <-"D:/GISbases/b065JaguarAlbersBuffer70/b065_Neotropic_Earthenv_dem90m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J65trk))
#crs(hf)   
#compareCRS(hf,get_crs(J65trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J65/UTM_b065_Neotropic_Earthenv_dem90m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J65/UTM_b065_Neotropic_Earthenv_dem90m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J66
#' b066_Neotropic_Earthenv_dem90m_tif_exp.
hf <-"D:/GISbases/b066JaguarAlbersBuffer70/b066_Neotropic_Earthenv_dem90m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J66trk))
#crs(hf)   
#compareCRS(hf,get_crs(J66trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J66/UTM_b066_Neotropic_Earthenv_dem90m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J66/UTM_b066_Neotropic_Earthenv_dem90m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J67
#' b067_Neotropic_Earthenv_dem90m_tif_exp.
hf <-"D:/GISbases/b067JaguarAlbersBuffer70/b067_Neotropic_Earthenv_dem90m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J67trk))
#crs(hf)   
#compareCRS(hf,get_crs(J67trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J67/UTM_b067_Neotropic_Earthenv_dem90m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J67/UTM_b067_Neotropic_Earthenv_dem90m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J68
#' b068_Neotropic_Earthenv_dem90m_tif_exp.
hf <-"D:/GISbases/b068JaguarAlbersBuffer70/b068_Neotropic_Earthenv_dem90m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J68trk))
#crs(hf)   
#compareCRS(hf,get_crs(J68trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J68/UTM_b068_Neotropic_Earthenv_dem90m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J68/UTM_b068_Neotropic_Earthenv_dem90m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J69
#' b069_Neotropic_Earthenv_dem90m_tif_exp.
hf <-"D:/GISbases/b069JaguarAlbersBuffer70/b069_Neotropic_Earthenv_dem90m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J69trk))
#crs(hf)   
#compareCRS(hf,get_crs(J69trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J69/UTM_b069_Neotropic_Earthenv_dem90m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J69/UTM_b069_Neotropic_Earthenv_dem90m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J70
#' b070_Neotropic_Earthenv_dem90m_tif_exp.
hf <-"D:/GISbases/b070JaguarAlbersBuffer70/b070_Neotropic_Earthenv_dem90m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J70trk))
#crs(hf)   
#compareCRS(hf,get_crs(J70trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J70/UTM_b070_Neotropic_Earthenv_dem90m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J70/UTM_b070_Neotropic_Earthenv_dem90m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J71
#' b071_Neotropic_Earthenv_dem90m_tif_exp.
hf <-"D:/GISbases/b071JaguarAlbersBuffer70/b071_Neotropic_Earthenv_dem90m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J71trk))
#crs(hf)   
#compareCRS(hf,get_crs(J71trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J71/UTM_b071_Neotropic_Earthenv_dem90m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J71/UTM_b071_Neotropic_Earthenv_dem90m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J72
#' b072_Neotropic_Earthenv_dem90m_tif_exp.
hf <-"D:/GISbases/b072JaguarAlbersBuffer70/b072_Neotropic_Earthenv_dem90m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J72trk))
#crs(hf)   
#compareCRS(hf,get_crs(J72trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J72/UTM_b072_Neotropic_Earthenv_dem90m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J72/UTM_b072_Neotropic_Earthenv_dem90m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J73
#' b073_Neotropic_Earthenv_dem90m_tif_exp.
hf <-"D:/GISbases/b073JaguarAlbersBuffer70/b073_Neotropic_Earthenv_dem90m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J73trk))
#crs(hf)   
#compareCRS(hf,get_crs(J73trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J73/UTM_b073_Neotropic_Earthenv_dem90m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J73/UTM_b073_Neotropic_Earthenv_dem90m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J74
#' b074_Neotropic_Earthenv_dem90m_tif_exp.
hf <-"D:/GISbases/b074JaguarAlbersBuffer70/b074_Neotropic_Earthenv_dem90m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J74trk))
#crs(hf)   
#compareCRS(hf,get_crs(J74trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J74/UTM_b074_Neotropic_Earthenv_dem90m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J74/UTM_b074_Neotropic_Earthenv_dem90m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J75
#' b075_Neotropic_Earthenv_dem90m_tif_exp.
hf <-"D:/GISbases/b075JaguarAlbersBuffer70/b075_Neotropic_Earthenv_dem90m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J75trk))
#crs(hf)   
#compareCRS(hf,get_crs(J75trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J75/UTM_b075_Neotropic_Earthenv_dem90m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J75/UTM_b075_Neotropic_Earthenv_dem90m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J76
#' b076_Neotropic_Earthenv_dem90m_tif_exp.
hf <-"D:/GISbases/b076JaguarAlbersBuffer70/b076_Neotropic_Earthenv_dem90m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J76trk))
#crs(hf)   
#compareCRS(hf,get_crs(J76trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J76/UTM_b076_Neotropic_Earthenv_dem90m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J76/UTM_b076_Neotropic_Earthenv_dem90m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J77
#' b077_Neotropic_Earthenv_dem90m_tif_exp.
hf <-"D:/GISbases/b077JaguarAlbersBuffer70/b077_Neotropic_Earthenv_dem90m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J77trk))
#crs(hf)   
#compareCRS(hf,get_crs(J77trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J77/UTM_b077_Neotropic_Earthenv_dem90m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J77/UTM_b077_Neotropic_Earthenv_dem90m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J78
#' b078_Neotropic_Earthenv_dem90m_tif_exp.
hf <-"D:/GISbases/b078JaguarAlbersBuffer70/b078_Neotropic_Earthenv_dem90m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J78trk))
#crs(hf)   
#compareCRS(hf,get_crs(J78trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J78/UTM_b078_Neotropic_Earthenv_dem90m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J78/UTM_b078_Neotropic_Earthenv_dem90m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J79
#' b079_Neotropic_Earthenv_dem90m_tif_exp.
hf <-"D:/GISbases/b079JaguarAlbersBuffer70/b079_Neotropic_Earthenv_dem90m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J79trk))
#crs(hf)   
#compareCRS(hf,get_crs(J79trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J79/UTM_b079_Neotropic_Earthenv_dem90m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J79/UTM_b079_Neotropic_Earthenv_dem90m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J80
#' b080_Neotropic_Earthenv_dem90m_tif_exp.
hf <-"D:/GISbases/b080JaguarAlbersBuffer70/b080_Neotropic_Earthenv_dem90m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J80trk))
#crs(hf)   
#compareCRS(hf,get_crs(J80trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J80/UTM_b080_Neotropic_Earthenv_dem90m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J80/UTM_b080_Neotropic_Earthenv_dem90m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J81
#' b081_Neotropic_Earthenv_dem90m_tif_exp.
hf <-"D:/GISbases/b081JaguarAlbersBuffer70/b081_Neotropic_Earthenv_dem90m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J81trk))
#crs(hf)   
#compareCRS(hf,get_crs(J81trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J81/UTM_b081_Neotropic_Earthenv_dem90m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J81/UTM_b081_Neotropic_Earthenv_dem90m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J82
#' b082_Neotropic_Earthenv_dem90m_tif_exp.
hf <-"D:/GISbases/b082JaguarAlbersBuffer70/b082_Neotropic_Earthenv_dem90m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J82trk))
#crs(hf)   
#compareCRS(hf,get_crs(J82trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J82/UTM_b082_Neotropic_Earthenv_dem90m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J82/UTM_b082_Neotropic_Earthenv_dem90m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J83
#' b083_Neotropic_Earthenv_dem90m_tif_exp.
hf <-"D:/GISbases/b083JaguarAlbersBuffer70/b083_Neotropic_Earthenv_dem90m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J83trk))
#crs(hf)   
#compareCRS(hf,get_crs(J83trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J83/UTM_b083_Neotropic_Earthenv_dem90m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J83/UTM_b083_Neotropic_Earthenv_dem90m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J84
#' b084_Neotropic_Earthenv_dem90m_tif_exp.
hf <-"D:/GISbases/b084JaguarAlbersBuffer70/b084_Neotropic_Earthenv_dem90m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J84trk))
#crs(hf)   
#compareCRS(hf,get_crs(J84trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J84/UTM_b084_Neotropic_Earthenv_dem90m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J84/UTM_b084_Neotropic_Earthenv_dem90m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J85
#' b085_Neotropic_Earthenv_dem90m_tif_exp.
hf <-"D:/GISbases/b085JaguarAlbersBuffer70/b085_Neotropic_Earthenv_dem90m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J85trk))
#crs(hf)   
#compareCRS(hf,get_crs(J85trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J85/UTM_b085_Neotropic_Earthenv_dem90m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J85/UTM_b085_Neotropic_Earthenv_dem90m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J86
#' b086_Neotropic_Earthenv_dem90m_tif_exp.
hf <-"D:/GISbases/b086JaguarAlbersBuffer70/b086_Neotropic_Earthenv_dem90m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J86trk))
#crs(hf)   
#compareCRS(hf,get_crs(J86trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J86/UTM_b086_Neotropic_Earthenv_dem90m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J86/UTM_b086_Neotropic_Earthenv_dem90m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J87
#' b087_Neotropic_Earthenv_dem90m_tif_exp.
hf <-"D:/GISbases/b087JaguarAlbersBuffer70/b087_Neotropic_Earthenv_dem90m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J87trk))
#crs(hf)   
#compareCRS(hf,get_crs(J87trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J87/UTM_b087_Neotropic_Earthenv_dem90m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J87/UTM_b087_Neotropic_Earthenv_dem90m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J88
#' b088_Neotropic_Earthenv_dem90m_tif_exp.
hf <-"D:/GISbases/b088JaguarAlbersBuffer70/b088_Neotropic_Earthenv_dem90m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J88trk))
#crs(hf)   
#compareCRS(hf,get_crs(J88trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J88/UTM_b088_Neotropic_Earthenv_dem90m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J88/UTM_b088_Neotropic_Earthenv_dem90m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J89
#' b089_Neotropic_Earthenv_dem90m_tif_exp.
hf <-"D:/GISbases/b089JaguarAlbersBuffer70/b089_Neotropic_Earthenv_dem90m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J89trk))
#crs(hf)   
#compareCRS(hf,get_crs(J89trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J89/UTM_b089_Neotropic_Earthenv_dem90m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J89/UTM_b089_Neotropic_Earthenv_dem90m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J90
#' b090_Neotropic_Earthenv_dem90m_tif_exp.
hf <-"D:/GISbases/b090JaguarAlbersBuffer70/b090_Neotropic_Earthenv_dem90m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J90trk))
#crs(hf)   
#compareCRS(hf,get_crs(J90trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J90/UTM_b090_Neotropic_Earthenv_dem90m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J90/UTM_b090_Neotropic_Earthenv_dem90m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J91
#' b091_Neotropic_Earthenv_dem90m_tif_exp.
hf <-"D:/GISbases/b091JaguarAlbersBuffer70/b091_Neotropic_Earthenv_dem90m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J91trk))
#crs(hf)   
#compareCRS(hf,get_crs(J91trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J91/UTM_b091_Neotropic_Earthenv_dem90m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J91/UTM_b091_Neotropic_Earthenv_dem90m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J92
#' b092_Neotropic_Earthenv_dem90m_tif_exp.
hf <-"D:/GISbases/b092JaguarAlbersBuffer70/b092_Neotropic_Earthenv_dem90m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J92trk))
#crs(hf)   
#compareCRS(hf,get_crs(J92trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J92/UTM_b092_Neotropic_Earthenv_dem90m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J92/UTM_b092_Neotropic_Earthenv_dem90m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J93
#' b093_Neotropic_Earthenv_dem90m_tif_exp.
hf <-"D:/GISbases/b093JaguarAlbersBuffer70/b093_Neotropic_Earthenv_dem90m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J93trk))
#crs(hf)   
#compareCRS(hf,get_crs(J93trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J93/UTM_b093_Neotropic_Earthenv_dem90m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J93/UTM_b093_Neotropic_Earthenv_dem90m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J94
#' b094_Neotropic_Earthenv_dem90m_tif_exp.
hf <-"D:/GISbases/b094JaguarAlbersBuffer70/b094_Neotropic_Earthenv_dem90m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J94trk))
#crs(hf)   
#compareCRS(hf,get_crs(J94trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J94/UTM_b094_Neotropic_Earthenv_dem90m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J94/UTM_b094_Neotropic_Earthenv_dem90m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J95
#' b095_Neotropic_Earthenv_dem90m_tif_exp.
hf <-"D:/GISbases/b095JaguarAlbersBuffer70/b095_Neotropic_Earthenv_dem90m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J95trk))
#crs(hf)   
#compareCRS(hf,get_crs(J95trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J95/UTM_b095_Neotropic_Earthenv_dem90m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J95/UTM_b095_Neotropic_Earthenv_dem90m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J96
#' b096_Neotropic_Earthenv_dem90m_tif_exp.
hf <-"D:/GISbases/b096JaguarAlbersBuffer70/b096_Neotropic_Earthenv_dem90m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J96trk))
#crs(hf)   
#compareCRS(hf,get_crs(J96trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J96/UTM_b096_Neotropic_Earthenv_dem90m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J96/UTM_b096_Neotropic_Earthenv_dem90m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J97
#' b097_Neotropic_Earthenv_dem90m_tif_exp.
hf <-"D:/GISbases/b097JaguarAlbersBuffer70/b097_Neotropic_Earthenv_dem90m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J97trk))
#crs(hf)   
#compareCRS(hf,get_crs(J97trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J97/UTM_b097_Neotropic_Earthenv_dem90m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J97/UTM_b097_Neotropic_Earthenv_dem90m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J98
#' b098_Neotropic_Earthenv_dem90m_tif_exp.
hf <-"D:/GISbases/b098JaguarAlbersBuffer70/b098_Neotropic_Earthenv_dem90m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J98trk))
#crs(hf)   
#compareCRS(hf,get_crs(J98trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J98/UTM_b098_Neotropic_Earthenv_dem90m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J98/UTM_b098_Neotropic_Earthenv_dem90m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J99
#' b099_Neotropic_Earthenv_dem90m_tif_exp.
hf <-"D:/GISbases/b099JaguarAlbersBuffer70/b099_Neotropic_Earthenv_dem90m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J99trk))
#crs(hf)   
#compareCRS(hf,get_crs(J99trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J99/UTM_b099_Neotropic_Earthenv_dem90m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J99/UTM_b099_Neotropic_Earthenv_dem90m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J100
#' b100_Neotropic_Earthenv_dem90m_tif_exp.
hf <-"D:/GISbases/b100JaguarAlbersBuffer70/b100_Neotropic_Earthenv_dem90m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J100trk))
#crs(hf)   
#compareCRS(hf,get_crs(J100trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J100/UTM_b100_Neotropic_Earthenv_dem90m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J100/UTM_b010_Neotropic_Earthenv_dem90m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J101
#' b101_Neotropic_Earthenv_dem90m_tif_exp.
hf <-"D:/GISbases/b101JaguarAlbersBuffer70/b101_Neotropic_Earthenv_dem90m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J101trk))
#crs(hf)   
#compareCRS(hf,get_crs(J101trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J101/UTM_b101_Neotropic_Earthenv_dem90m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J101/UTM_b101_Neotropic_Earthenv_dem90m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J102
#' b102_Neotropic_Earthenv_dem90m_tif_exp.
hf <-"D:/GISbases/b102JaguarAlbersBuffer70/b102_Neotropic_Earthenv_dem90m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J102trk))
#crs(hf)   
#compareCRS(hf,get_crs(J102trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J102/UTM_b102_Neotropic_Earthenv_dem90m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J102/UTM_b102_Neotropic_Earthenv_dem90m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J103
#' b103_Neotropic_Earthenv_dem90m_tif_exp.
hf <-"D:/GISbases/b103JaguarAlbersBuffer70/b103_Neotropic_Earthenv_dem90m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J103trk))
#crs(hf)   
#compareCRS(hf,get_crs(J103trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J103/UTM_b103_Neotropic_Earthenv_dem90m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J103/UTM_b103_Neotropic_Earthenv_dem90m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J104
#' b104_Neotropic_Earthenv_dem90m_tif_exp.
hf <-"D:/GISbases/b104JaguarAlbersBuffer70/b104_Neotropic_Earthenv_dem90m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J104trk))
#crs(hf)   
#compareCRS(hf,get_crs(J104trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J104/UTM_b104_Neotropic_Earthenv_dem90m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J104/UTM_b104_Neotropic_Earthenv_dem90m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J105
#' b105_Neotropic_Earthenv_dem90m_tif_exp.
hf <-"D:/GISbases/b105JaguarAlbersBuffer70/b105_Neotropic_Earthenv_dem90m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J105trk))
#crs(hf)   
#compareCRS(hf,get_crs(J105trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J105/UTM_b105_Neotropic_Earthenv_dem90m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J105/UTM_b105_Neotropic_Earthenv_dem90m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J106
#' b106_Neotropic_Earthenv_dem90m_tif_exp.
hf <-"D:/GISbases/b106JaguarAlbersBuffer70/b106_Neotropic_Earthenv_dem90m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J106trk))
#crs(hf)   
#compareCRS(hf,get_crs(J106trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J106/UTM_b106_Neotropic_Earthenv_dem90m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J106/UTM_b106_Neotropic_Earthenv_dem90m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J107
#' b107_Neotropic_Earthenv_dem90m_tif_exp.
hf <-"D:/GISbases/b107JaguarAlbersBuffer70/b107_Neotropic_Earthenv_dem90m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J107trk))
#crs(hf)   
#compareCRS(hf,get_crs(J107trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J107/UTM_b107_Neotropic_Earthenv_dem90m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J107/UTM_b107_Neotropic_Earthenv_dem90m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J108
#' b108_Neotropic_Earthenv_dem90m_tif_exp.
hf <-"D:/GISbases/b108JaguarAlbersBuffer70/b108_Neotropic_Earthenv_dem90m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J108trk))
#crs(hf)   
#compareCRS(hf,get_crs(J108trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J108/UTM_b108_Neotropic_Earthenv_dem90m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J108/UTM_b108_Neotropic_Earthenv_dem90m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J109
#' b109_Neotropic_Earthenv_dem90m_tif_exp.
hf <-"D:/GISbases/b109JaguarAlbersBuffer70/b109_Neotropic_Earthenv_dem90m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J109trk))
#crs(hf)   
#compareCRS(hf,get_crs(J109trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J109/UTM_b109_Neotropic_Earthenv_dem90m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J109/UTM_b109_Neotropic_Earthenv_dem90m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J110
#' b110_Neotropic_Earthenv_dem90m_tif_exp.
hf <-"D:/GISbases/b110JaguarAlbersBuffer70/b110_Neotropic_Earthenv_dem90m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J110trk))
#crs(hf)   
#compareCRS(hf,get_crs(J110trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J110/UTM_b110_Neotropic_Earthenv_dem90m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J110/UTM_b110_Neotropic_Earthenv_dem90m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J111
#' b111_Neotropic_Earthenv_dem90m_tif_exp.
hf <-"D:/GISbases/b111JaguarAlbersBuffer70/b111_Neotropic_Earthenv_dem90m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J111trk))
#crs(hf)   
#compareCRS(hf,get_crs(J111trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J111/UTM_b111_Neotropic_Earthenv_dem90m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J111/UTM_b111_Neotropic_Earthenv_dem90m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J112
#' b112_Neotropic_Earthenv_dem90m_tif_exp.
hf <-"D:/GISbases/b112JaguarAlbersBuffer70/b112_Neotropic_Earthenv_dem90m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J112trk))
#crs(hf)   
#compareCRS(hf,get_crs(J112trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J112/UTM_b112_Neotropic_Earthenv_dem90m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J112/UTM_b112_Neotropic_Earthenv_dem90m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J113
#' b113_Neotropic_Earthenv_dem90m_tif_exp.
hf <-"D:/GISbases/b113JaguarAlbersBuffer70/b113_Neotropic_Earthenv_dem90m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J113trk))
#crs(hf)   
#compareCRS(hf,get_crs(J113trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J113/UTM_b113_Neotropic_Earthenv_dem90m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J113/UTM_b113_Neotropic_Earthenv_dem90m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J114
#' b114_Neotropic_Earthenv_dem90m_tif_exp.
hf <-"D:/GISbases/b114JaguarAlbersBuffer70/b114_Neotropic_Earthenv_dem90m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J114trk))
#crs(hf)   
#compareCRS(hf,get_crs(J114trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J114/UTM_b114_Neotropic_Earthenv_dem90m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J114/UTM_b114_Neotropic_Earthenv_dem90m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J115
#' b115_Neotropic_Earthenv_dem90m_tif_exp.
hf <-"D:/GISbases/b115JaguarAlbersBuffer70/b115_Neotropic_Earthenv_dem90m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J115trk))
#crs(hf)   
#compareCRS(hf,get_crs(J115trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J115/UTM_b115_Neotropic_Earthenv_dem90m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J115/UTM_b115_Neotropic_Earthenv_dem90m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J116
#' b116_Neotropic_Earthenv_dem90m_tif_exp.
hf <-"D:/GISbases/b116JaguarAlbersBuffer70/b116_Neotropic_Earthenv_dem90m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J116trk))
#crs(hf)   
#compareCRS(hf,get_crs(J116trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J116/UTM_b116_Neotropic_Earthenv_dem90m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J116/UTM_b116_Neotropic_Earthenv_dem90m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J117
#' b117_Neotropic_Earthenv_dem90m_tif_exp.
hf <-"D:/GISbases/b117JaguarAlbersBuffer70/b117_Neotropic_Earthenv_dem90m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J117trk))
#crs(hf)   
#compareCRS(hf,get_crs(J117trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J117/UTM_b117_Neotropic_Earthenv_dem90m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J117/UTM_b117_Neotropic_Earthenv_dem90m_tif_exp.tif";(test=raster(test)); x11(); plot(test)


####################################################################################################################################
# water_frequency_2010_30m_tif_exp.tif
# water_frequency_2010_30m_tif_exp.

#' J1
#' water_frequency_2010_30m_tif_exp.tif
hf <-"D:/GISbases/b001JaguarAlbersBuffer70/b001_water_frequency_2010_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J1trk))
#crs(hf)   
#compareCRS(hf,get_crs(J1trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J1/UTM_b001_water_frequency_2010_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J1/UTM_b001_water_frequency_2010_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J2
#' b002_water_frequency_2010_30m_tif_exp.
hf <-"D:/GISbases/b002JaguarAlbersBuffer70/b002_water_frequency_2010_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J2trk))
#crs(hf)   
#compareCRS(hf,get_crs(J2trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J2/UTM_b002_water_frequency_2010_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J2/UTM_b002_water_frequency_2010_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J3
#' b003_water_frequency_2010_30m_tif_exp.
hf <-"D:/GISbases/b003JaguarAlbersBuffer70/b003_water_frequency_2010_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J3trk))
#crs(hf)   
#compareCRS(hf,get_crs(J3trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J3/UTM_b003_water_frequency_2010_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J3/UTM_b003_water_frequency_2010_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J4
#' b004_water_frequency_2010_30m_tif_exp.
hf <-"D:/GISbases/b004JaguarAlbersBuffer70/b004_water_frequency_2010_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J4trk))
#crs(hf)   
#compareCRS(hf,get_crs(J4trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J4/UTM_b004_water_frequency_2010_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J4/UTM_b004_water_frequency_2010_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J5
#' b005_water_frequency_2010_30m_tif_exp.
hf <-"D:/GISbases/b005JaguarAlbersBuffer70/b005_water_frequency_2010_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J5trk))
#crs(hf)   
#compareCRS(hf,get_crs(J5trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J5/UTM_b005_water_frequency_2010_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J5/UTM_b005_water_frequency_2010_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J6
#' b006_water_frequency_2010_30m_tif_exp.
hf <-"D:/GISbases/b006JaguarAlbersBuffer70/b006_water_frequency_2010_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J6trk))
#crs(hf)   
#compareCRS(hf,get_crs(J6trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J6/UTM_b006_water_frequency_2010_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J6/UTM_b006_water_frequency_2010_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J7
#' b007_water_frequency_2010_30m_tif_exp.
hf <-"D:/GISbases/b007JaguarAlbersBuffer70/b007_water_frequency_2010_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J7trk))
#crs(hf)   
#compareCRS(hf,get_crs(J7trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J7/UTM_b007_water_frequency_2010_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J7/UTM_b007_water_frequency_2010_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J8
#' b008_water_frequency_2010_30m_tif_exp.
hf <-"D:/GISbases/b008JaguarAlbersBuffer70/b008_water_frequency_2010_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J8trk))
#crs(hf)   
#compareCRS(hf,get_crs(J8trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J8/UTM_b008_water_frequency_2010_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J8/UTM_b008_water_frequency_2010_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J9
#' b009_water_frequency_2010_30m_tif_exp.
hf <-"D:/GISbases/b009JaguarAlbersBuffer70/b009_water_frequency_2010_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J9trk))
#crs(hf)   
#compareCRS(hf,get_crs(J9trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J9/UTM_b009_water_frequency_2010_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J9/UTM_b009_water_frequency_2010_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J10
#' b010_water_frequency_2010_30m_tif_exp.
hf <-"D:/GISbases/b010JaguarAlbersBuffer70/b010_water_frequency_2010_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J10trk))
#crs(hf)   
#compareCRS(hf,get_crs(J10trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J10/UTM_b010_water_frequency_2010_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J10/UTM_b010_water_frequency_2010_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J11
#' b011_water_frequency_2010_30m_tif_exp.
hf <-"D:/GISbases/b011JaguarAlbersBuffer70/b011_water_frequency_2010_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J11trk))
#crs(hf)   
#compareCRS(hf,get_crs(J11trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J11/UTM_b011_water_frequency_2010_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J11/UTM_b011_water_frequency_2010_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J12
#' b012_water_frequency_2010_30m_tif_exp.
hf <-"D:/GISbases/b012JaguarAlbersBuffer70/b012_water_frequency_2010_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J12trk))
#crs(hf)   
#compareCRS(hf,get_crs(J12trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J12/UTM_b012_water_frequency_2010_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J12/UTM_b012_water_frequency_2010_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J13
#' b013_water_frequency_2010_30m_tif_exp.
hf <-"D:/GISbases/b013JaguarAlbersBuffer70/b013_water_frequency_2010_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J13trk))
#crs(hf)   
#compareCRS(hf,get_crs(J13trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J13/UTM_b013_water_frequency_2010_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J13/UTM_b013_water_frequency_2010_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J14
#' b014_water_frequency_2010_30m_tif_exp.
hf <-"D:/GISbases/b014JaguarAlbersBuffer70/b014_water_frequency_2010_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J14trk))
#crs(hf)   
#compareCRS(hf,get_crs(J14trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J14/UTM_b014_water_frequency_2010_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J14/UTM_b014_water_frequency_2010_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J15
#' b015_water_frequency_2010_30m_tif_exp.
hf <-"D:/GISbases/b015JaguarAlbersBuffer70/b015_water_frequency_2010_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J15trk))
#crs(hf)   
#compareCRS(hf,get_crs(J15trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J15/UTM_b015_water_frequency_2010_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J15/UTM_b015_water_frequency_2010_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J16
#' b016_water_frequency_2010_30m_tif_exp.
hf <-"D:/GISbases/b016JaguarAlbersBuffer70/b016_water_frequency_2010_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J16trk))
#crs(hf)   
#compareCRS(hf,get_crs(J16trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J16/UTM_b016_water_frequency_2010_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J16/UTM_b016_water_frequency_2010_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J17
#' b017_water_frequency_2010_30m_tif_exp.
hf <-"D:/GISbases/b017JaguarAlbersBuffer70/b017_water_frequency_2010_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J17trk))
#crs(hf)   
#compareCRS(hf,get_crs(J17trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J17/UTM_b017_water_frequency_2010_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J17/UTM_b017_water_frequency_2010_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J18
#' b018_water_frequency_2010_30m_tif_exp.
hf <-"D:/GISbases/b018JaguarAlbersBuffer70/b018_water_frequency_2010_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J18trk))
#crs(hf)   
#compareCRS(hf,get_crs(J18trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J18/UTM_b018_water_frequency_2010_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J18/UTM_b018_water_frequency_2010_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J19
#' b019_water_frequency_2010_30m_tif_exp.
hf <-"D:/GISbases/b019JaguarAlbersBuffer70/b019_water_frequency_2010_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J19trk))
#crs(hf)   
#compareCRS(hf,get_crs(J19trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J19/UTM_b019_water_frequency_2010_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J19/UTM_b019_water_frequency_2010_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J20
#' b020_water_frequency_2010_30m_tif_exp.
hf <-"D:/GISbases/b020JaguarAlbersBuffer70/b020_water_frequency_2010_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J20trk))
#crs(hf)   
#compareCRS(hf,get_crs(J20trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J20/UTM_b020_water_frequency_2010_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J20/UTM_b020_water_frequency_2010_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J21
#' b021_water_frequency_2010_30m_tif_exp.
hf <-"D:/GISbases/b021JaguarAlbersBuffer70/b021_water_frequency_2010_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J21trk))
#crs(hf)   
#compareCRS(hf,get_crs(J21trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J21/UTM_b021_water_frequency_2010_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J21/UTM_b021_water_frequency_2010_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J22
#' b022_water_frequency_2010_30m_tif_exp.
hf <-"D:/GISbases/b022JaguarAlbersBuffer70/b022_water_frequency_2010_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J22trk))
#crs(hf)   
#compareCRS(hf,get_crs(J22trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J22/UTM_b022_water_frequency_2010_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J22/UTM_b022_water_frequency_2010_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J23
#' b023_water_frequency_2010_30m_tif_exp.
hf <-"D:/GISbases/b023JaguarAlbersBuffer70/b023_water_frequency_2010_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J23trk))
#crs(hf)   
#compareCRS(hf,get_crs(J23trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J23/UTM_b023_water_frequency_2010_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J23/UTM_b023_water_frequency_2010_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J24
#' b024_water_frequency_2010_30m_tif_exp.
hf <-"D:/GISbases/b024JaguarAlbersBuffer70/b024_water_frequency_2010_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J24trk))
#crs(hf)   
#compareCRS(hf,get_crs(J24trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J24/UTM_b024_water_frequency_2010_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J24/UTM_b024_water_frequency_2010_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J25
#' b025_water_frequency_2010_30m_tif_exp.
hf <-"D:/GISbases/b025JaguarAlbersBuffer70/b025_water_frequency_2010_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J25trk))
#crs(hf)   
#compareCRS(hf,get_crs(J25trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J25/UTM_b025_water_frequency_2010_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J25/UTM_b025_water_frequency_2010_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J26
#' b026_water_frequency_2010_30m_tif_exp.
hf <-"D:/GISbases/b026JaguarAlbersBuffer70/b026_water_frequency_2010_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J26trk))
#crs(hf)   
#compareCRS(hf,get_crs(J26trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J26/UTM_b026_water_frequency_2010_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J26/UTM_b026_water_frequency_2010_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J27
#' b027_water_frequency_2010_30m_tif_exp.
hf <-"D:/GISbases/b027JaguarAlbersBuffer70/b027_water_frequency_2010_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J27trk))
#crs(hf)   
#compareCRS(hf,get_crs(J27trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J27/UTM_b027_water_frequency_2010_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J27/UTM_b027_water_frequency_2010_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J28
#' b028_water_frequency_2010_30m_tif_exp.
hf <-"D:/GISbases/b028JaguarAlbersBuffer70/b028_water_frequency_2010_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J28trk))
#crs(hf)   
#compareCRS(hf,get_crs(J28trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J28/UTM_b028_water_frequency_2010_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J28/UTM_b028_water_frequency_2010_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J29
#' b029_water_frequency_2010_30m_tif_exp.
hf <-"D:/GISbases/b029JaguarAlbersBuffer70/b029_water_frequency_2010_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J29trk))
#crs(hf)   
#compareCRS(hf,get_crs(J29trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J29/UTM_b029_water_frequency_2010_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J29/UTM_b029_water_frequency_2010_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J30
#' b030_water_frequency_2010_30m_tif_exp.
hf <-"D:/GISbases/b030JaguarAlbersBuffer70/b030_water_frequency_2010_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J30trk))
#crs(hf)   
#compareCRS(hf,get_crs(J30trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J30/UTM_b030_water_frequency_2010_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J30/UTM_b030_water_frequency_2010_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J30
#' b030_water_frequency_2010_30m_tif_exp.
hf <-"D:/GISbases/b030JaguarAlbersBuffer70/b030_water_frequency_2010_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J30trk))
#crs(hf)   
#compareCRS(hf,get_crs(J30trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J30/UTM_b030_water_frequency_2010_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J30/UTM_b030_water_frequency_2010_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J31
#' b031_water_frequency_2010_30m_tif_exp.
hf <-"D:/GISbases/b031JaguarAlbersBuffer70/b031_water_frequency_2010_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J31trk))
#crs(hf)   
#compareCRS(hf,get_crs(J31trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J31/UTM_b031_water_frequency_2010_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J31/UTM_b031_water_frequency_2010_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J32
#' b032_water_frequency_2010_30m_tif_exp.
hf <-"D:/GISbases/b032JaguarAlbersBuffer70/b032_water_frequency_2010_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J32trk))
#crs(hf)   
#compareCRS(hf,get_crs(J32trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J32/UTM_b032_water_frequency_2010_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J32/UTM_b032_water_frequency_2010_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J33
#' b033_water_frequency_2010_30m_tif_exp.
hf <-"D:/GISbases/b033JaguarAlbersBuffer70/b033_water_frequency_2010_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J33trk))
#crs(hf)   
#compareCRS(hf,get_crs(J33trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J33/UTM_b033_water_frequency_2010_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J33/UTM_b033_water_frequency_2010_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J34
#' b034_water_frequency_2010_30m_tif_exp.
hf <-"D:/GISbases/b034JaguarAlbersBuffer70/b034_water_frequency_2010_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J34trk))
#crs(hf)   
#compareCRS(hf,get_crs(J34trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J34/UTM_b034_water_frequency_2010_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J34/UTM_b034_water_frequency_2010_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J35
#' b035_water_frequency_2010_30m_tif_exp.
hf <-"D:/GISbases/b035JaguarAlbersBuffer70/b035_water_frequency_2010_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J35trk))
#crs(hf)   
#compareCRS(hf,get_crs(J35trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J35/UTM_b035_water_frequency_2010_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J35/UTM_b035_water_frequency_2010_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J36
#' b036_water_frequency_2010_30m_tif_exp.
hf <-"D:/GISbases/b036JaguarAlbersBuffer70/b036_water_frequency_2010_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J36trk))
#crs(hf)   
#compareCRS(hf,get_crs(J36trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J36/UTM_b036_water_frequency_2010_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J36/UTM_b036_water_frequency_2010_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J37
#' b037_water_frequency_2010_30m_tif_exp.
hf <-"D:/GISbases/b037JaguarAlbersBuffer70/b037_water_frequency_2010_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J37trk))
#crs(hf)   
#compareCRS(hf,get_crs(J37trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J37/UTM_b037_water_frequency_2010_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J37/UTM_b037_water_frequency_2010_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J38
#' b038_water_frequency_2010_30m_tif_exp.
hf <-"D:/GISbases/b038JaguarAlbersBuffer70/b038_water_frequency_2010_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J38trk))
#crs(hf)   
#compareCRS(hf,get_crs(J38trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J38/UTM_b038_water_frequency_2010_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J38/UTM_b038_water_frequency_2010_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J39
#' b039_water_frequency_2010_30m_tif_exp.
hf <-"D:/GISbases/b039JaguarAlbersBuffer70/b039_water_frequency_2010_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J39trk))
#crs(hf)   
#compareCRS(hf,get_crs(J39trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J39/UTM_b039_water_frequency_2010_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J39/UTM_b039_water_frequency_2010_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J40
#' b040_water_frequency_2010_30m_tif_exp.
hf <-"D:/GISbases/b040JaguarAlbersBuffer70/b040_water_frequency_2010_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J40trk))
#crs(hf)   
#compareCRS(hf,get_crs(J40trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J40/UTM_b040_water_frequency_2010_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J40/UTM_b040_water_frequency_2010_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J41
#' b041_water_frequency_2010_30m_tif_exp.
hf <-"D:/GISbases/b041JaguarAlbersBuffer70/b041_water_frequency_2010_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J41trk))
#crs(hf)   
#compareCRS(hf,get_crs(J41trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J41/UTM_b041_water_frequency_2010_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J41/UTM_b041_water_frequency_2010_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J42
#' b042_water_frequency_2010_30m_tif_exp.
hf <-"D:/GISbases/b042JaguarAlbersBuffer70/b042_water_frequency_2010_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J42trk))
#crs(hf)   
#compareCRS(hf,get_crs(J42trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J42/UTM_b042_water_frequency_2010_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J42/UTM_b042_water_frequency_2010_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J43
#' b043_water_frequency_2010_30m_tif_exp.
hf <-"D:/GISbases/b043JaguarAlbersBuffer70/b043_water_frequency_2010_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J43trk))
#crs(hf)   
#compareCRS(hf,get_crs(J43trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J43/UTM_b043_water_frequency_2010_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J43/UTM_b043_water_frequency_2010_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J44
#' b044_water_frequency_2010_30m_tif_exp.
hf <-"D:/GISbases/b044JaguarAlbersBuffer70/b044_water_frequency_2010_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J44trk))
#crs(hf)   
#compareCRS(hf,get_crs(J44trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J44/UTM_b044_water_frequency_2010_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J44/UTM_b044_water_frequency_2010_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J45
#' b045_water_frequency_2010_30m_tif_exp.
hf <-"D:/GISbases/b045JaguarAlbersBuffer70/b045_water_frequency_2010_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J45trk))
#crs(hf)   
#compareCRS(hf,get_crs(J45trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J45/UTM_b045_water_frequency_2010_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J45/UTM_b045_water_frequency_2010_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J46
#' b046_water_frequency_2010_30m_tif_exp.
hf <-"D:/GISbases/b046JaguarAlbersBuffer70/b046_water_frequency_2010_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J46trk))
#crs(hf)   
#compareCRS(hf,get_crs(J46trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J46/UTM_b046_water_frequency_2010_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J46/UTM_b046_water_frequency_2010_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J47
#' b047_water_frequency_2010_30m_tif_exp.
hf <-"D:/GISbases/b047JaguarAlbersBuffer70/b047_water_frequency_2010_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J47trk))
#crs(hf)   
#compareCRS(hf,get_crs(J47trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J47/UTM_b047_water_frequency_2010_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J47/UTM_b047_water_frequency_2010_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J48
#' b048_water_frequency_2010_30m_tif_exp.
hf <-"D:/GISbases/b048JaguarAlbersBuffer70/b048_water_frequency_2010_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J48trk))
#crs(hf)   
#compareCRS(hf,get_crs(J48trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J48/UTM_b048_water_frequency_2010_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J48/UTM_b048_water_frequency_2010_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J49
#' b049_water_frequency_2010_30m_tif_exp.
hf <-"D:/GISbases/b049JaguarAlbersBuffer70/b049_water_frequency_2010_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J49trk))
#crs(hf)   
#compareCRS(hf,get_crs(J49trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J49/UTM_b049_water_frequency_2010_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J49/UTM_b049_water_frequency_2010_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J50
#' b050_water_frequency_2010_30m_tif_exp.
hf <-"D:/GISbases/b050JaguarAlbersBuffer70/b050_water_frequency_2010_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J50trk))
#crs(hf)   
#compareCRS(hf,get_crs(J50trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J50/UTM_b050_water_frequency_2010_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J50/UTM_b050_water_frequency_2010_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J51
#' b051_water_frequency_2010_30m_tif_exp.
hf <-"D:/GISbases/b051JaguarAlbersBuffer70/b051_water_frequency_2010_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J51trk))
#crs(hf)   
#compareCRS(hf,get_crs(J51trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J51/UTM_b051_water_frequency_2010_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J51/UTM_b051_water_frequency_2010_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J52
#' b052_water_frequency_2010_30m_tif_exp.
hf <-"D:/GISbases/b052JaguarAlbersBuffer70/b052_water_frequency_2010_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J52trk))
#crs(hf)   
#compareCRS(hf,get_crs(J52trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J52/UTM_b052_water_frequency_2010_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J52/UTM_b052_water_frequency_2010_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J53
#' b053_water_frequency_2010_30m_tif_exp.
hf <-"D:/GISbases/b053JaguarAlbersBuffer70/b053_water_frequency_2010_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J53trk))
#crs(hf)   
#compareCRS(hf,get_crs(J53trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J53/UTM_b053_water_frequency_2010_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J53/UTM_b053_water_frequency_2010_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J54
#' b054_water_frequency_2010_30m_tif_exp.
hf <-"D:/GISbases/b054JaguarAlbersBuffer70/b054_water_frequency_2010_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J54trk))
#crs(hf)   
#compareCRS(hf,get_crs(J54trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J54/UTM_b054_water_frequency_2010_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J54/UTM_b054_water_frequency_2010_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J55
#' b055_water_frequency_2010_30m_tif_exp.
hf <-"D:/GISbases/b055JaguarAlbersBuffer70/b055_water_frequency_2010_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J55trk))
#crs(hf)   
#compareCRS(hf,get_crs(J55trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J55/UTM_b055_water_frequency_2010_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J55/UTM_b055_water_frequency_2010_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J56
#' b056_water_frequency_2010_30m_tif_exp.
hf <-"D:/GISbases/b056JaguarAlbersBuffer70/b056_water_frequency_2010_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J56trk))
#crs(hf)   
#compareCRS(hf,get_crs(J56trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J56/UTM_b056_water_frequency_2010_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J56/UTM_b056_water_frequency_2010_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J57
#' b057_water_frequency_2010_30m_tif_exp.
hf <-"D:/GISbases/b057JaguarAlbersBuffer70/b057_water_frequency_2010_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J57trk))
#crs(hf)   
#compareCRS(hf,get_crs(J57trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J57/UTM_b057_water_frequency_2010_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J57/UTM_b057_water_frequency_2010_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J58
#' b058_water_frequency_2010_30m_tif_exp.
hf <-"D:/GISbases/b058JaguarAlbersBuffer70/b058_water_frequency_2010_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J58trk))
#crs(hf)   
#compareCRS(hf,get_crs(J58trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J58/UTM_b058_water_frequency_2010_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J58/UTM_b058_water_frequency_2010_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J59
#' b059_water_frequency_2010_30m_tif_exp.
hf <-"D:/GISbases/b059JaguarAlbersBuffer70/b059_water_frequency_2010_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J59trk))
#crs(hf)   
#compareCRS(hf,get_crs(J59trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J59/UTM_b059_water_frequency_2010_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J59/UTM_b059_water_frequency_2010_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J60
#' b060_water_frequency_2010_30m_tif_exp.
hf <-"D:/GISbases/b060JaguarAlbersBuffer70/b060_water_frequency_2010_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J60trk))
#crs(hf)   
#compareCRS(hf,get_crs(J60trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J60/UTM_b060_water_frequency_2010_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J60/UTM_b060_water_frequency_2010_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J61
#' b061_water_frequency_2010_30m_tif_exp.
hf <-"D:/GISbases/b061JaguarAlbersBuffer70/b061_water_frequency_2010_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J61trk))
#crs(hf)   
#compareCRS(hf,get_crs(J61trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J61/UTM_b061_water_frequency_2010_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J61/UTM_b061_water_frequency_2010_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J62
#' b062_water_frequency_2010_30m_tif_exp.
hf <-"D:/GISbases/b062JaguarAlbersBuffer70/b062_water_frequency_2010_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J62trk))
#crs(hf)   
#compareCRS(hf,get_crs(J62trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J62/UTM_b062_water_frequency_2010_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J62/UTM_b062_water_frequency_2010_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J63
#' b063_water_frequency_2010_30m_tif_exp.
hf <-"D:/GISbases/b063JaguarAlbersBuffer70/b063_water_frequency_2010_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J63trk))
#crs(hf)   
#compareCRS(hf,get_crs(J63trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J63/UTM_b063_water_frequency_2010_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J63/UTM_b063_water_frequency_2010_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J64
#' b064_water_frequency_2010_30m_tif_exp.
hf <-"D:/GISbases/b064JaguarAlbersBuffer70/b064_water_frequency_2010_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J64trk))
#crs(hf)   
#compareCRS(hf,get_crs(J64trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J64/UTM_b064_water_frequency_2010_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J64/UTM_b064_water_frequency_2010_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J65
#' b065_water_frequency_2010_30m_tif_exp.
hf <-"D:/GISbases/b065JaguarAlbersBuffer70/b065_water_frequency_2010_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J65trk))
#crs(hf)   
#compareCRS(hf,get_crs(J65trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J65/UTM_b065_water_frequency_2010_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J65/UTM_b065_water_frequency_2010_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J66
#' b066_water_frequency_2010_30m_tif_exp.
hf <-"D:/GISbases/b066JaguarAlbersBuffer70/b066_water_frequency_2010_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J66trk))
#crs(hf)   
#compareCRS(hf,get_crs(J66trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J66/UTM_b066_water_frequency_2010_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J66/UTM_b066_water_frequency_2010_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J67
#' b067_water_frequency_2010_30m_tif_exp.
hf <-"D:/GISbases/b067JaguarAlbersBuffer70/b067_water_frequency_2010_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J67trk))
#crs(hf)   
#compareCRS(hf,get_crs(J67trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J67/UTM_b067_water_frequency_2010_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J67/UTM_b067_water_frequency_2010_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J68
#' b068_water_frequency_2010_30m_tif_exp.
hf <-"D:/GISbases/b068JaguarAlbersBuffer70/b068_water_frequency_2010_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J68trk))
#crs(hf)   
#compareCRS(hf,get_crs(J68trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J68/UTM_b068_water_frequency_2010_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J68/UTM_b068_water_frequency_2010_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J69
#' b069_water_frequency_2010_30m_tif_exp.
hf <-"D:/GISbases/b069JaguarAlbersBuffer70/b069_water_frequency_2010_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J69trk))
#crs(hf)   
#compareCRS(hf,get_crs(J69trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J69/UTM_b069_water_frequency_2010_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J69/UTM_b069_water_frequency_2010_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J70
#' b070_water_frequency_2010_30m_tif_exp.
hf <-"D:/GISbases/b070JaguarAlbersBuffer70/b070_water_frequency_2010_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J70trk))
#crs(hf)   
#compareCRS(hf,get_crs(J70trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J70/UTM_b070_water_frequency_2010_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J70/UTM_b070_water_frequency_2010_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J71
#' b071_water_frequency_2010_30m_tif_exp.
hf <-"D:/GISbases/b071JaguarAlbersBuffer70/b071_water_frequency_2010_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J71trk))
#crs(hf)   
#compareCRS(hf,get_crs(J71trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J71/UTM_b071_water_frequency_2010_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J71/UTM_b071_water_frequency_2010_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J72
#' b072_water_frequency_2010_30m_tif_exp.
hf <-"D:/GISbases/b072JaguarAlbersBuffer70/b072_water_frequency_2010_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J72trk))
#crs(hf)   
#compareCRS(hf,get_crs(J72trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J72/UTM_b072_water_frequency_2010_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J72/UTM_b072_water_frequency_2010_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J73
#' b073_water_frequency_2010_30m_tif_exp.
hf <-"D:/GISbases/b073JaguarAlbersBuffer70/b073_water_frequency_2010_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J73trk))
#crs(hf)   
#compareCRS(hf,get_crs(J73trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J73/UTM_b073_water_frequency_2010_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J73/UTM_b073_water_frequency_2010_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J74
#' b074_water_frequency_2010_30m_tif_exp.
hf <-"D:/GISbases/b074JaguarAlbersBuffer70/b074_water_frequency_2010_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J74trk))
#crs(hf)   
#compareCRS(hf,get_crs(J74trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J74/UTM_b074_water_frequency_2010_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J74/UTM_b074_water_frequency_2010_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J75
#' b075_water_frequency_2010_30m_tif_exp.
hf <-"D:/GISbases/b075JaguarAlbersBuffer70/b075_water_frequency_2010_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J75trk))
#crs(hf)   
#compareCRS(hf,get_crs(J75trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J75/UTM_b075_water_frequency_2010_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J75/UTM_b075_water_frequency_2010_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J76
#' b076_water_frequency_2010_30m_tif_exp.
hf <-"D:/GISbases/b076JaguarAlbersBuffer70/b076_water_frequency_2010_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J76trk))
#crs(hf)   
#compareCRS(hf,get_crs(J76trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J76/UTM_b076_water_frequency_2010_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J76/UTM_b076_water_frequency_2010_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J77
#' b077_water_frequency_2010_30m_tif_exp.
hf <-"D:/GISbases/b077JaguarAlbersBuffer70/b077_water_frequency_2010_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J77trk))
#crs(hf)   
#compareCRS(hf,get_crs(J77trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J77/UTM_b077_water_frequency_2010_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J77/UTM_b077_water_frequency_2010_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J78
#' b078_water_frequency_2010_30m_tif_exp.
hf <-"D:/GISbases/b078JaguarAlbersBuffer70/b078_water_frequency_2010_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J78trk))
#crs(hf)   
#compareCRS(hf,get_crs(J78trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J78/UTM_b078_water_frequency_2010_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J78/UTM_b078_water_frequency_2010_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J79
#' b079_water_frequency_2010_30m_tif_exp.
hf <-"D:/GISbases/b079JaguarAlbersBuffer70/b079_water_frequency_2010_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J79trk))
#crs(hf)   
#compareCRS(hf,get_crs(J79trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J79/UTM_b079_water_frequency_2010_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J79/UTM_b079_water_frequency_2010_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J80
#' b080_water_frequency_2010_30m_tif_exp.
hf <-"D:/GISbases/b080JaguarAlbersBuffer70/b080_water_frequency_2010_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J80trk))
#crs(hf)   
#compareCRS(hf,get_crs(J80trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J80/UTM_b080_water_frequency_2010_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J80/UTM_b080_water_frequency_2010_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J81
#' b081_water_frequency_2010_30m_tif_exp.
hf <-"D:/GISbases/b081JaguarAlbersBuffer70/b081_water_frequency_2010_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J81trk))
#crs(hf)   
#compareCRS(hf,get_crs(J81trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J81/UTM_b081_water_frequency_2010_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J81/UTM_b081_water_frequency_2010_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J82
#' b082_water_frequency_2010_30m_tif_exp.
hf <-"D:/GISbases/b082JaguarAlbersBuffer70/b082_water_frequency_2010_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J82trk))
#crs(hf)   
#compareCRS(hf,get_crs(J82trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J82/UTM_b082_water_frequency_2010_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J82/UTM_b082_water_frequency_2010_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J83
#' b083_water_frequency_2010_30m_tif_exp.
hf <-"D:/GISbases/b083JaguarAlbersBuffer70/b083_water_frequency_2010_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J83trk))
#crs(hf)   
#compareCRS(hf,get_crs(J83trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J83/UTM_b083_water_frequency_2010_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J83/UTM_b083_water_frequency_2010_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J84
#' b084_water_frequency_2010_30m_tif_exp.
hf <-"D:/GISbases/b084JaguarAlbersBuffer70/b084_water_frequency_2010_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J84trk))
#crs(hf)   
#compareCRS(hf,get_crs(J84trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J84/UTM_b084_water_frequency_2010_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J84/UTM_b084_water_frequency_2010_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J85
#' b085_water_frequency_2010_30m_tif_exp.
hf <-"D:/GISbases/b085JaguarAlbersBuffer70/b085_water_frequency_2010_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J85trk))
#crs(hf)   
#compareCRS(hf,get_crs(J85trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J85/UTM_b085_water_frequency_2010_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J85/UTM_b085_water_frequency_2010_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J86
#' b086_water_frequency_2010_30m_tif_exp.
hf <-"D:/GISbases/b086JaguarAlbersBuffer70/b086_water_frequency_2010_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J86trk))
#crs(hf)   
#compareCRS(hf,get_crs(J86trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J86/UTM_b086_water_frequency_2010_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J86/UTM_b086_water_frequency_2010_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J87
#' b087_water_frequency_2010_30m_tif_exp.
hf <-"D:/GISbases/b087JaguarAlbersBuffer70/b087_water_frequency_2010_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J87trk))
#crs(hf)   
#compareCRS(hf,get_crs(J87trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J87/UTM_b087_water_frequency_2010_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J87/UTM_b087_water_frequency_2010_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J88
#' b088_water_frequency_2010_30m_tif_exp.
hf <-"D:/GISbases/b088JaguarAlbersBuffer70/b088_water_frequency_2010_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J88trk))
#crs(hf)   
#compareCRS(hf,get_crs(J88trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J88/UTM_b088_water_frequency_2010_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J88/UTM_b088_water_frequency_2010_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J89
#' b089_water_frequency_2010_30m_tif_exp.
hf <-"D:/GISbases/b089JaguarAlbersBuffer70/b089_water_frequency_2010_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J89trk))
#crs(hf)   
#compareCRS(hf,get_crs(J89trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J89/UTM_b089_water_frequency_2010_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J89/UTM_b089_water_frequency_2010_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J90
#' b090_water_frequency_2010_30m_tif_exp.
hf <-"D:/GISbases/b090JaguarAlbersBuffer70/b090_water_frequency_2010_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J90trk))
#crs(hf)   
#compareCRS(hf,get_crs(J90trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J90/UTM_b090_water_frequency_2010_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J90/UTM_b090_water_frequency_2010_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J91
#' b091_water_frequency_2010_30m_tif_exp.
hf <-"D:/GISbases/b091JaguarAlbersBuffer70/b091_water_frequency_2010_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J91trk))
#crs(hf)   
#compareCRS(hf,get_crs(J91trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J91/UTM_b091_water_frequency_2010_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J91/UTM_b091_water_frequency_2010_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J92
#' b092_water_frequency_2010_30m_tif_exp.
hf <-"D:/GISbases/b092JaguarAlbersBuffer70/b092_water_frequency_2010_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J92trk))
#crs(hf)   
#compareCRS(hf,get_crs(J92trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J92/UTM_b092_water_frequency_2010_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J92/UTM_b092_water_frequency_2010_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J93
#' b093_water_frequency_2010_30m_tif_exp.
hf <-"D:/GISbases/b093JaguarAlbersBuffer70/b093_water_frequency_2010_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J93trk))
#crs(hf)   
#compareCRS(hf,get_crs(J93trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J93/UTM_b093_water_frequency_2010_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J93/UTM_b093_water_frequency_2010_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J94
#' b094_water_frequency_2010_30m_tif_exp.
hf <-"D:/GISbases/b094JaguarAlbersBuffer70/b094_water_frequency_2010_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J94trk))
#crs(hf)   
#compareCRS(hf,get_crs(J94trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J94/UTM_b094_water_frequency_2010_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J94/UTM_b094_water_frequency_2010_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J95
#' b095_water_frequency_2010_30m_tif_exp.
hf <-"D:/GISbases/b095JaguarAlbersBuffer70/b095_water_frequency_2010_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J95trk))
#crs(hf)   
#compareCRS(hf,get_crs(J95trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J95/UTM_b095_water_frequency_2010_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J95/UTM_b095_water_frequency_2010_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J96
#' b096_water_frequency_2010_30m_tif_exp.
hf <-"D:/GISbases/b096JaguarAlbersBuffer70/b096_water_frequency_2010_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J96trk))
#crs(hf)   
#compareCRS(hf,get_crs(J96trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J96/UTM_b096_water_frequency_2010_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J96/UTM_b096_water_frequency_2010_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J97
#' b097_water_frequency_2010_30m_tif_exp.
hf <-"D:/GISbases/b097JaguarAlbersBuffer70/b097_water_frequency_2010_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J97trk))
#crs(hf)   
#compareCRS(hf,get_crs(J97trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J97/UTM_b097_water_frequency_2010_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J97/UTM_b097_water_frequency_2010_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J98
#' b098_water_frequency_2010_30m_tif_exp.
hf <-"D:/GISbases/b098JaguarAlbersBuffer70/b098_water_frequency_2010_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J98trk))
#crs(hf)   
#compareCRS(hf,get_crs(J98trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J98/UTM_b098_water_frequency_2010_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J98/UTM_b098_water_frequency_2010_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J99
#' b099_water_frequency_2010_30m_tif_exp.
hf <-"D:/GISbases/b099JaguarAlbersBuffer70/b099_water_frequency_2010_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J99trk))
#crs(hf)   
#compareCRS(hf,get_crs(J99trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J99/UTM_b099_water_frequency_2010_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J99/UTM_b099_water_frequency_2010_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J100
#' b100_water_frequency_2010_30m_tif_exp.
hf <-"D:/GISbases/b100JaguarAlbersBuffer70/b100_water_frequency_2010_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J100trk))
#crs(hf)   
#compareCRS(hf,get_crs(J100trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J100/UTM_b100_water_frequency_2010_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J100/UTM_b010_water_frequency_2010_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J101
#' b101_water_frequency_2010_30m_tif_exp.
hf <-"D:/GISbases/b101JaguarAlbersBuffer70/b101_water_frequency_2010_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J101trk))
#crs(hf)   
#compareCRS(hf,get_crs(J101trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J101/UTM_b101_water_frequency_2010_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J101/UTM_b101_water_frequency_2010_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J102
#' b102_water_frequency_2010_30m_tif_exp.
hf <-"D:/GISbases/b102JaguarAlbersBuffer70/b102_water_frequency_2010_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J102trk))
#crs(hf)   
#compareCRS(hf,get_crs(J102trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J102/UTM_b102_water_frequency_2010_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J102/UTM_b102_water_frequency_2010_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J103
#' b103_water_frequency_2010_30m_tif_exp.
hf <-"D:/GISbases/b103JaguarAlbersBuffer70/b103_water_frequency_2010_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J103trk))
#crs(hf)   
#compareCRS(hf,get_crs(J103trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J103/UTM_b103_water_frequency_2010_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J103/UTM_b103_water_frequency_2010_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J104
#' b104_water_frequency_2010_30m_tif_exp.
hf <-"D:/GISbases/b104JaguarAlbersBuffer70/b104_water_frequency_2010_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J104trk))
#crs(hf)   
#compareCRS(hf,get_crs(J104trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J104/UTM_b104_water_frequency_2010_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J104/UTM_b104_water_frequency_2010_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J105
#' b105_water_frequency_2010_30m_tif_exp.
hf <-"D:/GISbases/b105JaguarAlbersBuffer70/b105_water_frequency_2010_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J105trk))
#crs(hf)   
#compareCRS(hf,get_crs(J105trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J105/UTM_b105_water_frequency_2010_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J105/UTM_b105_water_frequency_2010_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J106
#' b106_water_frequency_2010_30m_tif_exp.
hf <-"D:/GISbases/b106JaguarAlbersBuffer70/b106_water_frequency_2010_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J106trk))
#crs(hf)   
#compareCRS(hf,get_crs(J106trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J106/UTM_b106_water_frequency_2010_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J106/UTM_b106_water_frequency_2010_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J107
#' b107_water_frequency_2010_30m_tif_exp.
hf <-"D:/GISbases/b107JaguarAlbersBuffer70/b107_water_frequency_2010_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J107trk))
#crs(hf)   
#compareCRS(hf,get_crs(J107trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J107/UTM_b107_water_frequency_2010_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J107/UTM_b107_water_frequency_2010_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J108
#' b108_water_frequency_2010_30m_tif_exp.
hf <-"D:/GISbases/b108JaguarAlbersBuffer70/b108_water_frequency_2010_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J108trk))
#crs(hf)   
#compareCRS(hf,get_crs(J108trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J108/UTM_b108_water_frequency_2010_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J108/UTM_b108_water_frequency_2010_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J109
#' b109_water_frequency_2010_30m_tif_exp.
hf <-"D:/GISbases/b109JaguarAlbersBuffer70/b109_water_frequency_2010_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J109trk))
#crs(hf)   
#compareCRS(hf,get_crs(J109trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J109/UTM_b109_water_frequency_2010_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J109/UTM_b109_water_frequency_2010_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J110
#' b110_water_frequency_2010_30m_tif_exp.
hf <-"D:/GISbases/b110JaguarAlbersBuffer70/b110_water_frequency_2010_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J110trk))
#crs(hf)   
#compareCRS(hf,get_crs(J110trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J110/UTM_b110_water_frequency_2010_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J110/UTM_b110_water_frequency_2010_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J111
#' b111_water_frequency_2010_30m_tif_exp.
hf <-"D:/GISbases/b111JaguarAlbersBuffer70/b111_water_frequency_2010_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J111trk))
#crs(hf)   
#compareCRS(hf,get_crs(J111trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J111/UTM_b111_water_frequency_2010_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J111/UTM_b111_water_frequency_2010_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J112
#' b112_water_frequency_2010_30m_tif_exp.
hf <-"D:/GISbases/b112JaguarAlbersBuffer70/b112_water_frequency_2010_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J112trk))
#crs(hf)   
#compareCRS(hf,get_crs(J112trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J112/UTM_b112_water_frequency_2010_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J112/UTM_b112_water_frequency_2010_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J113
#' b113_water_frequency_2010_30m_tif_exp.
hf <-"D:/GISbases/b113JaguarAlbersBuffer70/b113_water_frequency_2010_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J113trk))
#crs(hf)   
#compareCRS(hf,get_crs(J113trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J113/UTM_b113_water_frequency_2010_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J113/UTM_b113_water_frequency_2010_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J114
#' b114_water_frequency_2010_30m_tif_exp.
hf <-"D:/GISbases/b114JaguarAlbersBuffer70/b114_water_frequency_2010_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J114trk))
#crs(hf)   
#compareCRS(hf,get_crs(J114trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J114/UTM_b114_water_frequency_2010_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J114/UTM_b114_water_frequency_2010_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J115
#' b115_water_frequency_2010_30m_tif_exp.
hf <-"D:/GISbases/b115JaguarAlbersBuffer70/b115_water_frequency_2010_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J115trk))
#crs(hf)   
#compareCRS(hf,get_crs(J115trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J115/UTM_b115_water_frequency_2010_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J115/UTM_b115_water_frequency_2010_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J116
#' b116_water_frequency_2010_30m_tif_exp.
hf <-"D:/GISbases/b116JaguarAlbersBuffer70/b116_water_frequency_2010_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J116trk))
#crs(hf)   
#compareCRS(hf,get_crs(J116trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J116/UTM_b116_water_frequency_2010_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J116/UTM_b116_water_frequency_2010_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J117
#' b117_water_frequency_2010_30m_tif_exp.
hf <-"D:/GISbases/b117JaguarAlbersBuffer70/b117_water_frequency_2010_30m_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J117trk))
#crs(hf)   
#compareCRS(hf,get_crs(J117trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J117/UTM_b117_water_frequency_2010_30m_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J117/UTM_b117_water_frequency_2010_30m_tif_exp.tif";(test=raster(test)); x11(); plot(test)




################################################################################################################################
# dist2drainage_exp.

#' J1
#' dist2drainage_exp.tif
hf <-"D:/GISbases/b001JaguarAlbersBuffer70/b001_dist2drainage_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J1trk))
#crs(hf)   
#compareCRS(hf,get_crs(J1trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J1/UTM_b001_dist2drainage_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J1/UTM_b001_dist2drainage_exp.tif";(test=raster(test)); x11(); plot(test)


#' J2
#' b002_dist2drainage_exp.
hf <-"D:/GISbases/b002JaguarAlbersBuffer70/b002_dist2drainage_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J2trk))
#crs(hf)   
#compareCRS(hf,get_crs(J2trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J2/UTM_b002_dist2drainage_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J2/UTM_b002_dist2drainage_exp.tif";(test=raster(test)); x11(); plot(test)

#' J3
#' b003_dist2drainage_exp.
hf <-"D:/GISbases/b003JaguarAlbersBuffer70/b003_dist2drainage_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J3trk))
#crs(hf)   
#compareCRS(hf,get_crs(J3trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J3/UTM_b003_dist2drainage_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J3/UTM_b003_dist2drainage_exp.tif";(test=raster(test)); x11(); plot(test)

#' J4
#' b004_dist2drainage_exp.
hf <-"D:/GISbases/b004JaguarAlbersBuffer70/b004_dist2drainage_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J4trk))
#crs(hf)   
#compareCRS(hf,get_crs(J4trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J4/UTM_b004_dist2drainage_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J4/UTM_b004_dist2drainage_exp.tif";(test=raster(test)); x11(); plot(test)

#' J5
#' b005_dist2drainage_exp.
hf <-"D:/GISbases/b005JaguarAlbersBuffer70/b005_dist2drainage_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J5trk))
#crs(hf)   
#compareCRS(hf,get_crs(J5trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J5/UTM_b005_dist2drainage_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J5/UTM_b005_dist2drainage_exp.tif";(test=raster(test)); x11(); plot(test)

#' J6
#' b006_dist2drainage_exp.
hf <-"D:/GISbases/b006JaguarAlbersBuffer70/b006_dist2drainage_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J6trk))
#crs(hf)   
#compareCRS(hf,get_crs(J6trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J6/UTM_b006_dist2drainage_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J6/UTM_b006_dist2drainage_exp.tif";(test=raster(test)); x11(); plot(test)

#' J7
#' b007_dist2drainage_exp.
hf <-"D:/GISbases/b007JaguarAlbersBuffer70/b007_dist2drainage_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J7trk))
#crs(hf)   
#compareCRS(hf,get_crs(J7trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J7/UTM_b007_dist2drainage_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J7/UTM_b007_dist2drainage_exp.tif";(test=raster(test)); x11(); plot(test)

#' J8
#' b008_dist2drainage_exp.
hf <-"D:/GISbases/b008JaguarAlbersBuffer70/b008_dist2drainage_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J8trk))
#crs(hf)   
#compareCRS(hf,get_crs(J8trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J8/UTM_b008_dist2drainage_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J8/UTM_b008_dist2drainage_exp.tif";(test=raster(test)); x11(); plot(test)

#' J9
#' b009_dist2drainage_exp.
hf <-"D:/GISbases/b009JaguarAlbersBuffer70/b009_dist2drainage_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J9trk))
#crs(hf)   
#compareCRS(hf,get_crs(J9trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J9/UTM_b009_dist2drainage_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J9/UTM_b009_dist2drainage_exp.tif";(test=raster(test)); x11(); plot(test)

#' J10
#' b010_dist2drainage_exp.
hf <-"D:/GISbases/b010JaguarAlbersBuffer70/b010_dist2drainage_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J10trk))
#crs(hf)   
#compareCRS(hf,get_crs(J10trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J10/UTM_b010_dist2drainage_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J10/UTM_b010_dist2drainage_exp.tif";(test=raster(test)); x11(); plot(test)

gc()
gc()
gc()


#' J11
#' b011_dist2drainage_exp.
hf <-"D:/GISbases/b011JaguarAlbersBuffer70/b011_dist2drainage_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J11trk))
#crs(hf)   
#compareCRS(hf,get_crs(J11trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J11/UTM_b011_dist2drainage_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J11/UTM_b011_dist2drainage_exp.tif";(test=raster(test)); x11(); plot(test)

#' J12
#' b012_dist2drainage_exp.
hf <-"D:/GISbases/b012JaguarAlbersBuffer70/b012_dist2drainage_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J12trk))
#crs(hf)   
#compareCRS(hf,get_crs(J12trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J12/UTM_b012_dist2drainage_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J12/UTM_b012_dist2drainage_exp.tif";(test=raster(test)); x11(); plot(test)

#' J13
#' b013_dist2drainage_exp.
hf <-"D:/GISbases/b013JaguarAlbersBuffer70/b013_dist2drainage_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J13trk))
#crs(hf)   
#compareCRS(hf,get_crs(J13trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J13/UTM_b013_dist2drainage_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J13/UTM_b013_dist2drainage_exp.tif";(test=raster(test)); x11(); plot(test)

#' J14
#' b014_dist2drainage_exp.
hf <-"D:/GISbases/b014JaguarAlbersBuffer70/b014_dist2drainage_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J14trk))
#crs(hf)   
#compareCRS(hf,get_crs(J14trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J14/UTM_b014_dist2drainage_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J14/UTM_b014_dist2drainage_exp.tif";(test=raster(test)); x11(); plot(test)

#' J15
#' b015_dist2drainage_exp.
hf <-"D:/GISbases/b015JaguarAlbersBuffer70/b015_dist2drainage_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J15trk))
#crs(hf)   
#compareCRS(hf,get_crs(J15trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J15/UTM_b015_dist2drainage_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J15/UTM_b015_dist2drainage_exp.tif";(test=raster(test)); x11(); plot(test)

#' J16
#' b016_dist2drainage_exp.
hf <-"D:/GISbases/b016JaguarAlbersBuffer70/b016_dist2drainage_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J16trk))
#crs(hf)   
#compareCRS(hf,get_crs(J16trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J16/UTM_b016_dist2drainage_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J16/UTM_b016_dist2drainage_exp.tif";(test=raster(test)); x11(); plot(test)

#' J17
#' b017_dist2drainage_exp.
hf <-"D:/GISbases/b017JaguarAlbersBuffer70/b017_dist2drainage_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J17trk))
#crs(hf)   
#compareCRS(hf,get_crs(J17trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J17/UTM_b017_dist2drainage_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J17/UTM_b017_dist2drainage_exp.tif";(test=raster(test)); x11(); plot(test)

#' J18
#' b018_dist2drainage_exp.
hf <-"D:/GISbases/b018JaguarAlbersBuffer70/b018_dist2drainage_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J18trk))
#crs(hf)   
#compareCRS(hf,get_crs(J18trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J18/UTM_b018_dist2drainage_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J18/UTM_b018_dist2drainage_exp.tif";(test=raster(test)); x11(); plot(test)

#' J19
#' b019_dist2drainage_exp.
hf <-"D:/GISbases/b019JaguarAlbersBuffer70/b019_dist2drainage_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J19trk))
#crs(hf)   
#compareCRS(hf,get_crs(J19trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J19/UTM_b019_dist2drainage_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J19/UTM_b019_dist2drainage_exp.tif";(test=raster(test)); x11(); plot(test)

#' J20
#' b020_dist2drainage_exp.
hf <-"D:/GISbases/b020JaguarAlbersBuffer70/b020_dist2drainage_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J20trk))
#crs(hf)   
#compareCRS(hf,get_crs(J20trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J20/UTM_b020_dist2drainage_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J20/UTM_b020_dist2drainage_exp.tif";(test=raster(test)); x11(); plot(test)

#' J21
#' b021_dist2drainage_exp.
hf <-"D:/GISbases/b021JaguarAlbersBuffer70/b021_dist2drainage_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J21trk))
#crs(hf)   
#compareCRS(hf,get_crs(J21trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J21/UTM_b021_dist2drainage_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J21/UTM_b021_dist2drainage_exp.tif";(test=raster(test)); x11(); plot(test)

#' J22
#' b022_dist2drainage_exp.
hf <-"D:/GISbases/b022JaguarAlbersBuffer70/b022_dist2drainage_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J22trk))
#crs(hf)   
#compareCRS(hf,get_crs(J22trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J22/UTM_b022_dist2drainage_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J22/UTM_b022_dist2drainage_exp.tif";(test=raster(test)); x11(); plot(test)

#' J23
#' b023_dist2drainage_exp.
hf <-"D:/GISbases/b023JaguarAlbersBuffer70/b023_dist2drainage_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J23trk))
#crs(hf)   
#compareCRS(hf,get_crs(J23trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J23/UTM_b023_dist2drainage_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J23/UTM_b023_dist2drainage_exp.tif";(test=raster(test)); x11(); plot(test)

#' J24
#' b024_dist2drainage_exp.
hf <-"D:/GISbases/b024JaguarAlbersBuffer70/b024_dist2drainage_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J24trk))
#crs(hf)   
#compareCRS(hf,get_crs(J24trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J24/UTM_b024_dist2drainage_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J24/UTM_b024_dist2drainage_exp.tif";(test=raster(test)); x11(); plot(test)

#' J25
#' b025_dist2drainage_exp.
hf <-"D:/GISbases/b025JaguarAlbersBuffer70/b025_dist2drainage_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J25trk))
#crs(hf)   
#compareCRS(hf,get_crs(J25trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J25/UTM_b025_dist2drainage_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J25/UTM_b025_dist2drainage_exp.tif";(test=raster(test)); x11(); plot(test)

#' J26
#' b026_dist2drainage_exp.
hf <-"D:/GISbases/b026JaguarAlbersBuffer70/b026_dist2drainage_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J26trk))
#crs(hf)   
#compareCRS(hf,get_crs(J26trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J26/UTM_b026_dist2drainage_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J26/UTM_b026_dist2drainage_exp.tif";(test=raster(test)); x11(); plot(test)

#' J27
#' b027_dist2drainage_exp.
hf <-"D:/GISbases/b027JaguarAlbersBuffer70/b027_dist2drainage_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J27trk))
#crs(hf)   
#compareCRS(hf,get_crs(J27trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J27/UTM_b027_dist2drainage_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J27/UTM_b027_dist2drainage_exp.tif";(test=raster(test)); x11(); plot(test)

#' J28
#' b028_dist2drainage_exp.
hf <-"D:/GISbases/b028JaguarAlbersBuffer70/b028_dist2drainage_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J28trk))
#crs(hf)   
#compareCRS(hf,get_crs(J28trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J28/UTM_b028_dist2drainage_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J28/UTM_b028_dist2drainage_exp.tif";(test=raster(test)); x11(); plot(test)

#' J29
#' b029_dist2drainage_exp.
hf <-"D:/GISbases/b029JaguarAlbersBuffer70/b029_dist2drainage_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J29trk))
#crs(hf)   
#compareCRS(hf,get_crs(J29trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J29/UTM_b029_dist2drainage_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J29/UTM_b029_dist2drainage_exp.tif";(test=raster(test)); x11(); plot(test)

#' J30
#' b030_dist2drainage_exp.
hf <-"D:/GISbases/b030JaguarAlbersBuffer70/b030_dist2drainage_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J30trk))
#crs(hf)   
#compareCRS(hf,get_crs(J30trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J30/UTM_b030_dist2drainage_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J30/UTM_b030_dist2drainage_exp.tif";(test=raster(test)); x11(); plot(test)

#' J30
#' b030_dist2drainage_exp.
hf <-"D:/GISbases/b030JaguarAlbersBuffer70/b030_dist2drainage_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J30trk))
#crs(hf)   
#compareCRS(hf,get_crs(J30trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J30/UTM_b030_dist2drainage_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J30/UTM_b030_dist2drainage_exp.tif";(test=raster(test)); x11(); plot(test)

#' J31
#' b031_dist2drainage_exp.
hf <-"D:/GISbases/b031JaguarAlbersBuffer70/b031_dist2drainage_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J31trk))
#crs(hf)   
#compareCRS(hf,get_crs(J31trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J31/UTM_b031_dist2drainage_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J31/UTM_b031_dist2drainage_exp.tif";(test=raster(test)); x11(); plot(test)

#' J32
#' b032_dist2drainage_exp.
hf <-"D:/GISbases/b032JaguarAlbersBuffer70/b032_dist2drainage_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J32trk))
#crs(hf)   
#compareCRS(hf,get_crs(J32trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J32/UTM_b032_dist2drainage_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J32/UTM_b032_dist2drainage_exp.tif";(test=raster(test)); x11(); plot(test)

#' J33
#' b033_dist2drainage_exp.
hf <-"D:/GISbases/b033JaguarAlbersBuffer70/b033_dist2drainage_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J33trk))
#crs(hf)   
#compareCRS(hf,get_crs(J33trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J33/UTM_b033_dist2drainage_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J33/UTM_b033_dist2drainage_exp.tif";(test=raster(test)); x11(); plot(test)

#' J34
#' b034_dist2drainage_exp.
hf <-"D:/GISbases/b034JaguarAlbersBuffer70/b034_dist2drainage_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J34trk))
#crs(hf)   
#compareCRS(hf,get_crs(J34trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J34/UTM_b034_dist2drainage_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J34/UTM_b034_dist2drainage_exp.tif";(test=raster(test)); x11(); plot(test)

#' J35
#' b035_dist2drainage_exp.
hf <-"D:/GISbases/b035JaguarAlbersBuffer70/b035_dist2drainage_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J35trk))
#crs(hf)   
#compareCRS(hf,get_crs(J35trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J35/UTM_b035_dist2drainage_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J35/UTM_b035_dist2drainage_exp.tif";(test=raster(test)); x11(); plot(test)

#' J36
#' b036_dist2drainage_exp.
hf <-"D:/GISbases/b036JaguarAlbersBuffer70/b036_dist2drainage_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J36trk))
#crs(hf)   
#compareCRS(hf,get_crs(J36trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J36/UTM_b036_dist2drainage_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J36/UTM_b036_dist2drainage_exp.tif";(test=raster(test)); x11(); plot(test)

#' J37
#' b037_dist2drainage_exp.
hf <-"D:/GISbases/b037JaguarAlbersBuffer70/b037_dist2drainage_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J37trk))
#crs(hf)   
#compareCRS(hf,get_crs(J37trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J37/UTM_b037_dist2drainage_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J37/UTM_b037_dist2drainage_exp.tif";(test=raster(test)); x11(); plot(test)

#' J38
#' b038_dist2drainage_exp.
hf <-"D:/GISbases/b038JaguarAlbersBuffer70/b038_dist2drainage_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J38trk))
#crs(hf)   
#compareCRS(hf,get_crs(J38trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J38/UTM_b038_dist2drainage_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J38/UTM_b038_dist2drainage_exp.tif";(test=raster(test)); x11(); plot(test)

#' J39
#' b039_dist2drainage_exp.
hf <-"D:/GISbases/b039JaguarAlbersBuffer70/b039_dist2drainage_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J39trk))
#crs(hf)   
#compareCRS(hf,get_crs(J39trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J39/UTM_b039_dist2drainage_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J39/UTM_b039_dist2drainage_exp.tif";(test=raster(test)); x11(); plot(test)

#' J40
#' b040_dist2drainage_exp.
hf <-"D:/GISbases/b040JaguarAlbersBuffer70/b040_dist2drainage_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J40trk))
#crs(hf)   
#compareCRS(hf,get_crs(J40trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J40/UTM_b040_dist2drainage_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J40/UTM_b040_dist2drainage_exp.tif";(test=raster(test)); x11(); plot(test)

#' J41
#' b041_dist2drainage_exp.
hf <-"D:/GISbases/b041JaguarAlbersBuffer70/b041_dist2drainage_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J41trk))
#crs(hf)   
#compareCRS(hf,get_crs(J41trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J41/UTM_b041_dist2drainage_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J41/UTM_b041_dist2drainage_exp.tif";(test=raster(test)); x11(); plot(test)

#' J42
#' b042_dist2drainage_exp.
hf <-"D:/GISbases/b042JaguarAlbersBuffer70/b042_dist2drainage_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J42trk))
#crs(hf)   
#compareCRS(hf,get_crs(J42trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J42/UTM_b042_dist2drainage_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J42/UTM_b042_dist2drainage_exp.tif";(test=raster(test)); x11(); plot(test)

gc()
gc()
gc()

#' J43
#' b043_dist2drainage_exp.
hf <-"D:/GISbases/b043JaguarAlbersBuffer70/b043_dist2drainage_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J43trk))
#crs(hf)   
#compareCRS(hf,get_crs(J43trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J43/UTM_b043_dist2drainage_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J43/UTM_b043_dist2drainage_exp.tif";(test=raster(test)); x11(); plot(test)

#' J44
#' b044_dist2drainage_exp.
hf <-"D:/GISbases/b044JaguarAlbersBuffer70/b044_dist2drainage_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J44trk))
#crs(hf)   
#compareCRS(hf,get_crs(J44trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J44/UTM_b044_dist2drainage_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J44/UTM_b044_dist2drainage_exp.tif";(test=raster(test)); x11(); plot(test)

#' J45
#' b045_dist2drainage_exp.
hf <-"D:/GISbases/b045JaguarAlbersBuffer70/b045_dist2drainage_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J45trk))
#crs(hf)   
#compareCRS(hf,get_crs(J45trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J45/UTM_b045_dist2drainage_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J45/UTM_b045_dist2drainage_exp.tif";(test=raster(test)); x11(); plot(test)

#' J46
#' b046_dist2drainage_exp.
hf <-"D:/GISbases/b046JaguarAlbersBuffer70/b046_dist2drainage_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J46trk))
#crs(hf)   
#compareCRS(hf,get_crs(J46trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J46/UTM_b046_dist2drainage_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J46/UTM_b046_dist2drainage_exp.tif";(test=raster(test)); x11(); plot(test)

#' J47
#' b047_dist2drainage_exp.
hf <-"D:/GISbases/b047JaguarAlbersBuffer70/b047_dist2drainage_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J47trk))
#crs(hf)   
#compareCRS(hf,get_crs(J47trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J47/UTM_b047_dist2drainage_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J47/UTM_b047_dist2drainage_exp.tif";(test=raster(test)); x11(); plot(test)

#' J48
#' b048_dist2drainage_exp.
hf <-"D:/GISbases/b048JaguarAlbersBuffer70/b048_dist2drainage_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J48trk))
#crs(hf)   
#compareCRS(hf,get_crs(J48trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J48/UTM_b048_dist2drainage_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J48/UTM_b048_dist2drainage_exp.tif";(test=raster(test)); x11(); plot(test)

#' J49
#' b049_dist2drainage_exp.
hf <-"D:/GISbases/b049JaguarAlbersBuffer70/b049_dist2drainage_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J49trk))
#crs(hf)   
#compareCRS(hf,get_crs(J49trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J49/UTM_b049_dist2drainage_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J49/UTM_b049_dist2drainage_exp.tif";(test=raster(test)); x11(); plot(test)

#' J50
#' b050_dist2drainage_exp.
hf <-"D:/GISbases/b050JaguarAlbersBuffer70/b050_dist2drainage_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J50trk))
#crs(hf)   
#compareCRS(hf,get_crs(J50trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J50/UTM_b050_dist2drainage_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J50/UTM_b050_dist2drainage_exp.tif";(test=raster(test)); x11(); plot(test)

#' J51
#' b051_dist2drainage_exp.
hf <-"D:/GISbases/b051JaguarAlbersBuffer70/b051_dist2drainage_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J51trk))
#crs(hf)   
#compareCRS(hf,get_crs(J51trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J51/UTM_b051_dist2drainage_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J51/UTM_b051_dist2drainage_exp.tif";(test=raster(test)); x11(); plot(test)

#' J52
#' b052_dist2drainage_exp.
hf <-"D:/GISbases/b052JaguarAlbersBuffer70/b052_dist2drainage_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J52trk))
#crs(hf)   
#compareCRS(hf,get_crs(J52trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J52/UTM_b052_dist2drainage_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J52/UTM_b052_dist2drainage_exp.tif";(test=raster(test)); x11(); plot(test)

#' J53
#' b053_dist2drainage_exp.
hf <-"D:/GISbases/b053JaguarAlbersBuffer70/b053_dist2drainage_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J53trk))
#crs(hf)   
#compareCRS(hf,get_crs(J53trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J53/UTM_b053_dist2drainage_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J53/UTM_b053_dist2drainage_exp.tif";(test=raster(test)); x11(); plot(test)

#' J54
#' b054_dist2drainage_exp.
hf <-"D:/GISbases/b054JaguarAlbersBuffer70/b054_dist2drainage_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J54trk))
#crs(hf)   
#compareCRS(hf,get_crs(J54trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J54/UTM_b054_dist2drainage_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J54/UTM_b054_dist2drainage_exp.tif";(test=raster(test)); x11(); plot(test)

#' J55
#' b055_dist2drainage_exp.
hf <-"D:/GISbases/b055JaguarAlbersBuffer70/b055_dist2drainage_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J55trk))
#crs(hf)   
#compareCRS(hf,get_crs(J55trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J55/UTM_b055_dist2drainage_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J55/UTM_b055_dist2drainage_exp.tif";(test=raster(test)); x11(); plot(test)

#' J56
#' b056_dist2drainage_exp.
hf <-"D:/GISbases/b056JaguarAlbersBuffer70/b056_dist2drainage_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J56trk))
#crs(hf)   
#compareCRS(hf,get_crs(J56trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J56/UTM_b056_dist2drainage_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J56/UTM_b056_dist2drainage_exp.tif";(test=raster(test)); x11(); plot(test)

#' J57
#' b057_dist2drainage_exp.
hf <-"D:/GISbases/b057JaguarAlbersBuffer70/b057_dist2drainage_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J57trk))
#crs(hf)   
#compareCRS(hf,get_crs(J57trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J57/UTM_b057_dist2drainage_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J57/UTM_b057_dist2drainage_exp.tif";(test=raster(test)); x11(); plot(test)

#' J58
#' b058_dist2drainage_exp.
hf <-"D:/GISbases/b058JaguarAlbersBuffer70/b058_dist2drainage_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J58trk))
#crs(hf)   
#compareCRS(hf,get_crs(J58trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J58/UTM_b058_dist2drainage_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J58/UTM_b058_dist2drainage_exp.tif";(test=raster(test)); x11(); plot(test)

#' J59
#' b059_dist2drainage_exp.
hf <-"D:/GISbases/b059JaguarAlbersBuffer70/b059_dist2drainage_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J59trk))
#crs(hf)   
#compareCRS(hf,get_crs(J59trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J59/UTM_b059_dist2drainage_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J59/UTM_b059_dist2drainage_exp.tif";(test=raster(test)); x11(); plot(test)

#' J60
#' b060_dist2drainage_exp.
hf <-"D:/GISbases/b060JaguarAlbersBuffer70/b060_dist2drainage_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J60trk))
#crs(hf)   
#compareCRS(hf,get_crs(J60trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J60/UTM_b060_dist2drainage_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J60/UTM_b060_dist2drainage_exp.tif";(test=raster(test)); x11(); plot(test)

#' J61
#' b061_dist2drainage_exp.
hf <-"D:/GISbases/b061JaguarAlbersBuffer70/b061_dist2drainage_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J61trk))
#crs(hf)   
#compareCRS(hf,get_crs(J61trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J61/UTM_b061_dist2drainage_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J61/UTM_b061_dist2drainage_exp.tif";(test=raster(test)); x11(); plot(test)

#' J62
#' b062_dist2drainage_exp.
hf <-"D:/GISbases/b062JaguarAlbersBuffer70/b062_dist2drainage_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J62trk))
#crs(hf)   
#compareCRS(hf,get_crs(J62trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J62/UTM_b062_dist2drainage_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J62/UTM_b062_dist2drainage_exp.tif";(test=raster(test)); x11(); plot(test)

#' J63
#' b063_dist2drainage_exp.
hf <-"D:/GISbases/b063JaguarAlbersBuffer70/b063_dist2drainage_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J63trk))
#crs(hf)   
#compareCRS(hf,get_crs(J63trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J63/UTM_b063_dist2drainage_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J63/UTM_b063_dist2drainage_exp.tif";(test=raster(test)); x11(); plot(test)

gc()
gc()
gc()


#' J64
#' b064_dist2drainage_exp.
hf <-"D:/GISbases/b064JaguarAlbersBuffer70/b064_dist2drainage_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J64trk))
#crs(hf)   
#compareCRS(hf,get_crs(J64trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J64/UTM_b064_dist2drainage_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J64/UTM_b064_dist2drainage_exp.tif";(test=raster(test)); x11(); plot(test)

#' J65
#' b065_dist2drainage_exp.
hf <-"D:/GISbases/b065JaguarAlbersBuffer70/b065_dist2drainage_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J65trk))
#crs(hf)   
#compareCRS(hf,get_crs(J65trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J65/UTM_b065_dist2drainage_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J65/UTM_b065_dist2drainage_exp.tif";(test=raster(test)); x11(); plot(test)

#' J66
#' b066_dist2drainage_exp.
hf <-"D:/GISbases/b066JaguarAlbersBuffer70/b066_dist2drainage_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J66trk))
#crs(hf)   
#compareCRS(hf,get_crs(J66trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J66/UTM_b066_dist2drainage_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J66/UTM_b066_dist2drainage_exp.tif";(test=raster(test)); x11(); plot(test)

#' J67
#' b067_dist2drainage_exp.
hf <-"D:/GISbases/b067JaguarAlbersBuffer70/b067_dist2drainage_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J67trk))
#crs(hf)   
#compareCRS(hf,get_crs(J67trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J67/UTM_b067_dist2drainage_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J67/UTM_b067_dist2drainage_exp.tif";(test=raster(test)); x11(); plot(test)

#' J68
#' b068_dist2drainage_exp.
hf <-"D:/GISbases/b068JaguarAlbersBuffer70/b068_dist2drainage_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J68trk))
#crs(hf)   
#compareCRS(hf,get_crs(J68trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J68/UTM_b068_dist2drainage_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J68/UTM_b068_dist2drainage_exp.tif";(test=raster(test)); x11(); plot(test)

#' J69
#' b069_dist2drainage_exp.
hf <-"D:/GISbases/b069JaguarAlbersBuffer70/b069_dist2drainage_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J69trk))
#crs(hf)   
#compareCRS(hf,get_crs(J69trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J69/UTM_b069_dist2drainage_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J69/UTM_b069_dist2drainage_exp.tif";(test=raster(test)); x11(); plot(test)

gc()
gc()
gc()


#' J70
#' b070_dist2drainage_exp.
hf <-"D:/GISbases/b070JaguarAlbersBuffer70/b070_dist2drainage_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J70trk))
#crs(hf)   
#compareCRS(hf,get_crs(J70trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J70/UTM_b070_dist2drainage_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J70/UTM_b070_dist2drainage_exp.tif";(test=raster(test)); x11(); plot(test)

#' J71
#' b071_dist2drainage_exp.
hf <-"D:/GISbases/b071JaguarAlbersBuffer70/b071_dist2drainage_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J71trk))
#crs(hf)   
#compareCRS(hf,get_crs(J71trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J71/UTM_b071_dist2drainage_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J71/UTM_b071_dist2drainage_exp.tif";(test=raster(test)); x11(); plot(test)

#' J72
#' b072_dist2drainage_exp.
hf <-"D:/GISbases/b072JaguarAlbersBuffer70/b072_dist2drainage_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J72trk))
#crs(hf)   
#compareCRS(hf,get_crs(J72trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J72/UTM_b072_dist2drainage_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J72/UTM_b072_dist2drainage_exp.tif";(test=raster(test)); x11(); plot(test)

#' J73
#' b073_dist2drainage_exp.
hf <-"D:/GISbases/b073JaguarAlbersBuffer70/b073_dist2drainage_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J73trk))
#crs(hf)   
#compareCRS(hf,get_crs(J73trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J73/UTM_b073_dist2drainage_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J73/UTM_b073_dist2drainage_exp.tif";(test=raster(test)); x11(); plot(test)

#' J74
#' b074_dist2drainage_exp.
hf <-"D:/GISbases/b074JaguarAlbersBuffer70/b074_dist2drainage_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J74trk))
#crs(hf)   
#compareCRS(hf,get_crs(J74trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J74/UTM_b074_dist2drainage_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J74/UTM_b074_dist2drainage_exp.tif";(test=raster(test)); x11(); plot(test)

#' J75
#' b075_dist2drainage_exp.
hf <-"D:/GISbases/b075JaguarAlbersBuffer70/b075_dist2drainage_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J75trk))
#crs(hf)   
#compareCRS(hf,get_crs(J75trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J75/UTM_b075_dist2drainage_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J75/UTM_b075_dist2drainage_exp.tif";(test=raster(test)); x11(); plot(test)

#' J76
#' b076_dist2drainage_exp.
hf <-"D:/GISbases/b076JaguarAlbersBuffer70/b076_dist2drainage_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J76trk))
#crs(hf)   
#compareCRS(hf,get_crs(J76trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J76/UTM_b076_dist2drainage_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J76/UTM_b076_dist2drainage_exp.tif";(test=raster(test)); x11(); plot(test)

#' J77
#' b077_dist2drainage_exp.
hf <-"D:/GISbases/b077JaguarAlbersBuffer70/b077_dist2drainage_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J77trk))
#crs(hf)   
#compareCRS(hf,get_crs(J77trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J77/UTM_b077_dist2drainage_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J77/UTM_b077_dist2drainage_exp.tif";(test=raster(test)); x11(); plot(test)

#' J78
#' b078_dist2drainage_exp.
hf <-"D:/GISbases/b078JaguarAlbersBuffer70/b078_dist2drainage_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J78trk))
#crs(hf)   
#compareCRS(hf,get_crs(J78trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J78/UTM_b078_dist2drainage_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J78/UTM_b078_dist2drainage_exp.tif";(test=raster(test)); x11(); plot(test)

#' J79
#' b079_dist2drainage_exp.
hf <-"D:/GISbases/b079JaguarAlbersBuffer70/b079_dist2drainage_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J79trk))
#crs(hf)   
#compareCRS(hf,get_crs(J79trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J79/UTM_b079_dist2drainage_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J79/UTM_b079_dist2drainage_exp.tif";(test=raster(test)); x11(); plot(test)

#' J80
#' b080_dist2drainage_exp.
hf <-"D:/GISbases/b080JaguarAlbersBuffer70/b080_dist2drainage_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J80trk))
#crs(hf)   
#compareCRS(hf,get_crs(J80trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J80/UTM_b080_dist2drainage_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J80/UTM_b080_dist2drainage_exp.tif";(test=raster(test)); x11(); plot(test)

#' J81
#' b081_dist2drainage_exp.
hf <-"D:/GISbases/b081JaguarAlbersBuffer70/b081_dist2drainage_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J81trk))
#crs(hf)   
#compareCRS(hf,get_crs(J81trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J81/UTM_b081_dist2drainage_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J81/UTM_b081_dist2drainage_exp.tif";(test=raster(test)); x11(); plot(test)

#' J82
#' b082_dist2drainage_exp.
hf <-"D:/GISbases/b082JaguarAlbersBuffer70/b082_dist2drainage_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J82trk))
#crs(hf)   
#compareCRS(hf,get_crs(J82trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J82/UTM_b082_dist2drainage_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J82/UTM_b082_dist2drainage_exp.tif";(test=raster(test)); x11(); plot(test)

#' J83
#' b083_dist2drainage_exp.
hf <-"D:/GISbases/b083JaguarAlbersBuffer70/b083_dist2drainage_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J83trk))
#crs(hf)   
#compareCRS(hf,get_crs(J83trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J83/UTM_b083_dist2drainage_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J83/UTM_b083_dist2drainage_exp.tif";(test=raster(test)); x11(); plot(test)

#' J84
#' b084_dist2drainage_exp.
hf <-"D:/GISbases/b084JaguarAlbersBuffer70/b084_dist2drainage_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J84trk))
#crs(hf)   
#compareCRS(hf,get_crs(J84trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J84/UTM_b084_dist2drainage_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J84/UTM_b084_dist2drainage_exp.tif";(test=raster(test)); x11(); plot(test)

#' J85
#' b085_dist2drainage_exp.
hf <-"D:/GISbases/b085JaguarAlbersBuffer70/b085_dist2drainage_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J85trk))
#crs(hf)   
#compareCRS(hf,get_crs(J85trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J85/UTM_b085_dist2drainage_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J85/UTM_b085_dist2drainage_exp.tif";(test=raster(test)); x11(); plot(test)

#' J86
#' b086_dist2drainage_exp.
hf <-"D:/GISbases/b086JaguarAlbersBuffer70/b086_dist2drainage_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J86trk))
#crs(hf)   
#compareCRS(hf,get_crs(J86trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J86/UTM_b086_dist2drainage_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J86/UTM_b086_dist2drainage_exp.tif";(test=raster(test)); x11(); plot(test)

#' J87
#' b087_dist2drainage_exp.
hf <-"D:/GISbases/b087JaguarAlbersBuffer70/b087_dist2drainage_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J87trk))
#crs(hf)   
#compareCRS(hf,get_crs(J87trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J87/UTM_b087_dist2drainage_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J87/UTM_b087_dist2drainage_exp.tif";(test=raster(test)); x11(); plot(test)

#' J88
#' b088_dist2drainage_exp.
hf <-"D:/GISbases/b088JaguarAlbersBuffer70/b088_dist2drainage_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J88trk))
#crs(hf)   
#compareCRS(hf,get_crs(J88trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J88/UTM_b088_dist2drainage_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J88/UTM_b088_dist2drainage_exp.tif";(test=raster(test)); x11(); plot(test)

#' J89
#' b089_dist2drainage_exp.
hf <-"D:/GISbases/b089JaguarAlbersBuffer70/b089_dist2drainage_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J89trk))
#crs(hf)   
#compareCRS(hf,get_crs(J89trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J89/UTM_b089_dist2drainage_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J89/UTM_b089_dist2drainage_exp.tif";(test=raster(test)); x11(); plot(test)

#' J90
#' b090_dist2drainage_exp.
hf <-"D:/GISbases/b090JaguarAlbersBuffer70/b090_dist2drainage_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J90trk))
#crs(hf)   
#compareCRS(hf,get_crs(J90trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J90/UTM_b090_dist2drainage_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J90/UTM_b090_dist2drainage_exp.tif";(test=raster(test)); x11(); plot(test)

#' J91
#' b091_dist2drainage_exp.
hf <-"D:/GISbases/b091JaguarAlbersBuffer70/b091_dist2drainage_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J91trk))
#crs(hf)   
#compareCRS(hf,get_crs(J91trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J91/UTM_b091_dist2drainage_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J91/UTM_b091_dist2drainage_exp.tif";(test=raster(test)); x11(); plot(test)

#' J92
#' b092_dist2drainage_exp.
hf <-"D:/GISbases/b092JaguarAlbersBuffer70/b092_dist2drainage_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J92trk))
#crs(hf)   
#compareCRS(hf,get_crs(J92trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J92/UTM_b092_dist2drainage_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J92/UTM_b092_dist2drainage_exp.tif";(test=raster(test)); x11(); plot(test)

#' J93
#' b093_dist2drainage_exp.
hf <-"D:/GISbases/b093JaguarAlbersBuffer70/b093_dist2drainage_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J93trk))
#crs(hf)   
#compareCRS(hf,get_crs(J93trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J93/UTM_b093_dist2drainage_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J93/UTM_b093_dist2drainage_exp.tif";(test=raster(test)); x11(); plot(test)

#' J94
#' b094_dist2drainage_exp.
hf <-"D:/GISbases/b094JaguarAlbersBuffer70/b094_dist2drainage_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J94trk))
#crs(hf)   
#compareCRS(hf,get_crs(J94trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J94/UTM_b094_dist2drainage_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J94/UTM_b094_dist2drainage_exp.tif";(test=raster(test)); x11(); plot(test)

#' J95
#' b095_dist2drainage_exp.
hf <-"D:/GISbases/b095JaguarAlbersBuffer70/b095_dist2drainage_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J95trk))
#crs(hf)   
#compareCRS(hf,get_crs(J95trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J95/UTM_b095_dist2drainage_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J95/UTM_b095_dist2drainage_exp.tif";(test=raster(test)); x11(); plot(test)

#' J96
#' b096_dist2drainage_exp.
hf <-"D:/GISbases/b096JaguarAlbersBuffer70/b096_dist2drainage_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J96trk))
#crs(hf)   
#compareCRS(hf,get_crs(J96trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J96/UTM_b096_dist2drainage_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J96/UTM_b096_dist2drainage_exp.tif";(test=raster(test)); x11(); plot(test)

#' J97
#' b097_dist2drainage_exp.
hf <-"D:/GISbases/b097JaguarAlbersBuffer70/b097_dist2drainage_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J97trk))
#crs(hf)   
#compareCRS(hf,get_crs(J97trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J97/UTM_b097_dist2drainage_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J97/UTM_b097_dist2drainage_exp.tif";(test=raster(test)); x11(); plot(test)

#' J98
#' b098_dist2drainage_exp.
hf <-"D:/GISbases/b098JaguarAlbersBuffer70/b098_dist2drainage_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J98trk))
#crs(hf)   
#compareCRS(hf,get_crs(J98trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J98/UTM_b098_dist2drainage_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J98/UTM_b098_dist2drainage_exp.tif";(test=raster(test)); x11(); plot(test)

#' J99
#' b099_dist2drainage_exp.
hf <-"D:/GISbases/b099JaguarAlbersBuffer70/b099_dist2drainage_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J99trk))
#crs(hf)   
#compareCRS(hf,get_crs(J99trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J99/UTM_b099_dist2drainage_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J99/UTM_b099_dist2drainage_exp.tif";(test=raster(test)); x11(); plot(test)

#' J100
#' b100_dist2drainage_exp.
hf <-"D:/GISbases/b100JaguarAlbersBuffer70/b100_dist2drainage_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J100trk))
#crs(hf)   
#compareCRS(hf,get_crs(J100trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J100/UTM_b100_dist2drainage_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J100/UTM_b010_dist2drainage_exp.tif";(test=raster(test)); x11(); plot(test)

#' J101
#' b101_dist2drainage_exp.
hf <-"D:/GISbases/b101JaguarAlbersBuffer70/b101_dist2drainage_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J101trk))
#crs(hf)   
#compareCRS(hf,get_crs(J101trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J101/UTM_b101_dist2drainage_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J101/UTM_b101_dist2drainage_exp.tif";(test=raster(test)); x11(); plot(test)

#' J102
#' b102_dist2drainage_exp.
hf <-"D:/GISbases/b102JaguarAlbersBuffer70/b102_dist2drainage_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J102trk))
#crs(hf)   
#compareCRS(hf,get_crs(J102trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J102/UTM_b102_dist2drainage_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J102/UTM_b102_dist2drainage_exp.tif";(test=raster(test)); x11(); plot(test)

#' J103
#' b103_dist2drainage_exp.
hf <-"D:/GISbases/b103JaguarAlbersBuffer70/b103_dist2drainage_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J103trk))
#crs(hf)   
#compareCRS(hf,get_crs(J103trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J103/UTM_b103_dist2drainage_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J103/UTM_b103_dist2drainage_exp.tif";(test=raster(test)); x11(); plot(test)

#' J104
#' b104_dist2drainage_exp.
hf <-"D:/GISbases/b104JaguarAlbersBuffer70/b104_dist2drainage_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J104trk))
#crs(hf)   
#compareCRS(hf,get_crs(J104trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J104/UTM_b104_dist2drainage_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J104/UTM_b104_dist2drainage_exp.tif";(test=raster(test)); x11(); plot(test)

#' J105
#' b105_dist2drainage_exp.
hf <-"D:/GISbases/b105JaguarAlbersBuffer70/b105_dist2drainage_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J105trk))
#crs(hf)   
#compareCRS(hf,get_crs(J105trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J105/UTM_b105_dist2drainage_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J105/UTM_b105_dist2drainage_exp.tif";(test=raster(test)); x11(); plot(test)

#' J106
#' b106_dist2drainage_exp.
hf <-"D:/GISbases/b106JaguarAlbersBuffer70/b106_dist2drainage_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J106trk))
#crs(hf)   
#compareCRS(hf,get_crs(J106trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J106/UTM_b106_dist2drainage_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J106/UTM_b106_dist2drainage_exp.tif";(test=raster(test)); x11(); plot(test)

#' J107
#' b107_dist2drainage_exp.
hf <-"D:/GISbases/b107JaguarAlbersBuffer70/b107_dist2drainage_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J107trk))
#crs(hf)   
#compareCRS(hf,get_crs(J107trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J107/UTM_b107_dist2drainage_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J107/UTM_b107_dist2drainage_exp.tif";(test=raster(test)); x11(); plot(test)

#' J108
#' b108_dist2drainage_exp.
hf <-"D:/GISbases/b108JaguarAlbersBuffer70/b108_dist2drainage_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J108trk))
#crs(hf)   
#compareCRS(hf,get_crs(J108trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J108/UTM_b108_dist2drainage_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J108/UTM_b108_dist2drainage_exp.tif";(test=raster(test)); x11(); plot(test)

#' J109
#' b109_dist2drainage_exp.
hf <-"D:/GISbases/b109JaguarAlbersBuffer70/b109_dist2drainage_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J109trk))
#crs(hf)   
#compareCRS(hf,get_crs(J109trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J109/UTM_b109_dist2drainage_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J109/UTM_b109_dist2drainage_exp.tif";(test=raster(test)); x11(); plot(test)

#' J110
#' b110_dist2drainage_exp.
hf <-"D:/GISbases/b110JaguarAlbersBuffer70/b110_dist2drainage_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J110trk))
#crs(hf)   
#compareCRS(hf,get_crs(J110trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J110/UTM_b110_dist2drainage_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J110/UTM_b110_dist2drainage_exp.tif";(test=raster(test)); x11(); plot(test)

#' J111
#' b111_dist2drainage_exp.
hf <-"D:/GISbases/b111JaguarAlbersBuffer70/b111_dist2drainage_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J111trk))
#crs(hf)   
#compareCRS(hf,get_crs(J111trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J111/UTM_b111_dist2drainage_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J111/UTM_b111_dist2drainage_exp.tif";(test=raster(test)); x11(); plot(test)

#' J112
#' b112_dist2drainage_exp.
hf <-"D:/GISbases/b112JaguarAlbersBuffer70/b112_dist2drainage_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J112trk))
#crs(hf)   
#compareCRS(hf,get_crs(J112trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J112/UTM_b112_dist2drainage_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J112/UTM_b112_dist2drainage_exp.tif";(test=raster(test)); x11(); plot(test)

#' J113
#' b113_dist2drainage_exp.
hf <-"D:/GISbases/b113JaguarAlbersBuffer70/b113_dist2drainage_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J113trk))
#crs(hf)   
#compareCRS(hf,get_crs(J113trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J113/UTM_b113_dist2drainage_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J113/UTM_b113_dist2drainage_exp.tif";(test=raster(test)); x11(); plot(test)

#' J114
#' b114_dist2drainage_exp.
hf <-"D:/GISbases/b114JaguarAlbersBuffer70/b114_dist2drainage_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J114trk))
#crs(hf)   
#compareCRS(hf,get_crs(J114trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J114/UTM_b114_dist2drainage_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J114/UTM_b114_dist2drainage_exp.tif";(test=raster(test)); x11(); plot(test)

#' J115
#' b115_dist2drainage_exp.
hf <-"D:/GISbases/b115JaguarAlbersBuffer70/b115_dist2drainage_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J115trk))
#crs(hf)   
#compareCRS(hf,get_crs(J115trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J115/UTM_b115_dist2drainage_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J115/UTM_b115_dist2drainage_exp.tif";(test=raster(test)); x11(); plot(test)

#' J116
#' b116_dist2drainage_exp.
hf <-"D:/GISbases/b116JaguarAlbersBuffer70/b116_dist2drainage_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J116trk))
#crs(hf)   
#compareCRS(hf,get_crs(J116trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J116/UTM_b116_dist2drainage_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J116/UTM_b116_dist2drainage_exp.tif";(test=raster(test)); x11(); plot(test)

#' J117
#' b117_dist2drainage_exp.
hf <-"D:/GISbases/b117JaguarAlbersBuffer70/b117_dist2drainage_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J117trk))
#crs(hf)   
#compareCRS(hf,get_crs(J117trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J117/UTM_b117_dist2drainage_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J117/UTM_b117_dist2drainage_exp.tif";(test=raster(test)); x11(); plot(test)

#####################################################################################################################################
# dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.

gc()
gc()
gc()

#' J1
#' b001_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif
hf <-"D:/GISbases/b001JaguarAlbersBuffer70/b001_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J1trk))
#crs(hf)   
#compareCRS(hf,get_crs(J1trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J1/UTM_b001_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J1/UTM_b001_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif";(test=raster(test)); x11(); plot(test)

#' J2
#' b002_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.
hf <-"D:/GISbases/b002JaguarAlbersBuffer70/b002_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J2trk))
#crs(hf)   
#compareCRS(hf,get_crs(J2trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J2/UTM_b002_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J2/UTM_b002_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif";(test=raster(test)); x11(); plot(test)

#' J3
#' b003_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.
hf <-"D:/GISbases/b003JaguarAlbersBuffer70/b003_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J3trk))
#crs(hf)   
#compareCRS(hf,get_crs(J3trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J3/UTM_b003_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J3/UTM_b003_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif";(test=raster(test)); x11(); plot(test)

#' J4
#' b004_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.
hf <-"D:/GISbases/b004JaguarAlbersBuffer70/b004_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J4trk))
#crs(hf)   
#compareCRS(hf,get_crs(J4trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J4/UTM_b004_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J4/UTM_b004_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif";(test=raster(test)); x11(); plot(test)

#' J5
#' b005_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.
hf <-"D:/GISbases/b005JaguarAlbersBuffer70/b005_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J5trk))
#crs(hf)   
#compareCRS(hf,get_crs(J5trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J5/UTM_b005_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J5/UTM_b005_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif";(test=raster(test)); x11(); plot(test)

#' J6
#' b006_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.
hf <-"D:/GISbases/b006JaguarAlbersBuffer70/b006_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J6trk))
#crs(hf)   
#compareCRS(hf,get_crs(J6trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J6/UTM_b006_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J6/UTM_b006_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif";(test=raster(test)); x11(); plot(test)

#' J7
#' b007_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.
hf <-"D:/GISbases/b007JaguarAlbersBuffer70/b007_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J7trk))
#crs(hf)   
#compareCRS(hf,get_crs(J7trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J7/UTM_b007_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J7/UTM_b007_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif";(test=raster(test)); x11(); plot(test)

#' J8
#' b008_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.
hf <-"D:/GISbases/b008JaguarAlbersBuffer70/b008_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J8trk))
#crs(hf)   
#compareCRS(hf,get_crs(J8trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J8/UTM_b008_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J8/UTM_b008_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif";(test=raster(test)); x11(); plot(test)

#' J9
#' b009_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.
hf <-"D:/GISbases/b009JaguarAlbersBuffer70/b009_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J9trk))
#crs(hf)   
#compareCRS(hf,get_crs(J9trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J9/UTM_b009_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J9/UTM_b009_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif";(test=raster(test)); x11(); plot(test)

#' J10
#' b010_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.
hf <-"D:/GISbases/b010JaguarAlbersBuffer70/b010_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J10trk))
#crs(hf)   
#compareCRS(hf,get_crs(J10trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J10/UTM_b010_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J10/UTM_b010_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif";(test=raster(test)); x11(); plot(test)

gc()
gc()
gc()


#' J11
#' b011_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.
hf <-"D:/GISbases/b011JaguarAlbersBuffer70/b011_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J11trk))
#crs(hf)   
#compareCRS(hf,get_crs(J11trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J11/UTM_b011_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J11/UTM_b011_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif";(test=raster(test)); x11(); plot(test)

#' J12
#' b012_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.
hf <-"D:/GISbases/b012JaguarAlbersBuffer70/b012_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J12trk))
#crs(hf)   
#compareCRS(hf,get_crs(J12trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J12/UTM_b012_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J12/UTM_b012_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif";(test=raster(test)); x11(); plot(test)

#' J13
#' b013_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.
hf <-"D:/GISbases/b013JaguarAlbersBuffer70/b013_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J13trk))
#crs(hf)   
#compareCRS(hf,get_crs(J13trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J13/UTM_b013_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J13/UTM_b013_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif";(test=raster(test)); x11(); plot(test)

#' J14
#' b014_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.
hf <-"D:/GISbases/b014JaguarAlbersBuffer70/b014_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J14trk))
#crs(hf)   
#compareCRS(hf,get_crs(J14trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J14/UTM_b014_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J14/UTM_b014_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif";(test=raster(test)); x11(); plot(test)

#' J15
#' b015_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.
hf <-"D:/GISbases/b015JaguarAlbersBuffer70/b015_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J15trk))
#crs(hf)   
#compareCRS(hf,get_crs(J15trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J15/UTM_b015_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J15/UTM_b015_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif";(test=raster(test)); x11(); plot(test)

#' J16
#' b016_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.
hf <-"D:/GISbases/b016JaguarAlbersBuffer70/b016_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J16trk))
#crs(hf)   
#compareCRS(hf,get_crs(J16trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J16/UTM_b016_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J16/UTM_b016_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif";(test=raster(test)); x11(); plot(test)

#' J17
#' b017_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.
hf <-"D:/GISbases/b017JaguarAlbersBuffer70/b017_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J17trk))
#crs(hf)   
#compareCRS(hf,get_crs(J17trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J17/UTM_b017_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J17/UTM_b017_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif";(test=raster(test)); x11(); plot(test)

#' J18
#' b018_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.
hf <-"D:/GISbases/b018JaguarAlbersBuffer70/b018_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J18trk))
#crs(hf)   
#compareCRS(hf,get_crs(J18trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J18/UTM_b018_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J18/UTM_b018_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif";(test=raster(test)); x11(); plot(test)

#' J19
#' b019_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.
hf <-"D:/GISbases/b019JaguarAlbersBuffer70/b019_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J19trk))
#crs(hf)   
#compareCRS(hf,get_crs(J19trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J19/UTM_b019_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J19/UTM_b019_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif";(test=raster(test)); x11(); plot(test)

#' J20
#' b020_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.
hf <-"D:/GISbases/b020JaguarAlbersBuffer70/b020_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J20trk))
#crs(hf)   
#compareCRS(hf,get_crs(J20trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J20/UTM_b020_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J20/UTM_b020_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif";(test=raster(test)); x11(); plot(test)

#' J21
#' b021_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.
hf <-"D:/GISbases/b021JaguarAlbersBuffer70/b021_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J21trk))
#crs(hf)   
#compareCRS(hf,get_crs(J21trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J21/UTM_b021_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J21/UTM_b021_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif";(test=raster(test)); x11(); plot(test)

#' J22
#' b022_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.
hf <-"D:/GISbases/b022JaguarAlbersBuffer70/b022_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J22trk))
#crs(hf)   
#compareCRS(hf,get_crs(J22trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J22/UTM_b022_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J22/UTM_b022_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif";(test=raster(test)); x11(); plot(test)

#' J23
#' b023_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.
hf <-"D:/GISbases/b023JaguarAlbersBuffer70/b023_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J23trk))
#crs(hf)   
#compareCRS(hf,get_crs(J23trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J23/UTM_b023_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J23/UTM_b023_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif";(test=raster(test)); x11(); plot(test)

#' J24
#' b024_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.
hf <-"D:/GISbases/b024JaguarAlbersBuffer70/b024_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J24trk))
#crs(hf)   
#compareCRS(hf,get_crs(J24trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J24/UTM_b024_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J24/UTM_b024_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif";(test=raster(test)); x11(); plot(test)

#' J25
#' b025_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.
hf <-"D:/GISbases/b025JaguarAlbersBuffer70/b025_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J25trk))
#crs(hf)   
#compareCRS(hf,get_crs(J25trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J25/UTM_b025_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J25/UTM_b025_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif";(test=raster(test)); x11(); plot(test)

#' J26
#' b026_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.
hf <-"D:/GISbases/b026JaguarAlbersBuffer70/b026_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J26trk))
#crs(hf)   
#compareCRS(hf,get_crs(J26trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J26/UTM_b026_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J26/UTM_b026_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif";(test=raster(test)); x11(); plot(test)

#' J27
#' b027_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.
hf <-"D:/GISbases/b027JaguarAlbersBuffer70/b027_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J27trk))
#crs(hf)   
#compareCRS(hf,get_crs(J27trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J27/UTM_b027_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J27/UTM_b027_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif";(test=raster(test)); x11(); plot(test)

#' J28
#' b028_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.
hf <-"D:/GISbases/b028JaguarAlbersBuffer70/b028_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J28trk))
#crs(hf)   
#compareCRS(hf,get_crs(J28trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J28/UTM_b028_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J28/UTM_b028_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif";(test=raster(test)); x11(); plot(test)

#' J29
#' b029_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.
hf <-"D:/GISbases/b029JaguarAlbersBuffer70/b029_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J29trk))
#crs(hf)   
#compareCRS(hf,get_crs(J29trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J29/UTM_b029_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J29/UTM_b029_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif";(test=raster(test)); x11(); plot(test)

#' J30
#' b030_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.
hf <-"D:/GISbases/b030JaguarAlbersBuffer70/b030_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J30trk))
#crs(hf)   
#compareCRS(hf,get_crs(J30trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J30/UTM_b030_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J30/UTM_b030_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif";(test=raster(test)); x11(); plot(test)

#' J30
#' b030_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.
hf <-"D:/GISbases/b030JaguarAlbersBuffer70/b030_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J30trk))
#crs(hf)   
#compareCRS(hf,get_crs(J30trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J30/UTM_b030_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J30/UTM_b030_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif";(test=raster(test)); x11(); plot(test)

#' J31
#' b031_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.
hf <-"D:/GISbases/b031JaguarAlbersBuffer70/b031_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J31trk))
#crs(hf)   
#compareCRS(hf,get_crs(J31trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J31/UTM_b031_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J31/UTM_b031_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif";(test=raster(test)); x11(); plot(test)

#' J32
#' b032_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.
hf <-"D:/GISbases/b032JaguarAlbersBuffer70/b032_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J32trk))
#crs(hf)   
#compareCRS(hf,get_crs(J32trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J32/UTM_b032_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J32/UTM_b032_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif";(test=raster(test)); x11(); plot(test)

#' J33
#' b033_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.
hf <-"D:/GISbases/b033JaguarAlbersBuffer70/b033_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J33trk))
#crs(hf)   
#compareCRS(hf,get_crs(J33trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J33/UTM_b033_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J33/UTM_b033_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif";(test=raster(test)); x11(); plot(test)

#' J34
#' b034_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.
hf <-"D:/GISbases/b034JaguarAlbersBuffer70/b034_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J34trk))
#crs(hf)   
#compareCRS(hf,get_crs(J34trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J34/UTM_b034_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J34/UTM_b034_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif";(test=raster(test)); x11(); plot(test)

#' J35
#' b035_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.
hf <-"D:/GISbases/b035JaguarAlbersBuffer70/b035_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J35trk))
#crs(hf)   
#compareCRS(hf,get_crs(J35trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J35/UTM_b035_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J35/UTM_b035_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif";(test=raster(test)); x11(); plot(test)

#' J36
#' b036_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.
hf <-"D:/GISbases/b036JaguarAlbersBuffer70/b036_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J36trk))
#crs(hf)   
#compareCRS(hf,get_crs(J36trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J36/UTM_b036_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J36/UTM_b036_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif";(test=raster(test)); x11(); plot(test)

#' J37
#' b037_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.
hf <-"D:/GISbases/b037JaguarAlbersBuffer70/b037_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J37trk))
#crs(hf)   
#compareCRS(hf,get_crs(J37trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J37/UTM_b037_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J37/UTM_b037_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif";(test=raster(test)); x11(); plot(test)

#' J38
#' b038_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.
hf <-"D:/GISbases/b038JaguarAlbersBuffer70/b038_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J38trk))
#crs(hf)   
#compareCRS(hf,get_crs(J38trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J38/UTM_b038_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J38/UTM_b038_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif";(test=raster(test)); x11(); plot(test)

#' J39
#' b039_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.
hf <-"D:/GISbases/b039JaguarAlbersBuffer70/b039_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J39trk))
#crs(hf)   
#compareCRS(hf,get_crs(J39trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J39/UTM_b039_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J39/UTM_b039_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif";(test=raster(test)); x11(); plot(test)

#' J40
#' b040_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.
hf <-"D:/GISbases/b040JaguarAlbersBuffer70/b040_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J40trk))
#crs(hf)   
#compareCRS(hf,get_crs(J40trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J40/UTM_b040_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J40/UTM_b040_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif";(test=raster(test)); x11(); plot(test)

#' J41
#' b041_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.
hf <-"D:/GISbases/b041JaguarAlbersBuffer70/b041_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J41trk))
#crs(hf)   
#compareCRS(hf,get_crs(J41trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J41/UTM_b041_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J41/UTM_b041_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif";(test=raster(test)); x11(); plot(test)

#' J42
#' b042_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.
hf <-"D:/GISbases/b042JaguarAlbersBuffer70/b042_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J42trk))
#crs(hf)   
#compareCRS(hf,get_crs(J42trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J42/UTM_b042_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J42/UTM_b042_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif";(test=raster(test)); x11(); plot(test)

gc()
gc()
gc()

#' J43
#' b043_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.
hf <-"D:/GISbases/b043JaguarAlbersBuffer70/b043_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J43trk))
#crs(hf)   
#compareCRS(hf,get_crs(J43trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J43/UTM_b043_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J43/UTM_b043_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif";(test=raster(test)); x11(); plot(test)

#' J44
#' b044_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.
hf <-"D:/GISbases/b044JaguarAlbersBuffer70/b044_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J44trk))
#crs(hf)   
#compareCRS(hf,get_crs(J44trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J44/UTM_b044_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J44/UTM_b044_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif";(test=raster(test)); x11(); plot(test)

#' J45
#' b045_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.
hf <-"D:/GISbases/b045JaguarAlbersBuffer70/b045_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J45trk))
#crs(hf)   
#compareCRS(hf,get_crs(J45trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J45/UTM_b045_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J45/UTM_b045_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif";(test=raster(test)); x11(); plot(test)

#' J46
#' b046_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.
hf <-"D:/GISbases/b046JaguarAlbersBuffer70/b046_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J46trk))
#crs(hf)   
#compareCRS(hf,get_crs(J46trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J46/UTM_b046_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J46/UTM_b046_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif";(test=raster(test)); x11(); plot(test)

#' J47
#' b047_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.
hf <-"D:/GISbases/b047JaguarAlbersBuffer70/b047_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J47trk))
#crs(hf)   
#compareCRS(hf,get_crs(J47trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J47/UTM_b047_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J47/UTM_b047_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif";(test=raster(test)); x11(); plot(test)

#' J48
#' b048_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.
hf <-"D:/GISbases/b048JaguarAlbersBuffer70/b048_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J48trk))
#crs(hf)   
#compareCRS(hf,get_crs(J48trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J48/UTM_b048_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J48/UTM_b048_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif";(test=raster(test)); x11(); plot(test)

#' J49
#' b049_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.
hf <-"D:/GISbases/b049JaguarAlbersBuffer70/b049_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J49trk))
#crs(hf)   
#compareCRS(hf,get_crs(J49trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J49/UTM_b049_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J49/UTM_b049_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif";(test=raster(test)); x11(); plot(test)

#' J50
#' b050_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.
hf <-"D:/GISbases/b050JaguarAlbersBuffer70/b050_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J50trk))
#crs(hf)   
#compareCRS(hf,get_crs(J50trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J50/UTM_b050_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J50/UTM_b050_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif";(test=raster(test)); x11(); plot(test)

#' J51
#' b051_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.
hf <-"D:/GISbases/b051JaguarAlbersBuffer70/b051_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J51trk))
#crs(hf)   
#compareCRS(hf,get_crs(J51trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J51/UTM_b051_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J51/UTM_b051_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif";(test=raster(test)); x11(); plot(test)

#' J52
#' b052_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.
hf <-"D:/GISbases/b052JaguarAlbersBuffer70/b052_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J52trk))
#crs(hf)   
#compareCRS(hf,get_crs(J52trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J52/UTM_b052_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J52/UTM_b052_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif";(test=raster(test)); x11(); plot(test)

#' J53
#' b053_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.
hf <-"D:/GISbases/b053JaguarAlbersBuffer70/b053_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J53trk))
#crs(hf)   
#compareCRS(hf,get_crs(J53trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J53/UTM_b053_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J53/UTM_b053_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif";(test=raster(test)); x11(); plot(test)

#' J54
#' b054_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.
hf <-"D:/GISbases/b054JaguarAlbersBuffer70/b054_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J54trk))
#crs(hf)   
#compareCRS(hf,get_crs(J54trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J54/UTM_b054_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J54/UTM_b054_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif";(test=raster(test)); x11(); plot(test)

#' J55
#' b055_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.
hf <-"D:/GISbases/b055JaguarAlbersBuffer70/b055_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J55trk))
#crs(hf)   
#compareCRS(hf,get_crs(J55trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J55/UTM_b055_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J55/UTM_b055_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif";(test=raster(test)); x11(); plot(test)

#' J56
#' b056_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.
hf <-"D:/GISbases/b056JaguarAlbersBuffer70/b056_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J56trk))
#crs(hf)   
#compareCRS(hf,get_crs(J56trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J56/UTM_b056_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J56/UTM_b056_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif";(test=raster(test)); x11(); plot(test)

#' J57
#' b057_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.
hf <-"D:/GISbases/b057JaguarAlbersBuffer70/b057_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J57trk))
#crs(hf)   
#compareCRS(hf,get_crs(J57trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J57/UTM_b057_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J57/UTM_b057_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif";(test=raster(test)); x11(); plot(test)

#' J58
#' b058_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.
hf <-"D:/GISbases/b058JaguarAlbersBuffer70/b058_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J58trk))
#crs(hf)   
#compareCRS(hf,get_crs(J58trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J58/UTM_b058_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J58/UTM_b058_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif";(test=raster(test)); x11(); plot(test)

#' J59
#' b059_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.
hf <-"D:/GISbases/b059JaguarAlbersBuffer70/b059_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J59trk))
#crs(hf)   
#compareCRS(hf,get_crs(J59trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J59/UTM_b059_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J59/UTM_b059_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif";(test=raster(test)); x11(); plot(test)

#' J60
#' b060_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.
hf <-"D:/GISbases/b060JaguarAlbersBuffer70/b060_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J60trk))
#crs(hf)   
#compareCRS(hf,get_crs(J60trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J60/UTM_b060_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J60/UTM_b060_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif";(test=raster(test)); x11(); plot(test)

#' J61
#' b061_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.
hf <-"D:/GISbases/b061JaguarAlbersBuffer70/b061_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J61trk))
#crs(hf)   
#compareCRS(hf,get_crs(J61trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J61/UTM_b061_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J61/UTM_b061_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif";(test=raster(test)); x11(); plot(test)

#' J62
#' b062_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.
hf <-"D:/GISbases/b062JaguarAlbersBuffer70/b062_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J62trk))
#crs(hf)   
#compareCRS(hf,get_crs(J62trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J62/UTM_b062_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J62/UTM_b062_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif";(test=raster(test)); x11(); plot(test)

#' J63
#' b063_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.
hf <-"D:/GISbases/b063JaguarAlbersBuffer70/b063_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J63trk))
#crs(hf)   
#compareCRS(hf,get_crs(J63trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J63/UTM_b063_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J63/UTM_b063_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif";(test=raster(test)); x11(); plot(test)

gc()
gc()
gc()


#' J64
#' b064_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.
hf <-"D:/GISbases/b064JaguarAlbersBuffer70/b064_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J64trk))
#crs(hf)   
#compareCRS(hf,get_crs(J64trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J64/UTM_b064_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J64/UTM_b064_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif";(test=raster(test)); x11(); plot(test)

#' J65
#' b065_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.
hf <-"D:/GISbases/b065JaguarAlbersBuffer70/b065_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J65trk))
#crs(hf)   
#compareCRS(hf,get_crs(J65trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J65/UTM_b065_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J65/UTM_b065_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif";(test=raster(test)); x11(); plot(test)

#' J66
#' b066_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.
hf <-"D:/GISbases/b066JaguarAlbersBuffer70/b066_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J66trk))
#crs(hf)   
#compareCRS(hf,get_crs(J66trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J66/UTM_b066_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J66/UTM_b066_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif";(test=raster(test)); x11(); plot(test)

#' J67
#' b067_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.
hf <-"D:/GISbases/b067JaguarAlbersBuffer70/b067_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J67trk))
#crs(hf)   
#compareCRS(hf,get_crs(J67trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J67/UTM_b067_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J67/UTM_b067_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif";(test=raster(test)); x11(); plot(test)

#' J68
#' b068_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.
hf <-"D:/GISbases/b068JaguarAlbersBuffer70/b068_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J68trk))
#crs(hf)   
#compareCRS(hf,get_crs(J68trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J68/UTM_b068_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J68/UTM_b068_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif";(test=raster(test)); x11(); plot(test)

#' J69
#' b069_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.
hf <-"D:/GISbases/b069JaguarAlbersBuffer70/b069_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J69trk))
#crs(hf)   
#compareCRS(hf,get_crs(J69trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J69/UTM_b069_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J69/UTM_b069_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif";(test=raster(test)); x11(); plot(test)

gc()
gc()
gc()


#' J70
#' b070_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.
hf <-"D:/GISbases/b070JaguarAlbersBuffer70/b070_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J70trk))
#crs(hf)   
#compareCRS(hf,get_crs(J70trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J70/UTM_b070_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J70/UTM_b070_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif";(test=raster(test)); x11(); plot(test)

#' J71
#' b071_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.
hf <-"D:/GISbases/b071JaguarAlbersBuffer70/b071_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J71trk))
#crs(hf)   
#compareCRS(hf,get_crs(J71trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J71/UTM_b071_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J71/UTM_b071_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif";(test=raster(test)); x11(); plot(test)

#' J72
#' b072_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.
hf <-"D:/GISbases/b072JaguarAlbersBuffer70/b072_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J72trk))
#crs(hf)   
#compareCRS(hf,get_crs(J72trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J72/UTM_b072_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J72/UTM_b072_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif";(test=raster(test)); x11(); plot(test)

#' J73
#' b073_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.
hf <-"D:/GISbases/b073JaguarAlbersBuffer70/b073_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J73trk))
#crs(hf)   
#compareCRS(hf,get_crs(J73trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J73/UTM_b073_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J73/UTM_b073_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif";(test=raster(test)); x11(); plot(test)

#' J74
#' b074_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.
hf <-"D:/GISbases/b074JaguarAlbersBuffer70/b074_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J74trk))
#crs(hf)   
#compareCRS(hf,get_crs(J74trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J74/UTM_b074_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J74/UTM_b074_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif";(test=raster(test)); x11(); plot(test)

#' J75
#' b075_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.
hf <-"D:/GISbases/b075JaguarAlbersBuffer70/b075_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J75trk))
#crs(hf)   
#compareCRS(hf,get_crs(J75trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J75/UTM_b075_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J75/UTM_b075_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif";(test=raster(test)); x11(); plot(test)

#' J76
#' b076_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.
hf <-"D:/GISbases/b076JaguarAlbersBuffer70/b076_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J76trk))
#crs(hf)   
#compareCRS(hf,get_crs(J76trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J76/UTM_b076_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J76/UTM_b076_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif";(test=raster(test)); x11(); plot(test)

#' J77
#' b077_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.
hf <-"D:/GISbases/b077JaguarAlbersBuffer70/b077_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J77trk))
#crs(hf)   
#compareCRS(hf,get_crs(J77trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J77/UTM_b077_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J77/UTM_b077_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif";(test=raster(test)); x11(); plot(test)

#' J78
#' b078_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.
hf <-"D:/GISbases/b078JaguarAlbersBuffer70/b078_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J78trk))
#crs(hf)   
#compareCRS(hf,get_crs(J78trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J78/UTM_b078_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J78/UTM_b078_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif";(test=raster(test)); x11(); plot(test)

#' J79
#' b079_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.
hf <-"D:/GISbases/b079JaguarAlbersBuffer70/b079_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J79trk))
#crs(hf)   
#compareCRS(hf,get_crs(J79trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J79/UTM_b079_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J79/UTM_b079_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif";(test=raster(test)); x11(); plot(test)

#' J80
#' b080_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.
hf <-"D:/GISbases/b080JaguarAlbersBuffer70/b080_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J80trk))
#crs(hf)   
#compareCRS(hf,get_crs(J80trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J80/UTM_b080_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J80/UTM_b080_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif";(test=raster(test)); x11(); plot(test)

#' J81
#' b081_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.
hf <-"D:/GISbases/b081JaguarAlbersBuffer70/b081_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J81trk))
#crs(hf)   
#compareCRS(hf,get_crs(J81trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J81/UTM_b081_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J81/UTM_b081_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif";(test=raster(test)); x11(); plot(test)

#' J82
#' b082_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.
hf <-"D:/GISbases/b082JaguarAlbersBuffer70/b082_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J82trk))
#crs(hf)   
#compareCRS(hf,get_crs(J82trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J82/UTM_b082_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J82/UTM_b082_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif";(test=raster(test)); x11(); plot(test)

#' J83
#' b083_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.
hf <-"D:/GISbases/b083JaguarAlbersBuffer70/b083_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J83trk))
#crs(hf)   
#compareCRS(hf,get_crs(J83trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J83/UTM_b083_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J83/UTM_b083_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif";(test=raster(test)); x11(); plot(test)

#' J84
#' b084_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.
hf <-"D:/GISbases/b084JaguarAlbersBuffer70/b084_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J84trk))
#crs(hf)   
#compareCRS(hf,get_crs(J84trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J84/UTM_b084_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J84/UTM_b084_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif";(test=raster(test)); x11(); plot(test)

#' J85
#' b085_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.
hf <-"D:/GISbases/b085JaguarAlbersBuffer70/b085_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J85trk))
#crs(hf)   
#compareCRS(hf,get_crs(J85trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J85/UTM_b085_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J85/UTM_b085_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif";(test=raster(test)); x11(); plot(test)

#' J86
#' b086_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.
hf <-"D:/GISbases/b086JaguarAlbersBuffer70/b086_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J86trk))
#crs(hf)   
#compareCRS(hf,get_crs(J86trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J86/UTM_b086_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J86/UTM_b086_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif";(test=raster(test)); x11(); plot(test)

#' J87
#' b087_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.
hf <-"D:/GISbases/b087JaguarAlbersBuffer70/b087_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J87trk))
#crs(hf)   
#compareCRS(hf,get_crs(J87trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J87/UTM_b087_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J87/UTM_b087_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif";(test=raster(test)); x11(); plot(test)

#' J88
#' b088_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.
hf <-"D:/GISbases/b088JaguarAlbersBuffer70/b088_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J88trk))
#crs(hf)   
#compareCRS(hf,get_crs(J88trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J88/UTM_b088_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J88/UTM_b088_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif";(test=raster(test)); x11(); plot(test)

#' J89
#' b089_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.
hf <-"D:/GISbases/b089JaguarAlbersBuffer70/b089_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J89trk))
#crs(hf)   
#compareCRS(hf,get_crs(J89trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J89/UTM_b089_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J89/UTM_b089_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif";(test=raster(test)); x11(); plot(test)

#' J90
#' b090_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.
hf <-"D:/GISbases/b090JaguarAlbersBuffer70/b090_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J90trk))
#crs(hf)   
#compareCRS(hf,get_crs(J90trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J90/UTM_b090_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J90/UTM_b090_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif";(test=raster(test)); x11(); plot(test)

#' J91
#' b091_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.
hf <-"D:/GISbases/b091JaguarAlbersBuffer70/b091_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J91trk))
#crs(hf)   
#compareCRS(hf,get_crs(J91trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J91/UTM_b091_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J91/UTM_b091_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif";(test=raster(test)); x11(); plot(test)

#' J92
#' b092_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.
hf <-"D:/GISbases/b092JaguarAlbersBuffer70/b092_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J92trk))
#crs(hf)   
#compareCRS(hf,get_crs(J92trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J92/UTM_b092_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J92/UTM_b092_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif";(test=raster(test)); x11(); plot(test)

#' J93
#' b093_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.
hf <-"D:/GISbases/b093JaguarAlbersBuffer70/b093_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J93trk))
#crs(hf)   
#compareCRS(hf,get_crs(J93trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J93/UTM_b093_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J93/UTM_b093_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif";(test=raster(test)); x11(); plot(test)

#' J94
#' b094_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.
hf <-"D:/GISbases/b094JaguarAlbersBuffer70/b094_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J94trk))
#crs(hf)   
#compareCRS(hf,get_crs(J94trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J94/UTM_b094_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J94/UTM_b094_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif";(test=raster(test)); x11(); plot(test)

#' J95
#' b095_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.
hf <-"D:/GISbases/b095JaguarAlbersBuffer70/b095_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J95trk))
#crs(hf)   
#compareCRS(hf,get_crs(J95trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J95/UTM_b095_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J95/UTM_b095_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif";(test=raster(test)); x11(); plot(test)

#' J96
#' b096_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.
hf <-"D:/GISbases/b096JaguarAlbersBuffer70/b096_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J96trk))
#crs(hf)   
#compareCRS(hf,get_crs(J96trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J96/UTM_b096_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J96/UTM_b096_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif";(test=raster(test)); x11(); plot(test)

#' J97
#' b097_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.
hf <-"D:/GISbases/b097JaguarAlbersBuffer70/b097_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J97trk))
#crs(hf)   
#compareCRS(hf,get_crs(J97trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J97/UTM_b097_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J97/UTM_b097_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif";(test=raster(test)); x11(); plot(test)

#' J98
#' b098_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.
hf <-"D:/GISbases/b098JaguarAlbersBuffer70/b098_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J98trk))
#crs(hf)   
#compareCRS(hf,get_crs(J98trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J98/UTM_b098_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J98/UTM_b098_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif";(test=raster(test)); x11(); plot(test)

#' J99
#' b099_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.
hf <-"D:/GISbases/b099JaguarAlbersBuffer70/b099_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J99trk))
#crs(hf)   
#compareCRS(hf,get_crs(J99trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J99/UTM_b099_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J99/UTM_b099_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif";(test=raster(test)); x11(); plot(test)

#' J100
#' b100_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.
hf <-"D:/GISbases/b100JaguarAlbersBuffer70/b100_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J100trk))
#crs(hf)   
#compareCRS(hf,get_crs(J100trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J100/UTM_b100_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J100/UTM_b010_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif";(test=raster(test)); x11(); plot(test)

#' J101
#' b101_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.
hf <-"D:/GISbases/b101JaguarAlbersBuffer70/b101_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J101trk))
#crs(hf)   
#compareCRS(hf,get_crs(J101trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J101/UTM_b101_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J101/UTM_b101_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif";(test=raster(test)); x11(); plot(test)

#' J102
#' b102_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.
hf <-"D:/GISbases/b102JaguarAlbersBuffer70/b102_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J102trk))
#crs(hf)   
#compareCRS(hf,get_crs(J102trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J102/UTM_b102_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J102/UTM_b102_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif";(test=raster(test)); x11(); plot(test)

#' J103
#' b103_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.
hf <-"D:/GISbases/b103JaguarAlbersBuffer70/b103_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J103trk))
#crs(hf)   
#compareCRS(hf,get_crs(J103trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J103/UTM_b103_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J103/UTM_b103_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif";(test=raster(test)); x11(); plot(test)

#' J104
#' b104_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.
hf <-"D:/GISbases/b104JaguarAlbersBuffer70/b104_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J104trk))
#crs(hf)   
#compareCRS(hf,get_crs(J104trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J104/UTM_b104_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J104/UTM_b104_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif";(test=raster(test)); x11(); plot(test)

#' J105
#' b105_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.
hf <-"D:/GISbases/b105JaguarAlbersBuffer70/b105_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J105trk))
#crs(hf)   
#compareCRS(hf,get_crs(J105trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J105/UTM_b105_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J105/UTM_b105_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif";(test=raster(test)); x11(); plot(test)

#' J106
#' b106_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.
hf <-"D:/GISbases/b106JaguarAlbersBuffer70/b106_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J106trk))
#crs(hf)   
#compareCRS(hf,get_crs(J106trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J106/UTM_b106_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J106/UTM_b106_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif";(test=raster(test)); x11(); plot(test)

#' J107
#' b107_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.
hf <-"D:/GISbases/b107JaguarAlbersBuffer70/b107_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J107trk))
#crs(hf)   
#compareCRS(hf,get_crs(J107trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J107/UTM_b107_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J107/UTM_b107_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif";(test=raster(test)); x11(); plot(test)

#' J108
#' b108_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.
hf <-"D:/GISbases/b108JaguarAlbersBuffer70/b108_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J108trk))
#crs(hf)   
#compareCRS(hf,get_crs(J108trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J108/UTM_b108_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J108/UTM_b108_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif";(test=raster(test)); x11(); plot(test)

#' J109
#' b109_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.
hf <-"D:/GISbases/b109JaguarAlbersBuffer70/b109_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J109trk))
#crs(hf)   
#compareCRS(hf,get_crs(J109trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J109/UTM_b109_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J109/UTM_b109_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif";(test=raster(test)); x11(); plot(test)

#' J110
#' b110_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.
hf <-"D:/GISbases/b110JaguarAlbersBuffer70/b110_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J110trk))
#crs(hf)   
#compareCRS(hf,get_crs(J110trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J110/UTM_b110_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J110/UTM_b110_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif";(test=raster(test)); x11(); plot(test)

#' J111
#' b111_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.
hf <-"D:/GISbases/b111JaguarAlbersBuffer70/b111_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J111trk))
#crs(hf)   
#compareCRS(hf,get_crs(J111trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J111/UTM_b111_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J111/UTM_b111_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif";(test=raster(test)); x11(); plot(test)

#' J112
#' b112_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.
hf <-"D:/GISbases/b112JaguarAlbersBuffer70/b112_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J112trk))
#crs(hf)   
#compareCRS(hf,get_crs(J112trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J112/UTM_b112_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J112/UTM_b112_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif";(test=raster(test)); x11(); plot(test)

#' J113
#' b113_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.
hf <-"D:/GISbases/b113JaguarAlbersBuffer70/b113_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J113trk))
#crs(hf)   
#compareCRS(hf,get_crs(J113trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J113/UTM_b113_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J113/UTM_b113_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif";(test=raster(test)); x11(); plot(test)

#' J114
#' b114_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.
hf <-"D:/GISbases/b114JaguarAlbersBuffer70/b114_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J114trk))
#crs(hf)   
#compareCRS(hf,get_crs(J114trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J114/UTM_b114_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J114/UTM_b114_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif";(test=raster(test)); x11(); plot(test)

#' J115
#' b115_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.
hf <-"D:/GISbases/b115JaguarAlbersBuffer70/b115_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J115trk))
#crs(hf)   
#compareCRS(hf,get_crs(J115trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J115/UTM_b115_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J115/UTM_b115_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif";(test=raster(test)); x11(); plot(test)

#' J116
#' b116_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.
hf <-"D:/GISbases/b116JaguarAlbersBuffer70/b116_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J116trk))
#crs(hf)   
#compareCRS(hf,get_crs(J116trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J116/UTM_b116_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J116/UTM_b116_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif";(test=raster(test)); x11(); plot(test)

#' J117
#' b117_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.
hf <-"D:/GISbases/b117JaguarAlbersBuffer70/b117_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J117trk))
#crs(hf)   
#compareCRS(hf,get_crs(J117trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J117/UTM_b117_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J117/UTM_b117_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif";(test=raster(test)); x11(); plot(test)

######################################################################################################################################
# Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.
gc()
gc()
gc()

#' J1
#' b001_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif
hf <-"D:/GISbases/b001JaguarAlbersBuffer70/b001_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J1trk))
#crs(hf)   
#compareCRS(hf,get_crs(J1trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J1/UTM_b001_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J1/UTM_b001_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J2
#' b002_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b002JaguarAlbersBuffer70/b002_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J2trk))
#crs(hf)   
#compareCRS(hf,get_crs(J2trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J2/UTM_b002_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J2/UTM_b002_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J3
#' b003_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b003JaguarAlbersBuffer70/b003_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J3trk))
#crs(hf)   
#compareCRS(hf,get_crs(J3trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J3/UTM_b003_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J3/UTM_b003_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J4
#' b004_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b004JaguarAlbersBuffer70/b004_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J4trk))
#crs(hf)   
#compareCRS(hf,get_crs(J4trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J4/UTM_b004_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J4/UTM_b004_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J5
#' b005_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b005JaguarAlbersBuffer70/b005_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J5trk))
#crs(hf)   
#compareCRS(hf,get_crs(J5trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J5/UTM_b005_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J5/UTM_b005_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J6
#' b006_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b006JaguarAlbersBuffer70/b006_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J6trk))
#crs(hf)   
#compareCRS(hf,get_crs(J6trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J6/UTM_b006_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J6/UTM_b006_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J7
#' b007_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b007JaguarAlbersBuffer70/b007_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J7trk))
#crs(hf)   
#compareCRS(hf,get_crs(J7trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J7/UTM_b007_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J7/UTM_b007_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J8
#' b008_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b008JaguarAlbersBuffer70/b008_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J8trk))
#crs(hf)   
#compareCRS(hf,get_crs(J8trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J8/UTM_b008_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J8/UTM_b008_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J9
#' b009_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b009JaguarAlbersBuffer70/b009_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J9trk))
#crs(hf)   
#compareCRS(hf,get_crs(J9trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J9/UTM_b009_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J9/UTM_b009_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J10
#' b010_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b010JaguarAlbersBuffer70/b010_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J10trk))
#crs(hf)   
#compareCRS(hf,get_crs(J10trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J10/UTM_b010_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J10/UTM_b010_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

gc()
gc()
gc()


#' J11
#' b011_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b011JaguarAlbersBuffer70/b011_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J11trk))
#crs(hf)   
#compareCRS(hf,get_crs(J11trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J11/UTM_b011_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J11/UTM_b011_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J12
#' b012_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b012JaguarAlbersBuffer70/b012_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J12trk))
#crs(hf)   
#compareCRS(hf,get_crs(J12trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J12/UTM_b012_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J12/UTM_b012_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J13
#' b013_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b013JaguarAlbersBuffer70/b013_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J13trk))
#crs(hf)   
#compareCRS(hf,get_crs(J13trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J13/UTM_b013_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J13/UTM_b013_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J14
#' b014_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b014JaguarAlbersBuffer70/b014_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J14trk))
#crs(hf)   
#compareCRS(hf,get_crs(J14trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J14/UTM_b014_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J14/UTM_b014_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J15
#' b015_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b015JaguarAlbersBuffer70/b015_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J15trk))
#crs(hf)   
#compareCRS(hf,get_crs(J15trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J15/UTM_b015_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J15/UTM_b015_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J16
#' b016_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b016JaguarAlbersBuffer70/b016_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J16trk))
#crs(hf)   
#compareCRS(hf,get_crs(J16trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J16/UTM_b016_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J16/UTM_b016_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J17
#' b017_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b017JaguarAlbersBuffer70/b017_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J17trk))
#crs(hf)   
#compareCRS(hf,get_crs(J17trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J17/UTM_b017_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J17/UTM_b017_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J18
#' b018_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b018JaguarAlbersBuffer70/b018_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J18trk))
#crs(hf)   
#compareCRS(hf,get_crs(J18trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J18/UTM_b018_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J18/UTM_b018_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J19
#' b019_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b019JaguarAlbersBuffer70/b019_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J19trk))
#crs(hf)   
#compareCRS(hf,get_crs(J19trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J19/UTM_b019_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J19/UTM_b019_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J20
#' b020_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b020JaguarAlbersBuffer70/b020_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J20trk))
#crs(hf)   
#compareCRS(hf,get_crs(J20trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J20/UTM_b020_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J20/UTM_b020_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J21
#' b021_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b021JaguarAlbersBuffer70/b021_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J21trk))
#crs(hf)   
#compareCRS(hf,get_crs(J21trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J21/UTM_b021_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J21/UTM_b021_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J22
#' b022_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b022JaguarAlbersBuffer70/b022_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J22trk))
#crs(hf)   
#compareCRS(hf,get_crs(J22trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J22/UTM_b022_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J22/UTM_b022_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J23
#' b023_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b023JaguarAlbersBuffer70/b023_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J23trk))
#crs(hf)   
#compareCRS(hf,get_crs(J23trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J23/UTM_b023_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J23/UTM_b023_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J24
#' b024_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b024JaguarAlbersBuffer70/b024_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J24trk))
#crs(hf)   
#compareCRS(hf,get_crs(J24trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J24/UTM_b024_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J24/UTM_b024_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J25
#' b025_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b025JaguarAlbersBuffer70/b025_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J25trk))
#crs(hf)   
#compareCRS(hf,get_crs(J25trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J25/UTM_b025_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J25/UTM_b025_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J26
#' b026_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b026JaguarAlbersBuffer70/b026_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J26trk))
#crs(hf)   
#compareCRS(hf,get_crs(J26trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J26/UTM_b026_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J26/UTM_b026_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J27
#' b027_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b027JaguarAlbersBuffer70/b027_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J27trk))
#crs(hf)   
#compareCRS(hf,get_crs(J27trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J27/UTM_b027_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J27/UTM_b027_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J28
#' b028_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b028JaguarAlbersBuffer70/b028_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J28trk))
#crs(hf)   
#compareCRS(hf,get_crs(J28trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J28/UTM_b028_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J28/UTM_b028_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J29
#' b029_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b029JaguarAlbersBuffer70/b029_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J29trk))
#crs(hf)   
#compareCRS(hf,get_crs(J29trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J29/UTM_b029_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J29/UTM_b029_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J30
#' b030_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b030JaguarAlbersBuffer70/b030_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J30trk))
#crs(hf)   
#compareCRS(hf,get_crs(J30trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J30/UTM_b030_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J30/UTM_b030_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J30
#' b030_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b030JaguarAlbersBuffer70/b030_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J30trk))
#crs(hf)   
#compareCRS(hf,get_crs(J30trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J30/UTM_b030_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J30/UTM_b030_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J31
#' b031_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b031JaguarAlbersBuffer70/b031_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J31trk))
#crs(hf)   
#compareCRS(hf,get_crs(J31trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J31/UTM_b031_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J31/UTM_b031_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J32
#' b032_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b032JaguarAlbersBuffer70/b032_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J32trk))
#crs(hf)   
#compareCRS(hf,get_crs(J32trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J32/UTM_b032_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J32/UTM_b032_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J33
#' b033_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b033JaguarAlbersBuffer70/b033_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J33trk))
#crs(hf)   
#compareCRS(hf,get_crs(J33trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J33/UTM_b033_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J33/UTM_b033_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J34
#' b034_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b034JaguarAlbersBuffer70/b034_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J34trk))
#crs(hf)   
#compareCRS(hf,get_crs(J34trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J34/UTM_b034_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J34/UTM_b034_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J35
#' b035_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b035JaguarAlbersBuffer70/b035_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J35trk))
#crs(hf)   
#compareCRS(hf,get_crs(J35trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J35/UTM_b035_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J35/UTM_b035_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J36
#' b036_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b036JaguarAlbersBuffer70/b036_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J36trk))
#crs(hf)   
#compareCRS(hf,get_crs(J36trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J36/UTM_b036_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J36/UTM_b036_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J37
#' b037_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b037JaguarAlbersBuffer70/b037_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J37trk))
#crs(hf)   
#compareCRS(hf,get_crs(J37trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J37/UTM_b037_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J37/UTM_b037_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J38
#' b038_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b038JaguarAlbersBuffer70/b038_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J38trk))
#crs(hf)   
#compareCRS(hf,get_crs(J38trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J38/UTM_b038_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J38/UTM_b038_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J39
#' b039_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b039JaguarAlbersBuffer70/b039_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J39trk))
#crs(hf)   
#compareCRS(hf,get_crs(J39trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J39/UTM_b039_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J39/UTM_b039_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J40
#' b040_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b040JaguarAlbersBuffer70/b040_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J40trk))
#crs(hf)   
#compareCRS(hf,get_crs(J40trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J40/UTM_b040_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J40/UTM_b040_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J41
#' b041_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b041JaguarAlbersBuffer70/b041_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J41trk))
#crs(hf)   
#compareCRS(hf,get_crs(J41trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J41/UTM_b041_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J41/UTM_b041_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J42
#' b042_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b042JaguarAlbersBuffer70/b042_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J42trk))
#crs(hf)   
#compareCRS(hf,get_crs(J42trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J42/UTM_b042_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J42/UTM_b042_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

gc()
gc()
gc()

#' J43
#' b043_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b043JaguarAlbersBuffer70/b043_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J43trk))
#crs(hf)   
#compareCRS(hf,get_crs(J43trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J43/UTM_b043_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J43/UTM_b043_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J44
#' b044_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b044JaguarAlbersBuffer70/b044_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J44trk))
#crs(hf)   
#compareCRS(hf,get_crs(J44trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J44/UTM_b044_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J44/UTM_b044_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J45
#' b045_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b045JaguarAlbersBuffer70/b045_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J45trk))
#crs(hf)   
#compareCRS(hf,get_crs(J45trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J45/UTM_b045_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J45/UTM_b045_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J46
#' b046_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b046JaguarAlbersBuffer70/b046_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J46trk))
#crs(hf)   
#compareCRS(hf,get_crs(J46trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J46/UTM_b046_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J46/UTM_b046_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J47
#' b047_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b047JaguarAlbersBuffer70/b047_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J47trk))
#crs(hf)   
#compareCRS(hf,get_crs(J47trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J47/UTM_b047_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J47/UTM_b047_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J48
#' b048_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b048JaguarAlbersBuffer70/b048_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J48trk))
#crs(hf)   
#compareCRS(hf,get_crs(J48trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J48/UTM_b048_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J48/UTM_b048_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J49
#' b049_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b049JaguarAlbersBuffer70/b049_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J49trk))
#crs(hf)   
#compareCRS(hf,get_crs(J49trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J49/UTM_b049_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J49/UTM_b049_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J50
#' b050_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b050JaguarAlbersBuffer70/b050_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J50trk))
#crs(hf)   
#compareCRS(hf,get_crs(J50trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J50/UTM_b050_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J50/UTM_b050_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J51
#' b051_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b051JaguarAlbersBuffer70/b051_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J51trk))
#crs(hf)   
#compareCRS(hf,get_crs(J51trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J51/UTM_b051_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J51/UTM_b051_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J52
#' b052_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b052JaguarAlbersBuffer70/b052_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J52trk))
#crs(hf)   
#compareCRS(hf,get_crs(J52trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J52/UTM_b052_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J52/UTM_b052_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J53
#' b053_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b053JaguarAlbersBuffer70/b053_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J53trk))
#crs(hf)   
#compareCRS(hf,get_crs(J53trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J53/UTM_b053_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J53/UTM_b053_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J54
#' b054_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b054JaguarAlbersBuffer70/b054_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J54trk))
#crs(hf)   
#compareCRS(hf,get_crs(J54trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J54/UTM_b054_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J54/UTM_b054_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J55
#' b055_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b055JaguarAlbersBuffer70/b055_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J55trk))
#crs(hf)   
#compareCRS(hf,get_crs(J55trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J55/UTM_b055_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J55/UTM_b055_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J56
#' b056_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b056JaguarAlbersBuffer70/b056_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J56trk))
#crs(hf)   
#compareCRS(hf,get_crs(J56trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J56/UTM_b056_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J56/UTM_b056_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J57
#' b057_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b057JaguarAlbersBuffer70/b057_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J57trk))
#crs(hf)   
#compareCRS(hf,get_crs(J57trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J57/UTM_b057_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J57/UTM_b057_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J58
#' b058_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b058JaguarAlbersBuffer70/b058_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J58trk))
#crs(hf)   
#compareCRS(hf,get_crs(J58trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J58/UTM_b058_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J58/UTM_b058_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J59
#' b059_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b059JaguarAlbersBuffer70/b059_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J59trk))
#crs(hf)   
#compareCRS(hf,get_crs(J59trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J59/UTM_b059_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J59/UTM_b059_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J60
#' b060_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b060JaguarAlbersBuffer70/b060_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J60trk))
#crs(hf)   
#compareCRS(hf,get_crs(J60trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J60/UTM_b060_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J60/UTM_b060_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J61
#' b061_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b061JaguarAlbersBuffer70/b061_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J61trk))
#crs(hf)   
#compareCRS(hf,get_crs(J61trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J61/UTM_b061_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J61/UTM_b061_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J62
#' b062_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b062JaguarAlbersBuffer70/b062_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J62trk))
#crs(hf)   
#compareCRS(hf,get_crs(J62trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J62/UTM_b062_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J62/UTM_b062_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J63
#' b063_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b063JaguarAlbersBuffer70/b063_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J63trk))
#crs(hf)   
#compareCRS(hf,get_crs(J63trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J63/UTM_b063_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J63/UTM_b063_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

gc()
gc()
gc()


#' J64
#' b064_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b064JaguarAlbersBuffer70/b064_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J64trk))
#crs(hf)   
#compareCRS(hf,get_crs(J64trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J64/UTM_b064_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J64/UTM_b064_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J65
#' b065_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b065JaguarAlbersBuffer70/b065_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J65trk))
#crs(hf)   
#compareCRS(hf,get_crs(J65trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J65/UTM_b065_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J65/UTM_b065_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J66
#' b066_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b066JaguarAlbersBuffer70/b066_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J66trk))
#crs(hf)   
#compareCRS(hf,get_crs(J66trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J66/UTM_b066_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J66/UTM_b066_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J67
#' b067_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b067JaguarAlbersBuffer70/b067_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J67trk))
#crs(hf)   
#compareCRS(hf,get_crs(J67trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J67/UTM_b067_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J67/UTM_b067_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J68
#' b068_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b068JaguarAlbersBuffer70/b068_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J68trk))
#crs(hf)   
#compareCRS(hf,get_crs(J68trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J68/UTM_b068_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J68/UTM_b068_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J69
#' b069_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b069JaguarAlbersBuffer70/b069_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J69trk))
#crs(hf)   
#compareCRS(hf,get_crs(J69trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J69/UTM_b069_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J69/UTM_b069_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

gc()
gc()
gc()


#' J70
#' b070_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b070JaguarAlbersBuffer70/b070_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J70trk))
#crs(hf)   
#compareCRS(hf,get_crs(J70trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J70/UTM_b070_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J70/UTM_b070_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J71
#' b071_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b071JaguarAlbersBuffer70/b071_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J71trk))
#crs(hf)   
#compareCRS(hf,get_crs(J71trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J71/UTM_b071_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J71/UTM_b071_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J72
#' b072_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b072JaguarAlbersBuffer70/b072_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J72trk))
#crs(hf)   
#compareCRS(hf,get_crs(J72trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J72/UTM_b072_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J72/UTM_b072_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J73
#' b073_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b073JaguarAlbersBuffer70/b073_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J73trk))
#crs(hf)   
#compareCRS(hf,get_crs(J73trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J73/UTM_b073_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J73/UTM_b073_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J74
#' b074_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b074JaguarAlbersBuffer70/b074_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J74trk))
#crs(hf)   
#compareCRS(hf,get_crs(J74trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J74/UTM_b074_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J74/UTM_b074_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J75
#' b075_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b075JaguarAlbersBuffer70/b075_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J75trk))
#crs(hf)   
#compareCRS(hf,get_crs(J75trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J75/UTM_b075_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J75/UTM_b075_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J76
#' b076_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b076JaguarAlbersBuffer70/b076_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J76trk))
#crs(hf)   
#compareCRS(hf,get_crs(J76trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J76/UTM_b076_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J76/UTM_b076_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J77
#' b077_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b077JaguarAlbersBuffer70/b077_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J77trk))
#crs(hf)   
#compareCRS(hf,get_crs(J77trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J77/UTM_b077_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J77/UTM_b077_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J78
#' b078_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b078JaguarAlbersBuffer70/b078_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J78trk))
#crs(hf)   
#compareCRS(hf,get_crs(J78trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J78/UTM_b078_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J78/UTM_b078_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J79
#' b079_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b079JaguarAlbersBuffer70/b079_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J79trk))
#crs(hf)   
#compareCRS(hf,get_crs(J79trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J79/UTM_b079_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J79/UTM_b079_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J80
#' b080_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b080JaguarAlbersBuffer70/b080_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J80trk))
#crs(hf)   
#compareCRS(hf,get_crs(J80trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J80/UTM_b080_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J80/UTM_b080_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J81
#' b081_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b081JaguarAlbersBuffer70/b081_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J81trk))
#crs(hf)   
#compareCRS(hf,get_crs(J81trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J81/UTM_b081_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J81/UTM_b081_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J82
#' b082_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b082JaguarAlbersBuffer70/b082_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J82trk))
#crs(hf)   
#compareCRS(hf,get_crs(J82trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J82/UTM_b082_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J82/UTM_b082_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J83
#' b083_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b083JaguarAlbersBuffer70/b083_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J83trk))
#crs(hf)   
#compareCRS(hf,get_crs(J83trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J83/UTM_b083_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J83/UTM_b083_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J84
#' b084_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b084JaguarAlbersBuffer70/b084_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J84trk))
#crs(hf)   
#compareCRS(hf,get_crs(J84trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J84/UTM_b084_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J84/UTM_b084_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J85
#' b085_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b085JaguarAlbersBuffer70/b085_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J85trk))
#crs(hf)   
#compareCRS(hf,get_crs(J85trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J85/UTM_b085_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J85/UTM_b085_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J86
#' b086_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b086JaguarAlbersBuffer70/b086_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J86trk))
#crs(hf)   
#compareCRS(hf,get_crs(J86trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J86/UTM_b086_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J86/UTM_b086_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J87
#' b087_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b087JaguarAlbersBuffer70/b087_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J87trk))
#crs(hf)   
#compareCRS(hf,get_crs(J87trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J87/UTM_b087_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J87/UTM_b087_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J88
#' b088_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b088JaguarAlbersBuffer70/b088_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J88trk))
#crs(hf)   
#compareCRS(hf,get_crs(J88trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J88/UTM_b088_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J88/UTM_b088_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J89
#' b089_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b089JaguarAlbersBuffer70/b089_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J89trk))
#crs(hf)   
#compareCRS(hf,get_crs(J89trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J89/UTM_b089_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J89/UTM_b089_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J90
#' b090_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b090JaguarAlbersBuffer70/b090_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J90trk))
#crs(hf)   
#compareCRS(hf,get_crs(J90trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J90/UTM_b090_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J90/UTM_b090_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J91
#' b091_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b091JaguarAlbersBuffer70/b091_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J91trk))
#crs(hf)   
#compareCRS(hf,get_crs(J91trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J91/UTM_b091_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J91/UTM_b091_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J92
#' b092_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b092JaguarAlbersBuffer70/b092_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J92trk))
#crs(hf)   
#compareCRS(hf,get_crs(J92trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J92/UTM_b092_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J92/UTM_b092_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J93
#' b093_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b093JaguarAlbersBuffer70/b093_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J93trk))
#crs(hf)   
#compareCRS(hf,get_crs(J93trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J93/UTM_b093_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J93/UTM_b093_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J94
#' b094_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b094JaguarAlbersBuffer70/b094_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J94trk))
#crs(hf)   
#compareCRS(hf,get_crs(J94trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J94/UTM_b094_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J94/UTM_b094_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J95
#' b095_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b095JaguarAlbersBuffer70/b095_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J95trk))
#crs(hf)   
#compareCRS(hf,get_crs(J95trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J95/UTM_b095_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J95/UTM_b095_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J96
#' b096_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b096JaguarAlbersBuffer70/b096_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J96trk))
#crs(hf)   
#compareCRS(hf,get_crs(J96trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J96/UTM_b096_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J96/UTM_b096_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J97
#' b097_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b097JaguarAlbersBuffer70/b097_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J97trk))
#crs(hf)   
#compareCRS(hf,get_crs(J97trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J97/UTM_b097_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J97/UTM_b097_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J98
#' b098_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b098JaguarAlbersBuffer70/b098_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J98trk))
#crs(hf)   
#compareCRS(hf,get_crs(J98trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J98/UTM_b098_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J98/UTM_b098_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J99
#' b099_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b099JaguarAlbersBuffer70/b099_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J99trk))
#crs(hf)   
#compareCRS(hf,get_crs(J99trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J99/UTM_b099_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J99/UTM_b099_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J100
#' b100_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b100JaguarAlbersBuffer70/b100_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J100trk))
#crs(hf)   
#compareCRS(hf,get_crs(J100trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J100/UTM_b100_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J100/UTM_b010_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J101
#' b101_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b101JaguarAlbersBuffer70/b101_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J101trk))
#crs(hf)   
#compareCRS(hf,get_crs(J101trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J101/UTM_b101_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J101/UTM_b101_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J102
#' b102_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b102JaguarAlbersBuffer70/b102_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J102trk))
#crs(hf)   
#compareCRS(hf,get_crs(J102trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J102/UTM_b102_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J102/UTM_b102_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J103
#' b103_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b103JaguarAlbersBuffer70/b103_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J103trk))
#crs(hf)   
#compareCRS(hf,get_crs(J103trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J103/UTM_b103_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J103/UTM_b103_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J104
#' b104_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b104JaguarAlbersBuffer70/b104_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J104trk))
#crs(hf)   
#compareCRS(hf,get_crs(J104trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J104/UTM_b104_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J104/UTM_b104_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J105
#' b105_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b105JaguarAlbersBuffer70/b105_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J105trk))
#crs(hf)   
#compareCRS(hf,get_crs(J105trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J105/UTM_b105_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J105/UTM_b105_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J106
#' b106_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b106JaguarAlbersBuffer70/b106_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J106trk))
#crs(hf)   
#compareCRS(hf,get_crs(J106trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J106/UTM_b106_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J106/UTM_b106_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J107
#' b107_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b107JaguarAlbersBuffer70/b107_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J107trk))
#crs(hf)   
#compareCRS(hf,get_crs(J107trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J107/UTM_b107_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J107/UTM_b107_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J108
#' b108_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b108JaguarAlbersBuffer70/b108_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J108trk))
#crs(hf)   
#compareCRS(hf,get_crs(J108trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J108/UTM_b108_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J108/UTM_b108_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J109
#' b109_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b109JaguarAlbersBuffer70/b109_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J109trk))
#crs(hf)   
#compareCRS(hf,get_crs(J109trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J109/UTM_b109_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J109/UTM_b109_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J110
#' b110_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b110JaguarAlbersBuffer70/b110_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J110trk))
#crs(hf)   
#compareCRS(hf,get_crs(J110trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J110/UTM_b110_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J110/UTM_b110_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J111
#' b111_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b111JaguarAlbersBuffer70/b111_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J111trk))
#crs(hf)   
#compareCRS(hf,get_crs(J111trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J111/UTM_b111_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J111/UTM_b111_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J112
#' b112_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b112JaguarAlbersBuffer70/b112_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J112trk))
#crs(hf)   
#compareCRS(hf,get_crs(J112trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J112/UTM_b112_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J112/UTM_b112_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J113
#' b113_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b113JaguarAlbersBuffer70/b113_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J113trk))
#crs(hf)   
#compareCRS(hf,get_crs(J113trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J113/UTM_b113_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J113/UTM_b113_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J114
#' b114_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b114JaguarAlbersBuffer70/b114_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J114trk))
#crs(hf)   
#compareCRS(hf,get_crs(J114trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J114/UTM_b114_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J114/UTM_b114_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J115
#' b115_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b115JaguarAlbersBuffer70/b115_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J115trk))
#crs(hf)   
#compareCRS(hf,get_crs(J115trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J115/UTM_b115_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J115/UTM_b115_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J116
#' b116_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b116JaguarAlbersBuffer70/b116_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J116trk))
#crs(hf)   
#compareCRS(hf,get_crs(J116trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J116/UTM_b116_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J116/UTM_b116_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J117
#' b117_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b117JaguarAlbersBuffer70/b117_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J117trk))
#crs(hf)   
#compareCRS(hf,get_crs(J117trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J117/UTM_b117_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J117/UTM_b117_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#####################################################################################################################################
# Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.

gc()
gc()
gc()

#' J1
#' b001_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif
hf <-"D:/GISbases/b001JaguarAlbersBuffer70/b001_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J1trk))
#crs(hf)   
#compareCRS(hf,get_crs(J1trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J1/UTM_b001_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J1/UTM_b001_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J2
#' b002_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b002JaguarAlbersBuffer70/b002_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J2trk))
#crs(hf)   
#compareCRS(hf,get_crs(J2trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J2/UTM_b002_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J2/UTM_b002_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J3
#' b003_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b003JaguarAlbersBuffer70/b003_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J3trk))
#crs(hf)   
#compareCRS(hf,get_crs(J3trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J3/UTM_b003_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J3/UTM_b003_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J4
#' b004_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b004JaguarAlbersBuffer70/b004_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J4trk))
#crs(hf)   
#compareCRS(hf,get_crs(J4trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J4/UTM_b004_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J4/UTM_b004_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J5
#' b005_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b005JaguarAlbersBuffer70/b005_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J5trk))
#crs(hf)   
#compareCRS(hf,get_crs(J5trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J5/UTM_b005_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J5/UTM_b005_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J6
#' b006_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b006JaguarAlbersBuffer70/b006_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J6trk))
#crs(hf)   
#compareCRS(hf,get_crs(J6trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J6/UTM_b006_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J6/UTM_b006_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J7
#' b007_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b007JaguarAlbersBuffer70/b007_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J7trk))
#crs(hf)   
#compareCRS(hf,get_crs(J7trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J7/UTM_b007_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J7/UTM_b007_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J8
#' b008_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b008JaguarAlbersBuffer70/b008_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J8trk))
#crs(hf)   
#compareCRS(hf,get_crs(J8trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J8/UTM_b008_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J8/UTM_b008_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J9
#' b009_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b009JaguarAlbersBuffer70/b009_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J9trk))
#crs(hf)   
#compareCRS(hf,get_crs(J9trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J9/UTM_b009_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J9/UTM_b009_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J10
#' b010_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b010JaguarAlbersBuffer70/b010_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J10trk))
#crs(hf)   
#compareCRS(hf,get_crs(J10trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J10/UTM_b010_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J10/UTM_b010_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

gc()
gc()
gc()


#' J11
#' b011_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b011JaguarAlbersBuffer70/b011_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J11trk))
#crs(hf)   
#compareCRS(hf,get_crs(J11trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J11/UTM_b011_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J11/UTM_b011_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J12
#' b012_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b012JaguarAlbersBuffer70/b012_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J12trk))
#crs(hf)   
#compareCRS(hf,get_crs(J12trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J12/UTM_b012_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J12/UTM_b012_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J13
#' b013_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b013JaguarAlbersBuffer70/b013_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J13trk))
#crs(hf)   
#compareCRS(hf,get_crs(J13trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J13/UTM_b013_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J13/UTM_b013_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J14
#' b014_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b014JaguarAlbersBuffer70/b014_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J14trk))
#crs(hf)   
#compareCRS(hf,get_crs(J14trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J14/UTM_b014_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J14/UTM_b014_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J15
#' b015_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b015JaguarAlbersBuffer70/b015_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J15trk))
#crs(hf)   
#compareCRS(hf,get_crs(J15trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J15/UTM_b015_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J15/UTM_b015_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J16
#' b016_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b016JaguarAlbersBuffer70/b016_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J16trk))
#crs(hf)   
#compareCRS(hf,get_crs(J16trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J16/UTM_b016_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J16/UTM_b016_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J17
#' b017_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b017JaguarAlbersBuffer70/b017_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J17trk))
#crs(hf)   
#compareCRS(hf,get_crs(J17trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J17/UTM_b017_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J17/UTM_b017_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J18
#' b018_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b018JaguarAlbersBuffer70/b018_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J18trk))
#crs(hf)   
#compareCRS(hf,get_crs(J18trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J18/UTM_b018_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J18/UTM_b018_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J19
#' b019_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b019JaguarAlbersBuffer70/b019_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J19trk))
#crs(hf)   
#compareCRS(hf,get_crs(J19trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J19/UTM_b019_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J19/UTM_b019_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J20
#' b020_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b020JaguarAlbersBuffer70/b020_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J20trk))
#crs(hf)   
#compareCRS(hf,get_crs(J20trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J20/UTM_b020_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J20/UTM_b020_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J21
#' b021_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b021JaguarAlbersBuffer70/b021_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J21trk))
#crs(hf)   
#compareCRS(hf,get_crs(J21trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J21/UTM_b021_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J21/UTM_b021_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J22
#' b022_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b022JaguarAlbersBuffer70/b022_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J22trk))
#crs(hf)   
#compareCRS(hf,get_crs(J22trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J22/UTM_b022_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J22/UTM_b022_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J23
#' b023_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b023JaguarAlbersBuffer70/b023_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J23trk))
#crs(hf)   
#compareCRS(hf,get_crs(J23trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J23/UTM_b023_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J23/UTM_b023_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J24
#' b024_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b024JaguarAlbersBuffer70/b024_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J24trk))
#crs(hf)   
#compareCRS(hf,get_crs(J24trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J24/UTM_b024_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J24/UTM_b024_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J25
#' b025_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b025JaguarAlbersBuffer70/b025_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J25trk))
#crs(hf)   
#compareCRS(hf,get_crs(J25trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J25/UTM_b025_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J25/UTM_b025_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J26
#' b026_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b026JaguarAlbersBuffer70/b026_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J26trk))
#crs(hf)   
#compareCRS(hf,get_crs(J26trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J26/UTM_b026_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J26/UTM_b026_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J27
#' b027_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b027JaguarAlbersBuffer70/b027_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J27trk))
#crs(hf)   
#compareCRS(hf,get_crs(J27trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J27/UTM_b027_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J27/UTM_b027_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J28
#' b028_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b028JaguarAlbersBuffer70/b028_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J28trk))
#crs(hf)   
#compareCRS(hf,get_crs(J28trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J28/UTM_b028_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J28/UTM_b028_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J29
#' b029_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b029JaguarAlbersBuffer70/b029_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J29trk))
#crs(hf)   
#compareCRS(hf,get_crs(J29trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J29/UTM_b029_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J29/UTM_b029_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J30
#' b030_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b030JaguarAlbersBuffer70/b030_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J30trk))
#crs(hf)   
#compareCRS(hf,get_crs(J30trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J30/UTM_b030_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J30/UTM_b030_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J30
#' b030_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b030JaguarAlbersBuffer70/b030_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J30trk))
#crs(hf)   
#compareCRS(hf,get_crs(J30trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J30/UTM_b030_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J30/UTM_b030_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J31
#' b031_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b031JaguarAlbersBuffer70/b031_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J31trk))
#crs(hf)   
#compareCRS(hf,get_crs(J31trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J31/UTM_b031_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J31/UTM_b031_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J32
#' b032_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b032JaguarAlbersBuffer70/b032_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J32trk))
#crs(hf)   
#compareCRS(hf,get_crs(J32trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J32/UTM_b032_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J32/UTM_b032_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J33
#' b033_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b033JaguarAlbersBuffer70/b033_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J33trk))
#crs(hf)   
#compareCRS(hf,get_crs(J33trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J33/UTM_b033_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J33/UTM_b033_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J34
#' b034_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b034JaguarAlbersBuffer70/b034_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J34trk))
#crs(hf)   
#compareCRS(hf,get_crs(J34trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J34/UTM_b034_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J34/UTM_b034_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J35
#' b035_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b035JaguarAlbersBuffer70/b035_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J35trk))
#crs(hf)   
#compareCRS(hf,get_crs(J35trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J35/UTM_b035_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J35/UTM_b035_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J36
#' b036_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b036JaguarAlbersBuffer70/b036_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J36trk))
#crs(hf)   
#compareCRS(hf,get_crs(J36trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J36/UTM_b036_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J36/UTM_b036_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J37
#' b037_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b037JaguarAlbersBuffer70/b037_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J37trk))
#crs(hf)   
#compareCRS(hf,get_crs(J37trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J37/UTM_b037_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J37/UTM_b037_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J38
#' b038_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b038JaguarAlbersBuffer70/b038_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J38trk))
#crs(hf)   
#compareCRS(hf,get_crs(J38trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J38/UTM_b038_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J38/UTM_b038_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J39
#' b039_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b039JaguarAlbersBuffer70/b039_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J39trk))
#crs(hf)   
#compareCRS(hf,get_crs(J39trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J39/UTM_b039_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J39/UTM_b039_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J40
#' b040_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b040JaguarAlbersBuffer70/b040_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J40trk))
#crs(hf)   
#compareCRS(hf,get_crs(J40trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J40/UTM_b040_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J40/UTM_b040_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J41
#' b041_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b041JaguarAlbersBuffer70/b041_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J41trk))
#crs(hf)   
#compareCRS(hf,get_crs(J41trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J41/UTM_b041_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J41/UTM_b041_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J42
#' b042_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b042JaguarAlbersBuffer70/b042_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J42trk))
#crs(hf)   
#compareCRS(hf,get_crs(J42trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J42/UTM_b042_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J42/UTM_b042_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

gc()
gc()
gc()

#' J43
#' b043_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b043JaguarAlbersBuffer70/b043_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J43trk))
#crs(hf)   
#compareCRS(hf,get_crs(J43trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J43/UTM_b043_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J43/UTM_b043_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J44
#' b044_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b044JaguarAlbersBuffer70/b044_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J44trk))
#crs(hf)   
#compareCRS(hf,get_crs(J44trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J44/UTM_b044_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J44/UTM_b044_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J45
#' b045_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b045JaguarAlbersBuffer70/b045_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J45trk))
#crs(hf)   
#compareCRS(hf,get_crs(J45trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J45/UTM_b045_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J45/UTM_b045_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J46
#' b046_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b046JaguarAlbersBuffer70/b046_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J46trk))
#crs(hf)   
#compareCRS(hf,get_crs(J46trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J46/UTM_b046_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J46/UTM_b046_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J47
#' b047_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b047JaguarAlbersBuffer70/b047_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J47trk))
#crs(hf)   
#compareCRS(hf,get_crs(J47trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J47/UTM_b047_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J47/UTM_b047_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J48
#' b048_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b048JaguarAlbersBuffer70/b048_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J48trk))
#crs(hf)   
#compareCRS(hf,get_crs(J48trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J48/UTM_b048_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J48/UTM_b048_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J49
#' b049_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b049JaguarAlbersBuffer70/b049_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J49trk))
#crs(hf)   
#compareCRS(hf,get_crs(J49trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J49/UTM_b049_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J49/UTM_b049_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J50
#' b050_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b050JaguarAlbersBuffer70/b050_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J50trk))
#crs(hf)   
#compareCRS(hf,get_crs(J50trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J50/UTM_b050_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J50/UTM_b050_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J51
#' b051_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b051JaguarAlbersBuffer70/b051_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J51trk))
#crs(hf)   
#compareCRS(hf,get_crs(J51trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J51/UTM_b051_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J51/UTM_b051_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J52
#' b052_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b052JaguarAlbersBuffer70/b052_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J52trk))
#crs(hf)   
#compareCRS(hf,get_crs(J52trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J52/UTM_b052_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J52/UTM_b052_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J53
#' b053_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b053JaguarAlbersBuffer70/b053_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J53trk))
#crs(hf)   
#compareCRS(hf,get_crs(J53trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J53/UTM_b053_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J53/UTM_b053_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J54
#' b054_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b054JaguarAlbersBuffer70/b054_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J54trk))
#crs(hf)   
#compareCRS(hf,get_crs(J54trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J54/UTM_b054_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J54/UTM_b054_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J55
#' b055_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b055JaguarAlbersBuffer70/b055_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J55trk))
#crs(hf)   
#compareCRS(hf,get_crs(J55trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J55/UTM_b055_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J55/UTM_b055_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J56
#' b056_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b056JaguarAlbersBuffer70/b056_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J56trk))
#crs(hf)   
#compareCRS(hf,get_crs(J56trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J56/UTM_b056_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J56/UTM_b056_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J57
#' b057_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b057JaguarAlbersBuffer70/b057_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J57trk))
#crs(hf)   
#compareCRS(hf,get_crs(J57trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J57/UTM_b057_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J57/UTM_b057_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J58
#' b058_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b058JaguarAlbersBuffer70/b058_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J58trk))
#crs(hf)   
#compareCRS(hf,get_crs(J58trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J58/UTM_b058_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J58/UTM_b058_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J59
#' b059_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b059JaguarAlbersBuffer70/b059_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J59trk))
#crs(hf)   
#compareCRS(hf,get_crs(J59trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J59/UTM_b059_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J59/UTM_b059_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J60
#' b060_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b060JaguarAlbersBuffer70/b060_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J60trk))
#crs(hf)   
#compareCRS(hf,get_crs(J60trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J60/UTM_b060_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J60/UTM_b060_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J61
#' b061_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b061JaguarAlbersBuffer70/b061_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J61trk))
#crs(hf)   
#compareCRS(hf,get_crs(J61trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J61/UTM_b061_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J61/UTM_b061_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J62
#' b062_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b062JaguarAlbersBuffer70/b062_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J62trk))
#crs(hf)   
#compareCRS(hf,get_crs(J62trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J62/UTM_b062_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J62/UTM_b062_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J63
#' b063_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b063JaguarAlbersBuffer70/b063_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J63trk))
#crs(hf)   
#compareCRS(hf,get_crs(J63trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J63/UTM_b063_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J63/UTM_b063_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

gc()
gc()
gc()


#' J64
#' b064_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b064JaguarAlbersBuffer70/b064_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J64trk))
#crs(hf)   
#compareCRS(hf,get_crs(J64trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J64/UTM_b064_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J64/UTM_b064_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J65
#' b065_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b065JaguarAlbersBuffer70/b065_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J65trk))
#crs(hf)   
#compareCRS(hf,get_crs(J65trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J65/UTM_b065_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J65/UTM_b065_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J66
#' b066_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b066JaguarAlbersBuffer70/b066_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J66trk))
#crs(hf)   
#compareCRS(hf,get_crs(J66trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J66/UTM_b066_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J66/UTM_b066_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J67
#' b067_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b067JaguarAlbersBuffer70/b067_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J67trk))
#crs(hf)   
#compareCRS(hf,get_crs(J67trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J67/UTM_b067_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J67/UTM_b067_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J68
#' b068_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b068JaguarAlbersBuffer70/b068_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J68trk))
#crs(hf)   
#compareCRS(hf,get_crs(J68trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J68/UTM_b068_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J68/UTM_b068_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J69
#' b069_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b069JaguarAlbersBuffer70/b069_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J69trk))
#crs(hf)   
#compareCRS(hf,get_crs(J69trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J69/UTM_b069_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J69/UTM_b069_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

gc()
gc()
gc()


#' J70
#' b070_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b070JaguarAlbersBuffer70/b070_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J70trk))
#crs(hf)   
#compareCRS(hf,get_crs(J70trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J70/UTM_b070_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J70/UTM_b070_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J71
#' b071_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b071JaguarAlbersBuffer70/b071_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J71trk))
#crs(hf)   
#compareCRS(hf,get_crs(J71trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J71/UTM_b071_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J71/UTM_b071_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J72
#' b072_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b072JaguarAlbersBuffer70/b072_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J72trk))
#crs(hf)   
#compareCRS(hf,get_crs(J72trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J72/UTM_b072_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J72/UTM_b072_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J73
#' b073_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b073JaguarAlbersBuffer70/b073_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J73trk))
#crs(hf)   
#compareCRS(hf,get_crs(J73trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J73/UTM_b073_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J73/UTM_b073_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J74
#' b074_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b074JaguarAlbersBuffer70/b074_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J74trk))
#crs(hf)   
#compareCRS(hf,get_crs(J74trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J74/UTM_b074_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J74/UTM_b074_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J75
#' b075_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b075JaguarAlbersBuffer70/b075_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J75trk))
#crs(hf)   
#compareCRS(hf,get_crs(J75trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J75/UTM_b075_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J75/UTM_b075_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J76
#' b076_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b076JaguarAlbersBuffer70/b076_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J76trk))
#crs(hf)   
#compareCRS(hf,get_crs(J76trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J76/UTM_b076_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J76/UTM_b076_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J77
#' b077_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b077JaguarAlbersBuffer70/b077_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J77trk))
#crs(hf)   
#compareCRS(hf,get_crs(J77trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J77/UTM_b077_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J77/UTM_b077_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J78
#' b078_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b078JaguarAlbersBuffer70/b078_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J78trk))
#crs(hf)   
#compareCRS(hf,get_crs(J78trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J78/UTM_b078_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J78/UTM_b078_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J79
#' b079_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b079JaguarAlbersBuffer70/b079_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J79trk))
#crs(hf)   
#compareCRS(hf,get_crs(J79trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J79/UTM_b079_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J79/UTM_b079_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J80
#' b080_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b080JaguarAlbersBuffer70/b080_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J80trk))
#crs(hf)   
#compareCRS(hf,get_crs(J80trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J80/UTM_b080_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J80/UTM_b080_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J81
#' b081_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b081JaguarAlbersBuffer70/b081_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J81trk))
#crs(hf)   
#compareCRS(hf,get_crs(J81trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J81/UTM_b081_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J81/UTM_b081_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J82
#' b082_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b082JaguarAlbersBuffer70/b082_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J82trk))
#crs(hf)   
#compareCRS(hf,get_crs(J82trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J82/UTM_b082_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J82/UTM_b082_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J83
#' b083_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b083JaguarAlbersBuffer70/b083_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J83trk))
#crs(hf)   
#compareCRS(hf,get_crs(J83trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J83/UTM_b083_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J83/UTM_b083_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J84
#' b084_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b084JaguarAlbersBuffer70/b084_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J84trk))
#crs(hf)   
#compareCRS(hf,get_crs(J84trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J84/UTM_b084_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J84/UTM_b084_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J85
#' b085_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b085JaguarAlbersBuffer70/b085_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J85trk))
#crs(hf)   
#compareCRS(hf,get_crs(J85trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J85/UTM_b085_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J85/UTM_b085_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J86
#' b086_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b086JaguarAlbersBuffer70/b086_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J86trk))
#crs(hf)   
#compareCRS(hf,get_crs(J86trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J86/UTM_b086_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J86/UTM_b086_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J87
#' b087_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b087JaguarAlbersBuffer70/b087_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J87trk))
#crs(hf)   
#compareCRS(hf,get_crs(J87trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J87/UTM_b087_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J87/UTM_b087_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J88
#' b088_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b088JaguarAlbersBuffer70/b088_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J88trk))
#crs(hf)   
#compareCRS(hf,get_crs(J88trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J88/UTM_b088_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J88/UTM_b088_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J89
#' b089_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b089JaguarAlbersBuffer70/b089_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J89trk))
#crs(hf)   
#compareCRS(hf,get_crs(J89trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J89/UTM_b089_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J89/UTM_b089_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J90
#' b090_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b090JaguarAlbersBuffer70/b090_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J90trk))
#crs(hf)   
#compareCRS(hf,get_crs(J90trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J90/UTM_b090_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J90/UTM_b090_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J91
#' b091_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b091JaguarAlbersBuffer70/b091_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J91trk))
#crs(hf)   
#compareCRS(hf,get_crs(J91trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J91/UTM_b091_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J91/UTM_b091_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J92
#' b092_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b092JaguarAlbersBuffer70/b092_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J92trk))
#crs(hf)   
#compareCRS(hf,get_crs(J92trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J92/UTM_b092_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J92/UTM_b092_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J93
#' b093_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b093JaguarAlbersBuffer70/b093_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J93trk))
#crs(hf)   
#compareCRS(hf,get_crs(J93trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J93/UTM_b093_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J93/UTM_b093_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J94
#' b094_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b094JaguarAlbersBuffer70/b094_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J94trk))
#crs(hf)   
#compareCRS(hf,get_crs(J94trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J94/UTM_b094_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J94/UTM_b094_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J95
#' b095_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b095JaguarAlbersBuffer70/b095_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J95trk))
#crs(hf)   
#compareCRS(hf,get_crs(J95trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J95/UTM_b095_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J95/UTM_b095_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J96
#' b096_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b096JaguarAlbersBuffer70/b096_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J96trk))
#crs(hf)   
#compareCRS(hf,get_crs(J96trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J96/UTM_b096_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J96/UTM_b096_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J97
#' b097_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b097JaguarAlbersBuffer70/b097_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J97trk))
#crs(hf)   
#compareCRS(hf,get_crs(J97trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J97/UTM_b097_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J97/UTM_b097_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J98
#' b098_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b098JaguarAlbersBuffer70/b098_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J98trk))
#crs(hf)   
#compareCRS(hf,get_crs(J98trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J98/UTM_b098_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J98/UTM_b098_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J99
#' b099_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b099JaguarAlbersBuffer70/b099_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J99trk))
#crs(hf)   
#compareCRS(hf,get_crs(J99trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J99/UTM_b099_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J99/UTM_b099_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J100
#' b100_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b100JaguarAlbersBuffer70/b100_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J100trk))
#crs(hf)   
#compareCRS(hf,get_crs(J100trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J100/UTM_b100_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J100/UTM_b010_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J101
#' b101_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b101JaguarAlbersBuffer70/b101_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J101trk))
#crs(hf)   
#compareCRS(hf,get_crs(J101trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J101/UTM_b101_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J101/UTM_b101_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J102
#' b102_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b102JaguarAlbersBuffer70/b102_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J102trk))
#crs(hf)   
#compareCRS(hf,get_crs(J102trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J102/UTM_b102_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J102/UTM_b102_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J103
#' b103_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b103JaguarAlbersBuffer70/b103_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J103trk))
#crs(hf)   
#compareCRS(hf,get_crs(J103trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J103/UTM_b103_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J103/UTM_b103_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J104
#' b104_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b104JaguarAlbersBuffer70/b104_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J104trk))
#crs(hf)   
#compareCRS(hf,get_crs(J104trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J104/UTM_b104_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J104/UTM_b104_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J105
#' b105_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b105JaguarAlbersBuffer70/b105_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J105trk))
#crs(hf)   
#compareCRS(hf,get_crs(J105trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J105/UTM_b105_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J105/UTM_b105_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J106
#' b106_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b106JaguarAlbersBuffer70/b106_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J106trk))
#crs(hf)   
#compareCRS(hf,get_crs(J106trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J106/UTM_b106_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J106/UTM_b106_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J107
#' b107_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b107JaguarAlbersBuffer70/b107_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J107trk))
#crs(hf)   
#compareCRS(hf,get_crs(J107trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J107/UTM_b107_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J107/UTM_b107_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J108
#' b108_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b108JaguarAlbersBuffer70/b108_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J108trk))
#crs(hf)   
#compareCRS(hf,get_crs(J108trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J108/UTM_b108_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J108/UTM_b108_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J109
#' b109_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b109JaguarAlbersBuffer70/b109_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J109trk))
#crs(hf)   
#compareCRS(hf,get_crs(J109trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J109/UTM_b109_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J109/UTM_b109_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J110
#' b110_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b110JaguarAlbersBuffer70/b110_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J110trk))
#crs(hf)   
#compareCRS(hf,get_crs(J110trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J110/UTM_b110_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J110/UTM_b110_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J111
#' b111_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b111JaguarAlbersBuffer70/b111_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J111trk))
#crs(hf)   
#compareCRS(hf,get_crs(J111trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J111/UTM_b111_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J111/UTM_b111_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J112
#' b112_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b112JaguarAlbersBuffer70/b112_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J112trk))
#crs(hf)   
#compareCRS(hf,get_crs(J112trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J112/UTM_b112_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J112/UTM_b112_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J113
#' b113_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b113JaguarAlbersBuffer70/b113_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J113trk))
#crs(hf)   
#compareCRS(hf,get_crs(J113trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J113/UTM_b113_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J113/UTM_b113_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J114
#' b114_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b114JaguarAlbersBuffer70/b114_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J114trk))
#crs(hf)   
#compareCRS(hf,get_crs(J114trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J114/UTM_b114_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J114/UTM_b114_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J115
#' b115_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b115JaguarAlbersBuffer70/b115_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J115trk))
#crs(hf)   
#compareCRS(hf,get_crs(J115trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J115/UTM_b115_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J115/UTM_b115_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J116
#' b116_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b116JaguarAlbersBuffer70/b116_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J116trk))
#crs(hf)   
#compareCRS(hf,get_crs(J116trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J116/UTM_b116_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J116/UTM_b116_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J117
#' b117_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b117JaguarAlbersBuffer70/b117_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J117trk))
#crs(hf)   
#compareCRS(hf,get_crs(J117trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J117/UTM_b117_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J117/UTM_b117_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

######################################################################################################################################
# dist2waterbodies_exp.
gc()
gc()
gc()

#' J1
#' b001_dist2waterbodies_exp.tif
hf <-"D:/GISbases/b001JaguarAlbersBuffer70/b001_dist2waterbodies_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J1trk))
#crs(hf)   
#compareCRS(hf,get_crs(J1trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J1/UTM_b001_dist2waterbodies_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J1/UTM_b001_dist2waterbodies_exp.tif";(test=raster(test)); x11(); plot(test)

#' J2
#' b002_dist2waterbodies_exp.
hf <-"D:/GISbases/b002JaguarAlbersBuffer70/b002_dist2waterbodies_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J2trk))
#crs(hf)   
#compareCRS(hf,get_crs(J2trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J2/UTM_b002_dist2waterbodies_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J2/UTM_b002_dist2waterbodies_exp.tif";(test=raster(test)); x11(); plot(test)

#' J3
#' b003_dist2waterbodies_exp.
hf <-"D:/GISbases/b003JaguarAlbersBuffer70/b003_dist2waterbodies_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J3trk))
#crs(hf)   
#compareCRS(hf,get_crs(J3trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J3/UTM_b003_dist2waterbodies_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J3/UTM_b003_dist2waterbodies_exp.tif";(test=raster(test)); x11(); plot(test)

#' J4
#' b004_dist2waterbodies_exp.
hf <-"D:/GISbases/b004JaguarAlbersBuffer70/b004_dist2waterbodies_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J4trk))
#crs(hf)   
#compareCRS(hf,get_crs(J4trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J4/UTM_b004_dist2waterbodies_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J4/UTM_b004_dist2waterbodies_exp.tif";(test=raster(test)); x11(); plot(test)

#' J5
#' b005_dist2waterbodies_exp.
hf <-"D:/GISbases/b005JaguarAlbersBuffer70/b005_dist2waterbodies_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J5trk))
#crs(hf)   
#compareCRS(hf,get_crs(J5trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J5/UTM_b005_dist2waterbodies_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J5/UTM_b005_dist2waterbodies_exp.tif";(test=raster(test)); x11(); plot(test)

#' J6
#' b006_dist2waterbodies_exp.
hf <-"D:/GISbases/b006JaguarAlbersBuffer70/b006_dist2waterbodies_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J6trk))
#crs(hf)   
#compareCRS(hf,get_crs(J6trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J6/UTM_b006_dist2waterbodies_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J6/UTM_b006_dist2waterbodies_exp.tif";(test=raster(test)); x11(); plot(test)

#' J7
#' b007_dist2waterbodies_exp.
hf <-"D:/GISbases/b007JaguarAlbersBuffer70/b007_dist2waterbodies_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J7trk))
#crs(hf)   
#compareCRS(hf,get_crs(J7trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J7/UTM_b007_dist2waterbodies_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J7/UTM_b007_dist2waterbodies_exp.tif";(test=raster(test)); x11(); plot(test)

#' J8
#' b008_dist2waterbodies_exp.
hf <-"D:/GISbases/b008JaguarAlbersBuffer70/b008_dist2waterbodies_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J8trk))
#crs(hf)   
#compareCRS(hf,get_crs(J8trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J8/UTM_b008_dist2waterbodies_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J8/UTM_b008_dist2waterbodies_exp.tif";(test=raster(test)); x11(); plot(test)

#' J9
#' b009_dist2waterbodies_exp.
hf <-"D:/GISbases/b009JaguarAlbersBuffer70/b009_dist2waterbodies_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J9trk))
#crs(hf)   
#compareCRS(hf,get_crs(J9trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J9/UTM_b009_dist2waterbodies_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J9/UTM_b009_dist2waterbodies_exp.tif";(test=raster(test)); x11(); plot(test)

#' J10
#' b010_dist2waterbodies_exp.
hf <-"D:/GISbases/b010JaguarAlbersBuffer70/b010_dist2waterbodies_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J10trk))
#crs(hf)   
#compareCRS(hf,get_crs(J10trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J10/UTM_b010_dist2waterbodies_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J10/UTM_b010_dist2waterbodies_exp.tif";(test=raster(test)); x11(); plot(test)

gc()
gc()
gc()


#' J11
#' b011_dist2waterbodies_exp.
hf <-"D:/GISbases/b011JaguarAlbersBuffer70/b011_dist2waterbodies_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J11trk))
#crs(hf)   
#compareCRS(hf,get_crs(J11trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J11/UTM_b011_dist2waterbodies_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J11/UTM_b011_dist2waterbodies_exp.tif";(test=raster(test)); x11(); plot(test)

#' J12
#' b012_dist2waterbodies_exp.
hf <-"D:/GISbases/b012JaguarAlbersBuffer70/b012_dist2waterbodies_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J12trk))
#crs(hf)   
#compareCRS(hf,get_crs(J12trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J12/UTM_b012_dist2waterbodies_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J12/UTM_b012_dist2waterbodies_exp.tif";(test=raster(test)); x11(); plot(test)

#' J13
#' b013_dist2waterbodies_exp.
hf <-"D:/GISbases/b013JaguarAlbersBuffer70/b013_dist2waterbodies_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J13trk))
#crs(hf)   
#compareCRS(hf,get_crs(J13trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J13/UTM_b013_dist2waterbodies_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J13/UTM_b013_dist2waterbodies_exp.tif";(test=raster(test)); x11(); plot(test)

#' J14
#' b014_dist2waterbodies_exp.
hf <-"D:/GISbases/b014JaguarAlbersBuffer70/b014_dist2waterbodies_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J14trk))
#crs(hf)   
#compareCRS(hf,get_crs(J14trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J14/UTM_b014_dist2waterbodies_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J14/UTM_b014_dist2waterbodies_exp.tif";(test=raster(test)); x11(); plot(test)

#' J15
#' b015_dist2waterbodies_exp.
hf <-"D:/GISbases/b015JaguarAlbersBuffer70/b015_dist2waterbodies_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J15trk))
#crs(hf)   
#compareCRS(hf,get_crs(J15trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J15/UTM_b015_dist2waterbodies_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J15/UTM_b015_dist2waterbodies_exp.tif";(test=raster(test)); x11(); plot(test)

#' J16
#' b016_dist2waterbodies_exp.
hf <-"D:/GISbases/b016JaguarAlbersBuffer70/b016_dist2waterbodies_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J16trk))
#crs(hf)   
#compareCRS(hf,get_crs(J16trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J16/UTM_b016_dist2waterbodies_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J16/UTM_b016_dist2waterbodies_exp.tif";(test=raster(test)); x11(); plot(test)

#' J17
#' b017_dist2waterbodies_exp.
hf <-"D:/GISbases/b017JaguarAlbersBuffer70/b017_dist2waterbodies_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J17trk))
#crs(hf)   
#compareCRS(hf,get_crs(J17trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J17/UTM_b017_dist2waterbodies_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J17/UTM_b017_dist2waterbodies_exp.tif";(test=raster(test)); x11(); plot(test)

#' J18
#' b018_dist2waterbodies_exp.
hf <-"D:/GISbases/b018JaguarAlbersBuffer70/b018_dist2waterbodies_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J18trk))
#crs(hf)   
#compareCRS(hf,get_crs(J18trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J18/UTM_b018_dist2waterbodies_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J18/UTM_b018_dist2waterbodies_exp.tif";(test=raster(test)); x11(); plot(test)

#' J19
#' b019_dist2waterbodies_exp.
hf <-"D:/GISbases/b019JaguarAlbersBuffer70/b019_dist2waterbodies_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J19trk))
#crs(hf)   
#compareCRS(hf,get_crs(J19trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J19/UTM_b019_dist2waterbodies_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J19/UTM_b019_dist2waterbodies_exp.tif";(test=raster(test)); x11(); plot(test)

#' J20
#' b020_dist2waterbodies_exp.
hf <-"D:/GISbases/b020JaguarAlbersBuffer70/b020_dist2waterbodies_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J20trk))
#crs(hf)   
#compareCRS(hf,get_crs(J20trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J20/UTM_b020_dist2waterbodies_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J20/UTM_b020_dist2waterbodies_exp.tif";(test=raster(test)); x11(); plot(test)

#' J21
#' b021_dist2waterbodies_exp.
hf <-"D:/GISbases/b021JaguarAlbersBuffer70/b021_dist2waterbodies_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J21trk))
#crs(hf)   
#compareCRS(hf,get_crs(J21trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J21/UTM_b021_dist2waterbodies_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J21/UTM_b021_dist2waterbodies_exp.tif";(test=raster(test)); x11(); plot(test)

#' J22
#' b022_dist2waterbodies_exp.
hf <-"D:/GISbases/b022JaguarAlbersBuffer70/b022_dist2waterbodies_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J22trk))
#crs(hf)   
#compareCRS(hf,get_crs(J22trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J22/UTM_b022_dist2waterbodies_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J22/UTM_b022_dist2waterbodies_exp.tif";(test=raster(test)); x11(); plot(test)

#' J23
#' b023_dist2waterbodies_exp.
hf <-"D:/GISbases/b023JaguarAlbersBuffer70/b023_dist2waterbodies_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J23trk))
#crs(hf)   
#compareCRS(hf,get_crs(J23trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J23/UTM_b023_dist2waterbodies_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J23/UTM_b023_dist2waterbodies_exp.tif";(test=raster(test)); x11(); plot(test)

#' J24
#' b024_dist2waterbodies_exp.
hf <-"D:/GISbases/b024JaguarAlbersBuffer70/b024_dist2waterbodies_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J24trk))
#crs(hf)   
#compareCRS(hf,get_crs(J24trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J24/UTM_b024_dist2waterbodies_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J24/UTM_b024_dist2waterbodies_exp.tif";(test=raster(test)); x11(); plot(test)

#' J25
#' b025_dist2waterbodies_exp.
hf <-"D:/GISbases/b025JaguarAlbersBuffer70/b025_dist2waterbodies_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J25trk))
#crs(hf)   
#compareCRS(hf,get_crs(J25trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J25/UTM_b025_dist2waterbodies_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J25/UTM_b025_dist2waterbodies_exp.tif";(test=raster(test)); x11(); plot(test)

#' J26
#' b026_dist2waterbodies_exp.
hf <-"D:/GISbases/b026JaguarAlbersBuffer70/b026_dist2waterbodies_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J26trk))
#crs(hf)   
#compareCRS(hf,get_crs(J26trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J26/UTM_b026_dist2waterbodies_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J26/UTM_b026_dist2waterbodies_exp.tif";(test=raster(test)); x11(); plot(test)

#' J27
#' b027_dist2waterbodies_exp.
hf <-"D:/GISbases/b027JaguarAlbersBuffer70/b027_dist2waterbodies_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J27trk))
#crs(hf)   
#compareCRS(hf,get_crs(J27trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J27/UTM_b027_dist2waterbodies_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J27/UTM_b027_dist2waterbodies_exp.tif";(test=raster(test)); x11(); plot(test)

#' J28
#' b028_dist2waterbodies_exp.
hf <-"D:/GISbases/b028JaguarAlbersBuffer70/b028_dist2waterbodies_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J28trk))
#crs(hf)   
#compareCRS(hf,get_crs(J28trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J28/UTM_b028_dist2waterbodies_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J28/UTM_b028_dist2waterbodies_exp.tif";(test=raster(test)); x11(); plot(test)

#' J29
#' b029_dist2waterbodies_exp.
hf <-"D:/GISbases/b029JaguarAlbersBuffer70/b029_dist2waterbodies_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J29trk))
#crs(hf)   
#compareCRS(hf,get_crs(J29trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J29/UTM_b029_dist2waterbodies_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J29/UTM_b029_dist2waterbodies_exp.tif";(test=raster(test)); x11(); plot(test)

#' J30
#' b030_dist2waterbodies_exp.
hf <-"D:/GISbases/b030JaguarAlbersBuffer70/b030_dist2waterbodies_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J30trk))
#crs(hf)   
#compareCRS(hf,get_crs(J30trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J30/UTM_b030_dist2waterbodies_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J30/UTM_b030_dist2waterbodies_exp.tif";(test=raster(test)); x11(); plot(test)

#' J30
#' b030_dist2waterbodies_exp.
hf <-"D:/GISbases/b030JaguarAlbersBuffer70/b030_dist2waterbodies_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J30trk))
#crs(hf)   
#compareCRS(hf,get_crs(J30trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J30/UTM_b030_dist2waterbodies_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J30/UTM_b030_dist2waterbodies_exp.tif";(test=raster(test)); x11(); plot(test)

#' J31
#' b031_dist2waterbodies_exp.
hf <-"D:/GISbases/b031JaguarAlbersBuffer70/b031_dist2waterbodies_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J31trk))
#crs(hf)   
#compareCRS(hf,get_crs(J31trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J31/UTM_b031_dist2waterbodies_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J31/UTM_b031_dist2waterbodies_exp.tif";(test=raster(test)); x11(); plot(test)

#' J32
#' b032_dist2waterbodies_exp.
hf <-"D:/GISbases/b032JaguarAlbersBuffer70/b032_dist2waterbodies_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J32trk))
#crs(hf)   
#compareCRS(hf,get_crs(J32trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J32/UTM_b032_dist2waterbodies_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J32/UTM_b032_dist2waterbodies_exp.tif";(test=raster(test)); x11(); plot(test)

#' J33
#' b033_dist2waterbodies_exp.
hf <-"D:/GISbases/b033JaguarAlbersBuffer70/b033_dist2waterbodies_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J33trk))
#crs(hf)   
#compareCRS(hf,get_crs(J33trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J33/UTM_b033_dist2waterbodies_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J33/UTM_b033_dist2waterbodies_exp.tif";(test=raster(test)); x11(); plot(test)

#' J34
#' b034_dist2waterbodies_exp.
hf <-"D:/GISbases/b034JaguarAlbersBuffer70/b034_dist2waterbodies_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J34trk))
#crs(hf)   
#compareCRS(hf,get_crs(J34trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J34/UTM_b034_dist2waterbodies_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J34/UTM_b034_dist2waterbodies_exp.tif";(test=raster(test)); x11(); plot(test)

#' J35
#' b035_dist2waterbodies_exp.
hf <-"D:/GISbases/b035JaguarAlbersBuffer70/b035_dist2waterbodies_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J35trk))
#crs(hf)   
#compareCRS(hf,get_crs(J35trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J35/UTM_b035_dist2waterbodies_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J35/UTM_b035_dist2waterbodies_exp.tif";(test=raster(test)); x11(); plot(test)

#' J36
#' b036_dist2waterbodies_exp.
hf <-"D:/GISbases/b036JaguarAlbersBuffer70/b036_dist2waterbodies_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J36trk))
#crs(hf)   
#compareCRS(hf,get_crs(J36trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J36/UTM_b036_dist2waterbodies_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J36/UTM_b036_dist2waterbodies_exp.tif";(test=raster(test)); x11(); plot(test)

#' J37
#' b037_dist2waterbodies_exp.
hf <-"D:/GISbases/b037JaguarAlbersBuffer70/b037_dist2waterbodies_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J37trk))
#crs(hf)   
#compareCRS(hf,get_crs(J37trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J37/UTM_b037_dist2waterbodies_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J37/UTM_b037_dist2waterbodies_exp.tif";(test=raster(test)); x11(); plot(test)

#' J38
#' b038_dist2waterbodies_exp.
hf <-"D:/GISbases/b038JaguarAlbersBuffer70/b038_dist2waterbodies_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J38trk))
#crs(hf)   
#compareCRS(hf,get_crs(J38trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J38/UTM_b038_dist2waterbodies_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J38/UTM_b038_dist2waterbodies_exp.tif";(test=raster(test)); x11(); plot(test)

#' J39
#' b039_dist2waterbodies_exp.
hf <-"D:/GISbases/b039JaguarAlbersBuffer70/b039_dist2waterbodies_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J39trk))
#crs(hf)   
#compareCRS(hf,get_crs(J39trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J39/UTM_b039_dist2waterbodies_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J39/UTM_b039_dist2waterbodies_exp.tif";(test=raster(test)); x11(); plot(test)

#' J40
#' b040_dist2waterbodies_exp.
hf <-"D:/GISbases/b040JaguarAlbersBuffer70/b040_dist2waterbodies_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J40trk))
#crs(hf)   
#compareCRS(hf,get_crs(J40trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J40/UTM_b040_dist2waterbodies_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J40/UTM_b040_dist2waterbodies_exp.tif";(test=raster(test)); x11(); plot(test)

#' J41
#' b041_dist2waterbodies_exp.
hf <-"D:/GISbases/b041JaguarAlbersBuffer70/b041_dist2waterbodies_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J41trk))
#crs(hf)   
#compareCRS(hf,get_crs(J41trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J41/UTM_b041_dist2waterbodies_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J41/UTM_b041_dist2waterbodies_exp.tif";(test=raster(test)); x11(); plot(test)

#' J42
#' b042_dist2waterbodies_exp.
hf <-"D:/GISbases/b042JaguarAlbersBuffer70/b042_dist2waterbodies_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J42trk))
#crs(hf)   
#compareCRS(hf,get_crs(J42trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J42/UTM_b042_dist2waterbodies_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J42/UTM_b042_dist2waterbodies_exp.tif";(test=raster(test)); x11(); plot(test)

gc()
gc()
gc()

#' J43
#' b043_dist2waterbodies_exp.
hf <-"D:/GISbases/b043JaguarAlbersBuffer70/b043_dist2waterbodies_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J43trk))
#crs(hf)   
#compareCRS(hf,get_crs(J43trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J43/UTM_b043_dist2waterbodies_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J43/UTM_b043_dist2waterbodies_exp.tif";(test=raster(test)); x11(); plot(test)

#' J44
#' b044_dist2waterbodies_exp.
hf <-"D:/GISbases/b044JaguarAlbersBuffer70/b044_dist2waterbodies_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J44trk))
#crs(hf)   
#compareCRS(hf,get_crs(J44trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J44/UTM_b044_dist2waterbodies_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J44/UTM_b044_dist2waterbodies_exp.tif";(test=raster(test)); x11(); plot(test)

#' J45
#' b045_dist2waterbodies_exp.
hf <-"D:/GISbases/b045JaguarAlbersBuffer70/b045_dist2waterbodies_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J45trk))
#crs(hf)   
#compareCRS(hf,get_crs(J45trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J45/UTM_b045_dist2waterbodies_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J45/UTM_b045_dist2waterbodies_exp.tif";(test=raster(test)); x11(); plot(test)

#' J46
#' b046_dist2waterbodies_exp.
hf <-"D:/GISbases/b046JaguarAlbersBuffer70/b046_dist2waterbodies_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J46trk))
#crs(hf)   
#compareCRS(hf,get_crs(J46trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J46/UTM_b046_dist2waterbodies_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J46/UTM_b046_dist2waterbodies_exp.tif";(test=raster(test)); x11(); plot(test)

#' J47
#' b047_dist2waterbodies_exp.
hf <-"D:/GISbases/b047JaguarAlbersBuffer70/b047_dist2waterbodies_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J47trk))
#crs(hf)   
#compareCRS(hf,get_crs(J47trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J47/UTM_b047_dist2waterbodies_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J47/UTM_b047_dist2waterbodies_exp.tif";(test=raster(test)); x11(); plot(test)

#' J48
#' b048_dist2waterbodies_exp.
hf <-"D:/GISbases/b048JaguarAlbersBuffer70/b048_dist2waterbodies_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J48trk))
#crs(hf)   
#compareCRS(hf,get_crs(J48trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J48/UTM_b048_dist2waterbodies_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J48/UTM_b048_dist2waterbodies_exp.tif";(test=raster(test)); x11(); plot(test)

#' J49
#' b049_dist2waterbodies_exp.
hf <-"D:/GISbases/b049JaguarAlbersBuffer70/b049_dist2waterbodies_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J49trk))
#crs(hf)   
#compareCRS(hf,get_crs(J49trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J49/UTM_b049_dist2waterbodies_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J49/UTM_b049_dist2waterbodies_exp.tif";(test=raster(test)); x11(); plot(test)

#' J50
#' b050_dist2waterbodies_exp.
hf <-"D:/GISbases/b050JaguarAlbersBuffer70/b050_dist2waterbodies_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J50trk))
#crs(hf)   
#compareCRS(hf,get_crs(J50trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J50/UTM_b050_dist2waterbodies_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J50/UTM_b050_dist2waterbodies_exp.tif";(test=raster(test)); x11(); plot(test)

#' J51
#' b051_dist2waterbodies_exp.
hf <-"D:/GISbases/b051JaguarAlbersBuffer70/b051_dist2waterbodies_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J51trk))
#crs(hf)   
#compareCRS(hf,get_crs(J51trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J51/UTM_b051_dist2waterbodies_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J51/UTM_b051_dist2waterbodies_exp.tif";(test=raster(test)); x11(); plot(test)

#' J52
#' b052_dist2waterbodies_exp.
hf <-"D:/GISbases/b052JaguarAlbersBuffer70/b052_dist2waterbodies_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J52trk))
#crs(hf)   
#compareCRS(hf,get_crs(J52trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J52/UTM_b052_dist2waterbodies_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J52/UTM_b052_dist2waterbodies_exp.tif";(test=raster(test)); x11(); plot(test)

#' J53
#' b053_dist2waterbodies_exp.
hf <-"D:/GISbases/b053JaguarAlbersBuffer70/b053_dist2waterbodies_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J53trk))
#crs(hf)   
#compareCRS(hf,get_crs(J53trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J53/UTM_b053_dist2waterbodies_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J53/UTM_b053_dist2waterbodies_exp.tif";(test=raster(test)); x11(); plot(test)

#' J54
#' b054_dist2waterbodies_exp.
hf <-"D:/GISbases/b054JaguarAlbersBuffer70/b054_dist2waterbodies_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J54trk))
#crs(hf)   
#compareCRS(hf,get_crs(J54trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J54/UTM_b054_dist2waterbodies_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J54/UTM_b054_dist2waterbodies_exp.tif";(test=raster(test)); x11(); plot(test)

#' J55
#' b055_dist2waterbodies_exp.
hf <-"D:/GISbases/b055JaguarAlbersBuffer70/b055_dist2waterbodies_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J55trk))
#crs(hf)   
#compareCRS(hf,get_crs(J55trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J55/UTM_b055_dist2waterbodies_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J55/UTM_b055_dist2waterbodies_exp.tif";(test=raster(test)); x11(); plot(test)

#' J56
#' b056_dist2waterbodies_exp.
hf <-"D:/GISbases/b056JaguarAlbersBuffer70/b056_dist2waterbodies_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J56trk))
#crs(hf)   
#compareCRS(hf,get_crs(J56trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J56/UTM_b056_dist2waterbodies_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J56/UTM_b056_dist2waterbodies_exp.tif";(test=raster(test)); x11(); plot(test)

#' J57
#' b057_dist2waterbodies_exp.
hf <-"D:/GISbases/b057JaguarAlbersBuffer70/b057_dist2waterbodies_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J57trk))
#crs(hf)   
#compareCRS(hf,get_crs(J57trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J57/UTM_b057_dist2waterbodies_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J57/UTM_b057_dist2waterbodies_exp.tif";(test=raster(test)); x11(); plot(test)

#' J58
#' b058_dist2waterbodies_exp.
hf <-"D:/GISbases/b058JaguarAlbersBuffer70/b058_dist2waterbodies_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J58trk))
#crs(hf)   
#compareCRS(hf,get_crs(J58trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J58/UTM_b058_dist2waterbodies_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J58/UTM_b058_dist2waterbodies_exp.tif";(test=raster(test)); x11(); plot(test)

#' J59
#' b059_dist2waterbodies_exp.
hf <-"D:/GISbases/b059JaguarAlbersBuffer70/b059_dist2waterbodies_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J59trk))
#crs(hf)   
#compareCRS(hf,get_crs(J59trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J59/UTM_b059_dist2waterbodies_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J59/UTM_b059_dist2waterbodies_exp.tif";(test=raster(test)); x11(); plot(test)

#' J60
#' b060_dist2waterbodies_exp.
hf <-"D:/GISbases/b060JaguarAlbersBuffer70/b060_dist2waterbodies_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J60trk))
#crs(hf)   
#compareCRS(hf,get_crs(J60trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J60/UTM_b060_dist2waterbodies_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J60/UTM_b060_dist2waterbodies_exp.tif";(test=raster(test)); x11(); plot(test)

#' J61
#' b061_dist2waterbodies_exp.
hf <-"D:/GISbases/b061JaguarAlbersBuffer70/b061_dist2waterbodies_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J61trk))
#crs(hf)   
#compareCRS(hf,get_crs(J61trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J61/UTM_b061_dist2waterbodies_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J61/UTM_b061_dist2waterbodies_exp.tif";(test=raster(test)); x11(); plot(test)

#' J62
#' b062_dist2waterbodies_exp.
hf <-"D:/GISbases/b062JaguarAlbersBuffer70/b062_dist2waterbodies_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J62trk))
#crs(hf)   
#compareCRS(hf,get_crs(J62trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J62/UTM_b062_dist2waterbodies_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J62/UTM_b062_dist2waterbodies_exp.tif";(test=raster(test)); x11(); plot(test)

#' J63
#' b063_dist2waterbodies_exp.
hf <-"D:/GISbases/b063JaguarAlbersBuffer70/b063_dist2waterbodies_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J63trk))
#crs(hf)   
#compareCRS(hf,get_crs(J63trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J63/UTM_b063_dist2waterbodies_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J63/UTM_b063_dist2waterbodies_exp.tif";(test=raster(test)); x11(); plot(test)

gc()
gc()
gc()


#' J64
#' b064_dist2waterbodies_exp.
hf <-"D:/GISbases/b064JaguarAlbersBuffer70/b064_dist2waterbodies_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J64trk))
#crs(hf)   
#compareCRS(hf,get_crs(J64trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J64/UTM_b064_dist2waterbodies_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J64/UTM_b064_dist2waterbodies_exp.tif";(test=raster(test)); x11(); plot(test)

#' J65
#' b065_dist2waterbodies_exp.
hf <-"D:/GISbases/b065JaguarAlbersBuffer70/b065_dist2waterbodies_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J65trk))
#crs(hf)   
#compareCRS(hf,get_crs(J65trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J65/UTM_b065_dist2waterbodies_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J65/UTM_b065_dist2waterbodies_exp.tif";(test=raster(test)); x11(); plot(test)

#' J66
#' b066_dist2waterbodies_exp.
hf <-"D:/GISbases/b066JaguarAlbersBuffer70/b066_dist2waterbodies_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J66trk))
#crs(hf)   
#compareCRS(hf,get_crs(J66trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J66/UTM_b066_dist2waterbodies_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J66/UTM_b066_dist2waterbodies_exp.tif";(test=raster(test)); x11(); plot(test)

#' J67
#' b067_dist2waterbodies_exp.
hf <-"D:/GISbases/b067JaguarAlbersBuffer70/b067_dist2waterbodies_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J67trk))
#crs(hf)   
#compareCRS(hf,get_crs(J67trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J67/UTM_b067_dist2waterbodies_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J67/UTM_b067_dist2waterbodies_exp.tif";(test=raster(test)); x11(); plot(test)

#' J68
#' b068_dist2waterbodies_exp.
hf <-"D:/GISbases/b068JaguarAlbersBuffer70/b068_dist2waterbodies_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J68trk))
#crs(hf)   
#compareCRS(hf,get_crs(J68trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J68/UTM_b068_dist2waterbodies_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J68/UTM_b068_dist2waterbodies_exp.tif";(test=raster(test)); x11(); plot(test)

#' J69
#' b069_dist2waterbodies_exp.
hf <-"D:/GISbases/b069JaguarAlbersBuffer70/b069_dist2waterbodies_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J69trk))
#crs(hf)   
#compareCRS(hf,get_crs(J69trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J69/UTM_b069_dist2waterbodies_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J69/UTM_b069_dist2waterbodies_exp.tif";(test=raster(test)); x11(); plot(test)

gc()
gc()
gc()


#' J70
#' b070_dist2waterbodies_exp.
hf <-"D:/GISbases/b070JaguarAlbersBuffer70/b070_dist2waterbodies_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J70trk))
#crs(hf)   
#compareCRS(hf,get_crs(J70trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J70/UTM_b070_dist2waterbodies_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J70/UTM_b070_dist2waterbodies_exp.tif";(test=raster(test)); x11(); plot(test)

#' J71
#' b071_dist2waterbodies_exp.
hf <-"D:/GISbases/b071JaguarAlbersBuffer70/b071_dist2waterbodies_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J71trk))
#crs(hf)   
#compareCRS(hf,get_crs(J71trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J71/UTM_b071_dist2waterbodies_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J71/UTM_b071_dist2waterbodies_exp.tif";(test=raster(test)); x11(); plot(test)

#' J72
#' b072_dist2waterbodies_exp.
hf <-"D:/GISbases/b072JaguarAlbersBuffer70/b072_dist2waterbodies_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J72trk))
#crs(hf)   
#compareCRS(hf,get_crs(J72trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J72/UTM_b072_dist2waterbodies_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J72/UTM_b072_dist2waterbodies_exp.tif";(test=raster(test)); x11(); plot(test)

#' J73
#' b073_dist2waterbodies_exp.
hf <-"D:/GISbases/b073JaguarAlbersBuffer70/b073_dist2waterbodies_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J73trk))
#crs(hf)   
#compareCRS(hf,get_crs(J73trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J73/UTM_b073_dist2waterbodies_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J73/UTM_b073_dist2waterbodies_exp.tif";(test=raster(test)); x11(); plot(test)

#' J74
#' b074_dist2waterbodies_exp.
hf <-"D:/GISbases/b074JaguarAlbersBuffer70/b074_dist2waterbodies_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J74trk))
#crs(hf)   
#compareCRS(hf,get_crs(J74trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J74/UTM_b074_dist2waterbodies_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J74/UTM_b074_dist2waterbodies_exp.tif";(test=raster(test)); x11(); plot(test)

#' J75
#' b075_dist2waterbodies_exp.
hf <-"D:/GISbases/b075JaguarAlbersBuffer70/b075_dist2waterbodies_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J75trk))
#crs(hf)   
#compareCRS(hf,get_crs(J75trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J75/UTM_b075_dist2waterbodies_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J75/UTM_b075_dist2waterbodies_exp.tif";(test=raster(test)); x11(); plot(test)

#' J76
#' b076_dist2waterbodies_exp.
hf <-"D:/GISbases/b076JaguarAlbersBuffer70/b076_dist2waterbodies_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J76trk))
#crs(hf)   
#compareCRS(hf,get_crs(J76trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J76/UTM_b076_dist2waterbodies_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J76/UTM_b076_dist2waterbodies_exp.tif";(test=raster(test)); x11(); plot(test)

#' J77
#' b077_dist2waterbodies_exp.
hf <-"D:/GISbases/b077JaguarAlbersBuffer70/b077_dist2waterbodies_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J77trk))
#crs(hf)   
#compareCRS(hf,get_crs(J77trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J77/UTM_b077_dist2waterbodies_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J77/UTM_b077_dist2waterbodies_exp.tif";(test=raster(test)); x11(); plot(test)

#' J78
#' b078_dist2waterbodies_exp.
hf <-"D:/GISbases/b078JaguarAlbersBuffer70/b078_dist2waterbodies_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J78trk))
#crs(hf)   
#compareCRS(hf,get_crs(J78trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J78/UTM_b078_dist2waterbodies_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J78/UTM_b078_dist2waterbodies_exp.tif";(test=raster(test)); x11(); plot(test)

#' J79
#' b079_dist2waterbodies_exp.
hf <-"D:/GISbases/b079JaguarAlbersBuffer70/b079_dist2waterbodies_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J79trk))
#crs(hf)   
#compareCRS(hf,get_crs(J79trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J79/UTM_b079_dist2waterbodies_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J79/UTM_b079_dist2waterbodies_exp.tif";(test=raster(test)); x11(); plot(test)

#' J80
#' b080_dist2waterbodies_exp.
hf <-"D:/GISbases/b080JaguarAlbersBuffer70/b080_dist2waterbodies_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J80trk))
#crs(hf)   
#compareCRS(hf,get_crs(J80trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J80/UTM_b080_dist2waterbodies_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J80/UTM_b080_dist2waterbodies_exp.tif";(test=raster(test)); x11(); plot(test)

#' J81
#' b081_dist2waterbodies_exp.
hf <-"D:/GISbases/b081JaguarAlbersBuffer70/b081_dist2waterbodies_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J81trk))
#crs(hf)   
#compareCRS(hf,get_crs(J81trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J81/UTM_b081_dist2waterbodies_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J81/UTM_b081_dist2waterbodies_exp.tif";(test=raster(test)); x11(); plot(test)

#' J82
#' b082_dist2waterbodies_exp.
hf <-"D:/GISbases/b082JaguarAlbersBuffer70/b082_dist2waterbodies_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J82trk))
#crs(hf)   
#compareCRS(hf,get_crs(J82trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J82/UTM_b082_dist2waterbodies_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J82/UTM_b082_dist2waterbodies_exp.tif";(test=raster(test)); x11(); plot(test)

#' J83
#' b083_dist2waterbodies_exp.
hf <-"D:/GISbases/b083JaguarAlbersBuffer70/b083_dist2waterbodies_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J83trk))
#crs(hf)   
#compareCRS(hf,get_crs(J83trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J83/UTM_b083_dist2waterbodies_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J83/UTM_b083_dist2waterbodies_exp.tif";(test=raster(test)); x11(); plot(test)

#' J84
#' b084_dist2waterbodies_exp.
hf <-"D:/GISbases/b084JaguarAlbersBuffer70/b084_dist2waterbodies_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J84trk))
#crs(hf)   
#compareCRS(hf,get_crs(J84trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J84/UTM_b084_dist2waterbodies_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J84/UTM_b084_dist2waterbodies_exp.tif";(test=raster(test)); x11(); plot(test)

#' J85
#' b085_dist2waterbodies_exp.
hf <-"D:/GISbases/b085JaguarAlbersBuffer70/b085_dist2waterbodies_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J85trk))
#crs(hf)   
#compareCRS(hf,get_crs(J85trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J85/UTM_b085_dist2waterbodies_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J85/UTM_b085_dist2waterbodies_exp.tif";(test=raster(test)); x11(); plot(test)

#' J86
#' b086_dist2waterbodies_exp.
hf <-"D:/GISbases/b086JaguarAlbersBuffer70/b086_dist2waterbodies_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J86trk))
#crs(hf)   
#compareCRS(hf,get_crs(J86trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J86/UTM_b086_dist2waterbodies_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J86/UTM_b086_dist2waterbodies_exp.tif";(test=raster(test)); x11(); plot(test)

#' J87
#' b087_dist2waterbodies_exp.
hf <-"D:/GISbases/b087JaguarAlbersBuffer70/b087_dist2waterbodies_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J87trk))
#crs(hf)   
#compareCRS(hf,get_crs(J87trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J87/UTM_b087_dist2waterbodies_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J87/UTM_b087_dist2waterbodies_exp.tif";(test=raster(test)); x11(); plot(test)

#' J88
#' b088_dist2waterbodies_exp.
hf <-"D:/GISbases/b088JaguarAlbersBuffer70/b088_dist2waterbodies_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J88trk))
#crs(hf)   
#compareCRS(hf,get_crs(J88trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J88/UTM_b088_dist2waterbodies_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J88/UTM_b088_dist2waterbodies_exp.tif";(test=raster(test)); x11(); plot(test)

#' J89
#' b089_dist2waterbodies_exp.
hf <-"D:/GISbases/b089JaguarAlbersBuffer70/b089_dist2waterbodies_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J89trk))
#crs(hf)   
#compareCRS(hf,get_crs(J89trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J89/UTM_b089_dist2waterbodies_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J89/UTM_b089_dist2waterbodies_exp.tif";(test=raster(test)); x11(); plot(test)

#' J90
#' b090_dist2waterbodies_exp.
hf <-"D:/GISbases/b090JaguarAlbersBuffer70/b090_dist2waterbodies_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J90trk))
#crs(hf)   
#compareCRS(hf,get_crs(J90trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J90/UTM_b090_dist2waterbodies_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J90/UTM_b090_dist2waterbodies_exp.tif";(test=raster(test)); x11(); plot(test)

#' J91
#' b091_dist2waterbodies_exp.
hf <-"D:/GISbases/b091JaguarAlbersBuffer70/b091_dist2waterbodies_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J91trk))
#crs(hf)   
#compareCRS(hf,get_crs(J91trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J91/UTM_b091_dist2waterbodies_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J91/UTM_b091_dist2waterbodies_exp.tif";(test=raster(test)); x11(); plot(test)

#' J92
#' b092_dist2waterbodies_exp.
hf <-"D:/GISbases/b092JaguarAlbersBuffer70/b092_dist2waterbodies_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J92trk))
#crs(hf)   
#compareCRS(hf,get_crs(J92trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J92/UTM_b092_dist2waterbodies_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J92/UTM_b092_dist2waterbodies_exp.tif";(test=raster(test)); x11(); plot(test)

#' J93
#' b093_dist2waterbodies_exp.
hf <-"D:/GISbases/b093JaguarAlbersBuffer70/b093_dist2waterbodies_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J93trk))
#crs(hf)   
#compareCRS(hf,get_crs(J93trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J93/UTM_b093_dist2waterbodies_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J93/UTM_b093_dist2waterbodies_exp.tif";(test=raster(test)); x11(); plot(test)

#' J94
#' b094_dist2waterbodies_exp.
hf <-"D:/GISbases/b094JaguarAlbersBuffer70/b094_dist2waterbodies_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J94trk))
#crs(hf)   
#compareCRS(hf,get_crs(J94trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J94/UTM_b094_dist2waterbodies_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J94/UTM_b094_dist2waterbodies_exp.tif";(test=raster(test)); x11(); plot(test)

#' J95
#' b095_dist2waterbodies_exp.
hf <-"D:/GISbases/b095JaguarAlbersBuffer70/b095_dist2waterbodies_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J95trk))
#crs(hf)   
#compareCRS(hf,get_crs(J95trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J95/UTM_b095_dist2waterbodies_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J95/UTM_b095_dist2waterbodies_exp.tif";(test=raster(test)); x11(); plot(test)

#' J96
#' b096_dist2waterbodies_exp.
hf <-"D:/GISbases/b096JaguarAlbersBuffer70/b096_dist2waterbodies_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J96trk))
#crs(hf)   
#compareCRS(hf,get_crs(J96trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J96/UTM_b096_dist2waterbodies_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J96/UTM_b096_dist2waterbodies_exp.tif";(test=raster(test)); x11(); plot(test)

#' J97
#' b097_dist2waterbodies_exp.
hf <-"D:/GISbases/b097JaguarAlbersBuffer70/b097_dist2waterbodies_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J97trk))
#crs(hf)   
#compareCRS(hf,get_crs(J97trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J97/UTM_b097_dist2waterbodies_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J97/UTM_b097_dist2waterbodies_exp.tif";(test=raster(test)); x11(); plot(test)

#' J98
#' b098_dist2waterbodies_exp.
hf <-"D:/GISbases/b098JaguarAlbersBuffer70/b098_dist2waterbodies_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J98trk))
#crs(hf)   
#compareCRS(hf,get_crs(J98trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J98/UTM_b098_dist2waterbodies_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J98/UTM_b098_dist2waterbodies_exp.tif";(test=raster(test)); x11(); plot(test)

#' J99
#' b099_dist2waterbodies_exp.
hf <-"D:/GISbases/b099JaguarAlbersBuffer70/b099_dist2waterbodies_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J99trk))
#crs(hf)   
#compareCRS(hf,get_crs(J99trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J99/UTM_b099_dist2waterbodies_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J99/UTM_b099_dist2waterbodies_exp.tif";(test=raster(test)); x11(); plot(test)

#' J100
#' b100_dist2waterbodies_exp.
hf <-"D:/GISbases/b100JaguarAlbersBuffer70/b100_dist2waterbodies_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J100trk))
#crs(hf)   
#compareCRS(hf,get_crs(J100trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J100/UTM_b100_dist2waterbodies_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J100/UTM_b010_dist2waterbodies_exp.tif";(test=raster(test)); x11(); plot(test)

#' J101
#' b101_dist2waterbodies_exp.
hf <-"D:/GISbases/b101JaguarAlbersBuffer70/b101_dist2waterbodies_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J101trk))
#crs(hf)   
#compareCRS(hf,get_crs(J101trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J101/UTM_b101_dist2waterbodies_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J101/UTM_b101_dist2waterbodies_exp.tif";(test=raster(test)); x11(); plot(test)

#' J102
#' b102_dist2waterbodies_exp.
hf <-"D:/GISbases/b102JaguarAlbersBuffer70/b102_dist2waterbodies_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J102trk))
#crs(hf)   
#compareCRS(hf,get_crs(J102trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J102/UTM_b102_dist2waterbodies_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J102/UTM_b102_dist2waterbodies_exp.tif";(test=raster(test)); x11(); plot(test)

#' J103
#' b103_dist2waterbodies_exp.
hf <-"D:/GISbases/b103JaguarAlbersBuffer70/b103_dist2waterbodies_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J103trk))
#crs(hf)   
#compareCRS(hf,get_crs(J103trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J103/UTM_b103_dist2waterbodies_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J103/UTM_b103_dist2waterbodies_exp.tif";(test=raster(test)); x11(); plot(test)

#' J104
#' b104_dist2waterbodies_exp.
hf <-"D:/GISbases/b104JaguarAlbersBuffer70/b104_dist2waterbodies_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J104trk))
#crs(hf)   
#compareCRS(hf,get_crs(J104trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J104/UTM_b104_dist2waterbodies_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J104/UTM_b104_dist2waterbodies_exp.tif";(test=raster(test)); x11(); plot(test)

#' J105
#' b105_dist2waterbodies_exp.
hf <-"D:/GISbases/b105JaguarAlbersBuffer70/b105_dist2waterbodies_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J105trk))
#crs(hf)   
#compareCRS(hf,get_crs(J105trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J105/UTM_b105_dist2waterbodies_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J105/UTM_b105_dist2waterbodies_exp.tif";(test=raster(test)); x11(); plot(test)

#' J106
#' b106_dist2waterbodies_exp.
hf <-"D:/GISbases/b106JaguarAlbersBuffer70/b106_dist2waterbodies_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J106trk))
#crs(hf)   
#compareCRS(hf,get_crs(J106trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J106/UTM_b106_dist2waterbodies_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J106/UTM_b106_dist2waterbodies_exp.tif";(test=raster(test)); x11(); plot(test)

#' J107
#' b107_dist2waterbodies_exp.
hf <-"D:/GISbases/b107JaguarAlbersBuffer70/b107_dist2waterbodies_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J107trk))
#crs(hf)   
#compareCRS(hf,get_crs(J107trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J107/UTM_b107_dist2waterbodies_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J107/UTM_b107_dist2waterbodies_exp.tif";(test=raster(test)); x11(); plot(test)

#' J108
#' b108_dist2waterbodies_exp.
hf <-"D:/GISbases/b108JaguarAlbersBuffer70/b108_dist2waterbodies_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J108trk))
#crs(hf)   
#compareCRS(hf,get_crs(J108trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J108/UTM_b108_dist2waterbodies_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J108/UTM_b108_dist2waterbodies_exp.tif";(test=raster(test)); x11(); plot(test)

#' J109
#' b109_dist2waterbodies_exp.
hf <-"D:/GISbases/b109JaguarAlbersBuffer70/b109_dist2waterbodies_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J109trk))
#crs(hf)   
#compareCRS(hf,get_crs(J109trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J109/UTM_b109_dist2waterbodies_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J109/UTM_b109_dist2waterbodies_exp.tif";(test=raster(test)); x11(); plot(test)

#' J110
#' b110_dist2waterbodies_exp.
hf <-"D:/GISbases/b110JaguarAlbersBuffer70/b110_dist2waterbodies_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J110trk))
#crs(hf)   
#compareCRS(hf,get_crs(J110trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J110/UTM_b110_dist2waterbodies_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J110/UTM_b110_dist2waterbodies_exp.tif";(test=raster(test)); x11(); plot(test)

#' J111
#' b111_dist2waterbodies_exp.
hf <-"D:/GISbases/b111JaguarAlbersBuffer70/b111_dist2waterbodies_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J111trk))
#crs(hf)   
#compareCRS(hf,get_crs(J111trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J111/UTM_b111_dist2waterbodies_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J111/UTM_b111_dist2waterbodies_exp.tif";(test=raster(test)); x11(); plot(test)

#' J112
#' b112_dist2waterbodies_exp.
hf <-"D:/GISbases/b112JaguarAlbersBuffer70/b112_dist2waterbodies_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J112trk))
#crs(hf)   
#compareCRS(hf,get_crs(J112trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J112/UTM_b112_dist2waterbodies_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J112/UTM_b112_dist2waterbodies_exp.tif";(test=raster(test)); x11(); plot(test)

#' J113
#' b113_dist2waterbodies_exp.
hf <-"D:/GISbases/b113JaguarAlbersBuffer70/b113_dist2waterbodies_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J113trk))
#crs(hf)   
#compareCRS(hf,get_crs(J113trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J113/UTM_b113_dist2waterbodies_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J113/UTM_b113_dist2waterbodies_exp.tif";(test=raster(test)); x11(); plot(test)

#' J114
#' b114_dist2waterbodies_exp.
hf <-"D:/GISbases/b114JaguarAlbersBuffer70/b114_dist2waterbodies_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J114trk))
#crs(hf)   
#compareCRS(hf,get_crs(J114trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J114/UTM_b114_dist2waterbodies_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J114/UTM_b114_dist2waterbodies_exp.tif";(test=raster(test)); x11(); plot(test)

#' J115
#' b115_dist2waterbodies_exp.
hf <-"D:/GISbases/b115JaguarAlbersBuffer70/b115_dist2waterbodies_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J115trk))
#crs(hf)   
#compareCRS(hf,get_crs(J115trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J115/UTM_b115_dist2waterbodies_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J115/UTM_b115_dist2waterbodies_exp.tif";(test=raster(test)); x11(); plot(test)

#' J116
#' b116_dist2waterbodies_exp.
hf <-"D:/GISbases/b116JaguarAlbersBuffer70/b116_dist2waterbodies_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J116trk))
#crs(hf)   
#compareCRS(hf,get_crs(J116trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J116/UTM_b116_dist2waterbodies_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J116/UTM_b116_dist2waterbodies_exp.tif";(test=raster(test)); x11(); plot(test)

#' J117
#' b117_dist2waterbodies_exp.
hf <-"D:/GISbases/b117JaguarAlbersBuffer70/b117_dist2waterbodies_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J117trk))
#crs(hf)   
#compareCRS(hf,get_crs(J117trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J117/UTM_b117_dist2waterbodies_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J117/UTM_b117_dist2waterbodies_exp.tif";(test=raster(test)); x11(); plot(test)

####################################################################################################################################
# wc2.0_bio_30s_02_neotropic_albers_tif_exp.
# BIO2 = Mean Diurnal Range (Mean of monthly (max temp - min temp))
gc()
gc()
gc()

#' J1
#' b001_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif
hf <-"D:/GISbases/b001JaguarAlbersBuffer70/b001_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J1trk))
#crs(hf)   
#compareCRS(hf,get_crs(J1trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J1/UTM_b001_wc2.0_bio_30s_02_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J1/UTM_b001_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J2
#' b002_wc2.0_bio_30s_02_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b002JaguarAlbersBuffer70/b002_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J2trk))
#crs(hf)   
#compareCRS(hf,get_crs(J2trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J2/UTM_b002_wc2.0_bio_30s_02_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J2/UTM_b002_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J3
#' b003_wc2.0_bio_30s_02_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b003JaguarAlbersBuffer70/b003_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J3trk))
#crs(hf)   
#compareCRS(hf,get_crs(J3trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J3/UTM_b003_wc2.0_bio_30s_02_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J3/UTM_b003_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J4
#' b004_wc2.0_bio_30s_02_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b004JaguarAlbersBuffer70/b004_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J4trk))
#crs(hf)   
#compareCRS(hf,get_crs(J4trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J4/UTM_b004_wc2.0_bio_30s_02_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J4/UTM_b004_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J5
#' b005_wc2.0_bio_30s_02_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b005JaguarAlbersBuffer70/b005_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J5trk))
#crs(hf)   
#compareCRS(hf,get_crs(J5trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J5/UTM_b005_wc2.0_bio_30s_02_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J5/UTM_b005_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J6
#' b006_wc2.0_bio_30s_02_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b006JaguarAlbersBuffer70/b006_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J6trk))
#crs(hf)   
#compareCRS(hf,get_crs(J6trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J6/UTM_b006_wc2.0_bio_30s_02_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J6/UTM_b006_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J7
#' b007_wc2.0_bio_30s_02_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b007JaguarAlbersBuffer70/b007_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J7trk))
#crs(hf)   
#compareCRS(hf,get_crs(J7trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J7/UTM_b007_wc2.0_bio_30s_02_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J7/UTM_b007_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J8
#' b008_wc2.0_bio_30s_02_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b008JaguarAlbersBuffer70/b008_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J8trk))
#crs(hf)   
#compareCRS(hf,get_crs(J8trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J8/UTM_b008_wc2.0_bio_30s_02_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J8/UTM_b008_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J9
#' b009_wc2.0_bio_30s_02_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b009JaguarAlbersBuffer70/b009_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J9trk))
#crs(hf)   
#compareCRS(hf,get_crs(J9trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J9/UTM_b009_wc2.0_bio_30s_02_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J9/UTM_b009_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J10
#' b010_wc2.0_bio_30s_02_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b010JaguarAlbersBuffer70/b010_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J10trk))
#crs(hf)   
#compareCRS(hf,get_crs(J10trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J10/UTM_b010_wc2.0_bio_30s_02_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J10/UTM_b010_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

gc()
gc()
gc()


#' J11
#' b011_wc2.0_bio_30s_02_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b011JaguarAlbersBuffer70/b011_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J11trk))
#crs(hf)   
#compareCRS(hf,get_crs(J11trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J11/UTM_b011_wc2.0_bio_30s_02_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J11/UTM_b011_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J12
#' b012_wc2.0_bio_30s_02_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b012JaguarAlbersBuffer70/b012_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J12trk))
#crs(hf)   
#compareCRS(hf,get_crs(J12trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J12/UTM_b012_wc2.0_bio_30s_02_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J12/UTM_b012_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J13
#' b013_wc2.0_bio_30s_02_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b013JaguarAlbersBuffer70/b013_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J13trk))
#crs(hf)   
#compareCRS(hf,get_crs(J13trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J13/UTM_b013_wc2.0_bio_30s_02_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J13/UTM_b013_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J14
#' b014_wc2.0_bio_30s_02_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b014JaguarAlbersBuffer70/b014_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J14trk))
#crs(hf)   
#compareCRS(hf,get_crs(J14trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J14/UTM_b014_wc2.0_bio_30s_02_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J14/UTM_b014_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J15
#' b015_wc2.0_bio_30s_02_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b015JaguarAlbersBuffer70/b015_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J15trk))
#crs(hf)   
#compareCRS(hf,get_crs(J15trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J15/UTM_b015_wc2.0_bio_30s_02_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J15/UTM_b015_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J16
#' b016_wc2.0_bio_30s_02_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b016JaguarAlbersBuffer70/b016_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J16trk))
#crs(hf)   
#compareCRS(hf,get_crs(J16trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J16/UTM_b016_wc2.0_bio_30s_02_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J16/UTM_b016_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J17
#' b017_wc2.0_bio_30s_02_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b017JaguarAlbersBuffer70/b017_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J17trk))
#crs(hf)   
#compareCRS(hf,get_crs(J17trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J17/UTM_b017_wc2.0_bio_30s_02_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J17/UTM_b017_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J18
#' b018_wc2.0_bio_30s_02_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b018JaguarAlbersBuffer70/b018_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J18trk))
#crs(hf)   
#compareCRS(hf,get_crs(J18trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J18/UTM_b018_wc2.0_bio_30s_02_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J18/UTM_b018_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J19
#' b019_wc2.0_bio_30s_02_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b019JaguarAlbersBuffer70/b019_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J19trk))
#crs(hf)   
#compareCRS(hf,get_crs(J19trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J19/UTM_b019_wc2.0_bio_30s_02_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J19/UTM_b019_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J20
#' b020_wc2.0_bio_30s_02_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b020JaguarAlbersBuffer70/b020_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J20trk))
#crs(hf)   
#compareCRS(hf,get_crs(J20trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J20/UTM_b020_wc2.0_bio_30s_02_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J20/UTM_b020_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J21
#' b021_wc2.0_bio_30s_02_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b021JaguarAlbersBuffer70/b021_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J21trk))
#crs(hf)   
#compareCRS(hf,get_crs(J21trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J21/UTM_b021_wc2.0_bio_30s_02_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J21/UTM_b021_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J22
#' b022_wc2.0_bio_30s_02_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b022JaguarAlbersBuffer70/b022_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J22trk))
#crs(hf)   
#compareCRS(hf,get_crs(J22trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J22/UTM_b022_wc2.0_bio_30s_02_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J22/UTM_b022_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J23
#' b023_wc2.0_bio_30s_02_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b023JaguarAlbersBuffer70/b023_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J23trk))
#crs(hf)   
#compareCRS(hf,get_crs(J23trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J23/UTM_b023_wc2.0_bio_30s_02_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J23/UTM_b023_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J24
#' b024_wc2.0_bio_30s_02_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b024JaguarAlbersBuffer70/b024_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J24trk))
#crs(hf)   
#compareCRS(hf,get_crs(J24trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J24/UTM_b024_wc2.0_bio_30s_02_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J24/UTM_b024_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J25
#' b025_wc2.0_bio_30s_02_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b025JaguarAlbersBuffer70/b025_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J25trk))
#crs(hf)   
#compareCRS(hf,get_crs(J25trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J25/UTM_b025_wc2.0_bio_30s_02_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J25/UTM_b025_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J26
#' b026_wc2.0_bio_30s_02_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b026JaguarAlbersBuffer70/b026_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J26trk))
#crs(hf)   
#compareCRS(hf,get_crs(J26trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J26/UTM_b026_wc2.0_bio_30s_02_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J26/UTM_b026_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J27
#' b027_wc2.0_bio_30s_02_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b027JaguarAlbersBuffer70/b027_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J27trk))
#crs(hf)   
#compareCRS(hf,get_crs(J27trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J27/UTM_b027_wc2.0_bio_30s_02_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J27/UTM_b027_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J28
#' b028_wc2.0_bio_30s_02_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b028JaguarAlbersBuffer70/b028_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J28trk))
#crs(hf)   
#compareCRS(hf,get_crs(J28trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J28/UTM_b028_wc2.0_bio_30s_02_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J28/UTM_b028_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J29
#' b029_wc2.0_bio_30s_02_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b029JaguarAlbersBuffer70/b029_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J29trk))
#crs(hf)   
#compareCRS(hf,get_crs(J29trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J29/UTM_b029_wc2.0_bio_30s_02_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J29/UTM_b029_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J30
#' b030_wc2.0_bio_30s_02_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b030JaguarAlbersBuffer70/b030_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J30trk))
#crs(hf)   
#compareCRS(hf,get_crs(J30trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J30/UTM_b030_wc2.0_bio_30s_02_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J30/UTM_b030_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J30
#' b030_wc2.0_bio_30s_02_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b030JaguarAlbersBuffer70/b030_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J30trk))
#crs(hf)   
#compareCRS(hf,get_crs(J30trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J30/UTM_b030_wc2.0_bio_30s_02_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J30/UTM_b030_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J31
#' b031_wc2.0_bio_30s_02_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b031JaguarAlbersBuffer70/b031_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J31trk))
#crs(hf)   
#compareCRS(hf,get_crs(J31trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J31/UTM_b031_wc2.0_bio_30s_02_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J31/UTM_b031_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J32
#' b032_wc2.0_bio_30s_02_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b032JaguarAlbersBuffer70/b032_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J32trk))
#crs(hf)   
#compareCRS(hf,get_crs(J32trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J32/UTM_b032_wc2.0_bio_30s_02_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J32/UTM_b032_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J33
#' b033_wc2.0_bio_30s_02_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b033JaguarAlbersBuffer70/b033_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J33trk))
#crs(hf)   
#compareCRS(hf,get_crs(J33trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J33/UTM_b033_wc2.0_bio_30s_02_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J33/UTM_b033_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J34
#' b034_wc2.0_bio_30s_02_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b034JaguarAlbersBuffer70/b034_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J34trk))
#crs(hf)   
#compareCRS(hf,get_crs(J34trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J34/UTM_b034_wc2.0_bio_30s_02_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J34/UTM_b034_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J35
#' b035_wc2.0_bio_30s_02_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b035JaguarAlbersBuffer70/b035_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J35trk))
#crs(hf)   
#compareCRS(hf,get_crs(J35trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J35/UTM_b035_wc2.0_bio_30s_02_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J35/UTM_b035_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J36
#' b036_wc2.0_bio_30s_02_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b036JaguarAlbersBuffer70/b036_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J36trk))
#crs(hf)   
#compareCRS(hf,get_crs(J36trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J36/UTM_b036_wc2.0_bio_30s_02_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J36/UTM_b036_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J37
#' b037_wc2.0_bio_30s_02_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b037JaguarAlbersBuffer70/b037_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J37trk))
#crs(hf)   
#compareCRS(hf,get_crs(J37trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J37/UTM_b037_wc2.0_bio_30s_02_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J37/UTM_b037_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J38
#' b038_wc2.0_bio_30s_02_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b038JaguarAlbersBuffer70/b038_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J38trk))
#crs(hf)   
#compareCRS(hf,get_crs(J38trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J38/UTM_b038_wc2.0_bio_30s_02_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J38/UTM_b038_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J39
#' b039_wc2.0_bio_30s_02_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b039JaguarAlbersBuffer70/b039_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J39trk))
#crs(hf)   
#compareCRS(hf,get_crs(J39trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J39/UTM_b039_wc2.0_bio_30s_02_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J39/UTM_b039_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J40
#' b040_wc2.0_bio_30s_02_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b040JaguarAlbersBuffer70/b040_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J40trk))
#crs(hf)   
#compareCRS(hf,get_crs(J40trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J40/UTM_b040_wc2.0_bio_30s_02_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J40/UTM_b040_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J41
#' b041_wc2.0_bio_30s_02_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b041JaguarAlbersBuffer70/b041_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J41trk))
#crs(hf)   
#compareCRS(hf,get_crs(J41trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J41/UTM_b041_wc2.0_bio_30s_02_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J41/UTM_b041_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J42
#' b042_wc2.0_bio_30s_02_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b042JaguarAlbersBuffer70/b042_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J42trk))
#crs(hf)   
#compareCRS(hf,get_crs(J42trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J42/UTM_b042_wc2.0_bio_30s_02_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J42/UTM_b042_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

gc()
gc()
gc()

#' J43
#' b043_wc2.0_bio_30s_02_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b043JaguarAlbersBuffer70/b043_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J43trk))
#crs(hf)   
#compareCRS(hf,get_crs(J43trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J43/UTM_b043_wc2.0_bio_30s_02_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J43/UTM_b043_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J44
#' b044_wc2.0_bio_30s_02_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b044JaguarAlbersBuffer70/b044_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J44trk))
#crs(hf)   
#compareCRS(hf,get_crs(J44trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J44/UTM_b044_wc2.0_bio_30s_02_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J44/UTM_b044_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J45
#' b045_wc2.0_bio_30s_02_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b045JaguarAlbersBuffer70/b045_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J45trk))
#crs(hf)   
#compareCRS(hf,get_crs(J45trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J45/UTM_b045_wc2.0_bio_30s_02_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J45/UTM_b045_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J46
#' b046_wc2.0_bio_30s_02_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b046JaguarAlbersBuffer70/b046_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J46trk))
#crs(hf)   
#compareCRS(hf,get_crs(J46trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J46/UTM_b046_wc2.0_bio_30s_02_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J46/UTM_b046_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J47
#' b047_wc2.0_bio_30s_02_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b047JaguarAlbersBuffer70/b047_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J47trk))
#crs(hf)   
#compareCRS(hf,get_crs(J47trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J47/UTM_b047_wc2.0_bio_30s_02_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J47/UTM_b047_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J48
#' b048_wc2.0_bio_30s_02_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b048JaguarAlbersBuffer70/b048_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J48trk))
#crs(hf)   
#compareCRS(hf,get_crs(J48trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J48/UTM_b048_wc2.0_bio_30s_02_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J48/UTM_b048_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J49
#' b049_wc2.0_bio_30s_02_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b049JaguarAlbersBuffer70/b049_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J49trk))
#crs(hf)   
#compareCRS(hf,get_crs(J49trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J49/UTM_b049_wc2.0_bio_30s_02_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J49/UTM_b049_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J50
#' b050_wc2.0_bio_30s_02_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b050JaguarAlbersBuffer70/b050_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J50trk))
#crs(hf)   
#compareCRS(hf,get_crs(J50trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J50/UTM_b050_wc2.0_bio_30s_02_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J50/UTM_b050_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J51
#' b051_wc2.0_bio_30s_02_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b051JaguarAlbersBuffer70/b051_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J51trk))
#crs(hf)   
#compareCRS(hf,get_crs(J51trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J51/UTM_b051_wc2.0_bio_30s_02_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J51/UTM_b051_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J52
#' b052_wc2.0_bio_30s_02_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b052JaguarAlbersBuffer70/b052_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J52trk))
#crs(hf)   
#compareCRS(hf,get_crs(J52trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J52/UTM_b052_wc2.0_bio_30s_02_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J52/UTM_b052_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J53
#' b053_wc2.0_bio_30s_02_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b053JaguarAlbersBuffer70/b053_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J53trk))
#crs(hf)   
#compareCRS(hf,get_crs(J53trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J53/UTM_b053_wc2.0_bio_30s_02_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J53/UTM_b053_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J54
#' b054_wc2.0_bio_30s_02_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b054JaguarAlbersBuffer70/b054_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J54trk))
#crs(hf)   
#compareCRS(hf,get_crs(J54trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J54/UTM_b054_wc2.0_bio_30s_02_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J54/UTM_b054_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J55
#' b055_wc2.0_bio_30s_02_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b055JaguarAlbersBuffer70/b055_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J55trk))
#crs(hf)   
#compareCRS(hf,get_crs(J55trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J55/UTM_b055_wc2.0_bio_30s_02_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J55/UTM_b055_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J56
#' b056_wc2.0_bio_30s_02_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b056JaguarAlbersBuffer70/b056_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J56trk))
#crs(hf)   
#compareCRS(hf,get_crs(J56trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J56/UTM_b056_wc2.0_bio_30s_02_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J56/UTM_b056_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J57
#' b057_wc2.0_bio_30s_02_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b057JaguarAlbersBuffer70/b057_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J57trk))
#crs(hf)   
#compareCRS(hf,get_crs(J57trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J57/UTM_b057_wc2.0_bio_30s_02_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J57/UTM_b057_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J58
#' b058_wc2.0_bio_30s_02_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b058JaguarAlbersBuffer70/b058_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J58trk))
#crs(hf)   
#compareCRS(hf,get_crs(J58trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J58/UTM_b058_wc2.0_bio_30s_02_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J58/UTM_b058_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J59
#' b059_wc2.0_bio_30s_02_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b059JaguarAlbersBuffer70/b059_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J59trk))
#crs(hf)   
#compareCRS(hf,get_crs(J59trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J59/UTM_b059_wc2.0_bio_30s_02_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J59/UTM_b059_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J60
#' b060_wc2.0_bio_30s_02_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b060JaguarAlbersBuffer70/b060_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J60trk))
#crs(hf)   
#compareCRS(hf,get_crs(J60trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J60/UTM_b060_wc2.0_bio_30s_02_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J60/UTM_b060_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J61
#' b061_wc2.0_bio_30s_02_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b061JaguarAlbersBuffer70/b061_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J61trk))
#crs(hf)   
#compareCRS(hf,get_crs(J61trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J61/UTM_b061_wc2.0_bio_30s_02_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J61/UTM_b061_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J62
#' b062_wc2.0_bio_30s_02_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b062JaguarAlbersBuffer70/b062_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J62trk))
#crs(hf)   
#compareCRS(hf,get_crs(J62trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J62/UTM_b062_wc2.0_bio_30s_02_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J62/UTM_b062_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J63
#' b063_wc2.0_bio_30s_02_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b063JaguarAlbersBuffer70/b063_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J63trk))
#crs(hf)   
#compareCRS(hf,get_crs(J63trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J63/UTM_b063_wc2.0_bio_30s_02_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J63/UTM_b063_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

gc()
gc()
gc()


#' J64
#' b064_wc2.0_bio_30s_02_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b064JaguarAlbersBuffer70/b064_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J64trk))
#crs(hf)   
#compareCRS(hf,get_crs(J64trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J64/UTM_b064_wc2.0_bio_30s_02_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J64/UTM_b064_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J65
#' b065_wc2.0_bio_30s_02_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b065JaguarAlbersBuffer70/b065_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J65trk))
#crs(hf)   
#compareCRS(hf,get_crs(J65trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J65/UTM_b065_wc2.0_bio_30s_02_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J65/UTM_b065_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J66
#' b066_wc2.0_bio_30s_02_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b066JaguarAlbersBuffer70/b066_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J66trk))
#crs(hf)   
#compareCRS(hf,get_crs(J66trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J66/UTM_b066_wc2.0_bio_30s_02_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J66/UTM_b066_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J67
#' b067_wc2.0_bio_30s_02_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b067JaguarAlbersBuffer70/b067_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J67trk))
#crs(hf)   
#compareCRS(hf,get_crs(J67trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J67/UTM_b067_wc2.0_bio_30s_02_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J67/UTM_b067_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J68
#' b068_wc2.0_bio_30s_02_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b068JaguarAlbersBuffer70/b068_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J68trk))
#crs(hf)   
#compareCRS(hf,get_crs(J68trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J68/UTM_b068_wc2.0_bio_30s_02_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J68/UTM_b068_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J69
#' b069_wc2.0_bio_30s_02_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b069JaguarAlbersBuffer70/b069_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J69trk))
#crs(hf)   
#compareCRS(hf,get_crs(J69trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J69/UTM_b069_wc2.0_bio_30s_02_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J69/UTM_b069_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

gc()
gc()
gc()


#' J70
#' b070_wc2.0_bio_30s_02_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b070JaguarAlbersBuffer70/b070_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J70trk))
#crs(hf)   
#compareCRS(hf,get_crs(J70trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J70/UTM_b070_wc2.0_bio_30s_02_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J70/UTM_b070_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J71
#' b071_wc2.0_bio_30s_02_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b071JaguarAlbersBuffer70/b071_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J71trk))
#crs(hf)   
#compareCRS(hf,get_crs(J71trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J71/UTM_b071_wc2.0_bio_30s_02_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J71/UTM_b071_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J72
#' b072_wc2.0_bio_30s_02_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b072JaguarAlbersBuffer70/b072_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J72trk))
#crs(hf)   
#compareCRS(hf,get_crs(J72trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J72/UTM_b072_wc2.0_bio_30s_02_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J72/UTM_b072_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J73
#' b073_wc2.0_bio_30s_02_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b073JaguarAlbersBuffer70/b073_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J73trk))
#crs(hf)   
#compareCRS(hf,get_crs(J73trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J73/UTM_b073_wc2.0_bio_30s_02_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J73/UTM_b073_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J74
#' b074_wc2.0_bio_30s_02_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b074JaguarAlbersBuffer70/b074_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J74trk))
#crs(hf)   
#compareCRS(hf,get_crs(J74trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J74/UTM_b074_wc2.0_bio_30s_02_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J74/UTM_b074_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J75
#' b075_wc2.0_bio_30s_02_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b075JaguarAlbersBuffer70/b075_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J75trk))
#crs(hf)   
#compareCRS(hf,get_crs(J75trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J75/UTM_b075_wc2.0_bio_30s_02_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J75/UTM_b075_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J76
#' b076_wc2.0_bio_30s_02_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b076JaguarAlbersBuffer70/b076_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J76trk))
#crs(hf)   
#compareCRS(hf,get_crs(J76trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J76/UTM_b076_wc2.0_bio_30s_02_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J76/UTM_b076_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J77
#' b077_wc2.0_bio_30s_02_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b077JaguarAlbersBuffer70/b077_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J77trk))
#crs(hf)   
#compareCRS(hf,get_crs(J77trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J77/UTM_b077_wc2.0_bio_30s_02_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J77/UTM_b077_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J78
#' b078_wc2.0_bio_30s_02_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b078JaguarAlbersBuffer70/b078_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J78trk))
#crs(hf)   
#compareCRS(hf,get_crs(J78trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J78/UTM_b078_wc2.0_bio_30s_02_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J78/UTM_b078_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J79
#' b079_wc2.0_bio_30s_02_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b079JaguarAlbersBuffer70/b079_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J79trk))
#crs(hf)   
#compareCRS(hf,get_crs(J79trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J79/UTM_b079_wc2.0_bio_30s_02_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J79/UTM_b079_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J80
#' b080_wc2.0_bio_30s_02_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b080JaguarAlbersBuffer70/b080_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J80trk))
#crs(hf)   
#compareCRS(hf,get_crs(J80trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J80/UTM_b080_wc2.0_bio_30s_02_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J80/UTM_b080_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J81
#' b081_wc2.0_bio_30s_02_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b081JaguarAlbersBuffer70/b081_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J81trk))
#crs(hf)   
#compareCRS(hf,get_crs(J81trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J81/UTM_b081_wc2.0_bio_30s_02_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J81/UTM_b081_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J82
#' b082_wc2.0_bio_30s_02_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b082JaguarAlbersBuffer70/b082_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J82trk))
#crs(hf)   
#compareCRS(hf,get_crs(J82trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J82/UTM_b082_wc2.0_bio_30s_02_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J82/UTM_b082_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J83
#' b083_wc2.0_bio_30s_02_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b083JaguarAlbersBuffer70/b083_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J83trk))
#crs(hf)   
#compareCRS(hf,get_crs(J83trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J83/UTM_b083_wc2.0_bio_30s_02_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J83/UTM_b083_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J84
#' b084_wc2.0_bio_30s_02_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b084JaguarAlbersBuffer70/b084_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J84trk))
#crs(hf)   
#compareCRS(hf,get_crs(J84trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J84/UTM_b084_wc2.0_bio_30s_02_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J84/UTM_b084_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J85
#' b085_wc2.0_bio_30s_02_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b085JaguarAlbersBuffer70/b085_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J85trk))
#crs(hf)   
#compareCRS(hf,get_crs(J85trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J85/UTM_b085_wc2.0_bio_30s_02_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J85/UTM_b085_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J86
#' b086_wc2.0_bio_30s_02_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b086JaguarAlbersBuffer70/b086_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J86trk))
#crs(hf)   
#compareCRS(hf,get_crs(J86trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J86/UTM_b086_wc2.0_bio_30s_02_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J86/UTM_b086_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J87
#' b087_wc2.0_bio_30s_02_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b087JaguarAlbersBuffer70/b087_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J87trk))
#crs(hf)   
#compareCRS(hf,get_crs(J87trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J87/UTM_b087_wc2.0_bio_30s_02_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J87/UTM_b087_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J88
#' b088_wc2.0_bio_30s_02_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b088JaguarAlbersBuffer70/b088_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J88trk))
#crs(hf)   
#compareCRS(hf,get_crs(J88trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J88/UTM_b088_wc2.0_bio_30s_02_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J88/UTM_b088_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J89
#' b089_wc2.0_bio_30s_02_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b089JaguarAlbersBuffer70/b089_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J89trk))
#crs(hf)   
#compareCRS(hf,get_crs(J89trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J89/UTM_b089_wc2.0_bio_30s_02_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J89/UTM_b089_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J90
#' b090_wc2.0_bio_30s_02_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b090JaguarAlbersBuffer70/b090_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J90trk))
#crs(hf)   
#compareCRS(hf,get_crs(J90trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J90/UTM_b090_wc2.0_bio_30s_02_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J90/UTM_b090_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J91
#' b091_wc2.0_bio_30s_02_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b091JaguarAlbersBuffer70/b091_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J91trk))
#crs(hf)   
#compareCRS(hf,get_crs(J91trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J91/UTM_b091_wc2.0_bio_30s_02_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J91/UTM_b091_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J92
#' b092_wc2.0_bio_30s_02_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b092JaguarAlbersBuffer70/b092_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J92trk))
#crs(hf)   
#compareCRS(hf,get_crs(J92trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J92/UTM_b092_wc2.0_bio_30s_02_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J92/UTM_b092_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J93
#' b093_wc2.0_bio_30s_02_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b093JaguarAlbersBuffer70/b093_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J93trk))
#crs(hf)   
#compareCRS(hf,get_crs(J93trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J93/UTM_b093_wc2.0_bio_30s_02_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J93/UTM_b093_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J94
#' b094_wc2.0_bio_30s_02_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b094JaguarAlbersBuffer70/b094_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J94trk))
#crs(hf)   
#compareCRS(hf,get_crs(J94trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J94/UTM_b094_wc2.0_bio_30s_02_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J94/UTM_b094_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J95
#' b095_wc2.0_bio_30s_02_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b095JaguarAlbersBuffer70/b095_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J95trk))
#crs(hf)   
#compareCRS(hf,get_crs(J95trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J95/UTM_b095_wc2.0_bio_30s_02_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J95/UTM_b095_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J96
#' b096_wc2.0_bio_30s_02_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b096JaguarAlbersBuffer70/b096_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J96trk))
#crs(hf)   
#compareCRS(hf,get_crs(J96trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J96/UTM_b096_wc2.0_bio_30s_02_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J96/UTM_b096_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J97
#' b097_wc2.0_bio_30s_02_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b097JaguarAlbersBuffer70/b097_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J97trk))
#crs(hf)   
#compareCRS(hf,get_crs(J97trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J97/UTM_b097_wc2.0_bio_30s_02_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J97/UTM_b097_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J98
#' b098_wc2.0_bio_30s_02_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b098JaguarAlbersBuffer70/b098_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J98trk))
#crs(hf)   
#compareCRS(hf,get_crs(J98trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J98/UTM_b098_wc2.0_bio_30s_02_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J98/UTM_b098_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J99
#' b099_wc2.0_bio_30s_02_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b099JaguarAlbersBuffer70/b099_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J99trk))
#crs(hf)   
#compareCRS(hf,get_crs(J99trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J99/UTM_b099_wc2.0_bio_30s_02_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J99/UTM_b099_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J100
#' b100_wc2.0_bio_30s_02_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b100JaguarAlbersBuffer70/b100_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J100trk))
#crs(hf)   
#compareCRS(hf,get_crs(J100trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J100/UTM_b100_wc2.0_bio_30s_02_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J100/UTM_b010_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J101
#' b101_wc2.0_bio_30s_02_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b101JaguarAlbersBuffer70/b101_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J101trk))
#crs(hf)   
#compareCRS(hf,get_crs(J101trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J101/UTM_b101_wc2.0_bio_30s_02_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J101/UTM_b101_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J102
#' b102_wc2.0_bio_30s_02_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b102JaguarAlbersBuffer70/b102_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J102trk))
#crs(hf)   
#compareCRS(hf,get_crs(J102trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J102/UTM_b102_wc2.0_bio_30s_02_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J102/UTM_b102_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J103
#' b103_wc2.0_bio_30s_02_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b103JaguarAlbersBuffer70/b103_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J103trk))
#crs(hf)   
#compareCRS(hf,get_crs(J103trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J103/UTM_b103_wc2.0_bio_30s_02_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J103/UTM_b103_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J104
#' b104_wc2.0_bio_30s_02_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b104JaguarAlbersBuffer70/b104_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J104trk))
#crs(hf)   
#compareCRS(hf,get_crs(J104trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J104/UTM_b104_wc2.0_bio_30s_02_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J104/UTM_b104_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J105
#' b105_wc2.0_bio_30s_02_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b105JaguarAlbersBuffer70/b105_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J105trk))
#crs(hf)   
#compareCRS(hf,get_crs(J105trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J105/UTM_b105_wc2.0_bio_30s_02_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J105/UTM_b105_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J106
#' b106_wc2.0_bio_30s_02_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b106JaguarAlbersBuffer70/b106_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J106trk))
#crs(hf)   
#compareCRS(hf,get_crs(J106trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J106/UTM_b106_wc2.0_bio_30s_02_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J106/UTM_b106_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J107
#' b107_wc2.0_bio_30s_02_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b107JaguarAlbersBuffer70/b107_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J107trk))
#crs(hf)   
#compareCRS(hf,get_crs(J107trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J107/UTM_b107_wc2.0_bio_30s_02_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J107/UTM_b107_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J108
#' b108_wc2.0_bio_30s_02_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b108JaguarAlbersBuffer70/b108_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J108trk))
#crs(hf)   
#compareCRS(hf,get_crs(J108trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J108/UTM_b108_wc2.0_bio_30s_02_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J108/UTM_b108_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J109
#' b109_wc2.0_bio_30s_02_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b109JaguarAlbersBuffer70/b109_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J109trk))
#crs(hf)   
#compareCRS(hf,get_crs(J109trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J109/UTM_b109_wc2.0_bio_30s_02_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J109/UTM_b109_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J110
#' b110_wc2.0_bio_30s_02_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b110JaguarAlbersBuffer70/b110_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J110trk))
#crs(hf)   
#compareCRS(hf,get_crs(J110trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J110/UTM_b110_wc2.0_bio_30s_02_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J110/UTM_b110_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J111
#' b111_wc2.0_bio_30s_02_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b111JaguarAlbersBuffer70/b111_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J111trk))
#crs(hf)   
#compareCRS(hf,get_crs(J111trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J111/UTM_b111_wc2.0_bio_30s_02_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J111/UTM_b111_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J112
#' b112_wc2.0_bio_30s_02_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b112JaguarAlbersBuffer70/b112_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J112trk))
#crs(hf)   
#compareCRS(hf,get_crs(J112trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J112/UTM_b112_wc2.0_bio_30s_02_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J112/UTM_b112_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J113
#' b113_wc2.0_bio_30s_02_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b113JaguarAlbersBuffer70/b113_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J113trk))
#crs(hf)   
#compareCRS(hf,get_crs(J113trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J113/UTM_b113_wc2.0_bio_30s_02_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J113/UTM_b113_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J114
#' b114_wc2.0_bio_30s_02_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b114JaguarAlbersBuffer70/b114_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J114trk))
#crs(hf)   
#compareCRS(hf,get_crs(J114trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J114/UTM_b114_wc2.0_bio_30s_02_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J114/UTM_b114_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J115
#' b115_wc2.0_bio_30s_02_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b115JaguarAlbersBuffer70/b115_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J115trk))
#crs(hf)   
#compareCRS(hf,get_crs(J115trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J115/UTM_b115_wc2.0_bio_30s_02_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J115/UTM_b115_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J116
#' b116_wc2.0_bio_30s_02_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b116JaguarAlbersBuffer70/b116_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J116trk))
#crs(hf)   
#compareCRS(hf,get_crs(J116trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J116/UTM_b116_wc2.0_bio_30s_02_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J116/UTM_b116_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J117
#' b117_wc2.0_bio_30s_02_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b117JaguarAlbersBuffer70/b117_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J117trk))
#crs(hf)   
#compareCRS(hf,get_crs(J117trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J117/UTM_b117_wc2.0_bio_30s_02_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J117/UTM_b117_wc2.0_bio_30s_02_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#####################################################################################################################################
# wc2.0_bio_30s_04_neotropic_albers_tif_exp.
# BIO4 = Temperature Seasonality (standard deviation *100)
gc()
gc()
gc()

#' J1
#' b001_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif
hf <-"D:/GISbases/b001JaguarAlbersBuffer70/b001_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J1trk))
#crs(hf)   
#compareCRS(hf,get_crs(J1trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J1/UTM_b001_wc2.0_bio_30s_04_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J1/UTM_b001_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J2
#' b002_wc2.0_bio_30s_04_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b002JaguarAlbersBuffer70/b002_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J2trk))
#crs(hf)   
#compareCRS(hf,get_crs(J2trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J2/UTM_b002_wc2.0_bio_30s_04_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J2/UTM_b002_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J3
#' b003_wc2.0_bio_30s_04_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b003JaguarAlbersBuffer70/b003_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J3trk))
#crs(hf)   
#compareCRS(hf,get_crs(J3trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J3/UTM_b003_wc2.0_bio_30s_04_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J3/UTM_b003_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J4
#' b004_wc2.0_bio_30s_04_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b004JaguarAlbersBuffer70/b004_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J4trk))
#crs(hf)   
#compareCRS(hf,get_crs(J4trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J4/UTM_b004_wc2.0_bio_30s_04_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J4/UTM_b004_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J5
#' b005_wc2.0_bio_30s_04_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b005JaguarAlbersBuffer70/b005_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J5trk))
#crs(hf)   
#compareCRS(hf,get_crs(J5trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J5/UTM_b005_wc2.0_bio_30s_04_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J5/UTM_b005_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J6
#' b006_wc2.0_bio_30s_04_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b006JaguarAlbersBuffer70/b006_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J6trk))
#crs(hf)   
#compareCRS(hf,get_crs(J6trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J6/UTM_b006_wc2.0_bio_30s_04_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J6/UTM_b006_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J7
#' b007_wc2.0_bio_30s_04_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b007JaguarAlbersBuffer70/b007_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J7trk))
#crs(hf)   
#compareCRS(hf,get_crs(J7trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J7/UTM_b007_wc2.0_bio_30s_04_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J7/UTM_b007_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J8
#' b008_wc2.0_bio_30s_04_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b008JaguarAlbersBuffer70/b008_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J8trk))
#crs(hf)   
#compareCRS(hf,get_crs(J8trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J8/UTM_b008_wc2.0_bio_30s_04_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J8/UTM_b008_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J9
#' b009_wc2.0_bio_30s_04_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b009JaguarAlbersBuffer70/b009_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J9trk))
#crs(hf)   
#compareCRS(hf,get_crs(J9trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J9/UTM_b009_wc2.0_bio_30s_04_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J9/UTM_b009_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J10
#' b010_wc2.0_bio_30s_04_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b010JaguarAlbersBuffer70/b010_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J10trk))
#crs(hf)   
#compareCRS(hf,get_crs(J10trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J10/UTM_b010_wc2.0_bio_30s_04_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J10/UTM_b010_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

gc()
gc()
gc()


#' J11
#' b011_wc2.0_bio_30s_04_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b011JaguarAlbersBuffer70/b011_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J11trk))
#crs(hf)   
#compareCRS(hf,get_crs(J11trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J11/UTM_b011_wc2.0_bio_30s_04_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J11/UTM_b011_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J12
#' b012_wc2.0_bio_30s_04_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b012JaguarAlbersBuffer70/b012_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J12trk))
#crs(hf)   
#compareCRS(hf,get_crs(J12trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J12/UTM_b012_wc2.0_bio_30s_04_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J12/UTM_b012_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J13
#' b013_wc2.0_bio_30s_04_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b013JaguarAlbersBuffer70/b013_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J13trk))
#crs(hf)   
#compareCRS(hf,get_crs(J13trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J13/UTM_b013_wc2.0_bio_30s_04_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J13/UTM_b013_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J14
#' b014_wc2.0_bio_30s_04_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b014JaguarAlbersBuffer70/b014_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J14trk))
#crs(hf)   
#compareCRS(hf,get_crs(J14trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J14/UTM_b014_wc2.0_bio_30s_04_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J14/UTM_b014_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J15
#' b015_wc2.0_bio_30s_04_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b015JaguarAlbersBuffer70/b015_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J15trk))
#crs(hf)   
#compareCRS(hf,get_crs(J15trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J15/UTM_b015_wc2.0_bio_30s_04_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J15/UTM_b015_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J16
#' b016_wc2.0_bio_30s_04_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b016JaguarAlbersBuffer70/b016_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J16trk))
#crs(hf)   
#compareCRS(hf,get_crs(J16trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J16/UTM_b016_wc2.0_bio_30s_04_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J16/UTM_b016_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J17
#' b017_wc2.0_bio_30s_04_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b017JaguarAlbersBuffer70/b017_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J17trk))
#crs(hf)   
#compareCRS(hf,get_crs(J17trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J17/UTM_b017_wc2.0_bio_30s_04_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J17/UTM_b017_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J18
#' b018_wc2.0_bio_30s_04_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b018JaguarAlbersBuffer70/b018_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J18trk))
#crs(hf)   
#compareCRS(hf,get_crs(J18trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J18/UTM_b018_wc2.0_bio_30s_04_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J18/UTM_b018_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J19
#' b019_wc2.0_bio_30s_04_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b019JaguarAlbersBuffer70/b019_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J19trk))
#crs(hf)   
#compareCRS(hf,get_crs(J19trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J19/UTM_b019_wc2.0_bio_30s_04_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J19/UTM_b019_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J20
#' b020_wc2.0_bio_30s_04_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b020JaguarAlbersBuffer70/b020_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J20trk))
#crs(hf)   
#compareCRS(hf,get_crs(J20trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J20/UTM_b020_wc2.0_bio_30s_04_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J20/UTM_b020_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J21
#' b021_wc2.0_bio_30s_04_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b021JaguarAlbersBuffer70/b021_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J21trk))
#crs(hf)   
#compareCRS(hf,get_crs(J21trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J21/UTM_b021_wc2.0_bio_30s_04_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J21/UTM_b021_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J22
#' b022_wc2.0_bio_30s_04_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b022JaguarAlbersBuffer70/b022_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J22trk))
#crs(hf)   
#compareCRS(hf,get_crs(J22trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J22/UTM_b022_wc2.0_bio_30s_04_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J22/UTM_b022_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J23
#' b023_wc2.0_bio_30s_04_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b023JaguarAlbersBuffer70/b023_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J23trk))
#crs(hf)   
#compareCRS(hf,get_crs(J23trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J23/UTM_b023_wc2.0_bio_30s_04_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J23/UTM_b023_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J24
#' b024_wc2.0_bio_30s_04_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b024JaguarAlbersBuffer70/b024_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J24trk))
#crs(hf)   
#compareCRS(hf,get_crs(J24trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J24/UTM_b024_wc2.0_bio_30s_04_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J24/UTM_b024_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J25
#' b025_wc2.0_bio_30s_04_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b025JaguarAlbersBuffer70/b025_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J25trk))
#crs(hf)   
#compareCRS(hf,get_crs(J25trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J25/UTM_b025_wc2.0_bio_30s_04_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J25/UTM_b025_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J26
#' b026_wc2.0_bio_30s_04_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b026JaguarAlbersBuffer70/b026_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J26trk))
#crs(hf)   
#compareCRS(hf,get_crs(J26trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J26/UTM_b026_wc2.0_bio_30s_04_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J26/UTM_b026_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J27
#' b027_wc2.0_bio_30s_04_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b027JaguarAlbersBuffer70/b027_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J27trk))
#crs(hf)   
#compareCRS(hf,get_crs(J27trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J27/UTM_b027_wc2.0_bio_30s_04_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J27/UTM_b027_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J28
#' b028_wc2.0_bio_30s_04_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b028JaguarAlbersBuffer70/b028_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J28trk))
#crs(hf)   
#compareCRS(hf,get_crs(J28trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J28/UTM_b028_wc2.0_bio_30s_04_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J28/UTM_b028_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J29
#' b029_wc2.0_bio_30s_04_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b029JaguarAlbersBuffer70/b029_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J29trk))
#crs(hf)   
#compareCRS(hf,get_crs(J29trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J29/UTM_b029_wc2.0_bio_30s_04_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J29/UTM_b029_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J30
#' b030_wc2.0_bio_30s_04_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b030JaguarAlbersBuffer70/b030_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J30trk))
#crs(hf)   
#compareCRS(hf,get_crs(J30trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J30/UTM_b030_wc2.0_bio_30s_04_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J30/UTM_b030_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J30
#' b030_wc2.0_bio_30s_04_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b030JaguarAlbersBuffer70/b030_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J30trk))
#crs(hf)   
#compareCRS(hf,get_crs(J30trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J30/UTM_b030_wc2.0_bio_30s_04_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J30/UTM_b030_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J31
#' b031_wc2.0_bio_30s_04_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b031JaguarAlbersBuffer70/b031_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J31trk))
#crs(hf)   
#compareCRS(hf,get_crs(J31trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J31/UTM_b031_wc2.0_bio_30s_04_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J31/UTM_b031_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J32
#' b032_wc2.0_bio_30s_04_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b032JaguarAlbersBuffer70/b032_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J32trk))
#crs(hf)   
#compareCRS(hf,get_crs(J32trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J32/UTM_b032_wc2.0_bio_30s_04_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J32/UTM_b032_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J33
#' b033_wc2.0_bio_30s_04_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b033JaguarAlbersBuffer70/b033_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J33trk))
#crs(hf)   
#compareCRS(hf,get_crs(J33trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J33/UTM_b033_wc2.0_bio_30s_04_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J33/UTM_b033_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J34
#' b034_wc2.0_bio_30s_04_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b034JaguarAlbersBuffer70/b034_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J34trk))
#crs(hf)   
#compareCRS(hf,get_crs(J34trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J34/UTM_b034_wc2.0_bio_30s_04_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J34/UTM_b034_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J35
#' b035_wc2.0_bio_30s_04_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b035JaguarAlbersBuffer70/b035_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J35trk))
#crs(hf)   
#compareCRS(hf,get_crs(J35trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J35/UTM_b035_wc2.0_bio_30s_04_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J35/UTM_b035_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J36
#' b036_wc2.0_bio_30s_04_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b036JaguarAlbersBuffer70/b036_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J36trk))
#crs(hf)   
#compareCRS(hf,get_crs(J36trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J36/UTM_b036_wc2.0_bio_30s_04_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J36/UTM_b036_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J37
#' b037_wc2.0_bio_30s_04_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b037JaguarAlbersBuffer70/b037_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J37trk))
#crs(hf)   
#compareCRS(hf,get_crs(J37trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J37/UTM_b037_wc2.0_bio_30s_04_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J37/UTM_b037_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J38
#' b038_wc2.0_bio_30s_04_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b038JaguarAlbersBuffer70/b038_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J38trk))
#crs(hf)   
#compareCRS(hf,get_crs(J38trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J38/UTM_b038_wc2.0_bio_30s_04_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J38/UTM_b038_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J39
#' b039_wc2.0_bio_30s_04_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b039JaguarAlbersBuffer70/b039_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J39trk))
#crs(hf)   
#compareCRS(hf,get_crs(J39trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J39/UTM_b039_wc2.0_bio_30s_04_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J39/UTM_b039_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J40
#' b040_wc2.0_bio_30s_04_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b040JaguarAlbersBuffer70/b040_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J40trk))
#crs(hf)   
#compareCRS(hf,get_crs(J40trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J40/UTM_b040_wc2.0_bio_30s_04_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J40/UTM_b040_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J41
#' b041_wc2.0_bio_30s_04_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b041JaguarAlbersBuffer70/b041_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J41trk))
#crs(hf)   
#compareCRS(hf,get_crs(J41trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J41/UTM_b041_wc2.0_bio_30s_04_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J41/UTM_b041_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J42
#' b042_wc2.0_bio_30s_04_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b042JaguarAlbersBuffer70/b042_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J42trk))
#crs(hf)   
#compareCRS(hf,get_crs(J42trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J42/UTM_b042_wc2.0_bio_30s_04_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J42/UTM_b042_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

gc()
gc()
gc()

#' J43
#' b043_wc2.0_bio_30s_04_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b043JaguarAlbersBuffer70/b043_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J43trk))
#crs(hf)   
#compareCRS(hf,get_crs(J43trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J43/UTM_b043_wc2.0_bio_30s_04_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J43/UTM_b043_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J44
#' b044_wc2.0_bio_30s_04_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b044JaguarAlbersBuffer70/b044_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J44trk))
#crs(hf)   
#compareCRS(hf,get_crs(J44trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J44/UTM_b044_wc2.0_bio_30s_04_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J44/UTM_b044_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J45
#' b045_wc2.0_bio_30s_04_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b045JaguarAlbersBuffer70/b045_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J45trk))
#crs(hf)   
#compareCRS(hf,get_crs(J45trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J45/UTM_b045_wc2.0_bio_30s_04_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J45/UTM_b045_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J46
#' b046_wc2.0_bio_30s_04_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b046JaguarAlbersBuffer70/b046_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J46trk))
#crs(hf)   
#compareCRS(hf,get_crs(J46trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J46/UTM_b046_wc2.0_bio_30s_04_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J46/UTM_b046_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J47
#' b047_wc2.0_bio_30s_04_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b047JaguarAlbersBuffer70/b047_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J47trk))
#crs(hf)   
#compareCRS(hf,get_crs(J47trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J47/UTM_b047_wc2.0_bio_30s_04_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J47/UTM_b047_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J48
#' b048_wc2.0_bio_30s_04_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b048JaguarAlbersBuffer70/b048_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J48trk))
#crs(hf)   
#compareCRS(hf,get_crs(J48trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J48/UTM_b048_wc2.0_bio_30s_04_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J48/UTM_b048_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J49
#' b049_wc2.0_bio_30s_04_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b049JaguarAlbersBuffer70/b049_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J49trk))
#crs(hf)   
#compareCRS(hf,get_crs(J49trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J49/UTM_b049_wc2.0_bio_30s_04_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J49/UTM_b049_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J50
#' b050_wc2.0_bio_30s_04_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b050JaguarAlbersBuffer70/b050_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J50trk))
#crs(hf)   
#compareCRS(hf,get_crs(J50trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J50/UTM_b050_wc2.0_bio_30s_04_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J50/UTM_b050_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J51
#' b051_wc2.0_bio_30s_04_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b051JaguarAlbersBuffer70/b051_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J51trk))
#crs(hf)   
#compareCRS(hf,get_crs(J51trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J51/UTM_b051_wc2.0_bio_30s_04_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J51/UTM_b051_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J52
#' b052_wc2.0_bio_30s_04_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b052JaguarAlbersBuffer70/b052_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J52trk))
#crs(hf)   
#compareCRS(hf,get_crs(J52trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J52/UTM_b052_wc2.0_bio_30s_04_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J52/UTM_b052_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J53
#' b053_wc2.0_bio_30s_04_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b053JaguarAlbersBuffer70/b053_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J53trk))
#crs(hf)   
#compareCRS(hf,get_crs(J53trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J53/UTM_b053_wc2.0_bio_30s_04_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J53/UTM_b053_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J54
#' b054_wc2.0_bio_30s_04_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b054JaguarAlbersBuffer70/b054_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J54trk))
#crs(hf)   
#compareCRS(hf,get_crs(J54trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J54/UTM_b054_wc2.0_bio_30s_04_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J54/UTM_b054_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J55
#' b055_wc2.0_bio_30s_04_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b055JaguarAlbersBuffer70/b055_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J55trk))
#crs(hf)   
#compareCRS(hf,get_crs(J55trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J55/UTM_b055_wc2.0_bio_30s_04_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J55/UTM_b055_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J56
#' b056_wc2.0_bio_30s_04_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b056JaguarAlbersBuffer70/b056_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J56trk))
#crs(hf)   
#compareCRS(hf,get_crs(J56trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J56/UTM_b056_wc2.0_bio_30s_04_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J56/UTM_b056_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J57
#' b057_wc2.0_bio_30s_04_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b057JaguarAlbersBuffer70/b057_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J57trk))
#crs(hf)   
#compareCRS(hf,get_crs(J57trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J57/UTM_b057_wc2.0_bio_30s_04_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J57/UTM_b057_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J58
#' b058_wc2.0_bio_30s_04_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b058JaguarAlbersBuffer70/b058_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J58trk))
#crs(hf)   
#compareCRS(hf,get_crs(J58trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J58/UTM_b058_wc2.0_bio_30s_04_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J58/UTM_b058_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J59
#' b059_wc2.0_bio_30s_04_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b059JaguarAlbersBuffer70/b059_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J59trk))
#crs(hf)   
#compareCRS(hf,get_crs(J59trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J59/UTM_b059_wc2.0_bio_30s_04_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J59/UTM_b059_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J60
#' b060_wc2.0_bio_30s_04_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b060JaguarAlbersBuffer70/b060_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J60trk))
#crs(hf)   
#compareCRS(hf,get_crs(J60trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J60/UTM_b060_wc2.0_bio_30s_04_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J60/UTM_b060_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J61
#' b061_wc2.0_bio_30s_04_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b061JaguarAlbersBuffer70/b061_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J61trk))
#crs(hf)   
#compareCRS(hf,get_crs(J61trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J61/UTM_b061_wc2.0_bio_30s_04_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J61/UTM_b061_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J62
#' b062_wc2.0_bio_30s_04_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b062JaguarAlbersBuffer70/b062_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J62trk))
#crs(hf)   
#compareCRS(hf,get_crs(J62trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J62/UTM_b062_wc2.0_bio_30s_04_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J62/UTM_b062_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J63
#' b063_wc2.0_bio_30s_04_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b063JaguarAlbersBuffer70/b063_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J63trk))
#crs(hf)   
#compareCRS(hf,get_crs(J63trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J63/UTM_b063_wc2.0_bio_30s_04_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J63/UTM_b063_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

gc()
gc()
gc()


#' J64
#' b064_wc2.0_bio_30s_04_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b064JaguarAlbersBuffer70/b064_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J64trk))
#crs(hf)   
#compareCRS(hf,get_crs(J64trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J64/UTM_b064_wc2.0_bio_30s_04_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J64/UTM_b064_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J65
#' b065_wc2.0_bio_30s_04_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b065JaguarAlbersBuffer70/b065_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J65trk))
#crs(hf)   
#compareCRS(hf,get_crs(J65trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J65/UTM_b065_wc2.0_bio_30s_04_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J65/UTM_b065_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J66
#' b066_wc2.0_bio_30s_04_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b066JaguarAlbersBuffer70/b066_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J66trk))
#crs(hf)   
#compareCRS(hf,get_crs(J66trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J66/UTM_b066_wc2.0_bio_30s_04_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J66/UTM_b066_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J67
#' b067_wc2.0_bio_30s_04_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b067JaguarAlbersBuffer70/b067_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J67trk))
#crs(hf)   
#compareCRS(hf,get_crs(J67trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J67/UTM_b067_wc2.0_bio_30s_04_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J67/UTM_b067_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J68
#' b068_wc2.0_bio_30s_04_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b068JaguarAlbersBuffer70/b068_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J68trk))
#crs(hf)   
#compareCRS(hf,get_crs(J68trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J68/UTM_b068_wc2.0_bio_30s_04_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J68/UTM_b068_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J69
#' b069_wc2.0_bio_30s_04_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b069JaguarAlbersBuffer70/b069_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J69trk))
#crs(hf)   
#compareCRS(hf,get_crs(J69trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J69/UTM_b069_wc2.0_bio_30s_04_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J69/UTM_b069_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

gc()
gc()
gc()


#' J70
#' b070_wc2.0_bio_30s_04_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b070JaguarAlbersBuffer70/b070_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J70trk))
#crs(hf)   
#compareCRS(hf,get_crs(J70trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J70/UTM_b070_wc2.0_bio_30s_04_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J70/UTM_b070_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J71
#' b071_wc2.0_bio_30s_04_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b071JaguarAlbersBuffer70/b071_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J71trk))
#crs(hf)   
#compareCRS(hf,get_crs(J71trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J71/UTM_b071_wc2.0_bio_30s_04_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J71/UTM_b071_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J72
#' b072_wc2.0_bio_30s_04_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b072JaguarAlbersBuffer70/b072_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J72trk))
#crs(hf)   
#compareCRS(hf,get_crs(J72trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J72/UTM_b072_wc2.0_bio_30s_04_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J72/UTM_b072_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J73
#' b073_wc2.0_bio_30s_04_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b073JaguarAlbersBuffer70/b073_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J73trk))
#crs(hf)   
#compareCRS(hf,get_crs(J73trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J73/UTM_b073_wc2.0_bio_30s_04_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J73/UTM_b073_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J74
#' b074_wc2.0_bio_30s_04_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b074JaguarAlbersBuffer70/b074_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J74trk))
#crs(hf)   
#compareCRS(hf,get_crs(J74trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J74/UTM_b074_wc2.0_bio_30s_04_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J74/UTM_b074_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J75
#' b075_wc2.0_bio_30s_04_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b075JaguarAlbersBuffer70/b075_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J75trk))
#crs(hf)   
#compareCRS(hf,get_crs(J75trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J75/UTM_b075_wc2.0_bio_30s_04_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J75/UTM_b075_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J76
#' b076_wc2.0_bio_30s_04_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b076JaguarAlbersBuffer70/b076_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J76trk))
#crs(hf)   
#compareCRS(hf,get_crs(J76trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J76/UTM_b076_wc2.0_bio_30s_04_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J76/UTM_b076_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J77
#' b077_wc2.0_bio_30s_04_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b077JaguarAlbersBuffer70/b077_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J77trk))
#crs(hf)   
#compareCRS(hf,get_crs(J77trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J77/UTM_b077_wc2.0_bio_30s_04_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J77/UTM_b077_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J78
#' b078_wc2.0_bio_30s_04_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b078JaguarAlbersBuffer70/b078_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J78trk))
#crs(hf)   
#compareCRS(hf,get_crs(J78trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J78/UTM_b078_wc2.0_bio_30s_04_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J78/UTM_b078_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J79
#' b079_wc2.0_bio_30s_04_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b079JaguarAlbersBuffer70/b079_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J79trk))
#crs(hf)   
#compareCRS(hf,get_crs(J79trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J79/UTM_b079_wc2.0_bio_30s_04_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J79/UTM_b079_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J80
#' b080_wc2.0_bio_30s_04_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b080JaguarAlbersBuffer70/b080_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J80trk))
#crs(hf)   
#compareCRS(hf,get_crs(J80trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J80/UTM_b080_wc2.0_bio_30s_04_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J80/UTM_b080_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J81
#' b081_wc2.0_bio_30s_04_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b081JaguarAlbersBuffer70/b081_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J81trk))
#crs(hf)   
#compareCRS(hf,get_crs(J81trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J81/UTM_b081_wc2.0_bio_30s_04_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J81/UTM_b081_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J82
#' b082_wc2.0_bio_30s_04_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b082JaguarAlbersBuffer70/b082_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J82trk))
#crs(hf)   
#compareCRS(hf,get_crs(J82trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J82/UTM_b082_wc2.0_bio_30s_04_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J82/UTM_b082_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J83
#' b083_wc2.0_bio_30s_04_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b083JaguarAlbersBuffer70/b083_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J83trk))
#crs(hf)   
#compareCRS(hf,get_crs(J83trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J83/UTM_b083_wc2.0_bio_30s_04_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J83/UTM_b083_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J84
#' b084_wc2.0_bio_30s_04_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b084JaguarAlbersBuffer70/b084_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J84trk))
#crs(hf)   
#compareCRS(hf,get_crs(J84trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J84/UTM_b084_wc2.0_bio_30s_04_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J84/UTM_b084_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J85
#' b085_wc2.0_bio_30s_04_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b085JaguarAlbersBuffer70/b085_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J85trk))
#crs(hf)   
#compareCRS(hf,get_crs(J85trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J85/UTM_b085_wc2.0_bio_30s_04_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J85/UTM_b085_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J86
#' b086_wc2.0_bio_30s_04_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b086JaguarAlbersBuffer70/b086_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J86trk))
#crs(hf)   
#compareCRS(hf,get_crs(J86trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J86/UTM_b086_wc2.0_bio_30s_04_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J86/UTM_b086_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J87
#' b087_wc2.0_bio_30s_04_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b087JaguarAlbersBuffer70/b087_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J87trk))
#crs(hf)   
#compareCRS(hf,get_crs(J87trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J87/UTM_b087_wc2.0_bio_30s_04_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J87/UTM_b087_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J88
#' b088_wc2.0_bio_30s_04_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b088JaguarAlbersBuffer70/b088_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J88trk))
#crs(hf)   
#compareCRS(hf,get_crs(J88trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J88/UTM_b088_wc2.0_bio_30s_04_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J88/UTM_b088_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J89
#' b089_wc2.0_bio_30s_04_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b089JaguarAlbersBuffer70/b089_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J89trk))
#crs(hf)   
#compareCRS(hf,get_crs(J89trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J89/UTM_b089_wc2.0_bio_30s_04_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J89/UTM_b089_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J90
#' b090_wc2.0_bio_30s_04_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b090JaguarAlbersBuffer70/b090_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J90trk))
#crs(hf)   
#compareCRS(hf,get_crs(J90trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J90/UTM_b090_wc2.0_bio_30s_04_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J90/UTM_b090_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J91
#' b091_wc2.0_bio_30s_04_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b091JaguarAlbersBuffer70/b091_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J91trk))
#crs(hf)   
#compareCRS(hf,get_crs(J91trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J91/UTM_b091_wc2.0_bio_30s_04_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J91/UTM_b091_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J92
#' b092_wc2.0_bio_30s_04_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b092JaguarAlbersBuffer70/b092_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J92trk))
#crs(hf)   
#compareCRS(hf,get_crs(J92trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J92/UTM_b092_wc2.0_bio_30s_04_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J92/UTM_b092_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J93
#' b093_wc2.0_bio_30s_04_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b093JaguarAlbersBuffer70/b093_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J93trk))
#crs(hf)   
#compareCRS(hf,get_crs(J93trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J93/UTM_b093_wc2.0_bio_30s_04_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J93/UTM_b093_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J94
#' b094_wc2.0_bio_30s_04_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b094JaguarAlbersBuffer70/b094_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J94trk))
#crs(hf)   
#compareCRS(hf,get_crs(J94trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J94/UTM_b094_wc2.0_bio_30s_04_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J94/UTM_b094_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J95
#' b095_wc2.0_bio_30s_04_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b095JaguarAlbersBuffer70/b095_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J95trk))
#crs(hf)   
#compareCRS(hf,get_crs(J95trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J95/UTM_b095_wc2.0_bio_30s_04_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J95/UTM_b095_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J96
#' b096_wc2.0_bio_30s_04_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b096JaguarAlbersBuffer70/b096_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J96trk))
#crs(hf)   
#compareCRS(hf,get_crs(J96trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J96/UTM_b096_wc2.0_bio_30s_04_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J96/UTM_b096_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J97
#' b097_wc2.0_bio_30s_04_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b097JaguarAlbersBuffer70/b097_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J97trk))
#crs(hf)   
#compareCRS(hf,get_crs(J97trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J97/UTM_b097_wc2.0_bio_30s_04_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J97/UTM_b097_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J98
#' b098_wc2.0_bio_30s_04_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b098JaguarAlbersBuffer70/b098_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J98trk))
#crs(hf)   
#compareCRS(hf,get_crs(J98trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J98/UTM_b098_wc2.0_bio_30s_04_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J98/UTM_b098_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J99
#' b099_wc2.0_bio_30s_04_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b099JaguarAlbersBuffer70/b099_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J99trk))
#crs(hf)   
#compareCRS(hf,get_crs(J99trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J99/UTM_b099_wc2.0_bio_30s_04_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J99/UTM_b099_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J100
#' b100_wc2.0_bio_30s_04_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b100JaguarAlbersBuffer70/b100_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J100trk))
#crs(hf)   
#compareCRS(hf,get_crs(J100trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J100/UTM_b100_wc2.0_bio_30s_04_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J100/UTM_b010_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J101
#' b101_wc2.0_bio_30s_04_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b101JaguarAlbersBuffer70/b101_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J101trk))
#crs(hf)   
#compareCRS(hf,get_crs(J101trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J101/UTM_b101_wc2.0_bio_30s_04_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J101/UTM_b101_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J102
#' b102_wc2.0_bio_30s_04_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b102JaguarAlbersBuffer70/b102_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J102trk))
#crs(hf)   
#compareCRS(hf,get_crs(J102trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J102/UTM_b102_wc2.0_bio_30s_04_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J102/UTM_b102_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J103
#' b103_wc2.0_bio_30s_04_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b103JaguarAlbersBuffer70/b103_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J103trk))
#crs(hf)   
#compareCRS(hf,get_crs(J103trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J103/UTM_b103_wc2.0_bio_30s_04_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J103/UTM_b103_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J104
#' b104_wc2.0_bio_30s_04_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b104JaguarAlbersBuffer70/b104_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J104trk))
#crs(hf)   
#compareCRS(hf,get_crs(J104trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J104/UTM_b104_wc2.0_bio_30s_04_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J104/UTM_b104_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J105
#' b105_wc2.0_bio_30s_04_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b105JaguarAlbersBuffer70/b105_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J105trk))
#crs(hf)   
#compareCRS(hf,get_crs(J105trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J105/UTM_b105_wc2.0_bio_30s_04_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J105/UTM_b105_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J106
#' b106_wc2.0_bio_30s_04_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b106JaguarAlbersBuffer70/b106_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J106trk))
#crs(hf)   
#compareCRS(hf,get_crs(J106trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J106/UTM_b106_wc2.0_bio_30s_04_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J106/UTM_b106_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J107
#' b107_wc2.0_bio_30s_04_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b107JaguarAlbersBuffer70/b107_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J107trk))
#crs(hf)   
#compareCRS(hf,get_crs(J107trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J107/UTM_b107_wc2.0_bio_30s_04_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J107/UTM_b107_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J108
#' b108_wc2.0_bio_30s_04_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b108JaguarAlbersBuffer70/b108_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J108trk))
#crs(hf)   
#compareCRS(hf,get_crs(J108trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J108/UTM_b108_wc2.0_bio_30s_04_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J108/UTM_b108_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J109
#' b109_wc2.0_bio_30s_04_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b109JaguarAlbersBuffer70/b109_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J109trk))
#crs(hf)   
#compareCRS(hf,get_crs(J109trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J109/UTM_b109_wc2.0_bio_30s_04_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J109/UTM_b109_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J110
#' b110_wc2.0_bio_30s_04_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b110JaguarAlbersBuffer70/b110_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J110trk))
#crs(hf)   
#compareCRS(hf,get_crs(J110trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J110/UTM_b110_wc2.0_bio_30s_04_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J110/UTM_b110_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J111
#' b111_wc2.0_bio_30s_04_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b111JaguarAlbersBuffer70/b111_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J111trk))
#crs(hf)   
#compareCRS(hf,get_crs(J111trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J111/UTM_b111_wc2.0_bio_30s_04_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J111/UTM_b111_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J112
#' b112_wc2.0_bio_30s_04_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b112JaguarAlbersBuffer70/b112_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J112trk))
#crs(hf)   
#compareCRS(hf,get_crs(J112trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J112/UTM_b112_wc2.0_bio_30s_04_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J112/UTM_b112_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J113
#' b113_wc2.0_bio_30s_04_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b113JaguarAlbersBuffer70/b113_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J113trk))
#crs(hf)   
#compareCRS(hf,get_crs(J113trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J113/UTM_b113_wc2.0_bio_30s_04_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J113/UTM_b113_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J114
#' b114_wc2.0_bio_30s_04_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b114JaguarAlbersBuffer70/b114_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J114trk))
#crs(hf)   
#compareCRS(hf,get_crs(J114trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J114/UTM_b114_wc2.0_bio_30s_04_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J114/UTM_b114_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J115
#' b115_wc2.0_bio_30s_04_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b115JaguarAlbersBuffer70/b115_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J115trk))
#crs(hf)   
#compareCRS(hf,get_crs(J115trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J115/UTM_b115_wc2.0_bio_30s_04_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J115/UTM_b115_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J116
#' b116_wc2.0_bio_30s_04_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b116JaguarAlbersBuffer70/b116_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J116trk))
#crs(hf)   
#compareCRS(hf,get_crs(J116trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J116/UTM_b116_wc2.0_bio_30s_04_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J116/UTM_b116_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#' J117
#' b117_wc2.0_bio_30s_04_neotropic_albers_tif_exp.
hf <-"D:/GISbases/b117JaguarAlbersBuffer70/b117_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J117trk))
#crs(hf)   
#compareCRS(hf,get_crs(J117trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J117/UTM_b117_wc2.0_bio_30s_04_neotropic_albers_tif_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J117/UTM_b117_wc2.0_bio_30s_04_neotropic_albers_tif_exp.tif";(test=raster(test)); x11(); plot(test)

#####################################################################################################################################