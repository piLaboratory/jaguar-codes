#' #  **Transform GIS Bases**
#' Run JaguarDataPrep first !!! 
#' 
#' Load packages a few more packages                      
if(!require(install.load)) install.packages('install.load'); library(install.load)
install.load::install_load("raster","RCurl","rts","rvest","rgl","lubridate","lattice","rgdal","sp","stringr","methods",
                           "maptools", "vegan","sp","spatialEco","installr")

# Ecoregions2017_ecoregion_code_rast_exp.
# Ecoregions2017_ecoregion_code_rast_exp.

gc()
gc()
gc()

#' J1
#' b001_Ecoregions2017_ecoregion_code_rast_exp.tif
hf <-"D:/GISbases/b001JaguarAlbersBuffer70/b001_Ecoregions2017_ecoregion_code_rast_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J1trk))
#crs(hf)   
#compareCRS(hf,get_crs(J1trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J1/UTM_b001_Ecoregions2017_ecoregion_code_rast_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J1/UTM_b001_Ecoregions2017_ecoregion_code_rast_exp.tif";(test=raster(test)); x11(); plot(test)

#' J2
#' b002_Ecoregions2017_ecoregion_code_rast_exp.
hf <-"D:/GISbases/b002JaguarAlbersBuffer70/b002_Ecoregions2017_ecoregion_code_rast_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J2trk))
#crs(hf)   
#compareCRS(hf,get_crs(J2trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J2/UTM_b002_Ecoregions2017_ecoregion_code_rast_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J2/UTM_b002_Ecoregions2017_ecoregion_code_rast_exp.tif";(test=raster(test)); x11(); plot(test)

#' J3
#' b003_Ecoregions2017_ecoregion_code_rast_exp.
hf <-"D:/GISbases/b003JaguarAlbersBuffer70/b003_Ecoregions2017_ecoregion_code_rast_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J3trk))
#crs(hf)   
#compareCRS(hf,get_crs(J3trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J3/UTM_b003_Ecoregions2017_ecoregion_code_rast_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J3/UTM_b003_Ecoregions2017_ecoregion_code_rast_exp.tif";(test=raster(test)); x11(); plot(test)

#' J4
#' b004_Ecoregions2017_ecoregion_code_rast_exp.
hf <-"D:/GISbases/b004JaguarAlbersBuffer70/b004_Ecoregions2017_ecoregion_code_rast_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J4trk))
#crs(hf)   
#compareCRS(hf,get_crs(J4trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J4/UTM_b004_Ecoregions2017_ecoregion_code_rast_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J4/UTM_b004_Ecoregions2017_ecoregion_code_rast_exp.tif";(test=raster(test)); x11(); plot(test)

#' J5
#' b005_Ecoregions2017_ecoregion_code_rast_exp.
hf <-"D:/GISbases/b005JaguarAlbersBuffer70/b005_Ecoregions2017_ecoregion_code_rast_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J5trk))
#crs(hf)   
#compareCRS(hf,get_crs(J5trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J5/UTM_b005_Ecoregions2017_ecoregion_code_rast_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J5/UTM_b005_Ecoregions2017_ecoregion_code_rast_exp.tif";(test=raster(test)); x11(); plot(test)

#' J6
#' b006_Ecoregions2017_ecoregion_code_rast_exp.
hf <-"D:/GISbases/b006JaguarAlbersBuffer70/b006_Ecoregions2017_ecoregion_code_rast_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J6trk))
#crs(hf)   
#compareCRS(hf,get_crs(J6trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J6/UTM_b006_Ecoregions2017_ecoregion_code_rast_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J6/UTM_b006_Ecoregions2017_ecoregion_code_rast_exp.tif";(test=raster(test)); x11(); plot(test)

#' J7
#' b007_Ecoregions2017_ecoregion_code_rast_exp.
hf <-"D:/GISbases/b007JaguarAlbersBuffer70/b007_Ecoregions2017_ecoregion_code_rast_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J7trk))
#crs(hf)   
#compareCRS(hf,get_crs(J7trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J7/UTM_b007_Ecoregions2017_ecoregion_code_rast_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J7/UTM_b007_Ecoregions2017_ecoregion_code_rast_exp.tif";(test=raster(test)); x11(); plot(test)

#' J8
#' b008_Ecoregions2017_ecoregion_code_rast_exp.
hf <-"D:/GISbases/b008JaguarAlbersBuffer70/b008_Ecoregions2017_ecoregion_code_rast_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J8trk))
#crs(hf)   
#compareCRS(hf,get_crs(J8trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J8/UTM_b008_Ecoregions2017_ecoregion_code_rast_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J8/UTM_b008_Ecoregions2017_ecoregion_code_rast_exp.tif";(test=raster(test)); x11(); plot(test)

#' J9
#' b009_Ecoregions2017_ecoregion_code_rast_exp.
hf <-"D:/GISbases/b009JaguarAlbersBuffer70/b009_Ecoregions2017_ecoregion_code_rast_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J9trk))
#crs(hf)   
#compareCRS(hf,get_crs(J9trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J9/UTM_b009_Ecoregions2017_ecoregion_code_rast_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J9/UTM_b009_Ecoregions2017_ecoregion_code_rast_exp.tif";(test=raster(test)); x11(); plot(test)

#' J10
#' b010_Ecoregions2017_ecoregion_code_rast_exp.
hf <-"D:/GISbases/b010JaguarAlbersBuffer70/b010_Ecoregions2017_ecoregion_code_rast_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J10trk))
#crs(hf)   
#compareCRS(hf,get_crs(J10trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J10/UTM_b010_Ecoregions2017_ecoregion_code_rast_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J10/UTM_b010_Ecoregions2017_ecoregion_code_rast_exp.tif";(test=raster(test)); x11(); plot(test)

gc()
gc()
gc()


#' J11
#' b011_Ecoregions2017_ecoregion_code_rast_exp.
hf <-"D:/GISbases/b011JaguarAlbersBuffer70/b011_Ecoregions2017_ecoregion_code_rast_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J11trk))
#crs(hf)   
#compareCRS(hf,get_crs(J11trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J11/UTM_b011_Ecoregions2017_ecoregion_code_rast_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J11/UTM_b011_Ecoregions2017_ecoregion_code_rast_exp.tif";(test=raster(test)); x11(); plot(test)

#' J12
#' b012_Ecoregions2017_ecoregion_code_rast_exp.
hf <-"D:/GISbases/b012JaguarAlbersBuffer70/b012_Ecoregions2017_ecoregion_code_rast_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J12trk))
#crs(hf)   
#compareCRS(hf,get_crs(J12trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J12/UTM_b012_Ecoregions2017_ecoregion_code_rast_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J12/UTM_b012_Ecoregions2017_ecoregion_code_rast_exp.tif";(test=raster(test)); x11(); plot(test)

#' J13
#' b013_Ecoregions2017_ecoregion_code_rast_exp.
hf <-"D:/GISbases/b013JaguarAlbersBuffer70/b013_Ecoregions2017_ecoregion_code_rast_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J13trk))
#crs(hf)   
#compareCRS(hf,get_crs(J13trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J13/UTM_b013_Ecoregions2017_ecoregion_code_rast_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J13/UTM_b013_Ecoregions2017_ecoregion_code_rast_exp.tif";(test=raster(test)); x11(); plot(test)

#' J14
#' b014_Ecoregions2017_ecoregion_code_rast_exp.
hf <-"D:/GISbases/b014JaguarAlbersBuffer70/b014_Ecoregions2017_ecoregion_code_rast_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J14trk))
#crs(hf)   
#compareCRS(hf,get_crs(J14trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J14/UTM_b014_Ecoregions2017_ecoregion_code_rast_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J14/UTM_b014_Ecoregions2017_ecoregion_code_rast_exp.tif";(test=raster(test)); x11(); plot(test)

#' J15
#' b015_Ecoregions2017_ecoregion_code_rast_exp.
hf <-"D:/GISbases/b015JaguarAlbersBuffer70/b015_Ecoregions2017_ecoregion_code_rast_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J15trk))
#crs(hf)   
#compareCRS(hf,get_crs(J15trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J15/UTM_b015_Ecoregions2017_ecoregion_code_rast_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J15/UTM_b015_Ecoregions2017_ecoregion_code_rast_exp.tif";(test=raster(test)); x11(); plot(test)

#' J16
#' b016_Ecoregions2017_ecoregion_code_rast_exp.
hf <-"D:/GISbases/b016JaguarAlbersBuffer70/b016_Ecoregions2017_ecoregion_code_rast_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J16trk))
#crs(hf)   
#compareCRS(hf,get_crs(J16trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J16/UTM_b016_Ecoregions2017_ecoregion_code_rast_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J16/UTM_b016_Ecoregions2017_ecoregion_code_rast_exp.tif";(test=raster(test)); x11(); plot(test)

#' J17
#' b017_Ecoregions2017_ecoregion_code_rast_exp.
hf <-"D:/GISbases/b017JaguarAlbersBuffer70/b017_Ecoregions2017_ecoregion_code_rast_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J17trk))
#crs(hf)   
#compareCRS(hf,get_crs(J17trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J17/UTM_b017_Ecoregions2017_ecoregion_code_rast_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J17/UTM_b017_Ecoregions2017_ecoregion_code_rast_exp.tif";(test=raster(test)); x11(); plot(test)

#' J18
#' b018_Ecoregions2017_ecoregion_code_rast_exp.
hf <-"D:/GISbases/b018JaguarAlbersBuffer70/b018_Ecoregions2017_ecoregion_code_rast_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J18trk))
#crs(hf)   
#compareCRS(hf,get_crs(J18trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J18/UTM_b018_Ecoregions2017_ecoregion_code_rast_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J18/UTM_b018_Ecoregions2017_ecoregion_code_rast_exp.tif";(test=raster(test)); x11(); plot(test)

#' J19
#' b019_Ecoregions2017_ecoregion_code_rast_exp.
hf <-"D:/GISbases/b019JaguarAlbersBuffer70/b019_Ecoregions2017_ecoregion_code_rast_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J19trk))
#crs(hf)   
#compareCRS(hf,get_crs(J19trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J19/UTM_b019_Ecoregions2017_ecoregion_code_rast_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J19/UTM_b019_Ecoregions2017_ecoregion_code_rast_exp.tif";(test=raster(test)); x11(); plot(test)

#' J20
#' b020_Ecoregions2017_ecoregion_code_rast_exp.
hf <-"D:/GISbases/b020JaguarAlbersBuffer70/b020_Ecoregions2017_ecoregion_code_rast_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J20trk))
#crs(hf)   
#compareCRS(hf,get_crs(J20trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J20/UTM_b020_Ecoregions2017_ecoregion_code_rast_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J20/UTM_b020_Ecoregions2017_ecoregion_code_rast_exp.tif";(test=raster(test)); x11(); plot(test)

#' J21
#' b021_Ecoregions2017_ecoregion_code_rast_exp.
hf <-"D:/GISbases/b021JaguarAlbersBuffer70/b021_Ecoregions2017_ecoregion_code_rast_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J21trk))
#crs(hf)   
#compareCRS(hf,get_crs(J21trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J21/UTM_b021_Ecoregions2017_ecoregion_code_rast_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J21/UTM_b021_Ecoregions2017_ecoregion_code_rast_exp.tif";(test=raster(test)); x11(); plot(test)

#' J22
#' b022_Ecoregions2017_ecoregion_code_rast_exp.
hf <-"D:/GISbases/b022JaguarAlbersBuffer70/b022_Ecoregions2017_ecoregion_code_rast_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J22trk))
#crs(hf)   
#compareCRS(hf,get_crs(J22trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J22/UTM_b022_Ecoregions2017_ecoregion_code_rast_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J22/UTM_b022_Ecoregions2017_ecoregion_code_rast_exp.tif";(test=raster(test)); x11(); plot(test)

#' J23
#' b023_Ecoregions2017_ecoregion_code_rast_exp.
hf <-"D:/GISbases/b023JaguarAlbersBuffer70/b023_Ecoregions2017_ecoregion_code_rast_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J23trk))
#crs(hf)   
#compareCRS(hf,get_crs(J23trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J23/UTM_b023_Ecoregions2017_ecoregion_code_rast_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J23/UTM_b023_Ecoregions2017_ecoregion_code_rast_exp.tif";(test=raster(test)); x11(); plot(test)

#' J24
#' b024_Ecoregions2017_ecoregion_code_rast_exp.
hf <-"D:/GISbases/b024JaguarAlbersBuffer70/b024_Ecoregions2017_ecoregion_code_rast_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J24trk))
#crs(hf)   
#compareCRS(hf,get_crs(J24trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J24/UTM_b024_Ecoregions2017_ecoregion_code_rast_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J24/UTM_b024_Ecoregions2017_ecoregion_code_rast_exp.tif";(test=raster(test)); x11(); plot(test)

#' J25
#' b025_Ecoregions2017_ecoregion_code_rast_exp.
hf <-"D:/GISbases/b025JaguarAlbersBuffer70/b025_Ecoregions2017_ecoregion_code_rast_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J25trk))
#crs(hf)   
#compareCRS(hf,get_crs(J25trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J25/UTM_b025_Ecoregions2017_ecoregion_code_rast_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J25/UTM_b025_Ecoregions2017_ecoregion_code_rast_exp.tif";(test=raster(test)); x11(); plot(test)

#' J26
#' b026_Ecoregions2017_ecoregion_code_rast_exp.
hf <-"D:/GISbases/b026JaguarAlbersBuffer70/b026_Ecoregions2017_ecoregion_code_rast_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J26trk))
#crs(hf)   
#compareCRS(hf,get_crs(J26trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J26/UTM_b026_Ecoregions2017_ecoregion_code_rast_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J26/UTM_b026_Ecoregions2017_ecoregion_code_rast_exp.tif";(test=raster(test)); x11(); plot(test)

#' J27
#' b027_Ecoregions2017_ecoregion_code_rast_exp.
hf <-"D:/GISbases/b027JaguarAlbersBuffer70/b027_Ecoregions2017_ecoregion_code_rast_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J27trk))
#crs(hf)   
#compareCRS(hf,get_crs(J27trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J27/UTM_b027_Ecoregions2017_ecoregion_code_rast_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J27/UTM_b027_Ecoregions2017_ecoregion_code_rast_exp.tif";(test=raster(test)); x11(); plot(test)

#' J28
#' b028_Ecoregions2017_ecoregion_code_rast_exp.
hf <-"D:/GISbases/b028JaguarAlbersBuffer70/b028_Ecoregions2017_ecoregion_code_rast_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J28trk))
#crs(hf)   
#compareCRS(hf,get_crs(J28trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J28/UTM_b028_Ecoregions2017_ecoregion_code_rast_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J28/UTM_b028_Ecoregions2017_ecoregion_code_rast_exp.tif";(test=raster(test)); x11(); plot(test)

#' J29
#' b029_Ecoregions2017_ecoregion_code_rast_exp.
hf <-"D:/GISbases/b029JaguarAlbersBuffer70/b029_Ecoregions2017_ecoregion_code_rast_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J29trk))
#crs(hf)   
#compareCRS(hf,get_crs(J29trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J29/UTM_b029_Ecoregions2017_ecoregion_code_rast_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J29/UTM_b029_Ecoregions2017_ecoregion_code_rast_exp.tif";(test=raster(test)); x11(); plot(test)

#' J30
#' b030_Ecoregions2017_ecoregion_code_rast_exp.
hf <-"D:/GISbases/b030JaguarAlbersBuffer70/b030_Ecoregions2017_ecoregion_code_rast_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J30trk))
#crs(hf)   
#compareCRS(hf,get_crs(J30trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J30/UTM_b030_Ecoregions2017_ecoregion_code_rast_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J30/UTM_b030_Ecoregions2017_ecoregion_code_rast_exp.tif";(test=raster(test)); x11(); plot(test)

#' J30
#' b030_Ecoregions2017_ecoregion_code_rast_exp.
hf <-"D:/GISbases/b030JaguarAlbersBuffer70/b030_Ecoregions2017_ecoregion_code_rast_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J30trk))
#crs(hf)   
#compareCRS(hf,get_crs(J30trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J30/UTM_b030_Ecoregions2017_ecoregion_code_rast_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J30/UTM_b030_Ecoregions2017_ecoregion_code_rast_exp.tif";(test=raster(test)); x11(); plot(test)

#' J31
#' b031_Ecoregions2017_ecoregion_code_rast_exp.
hf <-"D:/GISbases/b031JaguarAlbersBuffer70/b031_Ecoregions2017_ecoregion_code_rast_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J31trk))
#crs(hf)   
#compareCRS(hf,get_crs(J31trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J31/UTM_b031_Ecoregions2017_ecoregion_code_rast_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J31/UTM_b031_Ecoregions2017_ecoregion_code_rast_exp.tif";(test=raster(test)); x11(); plot(test)

#' J32
#' b032_Ecoregions2017_ecoregion_code_rast_exp.
hf <-"D:/GISbases/b032JaguarAlbersBuffer70/b032_Ecoregions2017_ecoregion_code_rast_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J32trk))
#crs(hf)   
#compareCRS(hf,get_crs(J32trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J32/UTM_b032_Ecoregions2017_ecoregion_code_rast_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J32/UTM_b032_Ecoregions2017_ecoregion_code_rast_exp.tif";(test=raster(test)); x11(); plot(test)

#' J33
#' b033_Ecoregions2017_ecoregion_code_rast_exp.
hf <-"D:/GISbases/b033JaguarAlbersBuffer70/b033_Ecoregions2017_ecoregion_code_rast_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J33trk))
#crs(hf)   
#compareCRS(hf,get_crs(J33trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J33/UTM_b033_Ecoregions2017_ecoregion_code_rast_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J33/UTM_b033_Ecoregions2017_ecoregion_code_rast_exp.tif";(test=raster(test)); x11(); plot(test)

#' J34
#' b034_Ecoregions2017_ecoregion_code_rast_exp.
hf <-"D:/GISbases/b034JaguarAlbersBuffer70/b034_Ecoregions2017_ecoregion_code_rast_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J34trk))
#crs(hf)   
#compareCRS(hf,get_crs(J34trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J34/UTM_b034_Ecoregions2017_ecoregion_code_rast_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J34/UTM_b034_Ecoregions2017_ecoregion_code_rast_exp.tif";(test=raster(test)); x11(); plot(test)

#' J35
#' b035_Ecoregions2017_ecoregion_code_rast_exp.
hf <-"D:/GISbases/b035JaguarAlbersBuffer70/b035_Ecoregions2017_ecoregion_code_rast_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J35trk))
#crs(hf)   
#compareCRS(hf,get_crs(J35trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J35/UTM_b035_Ecoregions2017_ecoregion_code_rast_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J35/UTM_b035_Ecoregions2017_ecoregion_code_rast_exp.tif";(test=raster(test)); x11(); plot(test)

#' J36
#' b036_Ecoregions2017_ecoregion_code_rast_exp.
hf <-"D:/GISbases/b036JaguarAlbersBuffer70/b036_Ecoregions2017_ecoregion_code_rast_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J36trk))
#crs(hf)   
#compareCRS(hf,get_crs(J36trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J36/UTM_b036_Ecoregions2017_ecoregion_code_rast_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J36/UTM_b036_Ecoregions2017_ecoregion_code_rast_exp.tif";(test=raster(test)); x11(); plot(test)

#' J37
#' b037_Ecoregions2017_ecoregion_code_rast_exp.
hf <-"D:/GISbases/b037JaguarAlbersBuffer70/b037_Ecoregions2017_ecoregion_code_rast_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J37trk))
#crs(hf)   
#compareCRS(hf,get_crs(J37trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J37/UTM_b037_Ecoregions2017_ecoregion_code_rast_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J37/UTM_b037_Ecoregions2017_ecoregion_code_rast_exp.tif";(test=raster(test)); x11(); plot(test)

#' J38
#' b038_Ecoregions2017_ecoregion_code_rast_exp.
hf <-"D:/GISbases/b038JaguarAlbersBuffer70/b038_Ecoregions2017_ecoregion_code_rast_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J38trk))
#crs(hf)   
#compareCRS(hf,get_crs(J38trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J38/UTM_b038_Ecoregions2017_ecoregion_code_rast_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J38/UTM_b038_Ecoregions2017_ecoregion_code_rast_exp.tif";(test=raster(test)); x11(); plot(test)

#' J39
#' b039_Ecoregions2017_ecoregion_code_rast_exp.
hf <-"D:/GISbases/b039JaguarAlbersBuffer70/b039_Ecoregions2017_ecoregion_code_rast_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J39trk))
#crs(hf)   
#compareCRS(hf,get_crs(J39trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J39/UTM_b039_Ecoregions2017_ecoregion_code_rast_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J39/UTM_b039_Ecoregions2017_ecoregion_code_rast_exp.tif";(test=raster(test)); x11(); plot(test)

#' J40
#' b040_Ecoregions2017_ecoregion_code_rast_exp.
hf <-"D:/GISbases/b040JaguarAlbersBuffer70/b040_Ecoregions2017_ecoregion_code_rast_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J40trk))
#crs(hf)   
#compareCRS(hf,get_crs(J40trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J40/UTM_b040_Ecoregions2017_ecoregion_code_rast_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J40/UTM_b040_Ecoregions2017_ecoregion_code_rast_exp.tif";(test=raster(test)); x11(); plot(test)

#' J41
#' b041_Ecoregions2017_ecoregion_code_rast_exp.
hf <-"D:/GISbases/b041JaguarAlbersBuffer70/b041_Ecoregions2017_ecoregion_code_rast_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J41trk))
#crs(hf)   
#compareCRS(hf,get_crs(J41trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J41/UTM_b041_Ecoregions2017_ecoregion_code_rast_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J41/UTM_b041_Ecoregions2017_ecoregion_code_rast_exp.tif";(test=raster(test)); x11(); plot(test)

#' J42
#' b042_Ecoregions2017_ecoregion_code_rast_exp.
hf <-"D:/GISbases/b042JaguarAlbersBuffer70/b042_Ecoregions2017_ecoregion_code_rast_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J42trk))
#crs(hf)   
#compareCRS(hf,get_crs(J42trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J42/UTM_b042_Ecoregions2017_ecoregion_code_rast_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J42/UTM_b042_Ecoregions2017_ecoregion_code_rast_exp.tif";(test=raster(test)); x11(); plot(test)

gc()
gc()
gc()

#' J43
#' b043_Ecoregions2017_ecoregion_code_rast_exp.
hf <-"D:/GISbases/b043JaguarAlbersBuffer70/b043_Ecoregions2017_ecoregion_code_rast_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J43trk))
#crs(hf)   
#compareCRS(hf,get_crs(J43trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J43/UTM_b043_Ecoregions2017_ecoregion_code_rast_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J43/UTM_b043_Ecoregions2017_ecoregion_code_rast_exp.tif";(test=raster(test)); x11(); plot(test)

#' J44
#' b044_Ecoregions2017_ecoregion_code_rast_exp.
hf <-"D:/GISbases/b044JaguarAlbersBuffer70/b044_Ecoregions2017_ecoregion_code_rast_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J44trk))
#crs(hf)   
#compareCRS(hf,get_crs(J44trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J44/UTM_b044_Ecoregions2017_ecoregion_code_rast_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J44/UTM_b044_Ecoregions2017_ecoregion_code_rast_exp.tif";(test=raster(test)); x11(); plot(test)

#' J45
#' b045_Ecoregions2017_ecoregion_code_rast_exp.
hf <-"D:/GISbases/b045JaguarAlbersBuffer70/b045_Ecoregions2017_ecoregion_code_rast_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J45trk))
#crs(hf)   
#compareCRS(hf,get_crs(J45trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J45/UTM_b045_Ecoregions2017_ecoregion_code_rast_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J45/UTM_b045_Ecoregions2017_ecoregion_code_rast_exp.tif";(test=raster(test)); x11(); plot(test)

#' J46
#' b046_Ecoregions2017_ecoregion_code_rast_exp.
hf <-"D:/GISbases/b046JaguarAlbersBuffer70/b046_Ecoregions2017_ecoregion_code_rast_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J46trk))
#crs(hf)   
#compareCRS(hf,get_crs(J46trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J46/UTM_b046_Ecoregions2017_ecoregion_code_rast_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J46/UTM_b046_Ecoregions2017_ecoregion_code_rast_exp.tif";(test=raster(test)); x11(); plot(test)

#' J47
#' b047_Ecoregions2017_ecoregion_code_rast_exp.
hf <-"D:/GISbases/b047JaguarAlbersBuffer70/b047_Ecoregions2017_ecoregion_code_rast_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J47trk))
#crs(hf)   
#compareCRS(hf,get_crs(J47trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J47/UTM_b047_Ecoregions2017_ecoregion_code_rast_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J47/UTM_b047_Ecoregions2017_ecoregion_code_rast_exp.tif";(test=raster(test)); x11(); plot(test)

#' J48
#' b048_Ecoregions2017_ecoregion_code_rast_exp.
hf <-"D:/GISbases/b048JaguarAlbersBuffer70/b048_Ecoregions2017_ecoregion_code_rast_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J48trk))
#crs(hf)   
#compareCRS(hf,get_crs(J48trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J48/UTM_b048_Ecoregions2017_ecoregion_code_rast_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J48/UTM_b048_Ecoregions2017_ecoregion_code_rast_exp.tif";(test=raster(test)); x11(); plot(test)

#' J49
#' b049_Ecoregions2017_ecoregion_code_rast_exp.
hf <-"D:/GISbases/b049JaguarAlbersBuffer70/b049_Ecoregions2017_ecoregion_code_rast_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J49trk))
#crs(hf)   
#compareCRS(hf,get_crs(J49trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J49/UTM_b049_Ecoregions2017_ecoregion_code_rast_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J49/UTM_b049_Ecoregions2017_ecoregion_code_rast_exp.tif";(test=raster(test)); x11(); plot(test)

#' J50
#' b050_Ecoregions2017_ecoregion_code_rast_exp.
hf <-"D:/GISbases/b050JaguarAlbersBuffer70/b050_Ecoregions2017_ecoregion_code_rast_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J50trk))
#crs(hf)   
#compareCRS(hf,get_crs(J50trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J50/UTM_b050_Ecoregions2017_ecoregion_code_rast_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J50/UTM_b050_Ecoregions2017_ecoregion_code_rast_exp.tif";(test=raster(test)); x11(); plot(test)

#' J51
#' b051_Ecoregions2017_ecoregion_code_rast_exp.
hf <-"D:/GISbases/b051JaguarAlbersBuffer70/b051_Ecoregions2017_ecoregion_code_rast_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J51trk))
#crs(hf)   
#compareCRS(hf,get_crs(J51trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J51/UTM_b051_Ecoregions2017_ecoregion_code_rast_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J51/UTM_b051_Ecoregions2017_ecoregion_code_rast_exp.tif";(test=raster(test)); x11(); plot(test)

#' J52
#' b052_Ecoregions2017_ecoregion_code_rast_exp.
hf <-"D:/GISbases/b052JaguarAlbersBuffer70/b052_Ecoregions2017_ecoregion_code_rast_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J52trk))
#crs(hf)   
#compareCRS(hf,get_crs(J52trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J52/UTM_b052_Ecoregions2017_ecoregion_code_rast_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J52/UTM_b052_Ecoregions2017_ecoregion_code_rast_exp.tif";(test=raster(test)); x11(); plot(test)

#' J53
#' b053_Ecoregions2017_ecoregion_code_rast_exp.
hf <-"D:/GISbases/b053JaguarAlbersBuffer70/b053_Ecoregions2017_ecoregion_code_rast_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J53trk))
#crs(hf)   
#compareCRS(hf,get_crs(J53trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J53/UTM_b053_Ecoregions2017_ecoregion_code_rast_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J53/UTM_b053_Ecoregions2017_ecoregion_code_rast_exp.tif";(test=raster(test)); x11(); plot(test)

#' J54
#' b054_Ecoregions2017_ecoregion_code_rast_exp.
hf <-"D:/GISbases/b054JaguarAlbersBuffer70/b054_Ecoregions2017_ecoregion_code_rast_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J54trk))
#crs(hf)   
#compareCRS(hf,get_crs(J54trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J54/UTM_b054_Ecoregions2017_ecoregion_code_rast_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J54/UTM_b054_Ecoregions2017_ecoregion_code_rast_exp.tif";(test=raster(test)); x11(); plot(test)

#' J55
#' b055_Ecoregions2017_ecoregion_code_rast_exp.
hf <-"D:/GISbases/b055JaguarAlbersBuffer70/b055_Ecoregions2017_ecoregion_code_rast_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J55trk))
#crs(hf)   
#compareCRS(hf,get_crs(J55trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J55/UTM_b055_Ecoregions2017_ecoregion_code_rast_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J55/UTM_b055_Ecoregions2017_ecoregion_code_rast_exp.tif";(test=raster(test)); x11(); plot(test)

#' J56
#' b056_Ecoregions2017_ecoregion_code_rast_exp.
hf <-"D:/GISbases/b056JaguarAlbersBuffer70/b056_Ecoregions2017_ecoregion_code_rast_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J56trk))
#crs(hf)   
#compareCRS(hf,get_crs(J56trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J56/UTM_b056_Ecoregions2017_ecoregion_code_rast_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J56/UTM_b056_Ecoregions2017_ecoregion_code_rast_exp.tif";(test=raster(test)); x11(); plot(test)

#' J57
#' b057_Ecoregions2017_ecoregion_code_rast_exp.
hf <-"D:/GISbases/b057JaguarAlbersBuffer70/b057_Ecoregions2017_ecoregion_code_rast_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J57trk))
#crs(hf)   
#compareCRS(hf,get_crs(J57trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J57/UTM_b057_Ecoregions2017_ecoregion_code_rast_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J57/UTM_b057_Ecoregions2017_ecoregion_code_rast_exp.tif";(test=raster(test)); x11(); plot(test)

#' J58
#' b058_Ecoregions2017_ecoregion_code_rast_exp.
hf <-"D:/GISbases/b058JaguarAlbersBuffer70/b058_Ecoregions2017_ecoregion_code_rast_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J58trk))
#crs(hf)   
#compareCRS(hf,get_crs(J58trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J58/UTM_b058_Ecoregions2017_ecoregion_code_rast_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J58/UTM_b058_Ecoregions2017_ecoregion_code_rast_exp.tif";(test=raster(test)); x11(); plot(test)

#' J59
#' b059_Ecoregions2017_ecoregion_code_rast_exp.
hf <-"D:/GISbases/b059JaguarAlbersBuffer70/b059_Ecoregions2017_ecoregion_code_rast_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J59trk))
#crs(hf)   
#compareCRS(hf,get_crs(J59trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J59/UTM_b059_Ecoregions2017_ecoregion_code_rast_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J59/UTM_b059_Ecoregions2017_ecoregion_code_rast_exp.tif";(test=raster(test)); x11(); plot(test)

#' J60
#' b060_Ecoregions2017_ecoregion_code_rast_exp.
hf <-"D:/GISbases/b060JaguarAlbersBuffer70/b060_Ecoregions2017_ecoregion_code_rast_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J60trk))
#crs(hf)   
#compareCRS(hf,get_crs(J60trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J60/UTM_b060_Ecoregions2017_ecoregion_code_rast_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J60/UTM_b060_Ecoregions2017_ecoregion_code_rast_exp.tif";(test=raster(test)); x11(); plot(test)

#' J61
#' b061_Ecoregions2017_ecoregion_code_rast_exp.
hf <-"D:/GISbases/b061JaguarAlbersBuffer70/b061_Ecoregions2017_ecoregion_code_rast_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J61trk))
#crs(hf)   
#compareCRS(hf,get_crs(J61trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J61/UTM_b061_Ecoregions2017_ecoregion_code_rast_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J61/UTM_b061_Ecoregions2017_ecoregion_code_rast_exp.tif";(test=raster(test)); x11(); plot(test)

#' J62
#' b062_Ecoregions2017_ecoregion_code_rast_exp.
hf <-"D:/GISbases/b062JaguarAlbersBuffer70/b062_Ecoregions2017_ecoregion_code_rast_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J62trk))
#crs(hf)   
#compareCRS(hf,get_crs(J62trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J62/UTM_b062_Ecoregions2017_ecoregion_code_rast_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J62/UTM_b062_Ecoregions2017_ecoregion_code_rast_exp.tif";(test=raster(test)); x11(); plot(test)

#' J63
#' b063_Ecoregions2017_ecoregion_code_rast_exp.
hf <-"D:/GISbases/b063JaguarAlbersBuffer70/b063_Ecoregions2017_ecoregion_code_rast_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J63trk))
#crs(hf)   
#compareCRS(hf,get_crs(J63trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J63/UTM_b063_Ecoregions2017_ecoregion_code_rast_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J63/UTM_b063_Ecoregions2017_ecoregion_code_rast_exp.tif";(test=raster(test)); x11(); plot(test)

gc()
gc()
gc()


#' J64
#' b064_Ecoregions2017_ecoregion_code_rast_exp.
hf <-"D:/GISbases/b064JaguarAlbersBuffer70/b064_Ecoregions2017_ecoregion_code_rast_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J64trk))
#crs(hf)   
#compareCRS(hf,get_crs(J64trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J64/UTM_b064_Ecoregions2017_ecoregion_code_rast_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J64/UTM_b064_Ecoregions2017_ecoregion_code_rast_exp.tif";(test=raster(test)); x11(); plot(test)

#' J65
#' b065_Ecoregions2017_ecoregion_code_rast_exp.
hf <-"D:/GISbases/b065JaguarAlbersBuffer70/b065_Ecoregions2017_ecoregion_code_rast_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J65trk))
#crs(hf)   
#compareCRS(hf,get_crs(J65trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J65/UTM_b065_Ecoregions2017_ecoregion_code_rast_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J65/UTM_b065_Ecoregions2017_ecoregion_code_rast_exp.tif";(test=raster(test)); x11(); plot(test)

#' J66
#' b066_Ecoregions2017_ecoregion_code_rast_exp.
hf <-"D:/GISbases/b066JaguarAlbersBuffer70/b066_Ecoregions2017_ecoregion_code_rast_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J66trk))
#crs(hf)   
#compareCRS(hf,get_crs(J66trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J66/UTM_b066_Ecoregions2017_ecoregion_code_rast_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J66/UTM_b066_Ecoregions2017_ecoregion_code_rast_exp.tif";(test=raster(test)); x11(); plot(test)

#' J67
#' b067_Ecoregions2017_ecoregion_code_rast_exp.
hf <-"D:/GISbases/b067JaguarAlbersBuffer70/b067_Ecoregions2017_ecoregion_code_rast_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J67trk))
#crs(hf)   
#compareCRS(hf,get_crs(J67trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J67/UTM_b067_Ecoregions2017_ecoregion_code_rast_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J67/UTM_b067_Ecoregions2017_ecoregion_code_rast_exp.tif";(test=raster(test)); x11(); plot(test)

#' J68
#' b068_Ecoregions2017_ecoregion_code_rast_exp.
hf <-"D:/GISbases/b068JaguarAlbersBuffer70/b068_Ecoregions2017_ecoregion_code_rast_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J68trk))
#crs(hf)   
#compareCRS(hf,get_crs(J68trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J68/UTM_b068_Ecoregions2017_ecoregion_code_rast_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J68/UTM_b068_Ecoregions2017_ecoregion_code_rast_exp.tif";(test=raster(test)); x11(); plot(test)

#' J69
#' b069_Ecoregions2017_ecoregion_code_rast_exp.
hf <-"D:/GISbases/b069JaguarAlbersBuffer70/b069_Ecoregions2017_ecoregion_code_rast_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J69trk))
#crs(hf)   
#compareCRS(hf,get_crs(J69trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J69/UTM_b069_Ecoregions2017_ecoregion_code_rast_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J69/UTM_b069_Ecoregions2017_ecoregion_code_rast_exp.tif";(test=raster(test)); x11(); plot(test)

gc()
gc()
gc()


#' J70
#' b070_Ecoregions2017_ecoregion_code_rast_exp.
hf <-"D:/GISbases/b070JaguarAlbersBuffer70/b070_Ecoregions2017_ecoregion_code_rast_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J70trk))
#crs(hf)   
#compareCRS(hf,get_crs(J70trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J70/UTM_b070_Ecoregions2017_ecoregion_code_rast_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J70/UTM_b070_Ecoregions2017_ecoregion_code_rast_exp.tif";(test=raster(test)); x11(); plot(test)

#' J71
#' b071_Ecoregions2017_ecoregion_code_rast_exp.
hf <-"D:/GISbases/b071JaguarAlbersBuffer70/b071_Ecoregions2017_ecoregion_code_rast_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J71trk))
#crs(hf)   
#compareCRS(hf,get_crs(J71trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J71/UTM_b071_Ecoregions2017_ecoregion_code_rast_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J71/UTM_b071_Ecoregions2017_ecoregion_code_rast_exp.tif";(test=raster(test)); x11(); plot(test)

#' J72
#' b072_Ecoregions2017_ecoregion_code_rast_exp.
hf <-"D:/GISbases/b072JaguarAlbersBuffer70/b072_Ecoregions2017_ecoregion_code_rast_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J72trk))
#crs(hf)   
#compareCRS(hf,get_crs(J72trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J72/UTM_b072_Ecoregions2017_ecoregion_code_rast_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J72/UTM_b072_Ecoregions2017_ecoregion_code_rast_exp.tif";(test=raster(test)); x11(); plot(test)

#' J73
#' b073_Ecoregions2017_ecoregion_code_rast_exp.
hf <-"D:/GISbases/b073JaguarAlbersBuffer70/b073_Ecoregions2017_ecoregion_code_rast_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J73trk))
#crs(hf)   
#compareCRS(hf,get_crs(J73trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J73/UTM_b073_Ecoregions2017_ecoregion_code_rast_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J73/UTM_b073_Ecoregions2017_ecoregion_code_rast_exp.tif";(test=raster(test)); x11(); plot(test)

#' J74
#' b074_Ecoregions2017_ecoregion_code_rast_exp.
hf <-"D:/GISbases/b074JaguarAlbersBuffer70/b074_Ecoregions2017_ecoregion_code_rast_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J74trk))
#crs(hf)   
#compareCRS(hf,get_crs(J74trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J74/UTM_b074_Ecoregions2017_ecoregion_code_rast_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J74/UTM_b074_Ecoregions2017_ecoregion_code_rast_exp.tif";(test=raster(test)); x11(); plot(test)

#' J75
#' b075_Ecoregions2017_ecoregion_code_rast_exp.
hf <-"D:/GISbases/b075JaguarAlbersBuffer70/b075_Ecoregions2017_ecoregion_code_rast_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J75trk))
#crs(hf)   
#compareCRS(hf,get_crs(J75trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J75/UTM_b075_Ecoregions2017_ecoregion_code_rast_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J75/UTM_b075_Ecoregions2017_ecoregion_code_rast_exp.tif";(test=raster(test)); x11(); plot(test)

#' J76
#' b076_Ecoregions2017_ecoregion_code_rast_exp.
hf <-"D:/GISbases/b076JaguarAlbersBuffer70/b076_Ecoregions2017_ecoregion_code_rast_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J76trk))
#crs(hf)   
#compareCRS(hf,get_crs(J76trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J76/UTM_b076_Ecoregions2017_ecoregion_code_rast_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J76/UTM_b076_Ecoregions2017_ecoregion_code_rast_exp.tif";(test=raster(test)); x11(); plot(test)

#' J77
#' b077_Ecoregions2017_ecoregion_code_rast_exp.
hf <-"D:/GISbases/b077JaguarAlbersBuffer70/b077_Ecoregions2017_ecoregion_code_rast_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J77trk))
#crs(hf)   
#compareCRS(hf,get_crs(J77trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J77/UTM_b077_Ecoregions2017_ecoregion_code_rast_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J77/UTM_b077_Ecoregions2017_ecoregion_code_rast_exp.tif";(test=raster(test)); x11(); plot(test)

#' J78
#' b078_Ecoregions2017_ecoregion_code_rast_exp.
hf <-"D:/GISbases/b078JaguarAlbersBuffer70/b078_Ecoregions2017_ecoregion_code_rast_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J78trk))
#crs(hf)   
#compareCRS(hf,get_crs(J78trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J78/UTM_b078_Ecoregions2017_ecoregion_code_rast_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J78/UTM_b078_Ecoregions2017_ecoregion_code_rast_exp.tif";(test=raster(test)); x11(); plot(test)

#' J79
#' b079_Ecoregions2017_ecoregion_code_rast_exp.
hf <-"D:/GISbases/b079JaguarAlbersBuffer70/b079_Ecoregions2017_ecoregion_code_rast_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J79trk))
#crs(hf)   
#compareCRS(hf,get_crs(J79trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J79/UTM_b079_Ecoregions2017_ecoregion_code_rast_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J79/UTM_b079_Ecoregions2017_ecoregion_code_rast_exp.tif";(test=raster(test)); x11(); plot(test)

#' J80
#' b080_Ecoregions2017_ecoregion_code_rast_exp.
hf <-"D:/GISbases/b080JaguarAlbersBuffer70/b080_Ecoregions2017_ecoregion_code_rast_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J80trk))
#crs(hf)   
#compareCRS(hf,get_crs(J80trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J80/UTM_b080_Ecoregions2017_ecoregion_code_rast_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J80/UTM_b080_Ecoregions2017_ecoregion_code_rast_exp.tif";(test=raster(test)); x11(); plot(test)

#' J81
#' b081_Ecoregions2017_ecoregion_code_rast_exp.
hf <-"D:/GISbases/b081JaguarAlbersBuffer70/b081_Ecoregions2017_ecoregion_code_rast_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J81trk))
#crs(hf)   
#compareCRS(hf,get_crs(J81trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J81/UTM_b081_Ecoregions2017_ecoregion_code_rast_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J81/UTM_b081_Ecoregions2017_ecoregion_code_rast_exp.tif";(test=raster(test)); x11(); plot(test)

#' J82
#' b082_Ecoregions2017_ecoregion_code_rast_exp.
hf <-"D:/GISbases/b082JaguarAlbersBuffer70/b082_Ecoregions2017_ecoregion_code_rast_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J82trk))
#crs(hf)   
#compareCRS(hf,get_crs(J82trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J82/UTM_b082_Ecoregions2017_ecoregion_code_rast_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J82/UTM_b082_Ecoregions2017_ecoregion_code_rast_exp.tif";(test=raster(test)); x11(); plot(test)

#' J83
#' b083_Ecoregions2017_ecoregion_code_rast_exp.
hf <-"D:/GISbases/b083JaguarAlbersBuffer70/b083_Ecoregions2017_ecoregion_code_rast_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J83trk))
#crs(hf)   
#compareCRS(hf,get_crs(J83trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J83/UTM_b083_Ecoregions2017_ecoregion_code_rast_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J83/UTM_b083_Ecoregions2017_ecoregion_code_rast_exp.tif";(test=raster(test)); x11(); plot(test)

#' J84
#' b084_Ecoregions2017_ecoregion_code_rast_exp.
hf <-"D:/GISbases/b084JaguarAlbersBuffer70/b084_Ecoregions2017_ecoregion_code_rast_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J84trk))
#crs(hf)   
#compareCRS(hf,get_crs(J84trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J84/UTM_b084_Ecoregions2017_ecoregion_code_rast_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J84/UTM_b084_Ecoregions2017_ecoregion_code_rast_exp.tif";(test=raster(test)); x11(); plot(test)

#' J85
#' b085_Ecoregions2017_ecoregion_code_rast_exp.
hf <-"D:/GISbases/b085JaguarAlbersBuffer70/b085_Ecoregions2017_ecoregion_code_rast_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J85trk))
#crs(hf)   
#compareCRS(hf,get_crs(J85trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J85/UTM_b085_Ecoregions2017_ecoregion_code_rast_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J85/UTM_b085_Ecoregions2017_ecoregion_code_rast_exp.tif";(test=raster(test)); x11(); plot(test)

#' J86
#' b086_Ecoregions2017_ecoregion_code_rast_exp.
hf <-"D:/GISbases/b086JaguarAlbersBuffer70/b086_Ecoregions2017_ecoregion_code_rast_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J86trk))
#crs(hf)   
#compareCRS(hf,get_crs(J86trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J86/UTM_b086_Ecoregions2017_ecoregion_code_rast_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J86/UTM_b086_Ecoregions2017_ecoregion_code_rast_exp.tif";(test=raster(test)); x11(); plot(test)

#' J87
#' b087_Ecoregions2017_ecoregion_code_rast_exp.
hf <-"D:/GISbases/b087JaguarAlbersBuffer70/b087_Ecoregions2017_ecoregion_code_rast_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J87trk))
#crs(hf)   
#compareCRS(hf,get_crs(J87trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J87/UTM_b087_Ecoregions2017_ecoregion_code_rast_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J87/UTM_b087_Ecoregions2017_ecoregion_code_rast_exp.tif";(test=raster(test)); x11(); plot(test)

#' J88
#' b088_Ecoregions2017_ecoregion_code_rast_exp.
hf <-"D:/GISbases/b088JaguarAlbersBuffer70/b088_Ecoregions2017_ecoregion_code_rast_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J88trk))
#crs(hf)   
#compareCRS(hf,get_crs(J88trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J88/UTM_b088_Ecoregions2017_ecoregion_code_rast_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J88/UTM_b088_Ecoregions2017_ecoregion_code_rast_exp.tif";(test=raster(test)); x11(); plot(test)

#' J89
#' b089_Ecoregions2017_ecoregion_code_rast_exp.
hf <-"D:/GISbases/b089JaguarAlbersBuffer70/b089_Ecoregions2017_ecoregion_code_rast_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J89trk))
#crs(hf)   
#compareCRS(hf,get_crs(J89trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J89/UTM_b089_Ecoregions2017_ecoregion_code_rast_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J89/UTM_b089_Ecoregions2017_ecoregion_code_rast_exp.tif";(test=raster(test)); x11(); plot(test)

#' J90
#' b090_Ecoregions2017_ecoregion_code_rast_exp.
hf <-"D:/GISbases/b090JaguarAlbersBuffer70/b090_Ecoregions2017_ecoregion_code_rast_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J90trk))
#crs(hf)   
#compareCRS(hf,get_crs(J90trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J90/UTM_b090_Ecoregions2017_ecoregion_code_rast_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J90/UTM_b090_Ecoregions2017_ecoregion_code_rast_exp.tif";(test=raster(test)); x11(); plot(test)

#' J91
#' b091_Ecoregions2017_ecoregion_code_rast_exp.
hf <-"D:/GISbases/b091JaguarAlbersBuffer70/b091_Ecoregions2017_ecoregion_code_rast_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J91trk))
#crs(hf)   
#compareCRS(hf,get_crs(J91trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J91/UTM_b091_Ecoregions2017_ecoregion_code_rast_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J91/UTM_b091_Ecoregions2017_ecoregion_code_rast_exp.tif";(test=raster(test)); x11(); plot(test)

#' J92
#' b092_Ecoregions2017_ecoregion_code_rast_exp.
hf <-"D:/GISbases/b092JaguarAlbersBuffer70/b092_Ecoregions2017_ecoregion_code_rast_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J92trk))
#crs(hf)   
#compareCRS(hf,get_crs(J92trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J92/UTM_b092_Ecoregions2017_ecoregion_code_rast_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J92/UTM_b092_Ecoregions2017_ecoregion_code_rast_exp.tif";(test=raster(test)); x11(); plot(test)

#' J93
#' b093_Ecoregions2017_ecoregion_code_rast_exp.
hf <-"D:/GISbases/b093JaguarAlbersBuffer70/b093_Ecoregions2017_ecoregion_code_rast_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J93trk))
#crs(hf)   
#compareCRS(hf,get_crs(J93trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J93/UTM_b093_Ecoregions2017_ecoregion_code_rast_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J93/UTM_b093_Ecoregions2017_ecoregion_code_rast_exp.tif";(test=raster(test)); x11(); plot(test)

#' J94
#' b094_Ecoregions2017_ecoregion_code_rast_exp.
hf <-"D:/GISbases/b094JaguarAlbersBuffer70/b094_Ecoregions2017_ecoregion_code_rast_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J94trk))
#crs(hf)   
#compareCRS(hf,get_crs(J94trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J94/UTM_b094_Ecoregions2017_ecoregion_code_rast_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J94/UTM_b094_Ecoregions2017_ecoregion_code_rast_exp.tif";(test=raster(test)); x11(); plot(test)

#' J95
#' b095_Ecoregions2017_ecoregion_code_rast_exp.
hf <-"D:/GISbases/b095JaguarAlbersBuffer70/b095_Ecoregions2017_ecoregion_code_rast_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J95trk))
#crs(hf)   
#compareCRS(hf,get_crs(J95trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J95/UTM_b095_Ecoregions2017_ecoregion_code_rast_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J95/UTM_b095_Ecoregions2017_ecoregion_code_rast_exp.tif";(test=raster(test)); x11(); plot(test)

#' J96
#' b096_Ecoregions2017_ecoregion_code_rast_exp.
hf <-"D:/GISbases/b096JaguarAlbersBuffer70/b096_Ecoregions2017_ecoregion_code_rast_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J96trk))
#crs(hf)   
#compareCRS(hf,get_crs(J96trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J96/UTM_b096_Ecoregions2017_ecoregion_code_rast_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J96/UTM_b096_Ecoregions2017_ecoregion_code_rast_exp.tif";(test=raster(test)); x11(); plot(test)

#' J97
#' b097_Ecoregions2017_ecoregion_code_rast_exp.
hf <-"D:/GISbases/b097JaguarAlbersBuffer70/b097_Ecoregions2017_ecoregion_code_rast_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J97trk))
#crs(hf)   
#compareCRS(hf,get_crs(J97trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J97/UTM_b097_Ecoregions2017_ecoregion_code_rast_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J97/UTM_b097_Ecoregions2017_ecoregion_code_rast_exp.tif";(test=raster(test)); x11(); plot(test)

#' J98
#' b098_Ecoregions2017_ecoregion_code_rast_exp.
hf <-"D:/GISbases/b098JaguarAlbersBuffer70/b098_Ecoregions2017_ecoregion_code_rast_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J98trk))
#crs(hf)   
#compareCRS(hf,get_crs(J98trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J98/UTM_b098_Ecoregions2017_ecoregion_code_rast_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J98/UTM_b098_Ecoregions2017_ecoregion_code_rast_exp.tif";(test=raster(test)); x11(); plot(test)

#' J99
#' b099_Ecoregions2017_ecoregion_code_rast_exp.
hf <-"D:/GISbases/b099JaguarAlbersBuffer70/b099_Ecoregions2017_ecoregion_code_rast_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J99trk))
#crs(hf)   
#compareCRS(hf,get_crs(J99trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J99/UTM_b099_Ecoregions2017_ecoregion_code_rast_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J99/UTM_b099_Ecoregions2017_ecoregion_code_rast_exp.tif";(test=raster(test)); x11(); plot(test)

#' J100
#' b100_Ecoregions2017_ecoregion_code_rast_exp.
hf <-"D:/GISbases/b100JaguarAlbersBuffer70/b100_Ecoregions2017_ecoregion_code_rast_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J100trk))
#crs(hf)   
#compareCRS(hf,get_crs(J100trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J100/UTM_b100_Ecoregions2017_ecoregion_code_rast_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J100/UTM_b010_Ecoregions2017_ecoregion_code_rast_exp.tif";(test=raster(test)); x11(); plot(test)

#' J101
#' b101_Ecoregions2017_ecoregion_code_rast_exp.
hf <-"D:/GISbases/b101JaguarAlbersBuffer70/b101_Ecoregions2017_ecoregion_code_rast_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J101trk))
#crs(hf)   
#compareCRS(hf,get_crs(J101trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J101/UTM_b101_Ecoregions2017_ecoregion_code_rast_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J101/UTM_b101_Ecoregions2017_ecoregion_code_rast_exp.tif";(test=raster(test)); x11(); plot(test)

#' J102
#' b102_Ecoregions2017_ecoregion_code_rast_exp.
hf <-"D:/GISbases/b102JaguarAlbersBuffer70/b102_Ecoregions2017_ecoregion_code_rast_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J102trk))
#crs(hf)   
#compareCRS(hf,get_crs(J102trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J102/UTM_b102_Ecoregions2017_ecoregion_code_rast_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J102/UTM_b102_Ecoregions2017_ecoregion_code_rast_exp.tif";(test=raster(test)); x11(); plot(test)

#' J103
#' b103_Ecoregions2017_ecoregion_code_rast_exp.
hf <-"D:/GISbases/b103JaguarAlbersBuffer70/b103_Ecoregions2017_ecoregion_code_rast_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J103trk))
#crs(hf)   
#compareCRS(hf,get_crs(J103trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J103/UTM_b103_Ecoregions2017_ecoregion_code_rast_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J103/UTM_b103_Ecoregions2017_ecoregion_code_rast_exp.tif";(test=raster(test)); x11(); plot(test)

#' J104
#' b104_Ecoregions2017_ecoregion_code_rast_exp.
hf <-"D:/GISbases/b104JaguarAlbersBuffer70/b104_Ecoregions2017_ecoregion_code_rast_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J104trk))
#crs(hf)   
#compareCRS(hf,get_crs(J104trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J104/UTM_b104_Ecoregions2017_ecoregion_code_rast_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J104/UTM_b104_Ecoregions2017_ecoregion_code_rast_exp.tif";(test=raster(test)); x11(); plot(test)

#' J105
#' b105_Ecoregions2017_ecoregion_code_rast_exp.
hf <-"D:/GISbases/b105JaguarAlbersBuffer70/b105_Ecoregions2017_ecoregion_code_rast_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J105trk))
#crs(hf)   
#compareCRS(hf,get_crs(J105trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J105/UTM_b105_Ecoregions2017_ecoregion_code_rast_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J105/UTM_b105_Ecoregions2017_ecoregion_code_rast_exp.tif";(test=raster(test)); x11(); plot(test)

#' J106
#' b106_Ecoregions2017_ecoregion_code_rast_exp.
hf <-"D:/GISbases/b106JaguarAlbersBuffer70/b106_Ecoregions2017_ecoregion_code_rast_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J106trk))
#crs(hf)   
#compareCRS(hf,get_crs(J106trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J106/UTM_b106_Ecoregions2017_ecoregion_code_rast_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J106/UTM_b106_Ecoregions2017_ecoregion_code_rast_exp.tif";(test=raster(test)); x11(); plot(test)

#' J107
#' b107_Ecoregions2017_ecoregion_code_rast_exp.
hf <-"D:/GISbases/b107JaguarAlbersBuffer70/b107_Ecoregions2017_ecoregion_code_rast_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J107trk))
#crs(hf)   
#compareCRS(hf,get_crs(J107trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J107/UTM_b107_Ecoregions2017_ecoregion_code_rast_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J107/UTM_b107_Ecoregions2017_ecoregion_code_rast_exp.tif";(test=raster(test)); x11(); plot(test)

#' J108
#' b108_Ecoregions2017_ecoregion_code_rast_exp.
hf <-"D:/GISbases/b108JaguarAlbersBuffer70/b108_Ecoregions2017_ecoregion_code_rast_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J108trk))
#crs(hf)   
#compareCRS(hf,get_crs(J108trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J108/UTM_b108_Ecoregions2017_ecoregion_code_rast_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J108/UTM_b108_Ecoregions2017_ecoregion_code_rast_exp.tif";(test=raster(test)); x11(); plot(test)

#' J109
#' b109_Ecoregions2017_ecoregion_code_rast_exp.
hf <-"D:/GISbases/b109JaguarAlbersBuffer70/b109_Ecoregions2017_ecoregion_code_rast_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J109trk))
#crs(hf)   
#compareCRS(hf,get_crs(J109trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J109/UTM_b109_Ecoregions2017_ecoregion_code_rast_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J109/UTM_b109_Ecoregions2017_ecoregion_code_rast_exp.tif";(test=raster(test)); x11(); plot(test)

#' J110
#' b110_Ecoregions2017_ecoregion_code_rast_exp.
hf <-"D:/GISbases/b110JaguarAlbersBuffer70/b110_Ecoregions2017_ecoregion_code_rast_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J110trk))
#crs(hf)   
#compareCRS(hf,get_crs(J110trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J110/UTM_b110_Ecoregions2017_ecoregion_code_rast_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J110/UTM_b110_Ecoregions2017_ecoregion_code_rast_exp.tif";(test=raster(test)); x11(); plot(test)

#' J111
#' b111_Ecoregions2017_ecoregion_code_rast_exp.
hf <-"D:/GISbases/b111JaguarAlbersBuffer70/b111_Ecoregions2017_ecoregion_code_rast_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J111trk))
#crs(hf)   
#compareCRS(hf,get_crs(J111trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J111/UTM_b111_Ecoregions2017_ecoregion_code_rast_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J111/UTM_b111_Ecoregions2017_ecoregion_code_rast_exp.tif";(test=raster(test)); x11(); plot(test)

#' J112
#' b112_Ecoregions2017_ecoregion_code_rast_exp.
hf <-"D:/GISbases/b112JaguarAlbersBuffer70/b112_Ecoregions2017_ecoregion_code_rast_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J112trk))
#crs(hf)   
#compareCRS(hf,get_crs(J112trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J112/UTM_b112_Ecoregions2017_ecoregion_code_rast_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J112/UTM_b112_Ecoregions2017_ecoregion_code_rast_exp.tif";(test=raster(test)); x11(); plot(test)

#' J113
#' b113_Ecoregions2017_ecoregion_code_rast_exp.
hf <-"D:/GISbases/b113JaguarAlbersBuffer70/b113_Ecoregions2017_ecoregion_code_rast_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J113trk))
#crs(hf)   
#compareCRS(hf,get_crs(J113trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J113/UTM_b113_Ecoregions2017_ecoregion_code_rast_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J113/UTM_b113_Ecoregions2017_ecoregion_code_rast_exp.tif";(test=raster(test)); x11(); plot(test)

#' J114
#' b114_Ecoregions2017_ecoregion_code_rast_exp.
hf <-"D:/GISbases/b114JaguarAlbersBuffer70/b114_Ecoregions2017_ecoregion_code_rast_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J114trk))
#crs(hf)   
#compareCRS(hf,get_crs(J114trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J114/UTM_b114_Ecoregions2017_ecoregion_code_rast_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J114/UTM_b114_Ecoregions2017_ecoregion_code_rast_exp.tif";(test=raster(test)); x11(); plot(test)

#' J115
#' b115_Ecoregions2017_ecoregion_code_rast_exp.
hf <-"D:/GISbases/b115JaguarAlbersBuffer70/b115_Ecoregions2017_ecoregion_code_rast_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J115trk))
#crs(hf)   
#compareCRS(hf,get_crs(J115trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J115/UTM_b115_Ecoregions2017_ecoregion_code_rast_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J115/UTM_b115_Ecoregions2017_ecoregion_code_rast_exp.tif";(test=raster(test)); x11(); plot(test)

#' J116
#' b116_Ecoregions2017_ecoregion_code_rast_exp.
hf <-"D:/GISbases/b116JaguarAlbersBuffer70/b116_Ecoregions2017_ecoregion_code_rast_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J116trk))
#crs(hf)   
#compareCRS(hf,get_crs(J116trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J116/UTM_b116_Ecoregions2017_ecoregion_code_rast_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J116/UTM_b116_Ecoregions2017_ecoregion_code_rast_exp.tif";(test=raster(test)); x11(); plot(test)

#' J117
#' b117_Ecoregions2017_ecoregion_code_rast_exp.
hf <-"D:/GISbases/b117JaguarAlbersBuffer70/b117_Ecoregions2017_ecoregion_code_rast_exp.tif"
hf=raster(hf)
#crs(hf)   
#plot(hf)
hf <- projectRaster(hf, crs="+proj=longlat +datum=WGS84")
#crs(hf)   
hf <- projectRaster(hf, crs=get_crs(J117trk))
#crs(hf)   
#compareCRS(hf,get_crs(J117trk))
writeRaster(hf,filename=file.path("D:/GISUTM/J117/UTM_b117_Ecoregions2017_ecoregion_code_rast_exp."), format="GTiff", overwrite=TRUE)
# test <-"D:/GISUTM/J117/UTM_b117_Ecoregions2017_ecoregion_code_rast_exp.tif";(test=raster(test)); x11(); plot(test)

