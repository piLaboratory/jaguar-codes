###########################################################################################################
                            ############# Jaguar Data Preparation ##############   
###########################################################################################################						
######################################## Preparation ######################################################
rm(list= ls())                                                  ### For a fresh start
## Commented because is not generic
setwd("C:/RWorkDir/jaguardatapaper")                           ### Set directory
###########################################################################################################         																					
###Script modified from Bernardo Niebuhr & MoveBank ("John Fieberg")         																					
###Enter jaguar data from Morato et al. 2018
## Load packages
if(!require(install.load)) install.packages('install.load'); library(install.load)
install.load::install_load("ggmap","maptools",'move',"circular","RCurl","dplyr","readr","caTools","adehabitatLT","ctmm","ggsn","magick","rgl","tidyverse","htmltools","rglwidget","lubridate","raster","amt","tibble","knitr","leaflet","ezknitr","lattice","rgdal","sp")
### The require() can be used inside functions as it gives a warning message and returns a logical value. 
### FALSE if the requested package is not found and TRUE if the package is loaded.###(Movement data should be checked and cleaned (for outliers, duplicated timestamps, etc)first!!! 
### Load and adjust the data (already cleaned) and create a dataframe object
mov.data.org <- read.delim(file="c:/RWorkDir/jaguardatapaper/mov.data.org.txt")

mov.data.org <- dplyr::select(mov.data.org, -(individual.taxon.canonical.name:tag.local.identifier))
str(mov.data.org)

# Add Individual info
info <- read.delim(file="c:/RWorkDir/jaguardatapaper/info.txt")

#ind.info <- read.delim(file="c:/RWorkDir/jaguardatapaper/Jaguar_additional information.txt")
#info <- ind.info %>%
#   dplyr::select(ID, Sex, Estimated.Age, Weight, Collar.Type, Collar.Brand, Planned.Schedule)
#info <- info %>% rename(individual.local.identifier..ID.=ID)
#info <- info %>% rename(sex=Sex)
#info <- info %>% rename(age=Estimated.Age)
#info <- info %>% rename(weight=Weight)
#info <- info %>% rename(collar_type=Collar.Type)
#info <- info %>% rename(brand=Collar.Brand)
#info <- info %>% rename(schedule=Planned.Schedule)
### Movement Parameters (ctmm)
#movpar <- read.delim(file="c:/RWorkDir/jaguardatapaper/movementparameters.txt")
#str(movpar)
#movpar = movpar[-118,]   # there were an extra row of NA so I deleted that
#movpar <- movpar %>%
   #dplyr::select(Animal_ID, Model)
#movpar <- movpar %>% rename(individual.local.identifier..ID.=Animal_ID)
#movpar <- movpar %>% rename(model=Model)
#Merge movpar with id info and save txt
#info <- merge(info,movpar)   
#info <- info[with(info,order(individual.local.identifier..ID.)),]
#write.table(info,file="c:/RWorkDir/jaguardatapaper/info.txt",row.names = F,quote=F,col.names=T,sep="\t")

#Merge movement with individual info/parameters
merged<- merge(mov.data.org,info)
mov.data.org <- merged 
str(mov.data.org)
#write.table(info,file="c:/RWorkDir/jaguardatapaper/mov.data.org.txt",row.names = F,quote=F,col.names=T,sep="\t")
##########################

##########################
# Organize data
#mov.data.org 

# Add 2000 to years
get.year <- function(time.stamp) {
  init <- gregexpr('/', time.stamp, fixed = T)[[1]][2] + 1
  end <- gregexpr(' ', time.stamp, fixed = T)[[1]][1] - 1
  substr(time.stamp, init, end)
}

# Test
get.year(time.stamp = mov.data.org$timestamp[10000])
# All individuals
year <- as.numeric(sapply(mov.data.org$timestamp, get.year))
table(year)

# Add 1900/2000
new.year <- as.character(ifelse(year > 50, year + 1900, year + 2000))
table(new.year)

# New dates
set.year <- function(time.stamp, year) {
  init <- gregexpr('/', time.stamp, fixed = T)[[1]][2]
  end <- gregexpr(' ', time.stamp, fixed = T)[[1]][1]
  paste(substr(time.stamp, 1, init), year,
        substr(time.stamp, end, nchar(time.stamp)), sep = "")
}


# Test
set.year(time.stamp = as.character(mov.data.org$timestamp[10000]), year = '2013')
# All individuals
date.time <- as.character(mapply(set.year, as.character(mov.data.org$timestamp),
                                 new.year))
str(date.time)
#date.time
#########################################################

# All individuals
date.time <- as.character(mapply(set.year, as.character(mov.data.org$timestamp),
                                 new.year))
str(date.time)
#date.time
# Date/Time as POSIXct object
mov.data.org$timestamp.posix <- as.POSIXct(date.time, 
                                           format = "%m/%d/%Y %H:%M", tz = 'GMT')
str(mov.data.org)

##################################################################################################
###  Get local time!!!

# I included a column to represent the local timezone (already with the - signal) to them multiply the timestamp and get the difference:
mov.data.org$local_time <- mov.data.org$timestamp.posix + mov.data.org$timezone*60*60
mov.data.org$timestamp.posix.GMT <- mov.data.org$timestamp.posix
mov.data.org$timestamp.posix <- mov.data.org$local_time ### If we do that all the (timestamp.posix)'s calculations will be based on local time
str(mov.data.org)

################################################################

## adehabitatLT
# Transforms in ltraj object
coords <- data.frame(mov.data.org$location.long, mov.data.org$location.lat)
mov.traj <- as.ltraj(xy = coords, date=mov.data.org$timestamp.posix, 
                     id=mov.data.org$individual.local.identifier..ID., 
                     burst=mov.data.org$individual.local.identifier..ID., 
                     infolocs = mov.data.org[,-c(3:4, ncol(mov.data.org))])
mov.traj.df <- ld(mov.traj)
#mov.traj
#plot(mov.traj)

## move
# Organize data as a move package format
move.data <- move(x = mov.traj.df$x, y = mov.traj.df$y, 
                  time = mov.traj.df$date, 
                  proj = CRS('+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'),
                  data = mov.traj.df, animal = mov.traj.df$id, sensor = 'GPS')
				  
				  
#'+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs' "+proj=aea +lat_1=-5 +lat_2=-42 +lat_0=-32 +lon_0=-60 +x_0=0 +y_0=0 +ellps=aust_SA +towgs84=-57,1,-41,0,0,0,0 +units=m +no_defs"
move.data  ### moveStack
summary(move.data)

###Separate individual animals' trajectories 
unstacked <- split(move.data)
head(unstacked)
jaguar_df <- as(move.data, "data.frame")
#id <- as.integer(levels(jaguar_df$id))[jaguar_df$id]
age <- as.numeric(levels(jaguar_df$age))[jaguar_df$age]
weight <- as.numeric(levels(jaguar_df$weight))[jaguar_df$weight]
jaguar_df$id <- as.factor(jaguar_df$individual.local.identifier..ID.)  
jaguar_df$age <- age   ### converted to number
jaguar_df$weight <- weight   ### converted to number

head(jaguar_df)
str(jaguar_df)


#########################################################################################################							   							   									  
#' ######################## More Data cleaning ##########################################
###############################################################
#' Delete observations where missing lat or long or a timestamp.  There are no missing
#' observations in this data set, but it is still good practice to check.
ind<-complete.cases(jaguar_df[,c("y","x","date")])
jaguar_df<-jaguar_df[ind==TRUE,]

#' Check for duplicated observations (ones with same lat, long, timestamp,
#'  and individual identifier). There are no duplicate
#' observations in this data set, but it is still good practice to check.
ind2<-jaguar_df %>% select("date","x","y","id") %>% duplicated
sum(ind2) # no duplicates
jaguar_df<-jaguar_df[ind2!=TRUE,]
###or duplicatedLocs <- which(jaguar_df$date[1:(nrow(jaguar_df)-1)] == jaguar_df$date[2:(nrow(jaguar_df))])

###########   Clean suspectly close points!!!  Above 1200 secs or 20 min minimal interval ####################
excludes <- filter(jaguar_df, dt < 1200)
### see excludeds based on id
table(select(excludes, id))             
removed<- anti_join(jaguar_df, excludes)   ### test if were all removed
filter(removed, dt < 1200)
jaguar_df <- removed ###########   Clean suspectly close points!!!  Above 1200 secs or 20 min minimal interval #######################
str(jaguar_df)

library(dplyr)

# You don't need this line if you already have DateTime in proper format
df$DateTime <- as.POSIXct(df$DateTime)

# Add a date column (with whatever timezone you want)
df$date <- as.Date(df$DateTime, tz = 'EST')

# Following generates the sunrise and sunset times for the two example dates
sunRise <- c(as.POSIXct('2016-04-15 06:40:37'), as.POSIXct('2016-03-24 06:55:00'))
sunSet <- c(as.POSIXct('2016-04-15 18:40:37'), as.POSIXct('2016-03-24 18:25:00'))
sun <- data.frame(date = as.Date(sunRise, tz = 'EST'), sunRise = sunRise, sunSet = sunSet)
sun
        date             sunRise              sunSet
1 2016-04-15 2016-04-15 06:40:37 2016-04-15 18:40:37
2 2016-03-24 2016-03-24 06:55:00 2016-03-24 18:25:00

# Join the two tables and compute night/day
df <- inner_join(df, sun)
df$dayNight <- ifelse(df$DateTime > df$sunRise & df$DateTime < df$sunSet, 'day', 'night')




## Cleaning up columns which will be in excess due to repetition of analysis ########################################################################

jaguar_df$dx <- NULL
jaguar_df$dy <- NULL
jaguar_df$dist <- NULL
jaguar_df$dt <- NULL
jaguar_df$R2n <- NULL
jaguar_df$abs.angle <- NULL
jaguar_df$rel.angle  <- NULL
jaguar_df$location.lat <- NULL
jaguar_df$timestamps  <- NULL
jaguar_df$sensor <- NULL
jaguar_df$burst <- NULL
jaguar_df$optional <- NULL
jaguar_df$coords.x1 <- NULL
jaguar_df$coords.x2 <- NULL
jaguar_df$trackId <- NULL
jaguar_df$individual.local.identifier..ID.<- NULL
jaguar_df$study.name <- NULL
jaguar_df$collar_type <- NULL  
jaguar_df$brand <- NULL
jaguar_df$local_time <- NULL
jaguar_df$timezone <- NULL
head(jaguar_df) 
str(jaguar_df)


###############################################################################################################################
###  Add UTM  #################
###############################################################################################################################
###Separate individual animals' trajectories 
#unstacked <- split(move.data)
#head(unstacked)
#jaguar_df <- as(move.data, "data.frame")
table(jaguar_df$project_region)

#write.table(jaguar_df,file="c:/RWorkDir/jaguardatapaper/jaguar_df.txt",row.names = F,quote=F,col.names=T,sep="\t")
#jaguar_df <- read.delim(file="c:/RWorkDir/jaguardatapaper/mov.data.org.txt")
###############################################################################################################################
                                    ###   Atlantic Forest   ###
###############################################################################################################################
### # Atlantic Forest W1 
AFW1=subset(jaguar_df,project_region=='Atlantic Forest W1')
head(AFW1)
str(AFW1)
table(AFW1$id)
ft=ftable(AFW1$id)
f<-as.data.frame.table(ft)
f$Var1 <- NULL
f$Var2 <- NULL
subset(f,Freq>0)
range(AFW1$x);range(AFW1$y) 
coord.latlong = SpatialPoints(cbind(AFW1$x,AFW1$y), proj4string = CRS("+proj=longlat +datum=WGS84"))
coord.latlong
# Transforming coordinate to UTM 
coord.UTM <- spTransform(coord.latlong , CRS("+proj=utm +zone=22K +south +datum=WGS84 +units=m +no_defs"))
coord.UTM
coord.latlong.df <- as.data.frame(coord.latlong)
cootd.utm.df <- as.data.frame(coord.UTM)
head(coord.latlong.df)
head(cootd.utm.df)
par(mfrow = c(1, 2))
plot(coord.latlong.df, axes = TRUE, main = "Lat-Long Coordinates", cex.axis = 0.95)
plot(cootd.utm.df, axes = TRUE, main = "UTM Coordinates", col = "red", cex.axis = 0.95)
locsj_matx=cbind(coordinates(coord.UTM), coordinates(coord.latlong))
locsj_df<- as.data.frame(locsj_matx)
head(locsj_df)
head(AFW1)
point.names<-c("utm_x","utm_y","long_x","lat_y")
colnames(locsj_df)<-point.names
head(locsj_df)

AFW1=cbind(AFW1, locsj_df)
head(AFW1); str(AFW1)

#############################################
### # Atlantic Forest W2
AFW2=subset(jaguar_df,project_region=='Atlantic Forest W2')
head(AFW2)
str(AFW2)
table(AFW2$id)
ft=ftable(AFW2$id)
f<-as.data.frame.table(ft)
f$Var1 <- NULL
f$Var2 <- NULL
subset(f,Freq>0)
range(AFW2$x);range(AFW2$y) 
coord.latlong = SpatialPoints(cbind(AFW2$x,AFW2$y), proj4string = CRS("+proj=longlat +datum=WGS84"))
coord.latlong
# Transforming coordinate to UTM 
coord.UTM <- spTransform(coord.latlong , CRS("+proj=utm +zone=22K +south +datum=WGS84 +units=m +no_defs"))
coord.UTM
coord.latlong.df <- as.data.frame(coord.latlong)
cootd.utm.df <- as.data.frame(coord.UTM)
head(coord.latlong.df)
head(cootd.utm.df)
par(mfrow = c(1, 2))
plot(coord.latlong.df, axes = TRUE, main = "Lat-Long Coordinates", cex.axis = 0.95)
plot(cootd.utm.df, axes = TRUE, main = "UTM Coordinates", col = "red", cex.axis = 0.95)
locsj_matx=cbind(coordinates(coord.UTM), coordinates(coord.latlong))
locsj_df<- as.data.frame(locsj_matx)
head(locsj_df)
head(AFW2)
point.names<-c("utm_x","utm_y","long_x","lat_y")
colnames(locsj_df)<-point.names
head(locsj_df)

AFW2=cbind(AFW2, locsj_df)
head(AFW2)

AFW=rbind(AFW1,AFW2)
head(AFW)
ft=ftable(AFW$id)
f<-as.data.frame.table(ft)
f$Var1 <- NULL
f$Var2 <- NULL
subset(f,Freq>0)

#write.table(AFW,file="c:/RWorkDir/jaguardatapaper/AFW.txt",row.names = F,quote=F,col.names=T,sep="\t")
#AFW <- read.delim(file="c:/RWorkDir/jaguardatapaper/AFW.txt")
head(AFW); str(AFW)


###############################################################################################################################
                                    ###   CAATINGA  ###
###############################################################################################################################
### # Caatinga
###############################################################################################################################
Caatinga=subset(jaguar_df,project_region=='Caatinga')
head(Caatinga)
str(Caatinga)
table(Caatinga$id)
ft=ftable(Caatinga$id)
f<-as.data.frame.table(ft)
f$Var1 <- NULL
f$Var2 <- NULL
subset(f,Freq>0)
range(Caatinga$x);range(Caatinga$y) 
coord.latlong = SpatialPoints(cbind(Caatinga$x,Caatinga$y), proj4string = CRS("+proj=longlat +datum=WGS84"))
coord.latlong
# Transforming coordinate to UTM 
coord.UTM <- spTransform(coord.latlong , CRS("+proj=utm +zone=23L +south +datum=WGS84 +units=m +no_defs"))
coord.UTM
coord.latlong.df <- as.data.frame(coord.latlong)
cootd.utm.df <- as.data.frame(coord.UTM)
head(coord.latlong.df)
head(cootd.utm.df)
par(mfrow = c(1, 2))
plot(coord.latlong.df, axes = TRUE, main = "Lat-Long Coordinates", cex.axis = 0.95)
plot(cootd.utm.df, axes = TRUE, main = "UTM Coordinates", col = "red", cex.axis = 0.95)
locsj_matx=cbind(coordinates(coord.UTM ), coordinates(coord.latlong))
locsj_df<- as.data.frame(locsj_matx)
head(locsj_df)
head(Caatinga)
point.names<-c("utm_x","utm_y","long_x","lat_y")
colnames(locsj_df)<-point.names
head(locsj_df)

Caatinga=cbind(Caatinga, locsj_df)
head(Caatinga); str(Caatinga)
#write.table(Caatinga,file="c:/RWorkDir/jaguardatapaper/Caatinga.txt",row.names = F,quote=F,col.names=T,sep="\t")
#Caatinga <- read.delim(file="c:/RWorkDir/jaguardatapaper/Caatinga.txt")

###############################################################################################################################
                                    ###   CERRADO   ###
###############################################################################################################################
### # Cerrado1
###############################################################################################################################
Cerrado1=subset(jaguar_df,project_region=='Cerrado1')
head(Cerrado1)
str(Cerrado1)
table(Cerrado1$id)
ft=ftable(Cerrado1$id)
f<-as.data.frame.table(ft)
f$Var1 <- NULL
f$Var2 <- NULL
subset(f,Freq>0)
range(Cerrado1$x);range(Cerrado1$y) 
coord.latlong = SpatialPoints(cbind(Cerrado1$x,Cerrado1$y), proj4string = CRS("+proj=longlat +datum=WGS84"))
coord.latlong
# Transforming coordinate to UTM 
coord.UTM <- spTransform(coord.latlong , CRS("+proj=utm +zone=22K +south +datum=WGS84 +units=m +no_defs"))
coord.UTM
coord.latlong.df <- as.data.frame(coord.latlong)
cootd.utm.df <- as.data.frame(coord.UTM)
head(coord.latlong.df)
head(cootd.utm.df)
par(mfrow = c(1, 2))
plot(coord.latlong.df, axes = TRUE, main = "Lat-Long Coordinates", cex.axis = 0.95)
plot(cootd.utm.df, axes = TRUE, main = "UTM Coordinates", col = "red", cex.axis = 0.95)
locsj_matx=cbind(coordinates(coord.UTM ), coordinates(coord.latlong))
locsj_df<- as.data.frame(locsj_matx)
head(locsj_df)
head(Cerrado1)
point.names<-c("utm_x","utm_y","long_x","lat_y")
colnames(locsj_df)<-point.names
head(locsj_df)

Cerrado1=cbind(Cerrado1, locsj_df)
head(Cerrado1)

###############################################################################################################################
### # Cerrado2
###############################################################################################################################
Cerrado2=subset(jaguar_df,project_region=='Cerrado2')
#jfl89=subset(jaguar_df,id=='89')  # This animal occurs in 2 UTM zones 22K (40%) and 22L (60%)
head(Cerrado2)
str(Cerrado2)
table(Cerrado2$id)
ft=ftable(Cerrado2$id)
f<-as.data.frame.table(ft)
f$Var1 <- NULL
f$Var2 <- NULL
subset(f,Freq>0)
range(Cerrado2$x);range(Cerrado2$y) 
coord.latlong = SpatialPoints(cbind(Cerrado2$x,Cerrado2$y), proj4string = CRS("+proj=longlat +datum=WGS84"))
coord.latlong
# Transforming coordinate to UTM 
coord.UTM <- spTransform(coord.latlong , CRS("+proj=utm +zone=22L +south +datum=WGS84 +units=m +no_defs"))
coord.UTM
coord.latlong.df <- as.data.frame(coord.latlong)
cootd.utm.df <- as.data.frame(coord.UTM)
head(coord.latlong.df)
head(cootd.utm.df)
par(mfrow = c(1, 2))
plot(coord.latlong.df, axes = TRUE, main = "Lat-Long Coordinates", cex.axis = 0.95)
plot(cootd.utm.df, axes = TRUE, main = "UTM Coordinates", col = "red", cex.axis = 0.95)
locsj_matx=cbind(coordinates(coord.UTM ), coordinates(coord.latlong))
locsj_df<- as.data.frame(locsj_matx)
head(locsj_df)
head(Cerrado2)
point.names<-c("utm_x","utm_y","long_x","lat_y")
colnames(locsj_df)<-point.names
head(locsj_df)

Cerrado2=cbind(Cerrado2, locsj_df)
head(Cerrado2)

Cerrado=rbind(Cerrado1,Cerrado2)
head(Cerrado)
ft=ftable(Cerrado$id)
f<-as.data.frame.table(ft)
f$Var1 <- NULL
f$Var2 <- NULL
subset(f,Freq>0)

#write.table(Cerrado,file="c:/RWorkDir/jaguardatapaper/Cerrado.txt",row.names = F,quote=F,col.names=T,sep="\t")
#Cerrado <- read.delim(file="c:/RWorkDir/jaguardatapaper/Cerrado.txt")

###############################################################################################################################
                                    ###   COSTA RICA  ###
###############################################################################################################################
### # Costa Rica
###############################################################################################################################
CRica=subset(jaguar_df,project_region=='Costa Rica')
head(CRica)
str(CRica)
table(CRica$id)
ft=ftable(CRica$id)
f<-as.data.frame.table(ft)
f$Var1 <- NULL
f$Var2 <- NULL
subset(f,Freq>0)
range(CRica$x);range(CRica$y) 
coord.latlong = SpatialPoints(cbind(CRica$x,CRica$y), proj4string = CRS("+proj=longlat +datum=WGS84"))
coord.latlong
# Transforming coordinate to UTM 
coord.UTM <- spTransform(coord.latlong , CRS("+proj=utm +zone=16P +north +datum=WGS84 +units=m +no_defs"))
coord.UTM
coord.latlong.df <- as.data.frame(coord.latlong)
cootd.utm.df <- as.data.frame(coord.UTM)
head(coord.latlong.df)
head(cootd.utm.df)
par(mfrow = c(1, 2))
plot(coord.latlong.df, axes = TRUE, main = "Lat-Long Coordinates", cex.axis = 0.95)
plot(cootd.utm.df, axes = TRUE, main = "UTM Coordinates", col = "red", cex.axis = 0.95)
locsj_matx=cbind(coordinates(coord.UTM ), coordinates(coord.latlong))
locsj_df<- as.data.frame(locsj_matx)
head(locsj_df)
head(CRica)
point.names<-c("utm_x","utm_y","long_x","lat_y")
colnames(locsj_df)<-point.names
head(locsj_df)

CRica=cbind(CRica, locsj_df)
head(CRica)

#write.table(CRica,file="c:/RWorkDir/jaguardatapaper/CRica.txt",row.names = F,quote=F,col.names=T,sep="\t")
#CRica<- read.delim(file="c:/RWorkDir/jaguardatapaper/CRica.txt")

###############################################################################################################################
                                    ###   PANTANAL  ###
###############################################################################################################################
### # PantanalTotal Brazil&Paraguay
###############################################################################################################################

#12    13    14    15    18    19    22    23    25    27
#28    29    30    31    32    33    41    51    52    53    
#54    55    56    57    59    60    61    68    69    74    
#75    79    81    84    86    87    88    91    92   101   
#102   103   104   105   106   107   108   109   110  111   
#112   113   114   115   116   117


Pantanal =subset(jaguar_df,project_bioveg=='Pantanal')
head(Pantanal)
str(Pantanal)
table(Pantanal$id); max(table(Pantanal$id)); min(table(Pantanal$id))
ft=ftable(Pantanal$id)
f<-as.data.frame.table(ft)
f$Var1 <- NULL
f$Var2 <- NULL
subset(f,Freq>0)
head(table(Pantanal$id,Pantanal$project_region))
head(ifelse(table(Pantanal$id,Pantanal$project_region)>0,1,0))
Pant_projs=colSums(ifelse(table(Pantanal$id,Pantanal$project_region)>0,1,0))
Pant_projs
range(Pantanal$x);range(Pantanal$y) 
coord.latlong = SpatialPoints(cbind(Pantanal$x,Pantanal$y), proj4string = CRS("+proj=longlat +datum=WGS84"))
coord.latlong
# Transforming coordinate to UTM 
coord.UTM <- spTransform(coord.latlong , CRS("+proj=utm +zone=21K +south +datum=WGS84 +units=m +no_defs"))
coord.UTM
coord.latlong.df <- as.data.frame(coord.latlong)
cootd.utm.df <- as.data.frame(coord.UTM)
head(coord.latlong.df)
head(cootd.utm.df)
par(mfrow = c(1, 2))
plot(coord.latlong.df, axes = TRUE, main = "Lat-Long Coordinates", cex.axis = 0.95)
plot(cootd.utm.df, axes = TRUE, main = "UTM Coordinates", col = "red", cex.axis = 0.95)
locsj_matx=cbind(coordinates(coord.UTM ), coordinates(coord.latlong))
locsj_df<- as.data.frame(locsj_matx)
head(locsj_df)
head(Pantanal)
point.names<-c("utm_x","utm_y","long_x","lat_y")
colnames(locsj_df)<-point.names
head(locsj_df)

Pantanal=cbind(Pantanal, locsj_df)
head(Pantanal)

#write.table(Pantanal,file="c:/RWorkDir/jaguardatapaper/Pantanal.txt",row.names = F,quote=F,col.names=T,sep="\t")
#Pantanal <- read.delim(file="c:/RWorkDir/jaguardatapaper/Pantanal.txt")

###############################################################################################################################
                                    ###   DRY CHACO  ###
###############################################################################################################################
### # Dry Chaco
###############################################################################################################################
Drych=subset(jaguar_df,project_region=='Dry chaco')
head(Drych)
str(Drych)
table(Drych$id)
ft=ftable(Drych$id)
f<-as.data.frame.table(ft)
f$Var1 <- NULL
f$Var2 <- NULL
subset(f,Freq>0)
range(Drych$x);range(Drych$y) 
coord.latlong = SpatialPoints(cbind(Drych$x,Drych$y), proj4string = CRS("+proj=longlat +datum=WGS84"))
coord.latlong
# Transforming coordinate to UTM   # 16,70 => (most 20K), 76,77 => (20 K);      71,72,73 => (21 K) 
X16=subset(Drych,id=='16')
X70=subset(Drych,id=='70')
X71=subset(Drych,id=='71')
X72=subset(Drych,id=='72')
X73=subset(Drych,id=='73')
X76=subset(Drych,id=='76')
X77=subset(Drych,id=='77')

###  Drych1 and Drych2
Drych1=rbind(X16,X70,X76,X77)
Drych2=rbind(X71,X72,X73)

###  Drych1    ## 76,77  Zone 20 K; 16 & 70 (most is in Zone 20 K too, but a little bit in 21K)
range(Drych1$x);range(Drych1$y) 
coord.latlong = SpatialPoints(cbind(Drych1$x,Drych1$y), proj4string = CRS("+proj=longlat +datum=WGS84"))
coord.latlong
coord.UTM <- spTransform(coord.latlong , CRS("+proj=utm +zone=20K +south +datum=WGS84 +units=m +no_defs"))
#coord.UTM <- spTransform(coord.latlong , CRS("+proj=utm +zone=21K +south +datum=WGS84 +units=m +no_defs"))

coord.UTM
coord.latlong.df <- as.data.frame(coord.latlong)
cootd.utm.df <- as.data.frame(coord.UTM)
head(coord.latlong.df)
head(cootd.utm.df)
par(mfrow = c(1, 2))
plot(coord.latlong.df, axes = TRUE, main = "Lat-Long Coordinates", cex.axis = 0.95)
plot(cootd.utm.df, axes = TRUE, main = "UTM Coordinates", col = "red", cex.axis = 0.95)
locsj_matx=cbind(coordinates(coord.UTM ), coordinates(coord.latlong))
locsj_df<- as.data.frame(locsj_matx)
head(locsj_df)
head(Drych1)
point.names<-c("utm_x","utm_y","long_x","lat_y")
#point.names<-c("utm_x1","utm_y1","long_x","lat_y");colnames(locsj_df)<-point.names; locsj_df$long_x <- NULL; locsj_df$lat_y <-NULL
#Drych3 <- Drych1; Drych3=cbind(Drych3, locsj_df); # ;head(Drych3)
colnames(locsj_df)<-point.names
head(locsj_df)
Drych1=cbind(Drych1, locsj_df)
head(Drych1)
ft=ftable(Drych1$id)
f<-as.data.frame.table(ft)
f$Var1 <- NULL
f$Var2 <- NULL
subset(f,Freq>0)



###  Drych2    ##  
range(Drych2$x);range(Drych2$y) 
coord.latlong = SpatialPoints(cbind(Drych2$x,Drych2$y), proj4string = CRS("+proj=longlat +datum=WGS84"))
coord.latlong
coord.UTM <- spTransform(coord.latlong , CRS("+proj=utm +zone=21K +south +datum=WGS84 +units=m +no_defs"))
coord.UTM
coord.latlong.df <- as.data.frame(coord.latlong)
cootd.utm.df <- as.data.frame(coord.UTM)
head(coord.latlong.df)
head(cootd.utm.df)
par(mfrow = c(1, 2))
plot(coord.latlong.df, axes = TRUE, main = "Lat-Long Coordinates", cex.axis = 0.95)
plot(cootd.utm.df, axes = TRUE, main = "UTM Coordinates", col = "red", cex.axis = 0.95)
locsj_matx=cbind(coordinates(coord.UTM ), coordinates(coord.latlong))
locsj_df<- as.data.frame(locsj_matx)
head(locsj_df)
head(Drych2)
point.names<-c("utm_x","utm_y","long_x","lat_y")
colnames(locsj_df)<-point.names
head(locsj_df)
Drych2=cbind(Drych2, locsj_df)
head(Drych2)
ft=ftable(Drych2$id)
f<-as.data.frame.table(ft)
f$Var1 <- NULL
f$Var2 <- NULL
subset(f,Freq>0)


Drych=rbind(Drych1,Drych2)
head(Drych)
ft=ftable(Drych$id)
f<-as.data.frame.table(ft)
f$Var1 <- NULL
f$Var2 <- NULL
subset(f,Freq>0)
#write.table(Drych,file="c:/RWorkDir/jaguardatapaper/Drych.txt",row.names = F,quote=F,col.names=T,sep="\t")
#Drych <- read.delim(file="c:/RWorkDir/jaguardatapaper/Drych.txt")

# Time Posix class (2014-12-05 00:00:00 )
Drych$timestamp.posix <- as.character(Drych$timestamp.posix)
Drych$timestamp.posix<- as.POSIXct(Drych$timestamp.posix, format ="%Y-%m-%d %H:%M:%S", tz = 'GMT') 
Drych$date <- as.character(Drych$date)
Drych$date <- as.POSIXct(Drych$date, format ="%Y-%m-%d %H:%M:%S", tz = 'GMT') 
Drych$id <- as.factor(Drych$id)
#Drych$season <- as.factor(Drych$season) #only dry, flood, inter seasons
#Drych$dt <- as.numeric(Drych$dt)
#Drych$week <- as.numeric(strftime(as.POSIXlt(Drych$timestamp.posix),format="%W"))
#Drych$day <- as.numeric(strftime(as.POSIXlt(Drych$timestamp.posix),format="%j"))
#Drych$year <- as.numeric(strftime(as.POSIXlt(Drych$timestamp.posix),format="%Y"))
str(Drych)

###############################################################################################################################
### # Humid Chaco
###############################################################################################################################
Hch=subset(jaguar_df,project_region=='Humid chaco')
head(Hch)
str(Hch)
table(Hch$id)
ft=ftable(Hch$id)
f<-as.data.frame.table(ft)
f$Var1 <- NULL
f$Var2 <- NULL
subset(f,Freq>0)
range(Hch$x);range(Hch$y) 
coord.latlong = SpatialPoints(cbind(Hch$x,Hch$y), proj4string = CRS("+proj=longlat +datum=WGS84"))
coord.latlong
# Transforming coordinate to UTM 
coord.UTM <- spTransform(coord.latlong , CRS("+proj=utm +zone=21K +south +datum=WGS84 +units=m +no_defs"))
coord.UTM
coord.latlong.df <- as.data.frame(coord.latlong)
cootd.utm.df <- as.data.frame(coord.UTM)
head(coord.latlong.df)
head(cootd.utm.df)
par(mfrow = c(1, 2))
plot(coord.latlong.df, axes = TRUE, main = "Lat-Long Coordinates", cex.axis = 0.95)
plot(cootd.utm.df, axes = TRUE, main = "UTM Coordinates", col = "red", cex.axis = 0.95)
locsj_matx=cbind(coordinates(coord.UTM ), coordinates(coord.latlong))
locsj_df<- as.data.frame(locsj_matx)
head(locsj_df)
head(Hch)
point.names<-c("utm_x","utm_y","long_x","lat_y")
colnames(locsj_df)<-point.names
head(locsj_df)

#write.table(Hch,file="c:/RWorkDir/jaguardatapaper/Hch.txt",row.names = F,quote=F,col.names=T,sep="\t")
#Hch <- read.delim(file="c:/RWorkDir/jaguardatapaper/Hch.txt")
Hch=cbind(Hch, locsj_df)
head(Hch); str(Hch)

###############################################################################################################################
###          # Forest Paraguay 
###############################################################################################################################
FPy=subset(jaguar_df,project_region=='Forest Paraguay')
head(FPy)
str(FPy)
table(FPy$id)
ft=ftable(FPy$id)
f<-as.data.frame.table(ft)
f$Var1 <- NULL
f$Var2 <- NULL
subset(f,Freq>0)
range(FPy$x);range(FPy$y) 
coord.latlong = SpatialPoints(cbind(FPy$x,FPy$y), proj4string = CRS("+proj=longlat +datum=WGS84"))
coord.latlong
# Transforming coordinate to UTM 
coord.UTM <- spTransform(coord.latlong , CRS("+proj=utm +zone=21J +south +datum=WGS84 +units=m +no_defs"))
coord.UTM
coord.latlong.df <- as.data.frame(coord.latlong)
cootd.utm.df <- as.data.frame(coord.UTM)
head(coord.latlong.df)
head(cootd.utm.df)
par(mfrow = c(1, 2))
plot(coord.latlong.df, axes = TRUE, main = "Lat-Long Coordinates", cex.axis = 0.95)
plot(cootd.utm.df, axes = TRUE, main = "UTM Coordinates", col = "red", cex.axis = 0.95)
locsj_matx=cbind(coordinates(coord.UTM ), coordinates(coord.latlong))
locsj_df<- as.data.frame(locsj_matx)
head(locsj_df)
head(FPy)
point.names<-c("utm_x","utm_y","long_x","lat_y")
colnames(locsj_df)<-point.names
head(locsj_df)

FPy=cbind(FPy, locsj_df)
#write.table(FPy,file="c:/RWorkDir/jaguardatapaper/FPy.txt",row.names = F,quote=F,col.names=T,sep="\t")
#FPy <- read.delim(file="c:/RWorkDir/jaguardatapaper/FPy.txt")
head(FPy); str(FPy)

###############################################################################################################################
###          # Iguazu 
###############################################################################################################################
Iguazu=subset(jaguar_df,project_region=='Iguazu')
head(Iguazu)
str(Iguazu)
table(Iguazu$id)
ft=ftable(Iguazu$id)
f<-as.data.frame.table(ft)
f$Var1 <- NULL
f$Var2 <- NULL
subset(f,Freq>0)
range(Iguazu$x);range(Iguazu$y) 
coord.latlong = SpatialPoints(cbind(Iguazu$x,Iguazu$y), proj4string = CRS("+proj=longlat +datum=WGS84"))
coord.latlong
# Transforming coordinate to UTM   # 42,66,80,90 => (21 J);      83 => (22 J) 
X42=subset(Iguazu,id=='42')
X66=subset(Iguazu,id=='66')
X80=subset(Iguazu,id=='80')
X90=subset(Iguazu,id=='90')

X83=subset(Iguazu,id=='83')

###  Iguazu1 and Iguazu2
Iguazu1=rbind(X42,X66,X80,X90)
Iguazu2=(X83)

###  Iguazu1    ## # 42,66,80,90 => ( Zone 21 J)
range(Iguazu1$x);range(Iguazu1$y) 
coord.latlong = SpatialPoints(cbind(Iguazu1$x,Iguazu1$y), proj4string = CRS("+proj=longlat +datum=WGS84"))
coord.latlong
coord.UTM <- spTransform(coord.latlong , CRS("+proj=utm +zone=21J +south +datum=WGS84 +units=m +no_defs"))
coord.UTM
coord.latlong.df <- as.data.frame(coord.latlong)
cootd.utm.df <- as.data.frame(coord.UTM)
head(coord.latlong.df)
head(cootd.utm.df)
par(mfrow = c(1, 2))
plot(coord.latlong.df, axes = TRUE, main = "Lat-Long Coordinates", cex.axis = 0.95)
plot(cootd.utm.df, axes = TRUE, main = "UTM Coordinates", col = "red", cex.axis = 0.95)
locsj_matx=cbind(coordinates(coord.UTM ), coordinates(coord.latlong))
locsj_df<- as.data.frame(locsj_matx)
head(locsj_df)
head(Iguazu1)
point.names<-c("utm_x","utm_y","long_x","lat_y")
colnames(locsj_df)<-point.names
head(locsj_df)
Iguazu1=cbind(Iguazu1, locsj_df)
head(Iguazu1)
ft=ftable(Iguazu1$id)
f<-as.data.frame.table(ft)
f$Var1 <- NULL
f$Var2 <- NULL
subset(f,Freq>0)

###  Iguazu2    ##  83 => (22 J) 
range(Iguazu2$x);range(Iguazu2$y) 
coord.latlong = SpatialPoints(cbind(Iguazu2$x,Iguazu2$y), proj4string = CRS("+proj=longlat +datum=WGS84"))
coord.latlong
coord.UTM <- spTransform(coord.latlong , CRS("+proj=utm +zone=22J +south +datum=WGS84 +units=m +no_defs"))
coord.UTM
coord.latlong.df <- as.data.frame(coord.latlong)
cootd.utm.df <- as.data.frame(coord.UTM)
head(coord.latlong.df)
head(cootd.utm.df)
par(mfrow = c(1, 2))
plot(coord.latlong.df, axes = TRUE, main = "Lat-Long Coordinates", cex.axis = 0.95)
plot(cootd.utm.df, axes = TRUE, main = "UTM Coordinates", col = "red", cex.axis = 0.95)
locsj_matx=cbind(coordinates(coord.UTM ), coordinates(coord.latlong))
locsj_df<- as.data.frame(locsj_matx)
head(locsj_df)
head(Iguazu2)
point.names<-c("utm_x","utm_y","long_x","lat_y")
colnames(locsj_df)<-point.names
head(locsj_df)
Iguazu2=cbind(Iguazu2, locsj_df)
head(Iguazu2)
ft=ftable(Iguazu2$id)
f<-as.data.frame.table(ft)
f$Var1 <- NULL
f$Var2 <- NULL
subset(f,Freq>0)

Iguazu=rbind(Iguazu1,Iguazu2)
head(Iguazu)
ft=ftable(Iguazu$id)
f<-as.data.frame.table(ft)
f$Var1 <- NULL
f$Var2 <- NULL
subset(f,Freq>0)

#write.table(Iguazu,file="c:/RWorkDir/jaguardatapaper/Iguazu.txt",row.names = F,quote=F,col.names=T,sep="\t")
#Iguazu <- read.delim(file="c:/RWorkDir/jaguardatapaper/Iguazu.txt")
str(Iguazu)

###################################################################################################################################
                                    ###   AMAZONIA  ###
###############################################################################################################################
 ### Amazonia Mamiraua (Brazil)    ####  Flooded   ########################
Mamiraua =subset(jaguar_df,project_region=='Mamiraua')
head(Mamiraua)
str(Mamiraua )
table(Mamiraua$id)
ft=ftable(Mamiraua$id)
f<-as.data.frame.table(ft)
f$Var1 <- NULL
f$Var2 <- NULL
subset(f,Freq>0)
range(Mamiraua$x);range(Mamiraua$y) 
coord.latlong = SpatialPoints(cbind(Mamiraua$x,Mamiraua$y), proj4string = CRS("+proj=longlat +datum=WGS84"))
coord.latlong

# Transforming coordinate from WGS=84 to UTM Zone = 20 M 
coord.UTM <- spTransform(coord.latlong , CRS("+proj=utm +zone=20M +south +datum=WGS84 +units=m +no_defs"))
coord.UTM
coord.latlong.df <- as.data.frame(coord.latlong)
cootd.utm.df <- as.data.frame(coord.UTM)
head(coord.latlong.df)
head(cootd.utm.df)
par(mfrow = c(1, 2))
plot(coord.latlong.df, axes = TRUE, main = "Lat-Long Coordinates", cex.axis = 0.95)
plot(cootd.utm.df, axes = TRUE, main = "UTM Coordinates", col = "red", cex.axis = 0.95)
locsj_matx=cbind(coordinates(coord.UTM ), coordinates(coord.latlong))
locsj_df<- as.data.frame(locsj_matx)
head(locsj_df)
head(Mamiraua)
point.names<-c("utm_x","utm_y","long_x","lat_y")
colnames(locsj_df)<-point.names
head(locsj_df)

Mamiraua=cbind(Mamiraua, locsj_df)
#write.table(Mamiraua,file="c:/RWorkDir/jaguardatapaper/Mamiraua.txt",row.names = F,quote=F,col.names=T,sep="\t")
#Mamiraua <- read.delim(file="c:/RWorkDir/jaguardatapaper/Mamiraua.txt")
head(Mamiraua); str(Mamiraua)

##########################################   Dry Amazonia, PA, translocated ################################################### 
#IOP Para Amazonia   ### Translocated
iopPA=subset(jaguar_df,id=='24')
head(iopPA)
str(iopPA )
table(iopPA$id)
ft=ftable(iopPA$id)
f<-as.data.frame.table(ft)
f$Var1 <- NULL
f$Var2 <- NULL
subset(f,Freq>0)
range(iopPA$x);range(iopPA$y) 
coord.latlong = SpatialPoints(cbind(iopPA$x,iopPA$y), proj4string = CRS("+proj=longlat +datum=WGS84"))
coord.latlong

# Transforming coordinate from WGS=84 to UTM Zone = 22 M 
coord.UTM <- spTransform(coord.latlong , CRS("+proj=utm +zone=22M +south +datum=WGS84 +units=m +no_defs"))
coord.UTM
coord.latlong.df <- as.data.frame(coord.latlong)
cootd.utm.df <- as.data.frame(coord.UTM)
head(coord.latlong.df)
head(cootd.utm.df)
par(mfrow = c(1, 2))
plot(coord.latlong.df, axes = TRUE, main = "Lat-Long Coordinates", cex.axis = 0.95)
plot(cootd.utm.df, axes = TRUE, main = "UTM Coordinates", col = "red", cex.axis = 0.95)
locsj_matx=cbind(coordinates(coord.UTM ), coordinates(coord.latlong))
locsj_df<- as.data.frame(locsj_matx)
head(locsj_df)
head(iopPA)
point.names<-c("utm_x","utm_y","long_x","lat_y")
colnames(locsj_df)<-point.names
head(locsj_df)

iopPA=cbind(iopPA, locsj_df)
#write.table(iopPA,file="c:/RWorkDir/jaguardatapaper/iopPA.txt",row.names = F,quote=F,col.names=T,sep="\t")
#iopPA <- read.delim(file="c:/RWorkDir/jaguardatapaper/iopPA.txt")
head(iopPA); str(iopPA)

#####  Greater Lacandona, Mexico  #######################################################################################

Lacandona =subset(jaguar_df,project_region=='Greater Lacandona')
head(Lacandona)
str(Lacandona )
table(Lacandona$id)
ft=ftable(Lacandona$id)
f<-as.data.frame.table(ft)
f$Var1 <- NULL
f$Var2 <- NULL
subset(f,Freq>0)
range(Lacandona$x);range(Lacandona$y) 
coord.latlong = SpatialPoints(cbind(Lacandona$x,Lacandona$y), proj4string = CRS("+proj=longlat +datum=WGS84"))
coord.latlong

# Transforming coordinate from WGS=84 to UTM Zone = 15 Q
coord.UTM <- spTransform(coord.latlong , CRS("+proj=utm +zone=15Q +south +datum=WGS84 +units=m +no_defs"))
coord.UTM
coord.latlong.df <- as.data.frame(coord.latlong)
cootd.utm.df <- as.data.frame(coord.UTM)
head(coord.latlong.df)
head(cootd.utm.df)
par(mfrow = c(1, 2))
plot(coord.latlong.df, axes = TRUE, main = "Lat-Long Coordinates", cex.axis = 0.95)
plot(cootd.utm.df, axes = TRUE, main = "UTM Coordinates", col = "red", cex.axis = 0.95)
locsj_matx=cbind(coordinates(coord.UTM ), coordinates(coord.latlong))
locsj_df<- as.data.frame(locsj_matx)
head(locsj_df)
head(Lacandona)
point.names<-c("utm_x","utm_y","long_x","lat_y")
colnames(locsj_df)<-point.names
head(locsj_df)
Lacandona=cbind(Lacandona, locsj_df)
#write.table(Lacandona,file="c:/RWorkDir/jaguardatapaper/Lacandona.txt",row.names = F,quote=F,col.names=T,sep="\t")
#Lacandona <- read.delim(file="c:/RWorkDir/jaguardatapaper/Lacandona.txt")
head(Lacandona); str(Lacandona)

#####   Mexico East  ####################################################################################################   
MexEast =subset(jaguar_df,project_region=='Mexico East')
head(MexEast)
str(MexEast )
table(MexEast$id)
ft=ftable(MexEast$id)
f<-as.data.frame.table(ft)
f$Var1 <- NULL
f$Var2 <- NULL
subset(f,Freq>0)
range(MexEast$x);range(MexEast$y) 
coord.latlong = SpatialPoints(cbind(MexEast$x,MexEast$y), proj4string = CRS("+proj=longlat +datum=WGS84"))
coord.latlong

# Transforming coordinate from WGS=84 to UTM Zone = 16 Q
coord.UTM <- spTransform(coord.latlong , CRS("+proj=utm +zone=16Q +south +datum=WGS84 +units=m +no_defs"))
coord.UTM
coord.latlong.df <- as.data.frame(coord.latlong)
cootd.utm.df <- as.data.frame(coord.UTM)
head(coord.latlong.df)
head(cootd.utm.df)
par(mfrow = c(1, 2))
plot(coord.latlong.df, axes = TRUE, main = "Lat-Long Coordinates", cex.axis = 0.95)
plot(cootd.utm.df, axes = TRUE, main = "UTM Coordinates", col = "red", cex.axis = 0.95)
locsj_matx=cbind(coordinates(coord.UTM ), coordinates(coord.latlong))
locsj_df<- as.data.frame(locsj_matx)
head(locsj_df)
head(MexEast)
point.names<-c("utm_x","utm_y","long_x","lat_y")
colnames(locsj_df)<-point.names
head(locsj_df)
MexEast=cbind(MexEast, locsj_df)
write.table(MexEast,file="c:/RWorkDir/jaguardatapaper/MexEast.txt",row.names = F,quote=F,col.names=T,sep="\t")
MexEast <- read.delim(file="c:/RWorkDir/jaguardatapaper/MexEast.txt")
head(MexEast); str(MexEast)

#####   Mexico Sonora  ##################################################################################################
Sonora =subset(jaguar_df,project_region=='Mexico Sonora')
head(Sonora)
str(Sonora )
table(Sonora$id)
ft=ftable(Sonora$id)
f<-as.data.frame.table(ft)
f$Var1 <- NULL
f$Var2 <- NULL
subset(f,Freq>0)
range(Sonora$x);range(Sonora$y) 
coord.latlong = SpatialPoints(cbind(Sonora$x,Sonora$y), proj4string = CRS("+proj=longlat +datum=WGS84"))
coord.latlong

# Transforming coordinate from WGS=84 to UTM Zone = 12 R
coord.UTM <- spTransform(coord.latlong , CRS("+proj=utm +zone=12R +south +datum=WGS84 +units=m +no_defs"))
coord.UTM
coord.latlong.df <- as.data.frame(coord.latlong)
cootd.utm.df <- as.data.frame(coord.UTM)
head(coord.latlong.df)
head(cootd.utm.df)
par(mfrow = c(1, 2))
plot(coord.latlong.df, axes = TRUE, main = "Lat-Long Coordinates", cex.axis = 0.95)
plot(cootd.utm.df, axes = TRUE, main = "UTM Coordinates", col = "red", cex.axis = 0.95)
locsj_matx=cbind(coordinates(coord.UTM ), coordinates(coord.latlong))
locsj_df<- as.data.frame(locsj_matx)
head(locsj_df)
head(Sonora)
point.names<-c("utm_x","utm_y","long_x","lat_y")
colnames(locsj_df)<-point.names
head(locsj_df)
Sonora=cbind(Sonora, locsj_df)
#write.table(Sonora,file="c:/RWorkDir/jaguardatapaper/Sonora.txt",row.names = F,quote=F,col.names=T,sep="\t")
#Sonora <- read.delim(file="c:/RWorkDir/jaguardatapaper/Sonora.txt")
head(Sonora); str(Sonora)

#############  ###  Mexico ###

Mex=rbind(Lacandona,MexEast,Sonora)
head(Mex)
ft=ftable(Mex$id)
f<-as.data.frame.table(ft)
f$Var1 <- NULL
f$Var2 <- NULL
subset(f,Freq>0)
head(Mex)
#write.table(Mex,file="c:/RWorkDir/jaguardatapaper/Mex.txt",row.names = F,quote=F,col.names=T,sep="\t")
#Mex <- read.delim(file="c:/RWorkDir/jaguardatapaper/Mex.txt")

#######################  Jaguar Dataframe with UTMs ########################################################################
head(AFW);head(Caatinga);head(Cerrado);head(CRica);head(Pantanal);str(Drych);head(Hch);head(FPy);head(Iguazu);
head(Mamiraua);head(iopPA);head(Lacandona);head(MexEast);head(Sonora)

jaguar=rbind(AFW,Caatinga,Cerrado,CRica,Pantanal,Drych,Hch,FPy,Iguazu,Mamiraua,iopPA,Lacandona,MexEast,Sonora)
head(jaguar); str(jaguar)
jaguar_df <- jaguar

ft=ftable(jaguar$id)
f<-as.data.frame.table(ft)
f$Var1 <- NULL
f$Var2 <- NULL
subset(f,Freq>0)

#write.table(jaguar,file="c:/RWorkDir/jaguardatapaper/jaguar.txt",row.names = F,quote=F,col.names=T,sep="\t")
#write.table(jaguar_df,file="c:/RWorkDir/jaguardatapaper/jaguar_df.txt",row.names = F,quote=F,col.names=T,sep="\t")
#jaguar_df <- read.delim(file="c:/RWorkDir/jaguardatapaper/jaguar_df.txt")

################################################################################################################
################   ADEHABITAT LT   ###  Update the mov.traj object  ####  
## adehabitatLT
# Transforms in ltraj object

# Transforms in ltraj object
coords <- data.frame(jaguar$utm_x, jaguar$utm_y)
jaguar.traj <- as.ltraj(xy = coords, date=jaguar$timestamp.posix, 
                     id=jaguar$id, 
                     burst=jaguar$id, 
                     infolocs = jaguar)
jaguar.traj.df <- ld(jaguar.traj)
head(jaguar.traj.df) 
plot(jaguar.traj)
hist(jaguar.traj.df$dist)
max(jaguar.traj.df$dist,na.rm=TRUE) 
range(jaguar.traj.df$dist,na.rm=TRUE) 

jaguar_df<- jaguar.traj.df
jaguar_df$x <- NULL
jaguar_df$y <- NULL
jaguar_df$dx <- NULL
jaguar_df$dy <- NULL
jaguar_df$date.1 <- NULL
jaguar_df$id.1 <- NULL
jaguar_df$burst <- NULL
jaguar_df <- rename(jaguar_df, c("x.1"="x", "y.1"="y"))

jaguar_df <- jaguar_df%>%select(Event_ID,x,y,date, everything())

#write.table(jaguar_df,file="c:/RWorkDir/jaguardatapaper/jaguar_df.txt",row.names = F,quote=F,col.names=T,sep="\t")
#jaguar_df <- read.delim(file="c:/RWorkDir/jaguardatapaper/jaguar_df.txt")
AFW1=subset(jaguar_df,project_region=='Atlantic Forest W1'); head(AFW1) ### 1)
AFW2=subset(jaguar_df,project_region=='Atlantic Forest W2'); head(AFW2) ### 2)
Caatinga=subset(jaguar_df,project_region=='Caatinga'); head(Caatinga) ### 3)
Cerrado1=subset(jaguar_df,project_region=='Cerrado1'); head(Cerrado1) ### 4)
Cerrado2=subset(jaguar_df,project_region=='Cerrado2'); head(Cerrado2) ### C
CRica=subset(jaguar_df,project_region=='Costa Rica'); head(CRica) ### D
Pantanal =subset(jaguar_df,project_bioveg=='Pantanal'); head(Pantanal) ### E
Drych=subset(jaguar_df,project_region=='Dry chaco');  head(Drych) ### F
Hch=subset(jaguar_df,project_region=='Humid chaco');  head(Hch) ### G
FPy=subset(jaguar_df,project_region=='Forest Paraguay');  head(Hch) ### H
Iguazu=subset(jaguar_df,project_region=='Iguazu');  head(Iguazu) ### I
Mamiraua =subset(jaguar_df,project_region=='Mamiraua');  head(Mamiraua) ### J
iopPA=subset(jaguar_df,id=='24');  head(iopPA) ### K
Lacandona =subset(jaguar_df,project_region=='Greater Lacandona');  head(Lacandona) ### L
MexEast =subset(jaguar_df,project_region=='Mexico East');  head(MexEast) ### M
Sonora =subset(jaguar_df,project_region=='Mexico Sonora');  head(Sonora) ### N


# Organize data   
# Time Posix class (2014-12-05 00:00:00 )
jaguar_df$timestamp.posix <- as.character(jaguar_df$timestamp.posix)
jaguar_df$timestamp.posix<- as.POSIXct(jaguar_df$timestamp.posix, format ="%Y-%m-%d %H:%M:%S", tz = 'GMT') 
jaguar_df$date <- as.character(jaguar_df$date)
jaguar_df$date <- as.POSIXct(jaguar_df$date, format ="%Y-%m-%d %H:%M:%S", tz = 'GMT') 
#jaguar_df$id <- as.factor(jaguar_df$individual.local.identifier..ID.)
#jaguar_df$season <- as.factor(jaguar_df$season) #only dry, flood, inter seasons
jaguar_df$dt <- as.numeric(jaguar_df$dt)
jaguar_df$week <- as.numeric(strftime(as.POSIXlt(jaguar_df$timestamp.posix),format="%W"))
jaguar_df$day <- as.numeric(strftime(as.POSIXlt(jaguar_df$timestamp.posix),format="%j"))
jaguar_df$year <- as.numeric(strftime(as.POSIXlt(jaguar_df$timestamp.posix),format="%Y"))
str(jaguar_df)


###############################################################################################################
###          Update the MoveStack object (move.data)to contain only the cleaned points
move.data <- move(x = jaguar_df$x, y = jaguar_df$y, 
                  time = jaguar_df$date, 
                  proj = CRS('+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'),
                  data = jaguar_df, animal = jaguar_df$id, sensor = 'GPS')
		
#"+proj=aea +lat_1=-5 +lat_2=-42 +lat_0=-32 +lon_0=-60 +x_0=0 +y_0=0 +ellps=aust_SA +towgs84=-57,1,-41,0,0,0,0 +units=m +no_defs"		
################################################################################################################

##############################################################################################################################
# now do for the residence time (method of Barraquand & Benhamou 2008)  ####  Require to be in UTM
# choose 1km radius and 1 hours as maximum time threshold that the animal is allowed to spend outside the patch 
# before we consider that the animal actually left the patch

###          ltraj object from Cerrado  ####  

#mov.traj<-move2ade(move.data)
# Transforms in ltraj object
coords <- data.frame(jaguar_df$x, jaguar_df$y)
mov.traj <- as.ltraj(xy = coords, date=jaguar_df$timestamp.posix, 
                     id=jaguar_df$individual.local.identifier..ID., 
                     burst=jaguar_df$individual.local.identifier..ID., 
                     infolocs = jaguar_df)
mov.traj.df <- ld(mov.traj)
#mov.traj
#plot(mov.traj)
hist(mov.traj.df$dist)


  cerrado1<- mov.traj[1]
  res1 <- residenceTime(cerrado1, radius = 1000, maxt=24, units="hour")
  plot(res1)

# maybe lets increase to 2km and 12 hours
  res4049 <- residenceTime(elk4049, radius = 2000, maxt=12, units="hour")
  plot(res4049)

# add this to the infolocs slot

  res4049 <- residenceTime(elk4049, radius = 2000, maxt=12, addinfo = TRUE, units="hour")
  res4049

# repeat for the other individuals
jtraj1<- mov.traj[1]
  res1 <- residenceTime(jtraj1, radius = 1, maxt=24, units="hour")
  plot(res1)



###############################################################################################################
# Prepare data in the ctmm format
mov.tel <- as.telemetry(move.data)
str(mov.tel)
###############################################################################################################

# it can come handy to add an identifier variable, identifying the first location of each individual - let us call this 'firstRec'
 foo <- which(jaguar_df$individual.local.identifier..ID.[1:(nrow(jaguar_df)-1)] != jaguar_df$individual.local.identifier..ID.[2:(nrow(jaguar_df))])
 jaguar_df$firstRec <- rep(0,nrow(jaguar_df))
 jaguar_df$firstRec[foo+1] <- 1
 jaguar_df$firstRec[1] <- 1

 # let us check if this is correct
  length(unique(jaguar_df$individual.local.identifier..ID.)) # count N individuals
  sum(jaguar_df$firstRec)	# sum = 117, looks fine!
  jaguar_df[sort(c(foo-1,foo,foo+1)),c('individual.local.identifier..ID.','date','firstRec')] # first records seem correctly identified
  rm(foo) 	# keep workspace clean

  # time between steps
# ?difftime
# syntax: difftime(then,now,units="secs")
  foo <- difftime(jaguar_df$date[2:nrow(jaguar_df)], jaguar_df$date[1:(nrow(jaguar_df)-1)], units = "secs")
  foo <- c(NA, foo)
  summary(as.numeric(foo))
  foo <- ifelse(jaguar_df$firstRec == 1, NA, foo)
  summary(as.numeric(foo))
  jaguar_df$dt <- foo
  rm(foo)  ### clean foo
  
# investigate distribution of step length and time lags between steps
  hist(jaguar_df$dt); summary(jaguar_df$dt)	
# think about the sampling regime -- 7200s = 2h; 3600s = 1h; 900s = 15min; likely a resampling and/or imputation will be required for analyses
  hist(jaguar_df[jaguar_df$dt <86400,'dt'])		
  hist(jaguar_df[jaguar_df$dt < 172800 & jaguar_df$dt > 1200,'dt'])

test <- filter(jaguar_df, dt > 100000) ###  2489 obs. for dt > 100000 approx 1 day
hist(jaguar_df[jaguar_df$dt > 100000& jaguar_df$dt <2600000,'dt']); summary(test$dt) ### approx between 1 day and a month
testb <- filter(jaguar_df, dt > 2600000); summary(testb$dt) ### only 54 obs. with more than a month
hist(jaguar_df[jaguar_df$dt > 2600000,'dt'])
testb2 <- testb %>% group_by(id)%>%nest(id, dt)%>% print(n = Inf)  ### we can see the frequency which animals had highest timelags
rm(test) ; rm(testb);rm(testb2)  ### clean tests

##############################################################################################################
#' ## Creating a track in amt (commom to both RSF and SSF)  
##############################################################################################################
#' Before we can use the amt package to calculate step lengths, turn angles, and bearings
#' for fisher data, we need to add a class (track) to the data. Then, we can summarize 
#' the data by individual, month, etc. First, create a track using utms and the timestamp.
#' 
#' If we have a data set with locations in utms, we could use:
#trk.temp <- make_track(fisher.dat, .x=utm.easting, .y=utm.northing, .t=timestamp, id = individual_local.identifier)
#trk.temp

#' Note: we did not need to explicitly specify x, y and t (but probably good to do so).
#' This would also have worked
#' trk <- make_track(fisher.dat, utm.easting, utm.northing, timestamp, id = local_identifier)

#' We can also use lat, long, which will allow us to determine
#' time of day 
trk <- mk_track(jaguar_df,.x=x, .y=y, .t=date, Event_ID=Event_ID, id=id, sex=sex,
age=age, weight=weight, status=status,dist=dist,project_region=project_region, crs = CRS("+init=epsg:4326"))    ### season=season, speed_kmh=speed_kmh,
trk <- trk %>% arrange(id)
trk

trk <- transform_coords(trk, sp::CRS("+proj=aea +lat_1=-5 +lat_2=-42 +lat_0=-32 +lon_0=-60 +x_0=0 +y_0=0 +ellps=aust_SA +units=m +no_defs"))
trk

trk <- trk %>% time_of_day()  
trk.class<-class(trk)

#' ## Movement Characteristics
nesttrk<-trk%>%nest(-id)
nesttrk      #%>% print(n = 117)

#' Each row contains data from an individual.  For example, we can access data
#' from the first individual using:
nesttrk$data[[1]]

#' We could calculate movement characteristics by individual using:
temp<-direction_rel(nesttrk$data[[1]])
head(temp)

#' or:
temp<-trk %>% filter(id=="1") %>% direction_rel
head(temp)

#' Or, we can add a columns to each nested column of data using purrr::map
trk<-trk %>% nest(-id) %>% 
  mutate(dir_abs = map(data, direction_abs,full_circle=TRUE, zero="N"), 
         dir_rel = map(data, direction_rel), 
         sl = map(data, step_lengths),
         nsd_=map(data, nsd))%>%unnest()
		 

trk
#' Now, calculate month, year, hour, week of each observation and append these to the dataset
#' Unlike the movement charactersitics, these calculations can be done all at once, 
#' since they do not utilize successive observations (like step lengths and turn angles do).
trk<-trk%>% 
  mutate(
    week=week(t_),
    month = month(t_, label=TRUE), 
    year=year(t_),
    hour = hour(t_)
  )

#' Now, we need to again tell R that this is a track (rather 
#' than just a data frame)
class(trk)
class(trk)<-trk.class

#' Lets take a look at what we created
trk

a <-trk$sl
b <-trk$dist
c <-a-b
head(c)
max(c, na.rm=TRUE)

trk<-trk%>% mutate(difdist=(trk$sl-trk$dist))

tapply(abs(trk$difdist),list(trk$id),na.rm=TRUE,max)
tapply(trk$difdist,list(trk$id),na.rm=TRUE,mean)
tapply(trk$difdist,list(trk$project_region,trk$id),na.rm=TRUE,max)
proj_dif=tapply(abs(trk$difdist),list(trk$id),na.rm=TRUE,max)
plot(proj_dif)

xyplot(difdist~id, data = trk, groups = id)
xyplot(difdist~id, data = trk, groups = project_region, auto.key=list(columns  = 2))
xyplot(difdist~project_region, data = trk, groups = project_region, auto.key=list(columns  = 2))


table(trk1$id)
ft=ftable(AFW1$id)
f<-as.data.frame.table(ft)
f$Var1 <- NULL
f$Var2 <- NULL
subset(f,Freq>0)


ftable(tapply(Frequency,list(ScenarioC,BiodUNEP,Biome),sum))

(difdist<-trk %>% nest(-id,-project_region) %>% mutate(sr = map(data, summarize_sampling_rate, time_unit = "hour")) %>%
    select(id,sex,age,weight,status,sr) %>% unnest)


nesttrk<-trk%>%nest(-id,-project_region)%>% max(trk$difdist, na.rm=TRUE)
nesttrk      #%>% print(n = 117)

########################################################################################################################################
#######################  Jaguar Dataframe with UTMs ########################################################################
head(AFW1);head(AFW2),head(Caatinga);head(Cerrado1);head(Cerrado2);head(CRica);head(Pantanal);head(Drych1);str(Drych2);head(Hch);
head(FPy);head(Iguazu1);head(Iguazu2);head(Mamiraua);head(iopPA);head(Lacandona);head(MexEast);head(Sonora)

jaguar=rbind(AFW1,AFW2,Caatinga,Cerrado1,Cerrado2,CRica,Pantanal,Drych1,Drych2,Hch,FPy,Iguazu1,Iguazu2,Mamiraua,iopPA,Lacandona,MexEast,Sonora)
#head(jaguar); str(jaguar)
#jaguar_df <- jaguar

#write.table(jaguar,file="c:/RWorkDir/jaguardatapaper/jaguar.txt",row.names = F,quote=F,col.names=T,sep="\t")
#write.table(jaguar_df,file="c:/RWorkDir/jaguardatapaper/jaguar_df.txt",row.names = F,quote=F,col.names=T,sep="\t")
#jaguar_df <- read.delim(file="c:/RWorkDir/jaguardatapaper/jaguar_df.txt")


##############################################################################################################
#' ## Creating a track in amt (commom to both RSF and SSF)  
##############################################################################################################
### (1) ############################################## 
### #      Atlantic Forest W1 
########################################################
#AFW1 <- read.delim(file="c:/RWorkDir/jaguardatapaper/AFW1.txt")  # <-- if use read table it will be required to adjust variables again!!!
#age <- as.numeric(levels(AFW1$age))[AFW1$age]
#weight <- as.numeric(levels(AFW1$weight))[AFW1$weight]
#AFW1$id <- as.factor(AFW1$individual.local.identifier..ID.)  
#AFW1$age <- age   ### converted to number
#AFW1$weight <- weight   ### converted to number
#AFW1 <- read.delim(file="c:/RWorkDir/jaguardatapaper/AFW1.txt")
#AFW1$timestamp.posix <- as.character(AFW1$timestamp.posix)
#AFW1$timestamp.posix<- as.POSIXct(AFW1$timestamp.posix, format ="%Y-%m-%d %H:%M:%S", tz = 'GMT') 
#AFW1$date <- as.character(AFW1$date)
#AFW1$date <- as.POSIXct(AFW1$date, format ="%Y-%m-%d %H:%M:%S", tz = 'GMT') 
head(AFW1);str(AFW1)

#" Checking order again
jaguar_ord <- AFW1[order(AFW1$id,AFW1$date),]
all.equal(AFW1,jaguar_ord)

trk <- mk_track(AFW1,.x=utm_x, .y=utm_y, .t=date, id=id,Event_ID=Event_ID, sex=sex, age=age, weight=weight, 
status=status, period=period,long_x=long_x, lat_y=lat_y, crs = CRS("+proj=utm +zone=22K +south +datum=WGS84 +units=m +no_defs"))   ###  WGS84
trk <- trk %>% arrange(id)
trk
#' ## Movement Characteristics
nesttrk<-trk%>%nest(-id)
nesttrk      #%>% print(n = 117)
#' We can add a columns to each nested column of data using purrr::map
trk<-trk %>% nest(-id) %>% 
  mutate(dir_abs = map(data, direction_abs,full_circle=TRUE, zero="N"), 
         dir_rel = map(data, direction_rel), 
         sl = map(data, step_lengths),
         nsd_=map(data, nsd))%>%unnest()
trk.class<-class(trk)
trk

#' Calculate month, year, hour, week of each observation and append these to the dataset
#' Unlike the movement charactersitics, these calculations can be done all at once, 
#' since they do not utilize successive observations (like step lengths and turn angles do).
trk<-trk%>% 
  mutate(
    week=week(t_),
    month = month(t_, label=TRUE), 
    year=year(t_),
    hour = hour(t_)
  )

#' Now, we need to again tell R that this is a track (rather than just a data frame)
class(trk)
class(trk)<-trk.class

#' Lets take a look at what we created
trk
AFW1trk <-trk
########################################

### (2) ############################################## 
### #      Atlantic Forest W2 
########################################################
#AFW2 <- read.delim(file="c:/RWorkDir/jaguardatapaper/AFW2.txt")  # <-- if use read table it will be required to adjust variables again!!!
#age <- as.numeric(levels(AFW2$age))[AFW2$age]
#weight <- as.numeric(levels(AFW2$weight))[AFW2$weight]
#AFW2$id <- as.factor(AFW2$individual.local.identifier..ID.)  
#AFW2$age <- age   ### converted to number
#AFW2$weight <- weight   ### converted to number
#AFW2 <- read.delim(file="c:/RWorkDir/jaguardatapaper/AFW2.txt")
#AFW2$timestamp.posix <- as.character(AFW2$timestamp.posix)
#AFW2$timestamp.posix<- as.POSIXct(AFW2$timestamp.posix, format ="%Y-%m-%d %H:%M:%S", tz = 'GMT') 
#AFW2$date <- as.character(AFW2$date)
#AFW2$date <- as.POSIXct(AFW2$date, format ="%Y-%m-%d %H:%M:%S", tz = 'GMT') 
head(AFW2);str(AFW2)

#" Checking order again
jaguar_ord <- AFW2[order(AFW2$id,AFW2$date),]
all.equal(AFW2,jaguar_ord)

trk <- mk_track(AFW2,.x=utm_x, .y=utm_y, .t=date, id=id,Event_ID=Event_ID, sex=sex, age=age, weight=weight, 
status=status, period=period,long_x=long_x, lat_y=lat_y, crs = CRS("+proj=utm +zone=22K +south +datum=WGS84 +units=m +no_defs"))   ###  WGS84
trk <- trk %>% arrange(id)
trk
#' ## Movement Characteristics
nesttrk<-trk%>%nest(-id)
nesttrk      #%>% print(n = 117)
#' We can add a columns to each nested column of data using purrr::map
trk<-trk %>% nest(-id) %>% 
  mutate(dir_abs = map(data, direction_abs,full_circle=TRUE, zero="N"), 
         dir_rel = map(data, direction_rel), 
         sl = map(data, step_lengths),
         nsd_=map(data, nsd))%>%unnest()
trk.class<-class(trk)
trk

#' Calculate month, year, hour, week of each observation and append these to the dataset
#' Unlike the movement charactersitics, these calculations can be done all at once, 
#' since they do not utilize successive observations (like step lengths and turn angles do).
trk<-trk%>% 
  mutate(
    week=week(t_),
    month = month(t_, label=TRUE), 
    year=year(t_),
    hour = hour(t_)
  )

#' Now, we need to again tell R that this is a track (rather than just a data frame)
class(trk)
class(trk)<-trk.class

#' Lets take a look at what we created
trk
AFW2trk <-trk
#################################################################################
      ###  Atlantic Forest West =>  AFWtrk   ####
   AFWtrk=rbind(AFW1trk,AFW2trk); AFWtrk
################################################################################

