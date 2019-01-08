#############################################################################################################################
                            ############# Jaguar Data Preparation ##############    ### jaguar data from Morato et al. 2018    
#############################################################################################################################						
### Alan E. de Barros/ adapted from Bernardo Niebuhr data preparation, Luca Borger lectures and amt John Fieberg's scripts (from MoveBank) 																	  
rm(list= ls())                                                  ### For a fresh start
setwd("C:/RWorkDir/jaguardatapaper")                            ### Set directory
###########################################################################################################         																					
## Load packages
if(!require(install.load)) install.packages('install.load'); library(install.load)
install.load::install_load("ggmap","maptools",'move',"circular","RCurl","dplyr","readr","caTools","adehabitatLT","rgl",
"lubridate","raster","amt","tibble","knitr","leaflet","ezknitr","lattice","rgdal","sp")

### Load and adjust the data and create a dataframe object
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
####################################################

### Organize data
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
##################################################################################################
# I included a column to represent the local timezone (already with the - signal) to them multiply the timestamp and get the difference:
mov.data.org$local_time <- mov.data.org$timestamp.posix + mov.data.org$timezone*60*60
mov.data.org$timestamp.posix.GMT <- mov.data.org$timestamp.posix
mov.data.org$timestamp.posix <- mov.data.org$local_time ### If we do that all the (timestamp.posix)'s calculations will be based on local time
str(mov.data.org)

################################################################
## adehabitatLT
################
# Transforms in ltraj object
coords <- data.frame(mov.data.org$location.long, mov.data.org$location.lat)
mov.traj <- as.ltraj(xy = coords, date=mov.data.org$timestamp.posix, 
                     id=mov.data.org$individual.local.identifier..ID., 
                     burst=mov.data.org$individual.local.identifier..ID., 
                     infolocs = mov.data.org[,-c(3:4, ncol(mov.data.org))])
mov.traj.df <- ld(mov.traj)
#mov.traj
#plot(mov.traj)
###############
## move
###############
# Organize data as a move package format
move.data <- move(x = mov.traj.df$x, y = mov.traj.df$y, 
                  time = mov.traj.df$date, 
                  proj = CRS('+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'),
                  data = mov.traj.df, animal = mov.traj.df$id, sensor = 'GPS')
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

jaguar_df$week <- as.numeric(strftime(as.POSIXlt(jaguar_df$timestamp.posix),format="%W"))
jaguar_df$day <- as.numeric(strftime(as.POSIXlt(jaguar_df$timestamp.posix),format="%j"))
jaguar_df$year <- as.numeric(strftime(as.POSIXlt(jaguar_df$timestamp.posix),format="%Y"))
jaguar_df$hour <- as.numeric(strftime(as.POSIXlt(jaguar_df$timestamp.posix),format="%H"))
jaguar_df$min <- as.numeric(strftime(as.POSIXlt(jaguar_df$timestamp.posix),format="%M"))
jaguar_df$time <- jaguar_df$hour + (jaguar_df$min)/60

# day  7 to 16  and  night 19 to 4 and  sun riseset 5,6,17,18
#jaguar_df$month=as.numeric(substr(jaguar_df$date,6,7))
jaguar_df$period=ifelse(jaguar_df$hour==7,"day",ifelse(jaguar_df$hour==8,"day",ifelse(jaguar_df$hour==9,"day",
ifelse(jaguar_df$hour==10,"day", ifelse(jaguar_df$hour==11,"day", ifelse(jaguar_df$hour==12,"day", ifelse(jaguar_df$hour==13,"day",
ifelse(jaguar_df$hour==14,"day",ifelse(jaguar_df$hour==15,"day", ifelse(jaguar_df$hour==16,"day",
ifelse(jaguar_df$hour==19,"night",ifelse(jaguar_df$hour==20,"night",ifelse(jaguar_df$hour==21,"night",ifelse(jaguar_df$hour==22,"night",
ifelse(jaguar_df$hour==23,"night",ifelse(jaguar_df$hour==0,"night",ifelse(jaguar_df$hour==1,"night",ifelse(jaguar_df$hour==2,"night",
ifelse(jaguar_df$hour==3,"night",ifelse(jaguar_df$hour==4,"night","riseset"))))))))))))))))))))


head(jaguar_df)
str(jaguar_df)

#########################################################################################################							   							   									  
#' ######################## More Data checking and cleaning ##########################################
###############################################################
#' Delete observations where missing lat or long or a timestamp. (There are no missing observations but it is a good practice)
ind<-complete.cases(jaguar_df[,c("y","x","date")])
jaguar_df<-jaguar_df[ind==TRUE,]

#" Check order (Dataset had been previously ordered)
jaguar_ord <- jaguar_df[order(jaguar_df$id,jaguar_df$date),]
#jaguar_ord <- jaguar_df[order(jaguar_df$Event_ID),] ### Strangely when ordered by Event_ID results are not true
all.equal(jaguar_df,jaguar_ord)

#' Check for duplicated observations (ones with same lat, long, timestamp,
#'  and individual identifier). No duplicate observations in this data set
ind2<-jaguar_df %>% select("date","x","y","id") %>% duplicated
sum(ind2) # no duplicates
jaguar_df<-jaguar_df[ind2!=TRUE,]
### or duplicatedLocs <- which(jaguar_df$date[1:(nrow(jaguar_df)-1)] == jaguar_df$date[2:(nrow(jaguar_df))])

###########   Clean suspectly close points!!!  Above 1200 secs or 20 min minimal interval ####################
### Further removals can be done later if required
excludes <- filter(jaguar_df, dt < 1200)
### see excludeds based on id
table(select(excludes, id))             
removed<- anti_join(jaguar_df, excludes)   ### test if were all removed
filter(removed, dt < 1200)
jaguar_df <- removed 
str(jaguar_df)

## Cleaning up columns which will be in excess due to repetition of analysis ########################################################################
jaguar_df$dx <- NULL
jaguar_df$dy <- NULL
jaguar_df$dist <- NULL
#jaguar_df$dt <- NULL
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
#jaguar_df$timezone <- NULL
head(jaguar_df) 
str(jaguar_df)

#write.table(jaguar_df,file="c:/RWorkDir/jaguardatapaper/jaguar_df.txt",row.names = F,quote=F,col.names=T,sep="\t")
#jaguar_df <- read.delim(file="c:/RWorkDir/jaguardatapaper/jaguar_df.txt")

###############################################################################################################################
###  Add UTMs  #################
###############################################################################################################################
### Looking at project regions and grouping then when they occur within the same UTM zone

table(jaguar_df$project_region)

###############################################################################################################################
                               ###  ATLANTIC FOREST WEST  ###    (Project regions 1 and 2)
###############################################################################################################################
### # 1) Atlantic Forest W1 
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
coord.utm.df <- as.data.frame(coord.UTM)
head(coord.latlong.df)
head(coord.utm.df)
par(mfrow = c(1, 2))
plot(coord.latlong.df, axes = TRUE, main = "Lat-Long Coordinates", cex.axis = 0.95)
plot(coord.utm.df, axes = TRUE, main = "UTM Coordinates", col = "red", cex.axis = 0.95)
locsj_matx=cbind(coordinates(coord.UTM), coordinates(coord.latlong))
locsj_df<- as.data.frame(locsj_matx)
head(locsj_df)
head(AFW1)
point.names<-c("utm_x","utm_y","long_x","lat_y")
colnames(locsj_df)<-point.names
head(locsj_df)

AFW1=cbind(AFW1, locsj_df)
head(AFW1); str(AFW1)
#write.table(AFW1,file="c:/RWorkDir/jaguardatapaper/AFW1.txt",row.names = F,quote=F,col.names=T,sep="\t")
#AFW1 <- read.delim(file="c:/RWorkDir/jaguardatapaper/AFW1.txt")
head(AFW1); str(AFW1)

#############################################
### # 2) Atlantic Forest W2
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
coord.utm.df <- as.data.frame(coord.UTM)
head(coord.latlong.df)
head(coord.utm.df)
par(mfrow = c(1, 2))
plot(coord.latlong.df, axes = TRUE, main = "Lat-Long Coordinates", cex.axis = 0.95)
plot(coord.utm.df, axes = TRUE, main = "UTM Coordinates", col = "red", cex.axis = 0.95)
locsj_matx=cbind(coordinates(coord.UTM), coordinates(coord.latlong))
locsj_df<- as.data.frame(locsj_matx)
head(locsj_df)
head(AFW2)
point.names<-c("utm_x","utm_y","long_x","lat_y")
colnames(locsj_df)<-point.names
head(locsj_df)

AFW2=cbind(AFW2, locsj_df)
head(AFW2)

#write.table(AFW2,file="c:/RWorkDir/jaguardatapaper/AFW2.txt",row.names = F,quote=F,col.names=T,sep="\t")
#AFW1 <- read.delim(file="c:/RWorkDir/jaguardatapaper/AFW1.txt")
head(AFW2); str(AFW2)


AFW=rbind(AFW1,AFW2)
head(AFW)
ft=ftable(AFW$id)
f<-as.data.frame.table(ft)
f$Var1 <- NULL
f$Var2 <- NULL
subset(f,Freq>0)

write.table(AFW,file="c:/RWorkDir/jaguardatapaper/AFW.txt",row.names = F,quote=F,col.names=T,sep="\t")
#AFW <- read.delim(file="c:/RWorkDir/jaguardatapaper/AFW.txt")
head(AFW); str(AFW)

###############################################################################################################################
                                    ###  CAATINGA  ###  
###############################################################################################################################
### # 3) Caatinga

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
coord.utm.df <- as.data.frame(coord.UTM)
head(coord.latlong.df)
head(coord.utm.df)
par(mfrow = c(1, 2))
plot(coord.latlong.df, axes = TRUE, main = "Lat-Long Coordinates", cex.axis = 0.95)
plot(coord.utm.df, axes = TRUE, main = "UTM Coordinates", col = "red", cex.axis = 0.95)
locsj_matx=cbind(coordinates(coord.UTM ), coordinates(coord.latlong))
locsj_df<- as.data.frame(locsj_matx)
head(locsj_df)
head(Caatinga)
point.names<-c("utm_x","utm_y","long_x","lat_y")
colnames(locsj_df)<-point.names
head(locsj_df)

Caatinga=cbind(Caatinga, locsj_df)
head(Caatinga); str(Caatinga)
write.table(Caatinga,file="c:/RWorkDir/jaguardatapaper/Caatinga.txt",row.names = F,quote=F,col.names=T,sep="\t")
#Caatinga <- read.delim(file="c:/RWorkDir/jaguardatapaper/Caatinga.txt")

###############################################################################################################################
                                    ###  CERRADO   ###    2 Projects (4 and 5)
###############################################################################################################################
### # 4)  Cerrado1

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
coord.utm.df <- as.data.frame(coord.UTM)
head(coord.latlong.df)
head(coord.utm.df)
par(mfrow = c(1, 2))
plot(coord.latlong.df, axes = TRUE, main = "Lat-Long Coordinates", cex.axis = 0.95)
plot(coord.utm.df, axes = TRUE, main = "UTM Coordinates", col = "red", cex.axis = 0.95)
locsj_matx=cbind(coordinates(coord.UTM ), coordinates(coord.latlong))
locsj_df<- as.data.frame(locsj_matx)
head(locsj_df)
head(Cerrado1)
point.names<-c("utm_x","utm_y","long_x","lat_y")
colnames(locsj_df)<-point.names
head(locsj_df)

Cerrado1=cbind(Cerrado1, locsj_df)
head(Cerrado1)
write.table(Cerrado1,file="c:/RWorkDir/jaguardatapaper/Cerrado1.txt",row.names = F,quote=F,col.names=T,sep="\t")
#Cerrado1 <- read.delim(file="c:/RWorkDir/jaguardatapaper/Cerrado1.txt")
head(Cerrado1); str(Cerrado1)


###############################################################################################################################
### # 5)  Cerrado2
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
coord.utm.df <- as.data.frame(coord.UTM)
head(coord.latlong.df)
head(coord.utm.df)
par(mfrow = c(1, 2))
plot(coord.latlong.df, axes = TRUE, main = "Lat-Long Coordinates", cex.axis = 0.95)
plot(coord.utm.df, axes = TRUE, main = "UTM Coordinates", col = "red", cex.axis = 0.95)
locsj_matx=cbind(coordinates(coord.UTM ), coordinates(coord.latlong))
locsj_df<- as.data.frame(locsj_matx)
head(locsj_df)
head(Cerrado2)
point.names<-c("utm_x","utm_y","long_x","lat_y")
colnames(locsj_df)<-point.names
head(locsj_df)

Cerrado2=cbind(Cerrado2, locsj_df)
head(Cerrado2)
write.table(Cerrado2,file="c:/RWorkDir/jaguardatapaper/Cerrado2.txt",row.names = F,quote=F,col.names=T,sep="\t")
#Cerrado2 <- read.delim(file="c:/RWorkDir/jaguardatapaper/Cerrado2.txt")
head(Cerrado2); str(Cerrado2)


Cerrado=rbind(Cerrado1,Cerrado2)
head(Cerrado)
ft=ftable(Cerrado$id)
f<-as.data.frame.table(ft)
f$Var1 <- NULL
f$Var2 <- NULL
subset(f,Freq>0)

write.table(Cerrado,file="c:/RWorkDir/jaguardatapaper/Cerrado.txt",row.names = F,quote=F,col.names=T,sep="\t")
#Cerrado <- read.delim(file="c:/RWorkDir/jaguardatapaper/Cerrado.txt")

###############################################################################################################################
                                    ###  COSTA RICA  ### 
###############################################################################################################################
### # 6) Costa Rica

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
coord.utm.df <- as.data.frame(coord.UTM)
head(coord.latlong.df)
head(coord.utm.df)
par(mfrow = c(1, 2))
plot(coord.latlong.df, axes = TRUE, main = "Lat-Long Coordinates", cex.axis = 0.95)
plot(coord.utm.df, axes = TRUE, main = "UTM Coordinates", col = "red", cex.axis = 0.95)
locsj_matx=cbind(coordinates(coord.UTM ), coordinates(coord.latlong))
locsj_df<- as.data.frame(locsj_matx)
head(locsj_df)
head(CRica)
point.names<-c("utm_x","utm_y","long_x","lat_y")
colnames(locsj_df)<-point.names
head(locsj_df)

CRica=cbind(CRica, locsj_df)
head(CRica)

write.table(CRica,file="c:/RWorkDir/jaguardatapaper/CRica.txt",row.names = F,quote=F,col.names=T,sep="\t")
#CRica<- read.delim(file="c:/RWorkDir/jaguardatapaper/CRica.txt")

###############################################################################################################################
                                    ###   DRY CHACO  ###   
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

### # 7) Drych1    ## 76,77  Zone 20 K; 16 & 70 (most is in Zone 20 K too, but a little bit in 21K)
range(Drych1$x);range(Drych1$y) 
coord.latlong = SpatialPoints(cbind(Drych1$x,Drych1$y), proj4string = CRS("+proj=longlat +datum=WGS84"))
coord.latlong
coord.UTM <- spTransform(coord.latlong , CRS("+proj=utm +zone=20K +south +datum=WGS84 +units=m +no_defs"))
#coord.UTM <- spTransform(coord.latlong , CRS("+proj=utm +zone=21K +south +datum=WGS84 +units=m +no_defs"))

coord.UTM
coord.latlong.df <- as.data.frame(coord.latlong)
coord.utm.df <- as.data.frame(coord.UTM)
head(coord.latlong.df)
head(coord.utm.df)
par(mfrow = c(1, 2))
plot(coord.latlong.df, axes = TRUE, main = "Lat-Long Coordinates", cex.axis = 0.95)
plot(coord.utm.df, axes = TRUE, main = "UTM Coordinates", col = "red", cex.axis = 0.95)
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
write.table(Drych1,file="c:/RWorkDir/jaguardatapaper/Drych1.txt",row.names = F,quote=F,col.names=T,sep="\t")
#Drych1 <- read.delim(file="c:/RWorkDir/jaguardatapaper/Drych1.txt")
head(Drych1); str(Drych1)



### # 8) Drych2    ##  
range(Drych2$x);range(Drych2$y) 
coord.latlong = SpatialPoints(cbind(Drych2$x,Drych2$y), proj4string = CRS("+proj=longlat +datum=WGS84"))
coord.latlong
coord.UTM <- spTransform(coord.latlong , CRS("+proj=utm +zone=21K +south +datum=WGS84 +units=m +no_defs"))
coord.UTM
coord.latlong.df <- as.data.frame(coord.latlong)
coord.utm.df <- as.data.frame(coord.UTM)
head(coord.latlong.df)
head(coord.utm.df)
par(mfrow = c(1, 2))
plot(coord.latlong.df, axes = TRUE, main = "Lat-Long Coordinates", cex.axis = 0.95)
plot(coord.utm.df, axes = TRUE, main = "UTM Coordinates", col = "red", cex.axis = 0.95)
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
write.table(Drych2,file="c:/RWorkDir/jaguardatapaper/Drych2.txt",row.names = F,quote=F,col.names=T,sep="\t")
#Drych2 <- read.delim(file="c:/RWorkDir/jaguardatapaper/Drych2.txt")
head(Drych2); str(Drych2)


Drych=rbind(Drych1,Drych2)
head(Drych)
ft=ftable(Drych$id)
f<-as.data.frame.table(ft)
f$Var1 <- NULL
f$Var2 <- NULL
subset(f,Freq>0)
str(Drych)
#write.table(Drych,file="c:/RWorkDir/jaguardatapaper/Drych.txt",row.names = F,quote=F,col.names=T,sep="\t")
#Drych <- read.delim(file="c:/RWorkDir/jaguardatapaper/Drych.txt")
###############################################################################################################################
                                 ###    HUMID CHACO  ###
###############################################################################################################################
### # 9) Humid Chaco       
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
coord.utm.df <- as.data.frame(coord.UTM)
head(coord.latlong.df)
head(coord.utm.df)
par(mfrow = c(1, 2))
plot(coord.latlong.df, axes = TRUE, main = "Lat-Long Coordinates", cex.axis = 0.95)
plot(coord.utm.df, axes = TRUE, main = "UTM Coordinates", col = "red", cex.axis = 0.95)
locsj_matx=cbind(coordinates(coord.UTM ), coordinates(coord.latlong))
locsj_df<- as.data.frame(locsj_matx)
head(locsj_df)
head(Hch)
point.names<-c("utm_x","utm_y","long_x","lat_y")
colnames(locsj_df)<-point.names
head(locsj_df)

write.table(Hch,file="c:/RWorkDir/jaguardatapaper/Hch.txt",row.names = F,quote=F,col.names=T,sep="\t")
#Hch <- read.delim(file="c:/RWorkDir/jaguardatapaper/Hch.txt")
Hch=cbind(Hch, locsj_df)
head(Hch); str(Hch)

###############################################################################################################################
                               #######   FOREST PARAGUAY ######
###############################################################################################################################
### # 10) Forest Paraguay   
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
coord.utm.df <- as.data.frame(coord.UTM)
head(coord.latlong.df)
head(coord.utm.df)
par(mfrow = c(1, 2))
plot(coord.latlong.df, axes = TRUE, main = "Lat-Long Coordinates", cex.axis = 0.95)
plot(coord.utm.df, axes = TRUE, main = "UTM Coordinates", col = "red", cex.axis = 0.95)
locsj_matx=cbind(coordinates(coord.UTM ), coordinates(coord.latlong))
locsj_df<- as.data.frame(locsj_matx)
head(locsj_df)
head(FPy)
point.names<-c("utm_x","utm_y","long_x","lat_y")
colnames(locsj_df)<-point.names
head(locsj_df)

FPy=cbind(FPy, locsj_df)

write.table(FPy,file="c:/RWorkDir/jaguardatapaper/FPy.txt",row.names = F,quote=F,col.names=T,sep="\t")
#FPy <- read.delim(file="c:/RWorkDir/jaguardatapaper/FPy.txt")
head(FPy); str(FPy)

###############################################################################################################################
                           ############# Iguazu  #################  
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

### #11) Iguazu1    ## # 42,66,80,90 => ( Zone 21 J)
range(Iguazu1$x);range(Iguazu1$y) 
coord.latlong = SpatialPoints(cbind(Iguazu1$x,Iguazu1$y), proj4string = CRS("+proj=longlat +datum=WGS84"))
coord.latlong
coord.UTM <- spTransform(coord.latlong , CRS("+proj=utm +zone=21J +south +datum=WGS84 +units=m +no_defs"))
coord.UTM
coord.latlong.df <- as.data.frame(coord.latlong)
coord.utm.df <- as.data.frame(coord.UTM)
head(coord.latlong.df)
head(coord.utm.df)
par(mfrow = c(1, 2))
plot(coord.latlong.df, axes = TRUE, main = "Lat-Long Coordinates", cex.axis = 0.95)
plot(coord.utm.df, axes = TRUE, main = "UTM Coordinates", col = "red", cex.axis = 0.95)
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

write.table(Iguazu1,file="c:/RWorkDir/jaguardatapaper/Iguazu1.txt",row.names = F,quote=F,col.names=T,sep="\t")
#Iguazu1 <- read.delim(file="c:/RWorkDir/jaguardatapaper/Iguazu1.txt")
head(Iguazu1); str(Iguazu1)


### #12) Iguazu2    ##  83 => (22 J) 
range(Iguazu2$x);range(Iguazu2$y) 
coord.latlong = SpatialPoints(cbind(Iguazu2$x,Iguazu2$y), proj4string = CRS("+proj=longlat +datum=WGS84"))
coord.latlong
coord.UTM <- spTransform(coord.latlong , CRS("+proj=utm +zone=22J +south +datum=WGS84 +units=m +no_defs"))
coord.UTM
coord.latlong.df <- as.data.frame(coord.latlong)
coord.utm.df <- as.data.frame(coord.UTM)
head(coord.latlong.df)
head(coord.utm.df)
par(mfrow = c(1, 2))
plot(coord.latlong.df, axes = TRUE, main = "Lat-Long Coordinates", cex.axis = 0.95)
plot(coord.utm.df, axes = TRUE, main = "UTM Coordinates", col = "red", cex.axis = 0.95)
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

write.table(Iguazu2,file="c:/RWorkDir/jaguardatapaper/Iguazu2.txt",row.names = F,quote=F,col.names=T,sep="\t")
#Iguazu2 <- read.delim(file="c:/RWorkDir/jaguardatapaper/Iguazu2.txt")
head(Iguazu2); str(Iguazu2)


Iguazu=rbind(Iguazu1,Iguazu2)
head(Iguazu)
ft=ftable(Iguazu$id)
f<-as.data.frame.table(ft)
f$Var1 <- NULL
f$Var2 <- NULL
subset(f,Freq>0)

write.table(Iguazu,file="c:/RWorkDir/jaguardatapaper/Iguazu.txt",row.names = F,quote=F,col.names=T,sep="\t")
#Iguazu <- read.delim(file="c:/RWorkDir/jaguardatapaper/Iguazu.txt")
str(Iguazu)

###################################################################################################################################
                                    ###   AMAZONIA  ###            
###############################################################################################################################
 ### # 13)  Amazonia Mamiraua (Brazil)    ####  Flooded  ########################
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
coord.utm.df <- as.data.frame(coord.UTM)
head(coord.latlong.df)
head(coord.utm.df)
par(mfrow = c(1, 2))
plot(coord.latlong.df, axes = TRUE, main = "Lat-Long Coordinates", cex.axis = 0.95)
plot(coord.utm.df, axes = TRUE, main = "UTM Coordinates", col = "red", cex.axis = 0.95)
locsj_matx=cbind(coordinates(coord.UTM ), coordinates(coord.latlong))
locsj_df<- as.data.frame(locsj_matx)
head(locsj_df)
head(Mamiraua)
point.names<-c("utm_x","utm_y","long_x","lat_y")
colnames(locsj_df)<-point.names
head(locsj_df)

Mamiraua=cbind(Mamiraua, locsj_df)
write.table(Mamiraua,file="c:/RWorkDir/jaguardatapaper/Mamiraua.txt",row.names = F,quote=F,col.names=T,sep="\t")
#Mamiraua <- read.delim(file="c:/RWorkDir/jaguardatapaper/Mamiraua.txt")
head(Mamiraua); str(Mamiraua)

##########################################   Dry Amazonia, PA, translocated ################################################### 
### # 14 IOP Para Amazonia   ### Translocated         
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
coord.utm.df <- as.data.frame(coord.UTM)
head(coord.latlong.df)
head(coord.utm.df)
par(mfrow = c(1, 2))
plot(coord.latlong.df, axes = TRUE, main = "Lat-Long Coordinates", cex.axis = 0.95)
plot(coord.utm.df, axes = TRUE, main = "UTM Coordinates", col = "red", cex.axis = 0.95)
locsj_matx=cbind(coordinates(coord.UTM ), coordinates(coord.latlong))
locsj_df<- as.data.frame(locsj_matx)
head(locsj_df)
head(iopPA)
point.names<-c("utm_x","utm_y","long_x","lat_y")
colnames(locsj_df)<-point.names
head(locsj_df)

iopPA=cbind(iopPA, locsj_df)
write.table(iopPA,file="c:/RWorkDir/jaguardatapaper/iopPA.txt",row.names = F,quote=F,col.names=T,sep="\t")
#iopPA <- read.delim(file="c:/RWorkDir/jaguardatapaper/iopPA.txt")
head(iopPA); str(iopPA)

###################################################################################################################################
                                    ###   Mexico  ###           
###############################################################################################################################
#####  Greater Lacandona # 15) , Mexico  ###############################

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
coord.utm.df <- as.data.frame(coord.UTM)
head(coord.latlong.df)
head(coord.utm.df)
par(mfrow = c(1, 2))
plot(coord.latlong.df, axes = TRUE, main = "Lat-Long Coordinates", cex.axis = 0.95)
plot(coord.utm.df, axes = TRUE, main = "UTM Coordinates", col = "red", cex.axis = 0.95)
locsj_matx=cbind(coordinates(coord.UTM ), coordinates(coord.latlong))
locsj_df<- as.data.frame(locsj_matx)
head(locsj_df)
head(Lacandona)
point.names<-c("utm_x","utm_y","long_x","lat_y")
colnames(locsj_df)<-point.names
head(locsj_df)
Lacandona=cbind(Lacandona, locsj_df)
write.table(Lacandona,file="c:/RWorkDir/jaguardatapaper/Lacandona.txt",row.names = F,quote=F,col.names=T,sep="\t")
#Lacandona <- read.delim(file="c:/RWorkDir/jaguardatapaper/Lacandona.txt")
head(Lacandona); str(Lacandona)

#####   Mexico East # 16)  ####################################################################################################   
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
coord.utm.df <- as.data.frame(coord.UTM)
head(coord.latlong.df)
head(coord.utm.df)
par(mfrow = c(1, 2))
plot(coord.latlong.df, axes = TRUE, main = "Lat-Long Coordinates", cex.axis = 0.95)
plot(coord.utm.df, axes = TRUE, main = "UTM Coordinates", col = "red", cex.axis = 0.95)
locsj_matx=cbind(coordinates(coord.UTM ), coordinates(coord.latlong))
locsj_df<- as.data.frame(locsj_matx)
head(locsj_df)
head(MexEast)
point.names<-c("utm_x","utm_y","long_x","lat_y")
colnames(locsj_df)<-point.names
head(locsj_df)
MexEast=cbind(MexEast, locsj_df)
write.table(MexEast,file="c:/RWorkDir/jaguardatapaper/MexEast.txt",row.names = F,quote=F,col.names=T,sep="\t")
#MexEast <- read.delim(file="c:/RWorkDir/jaguardatapaper/MexEast.txt")
head(MexEast); str(MexEast)

#####   Mexico Sonora    # 17)  ##################################################################################################
Sonora =subset(jaguar_df,project_region=='Mexico Sonora')
head(Sonora)
str(Sonora)
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
coord.utm.df <- as.data.frame(coord.UTM)
head(coord.latlong.df)
head(coord.utm.df)
par(mfrow = c(1, 2))
plot(coord.latlong.df, axes = TRUE, main = "Lat-Long Coordinates", cex.axis = 0.95)
plot(coord.utm.df, axes = TRUE, main = "UTM Coordinates", col = "red", cex.axis = 0.95)
locsj_matx=cbind(coordinates(coord.UTM ), coordinates(coord.latlong))
locsj_df<- as.data.frame(locsj_matx)
head(locsj_df)
head(Sonora)
point.names<-c("utm_x","utm_y","long_x","lat_y")
colnames(locsj_df)<-point.names
head(locsj_df)
Sonora=cbind(Sonora, locsj_df)
write.table(Sonora,file="c:/RWorkDir/jaguardatapaper/Sonora.txt",row.names = F,quote=F,col.names=T,sep="\t")
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
write.table(Mex,file="c:/RWorkDir/jaguardatapaper/Mex.txt",row.names = F,quote=F,col.names=T,sep="\t")
#Mex <- read.delim(file="c:/RWorkDir/jaguardatapaper/Mex.txt")

###############################################################################################################################
                                    ###   PANTANAL  ###  
###############################################################################################################################
### # 18) PantanalTotal Brazil&Paraguay     7 Projects in total, all in the same UTM zone
###############################################################################################################################

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
coord.utm.df <- as.data.frame(coord.UTM)
head(coord.latlong.df)
head(coord.utm.df)
par(mfrow = c(1, 2))
plot(coord.latlong.df, axes = TRUE, main = "Lat-Long Coordinates", cex.axis = 0.95)
plot(coord.utm.df, axes = TRUE, main = "UTM Coordinates", col = "red", cex.axis = 0.95)
locsj_matx=cbind(coordinates(coord.UTM ), coordinates(coord.latlong))
locsj_df<- as.data.frame(locsj_matx)
head(locsj_df)
head(Pantanal)
point.names<-c("utm_x","utm_y","long_x","lat_y")
colnames(locsj_df)<-point.names
head(locsj_df)

Pantanal=cbind(Pantanal, locsj_df)
head(Pantanal)
 
write.table(Pantanal,file="c:/RWorkDir/jaguardatapaper/Pantanal.txt",row.names = F,quote=F,col.names=T,sep="\t")
#Pantanal <- read.delim(file="c:/RWorkDir/jaguardatapaper/Pantanal.txt")

#######################  Jaguar Dataframe with UTMs ########################################################################
head(AFW1);head(AFW2),head(Caatinga);head(Cerrado1);head(Cerrado2);head(CRica);head(Pantanal);head(Drych1);str(Drych2);head(Hch);
head(FPy);head(Iguazu1);head(Iguazu2);head(Mamiraua);head(iopPA);head(Lacandona);head(MexEast);head(Sonora)

jaguar=rbind(AFW1,AFW2,Caatinga,Cerrado1,Cerrado2,CRica,Pantanal,Drych1,Drych2,Hch,FPy,Iguazu1,Iguazu2,Mamiraua,iopPA,Lacandona,MexEast,Sonora)
#head(jaguar); str(jaguar)
jaguar_df <- jaguar

write.table(jaguar,file="c:/RWorkDir/jaguardatapaper/jaguar.txt",row.names = F,quote=F,col.names=T,sep="\t")
#write.table(jaguar_df,file="c:/RWorkDir/jaguardatapaper/jaguar_df.txt",row.names = F,quote=F,col.names=T,sep="\t")
#jaguar <- read.delim(file="c:/RWorkDir/jaguardatapaper/jaguar.txt")



###############################################################################################################################
#'              ### Creating a track in amt (commom to both RSF and SSF) based on UTM and project regions 
###############################################################################################################################
 
                               ### ATLANTIC FOREST WEST  ###    (Project regions 1 and 2)

#AFW <- read.delim(file="c:/RWorkDir/jaguardatapaper/AFW.txt")  # <-- if use read table it will be required to adjust variables again!!!
#age <- as.numeric(levels(AFW$age))[AFW$age]
#weight <- as.numeric(levels(AFW$weight))[AFW$weight]
#AFW$id <- as.factor(AFW$individual.local.identifier..ID.)  
#AFW$age <- age   ### converted to number
#AFW$weight <- weight   ### converted to number
#AFW <- read.delim(file="c:/RWorkDir/jaguardatapaper/AFW.txt")
#AFW$timestamp.posix <- as.character(AFW$timestamp.posix)
#AFW$timestamp.posix<- as.POSIXct(AFW$timestamp.posix, format ="%Y-%m-%d %H:%M:%S", tz = 'GMT') 
#AFW$date <- as.character(AFW$date)
#AFW$date <- as.POSIXct(AFW$date, format ="%Y-%m-%d %H:%M:%S", tz = 'GMT') 
head(AFW);str(AFW)

trk <- mk_track(AFW,.x=utm_x, .y=utm_y, .t=date, id=id,Event_ID=Event_ID,project_region=project_region, sex=sex, age=age, weight=weight, 
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
      ###  Atlantic Forest West =>  AFWtrk   ####
AFWtrk <-trk; AFWtrk
########################################
 AFWtrk=rbind(AFW1trk,AFW2trk); AFWtrk
   ### Objects for AFW project regions
 (AFW1trk=AFWtrk %>% filter(project_region == "Atlantic Forest W1"))
 (AFW2trk=AFWtrk %>% filter(project_region == "Atlantic Forest W2"))

################################################################################
                                      ###  CAATINGA  ###  
### # 3) Caatinga

#Caatinga <- read.delim(file="c:/RWorkDir/jaguardatapaper/Caatinga.txt")  # <-- if use read table it will be required to adjust variables again!!!
#age <- as.numeric(levels(Caatinga$age))[Caatinga$age]
#weight <- as.numeric(levels(Caatinga$weight))[Caatinga$weight]
#Caatinga$id <- as.factor(Caatinga$individual.local.identifier..ID.)  
#Caatinga$age <- age   ### converted to number
#Caatinga$weight <- weight   ### converted to number
#Caatinga <- read.delim(file="c:/RWorkDir/jaguardatapaper/Caatinga.txt")
#Caatinga$timestamp.posix <- as.character(Caatinga$timestamp.posix)
#Caatinga$timestamp.posix<- as.POSIXct(Caatinga$timestamp.posix, format ="%Y-%m-%d %H:%M:%S", tz = 'GMT') 
#Caatinga$date <- as.character(Caatinga$date)
#Caatinga$date <- as.POSIXct(Caatinga$date, format ="%Y-%m-%d %H:%M:%S", tz = 'GMT') 
head(Caatinga);str(Caatinga)

trk <- mk_track(Caatinga,.x=utm_x, .y=utm_y, .t=date, id=id,Event_ID=Event_ID,project_region=project_region, sex=sex, age=age, weight=weight, 
status=status, period=period,long_x=long_x, lat_y=lat_y, crs = CRS("+proj=utm +zone=23L +south +datum=WGS84 +units=m +no_defs"))   ###  WGS84
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
Caatingatrk <-trk; Caatingatrk
###########################################################################################################

###############################################################################################################################
                               ###  CERRADO   ###       (2 Projects: 4 and 5)
									
### # 4)  Cerrado1
#Cerrado1 <- read.delim(file="c:/RWorkDir/jaguardatapaper/Cerrado1.txt")  # <-- if use read table it will be required to adjust variables again!!!
#age <- as.numeric(levels(Cerrado1$age))[Cerrado1$age]
#weight <- as.numeric(levels(Cerrado1$weight))[Cerrado1$weight]
#Cerrado1$id <- as.factor(Cerrado1$individual.local.identifier..ID.)  
#Cerrado1$age <- age   ### converted to number
#Cerrado1$weight <- weight   ### converted to number
#Cerrado1 <- read.delim(file="c:/RWorkDir/jaguardatapaper/Cerrado1.txt")
#Cerrado1$timestamp.posix <- as.character(Cerrado1$timestamp.posix)
#Cerrado1$timestamp.posix<- as.POSIXct(Cerrado1$timestamp.posix, format ="%Y-%m-%d %H:%M:%S", tz = 'GMT') 
#Cerrado1$date <- as.character(Cerrado1$date)
#Cerrado1$date <- as.POSIXct(Cerrado1$date, format ="%Y-%m-%d %H:%M:%S", tz = 'GMT') 
head(Cerrado1);str(Cerrado1)

trk <- mk_track(Cerrado1,.x=utm_x, .y=utm_y, .t=date, id=id,Event_ID=Event_ID, project_region=project_region,sex=sex, age=age, weight=weight, 
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
Cerrado1trk <-trk; Cerrado1trk
########################################

### # 5)   Cerrado2

#Cerrado2 <- read.delim(file="c:/RWorkDir/jaguardatapaper/Cerrado2.txt")  # <-- if use read table it will be required to adjust variables again!!!
#age <- as.numeric(levels(Cerrado2$age))[Cerrado2$age]
#weight <- as.numeric(levels(Cerrado2$weight))[Cerrado2$weight]
#Cerrado2$id <- as.factor(Cerrado2$individual.local.identifier..ID.)  
#Cerrado2$age <- age   ### converted to number
#Cerrado2$weight <- weight   ### converted to number
#Cerrado2 <- read.delim(file="c:/RWorkDir/jaguardatapaper/Cerrado2.txt")
#Cerrado2$timestamp.posix <- as.character(Cerrado2$timestamp.posix)
#Cerrado2$timestamp.posix<- as.POSIXct(Cerrado2$timestamp.posix, format ="%Y-%m-%d %H:%M:%S", tz = 'GMT') 
#Cerrado2$date <- as.character(Cerrado2$date)
#Cerrado2$date <- as.POSIXct(Cerrado2$date, format ="%Y-%m-%d %H:%M:%S", tz = 'GMT') 
head(Cerrado2);str(Cerrado2)

trk <- mk_track(Cerrado2,.x=utm_x, .y=utm_y, .t=date, id=id,Event_ID=Event_ID,project_region=project_region, sex=sex, age=age, weight=weight, 
status=status, period=period,long_x=long_x, lat_y=lat_y, crs = CRS("+proj=utm +zone=22L +south +datum=WGS84 +units=m +no_defs"))   ###  WGS84
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
Cerrado2trk <-trk; Cerrado2trk
########################################
###  Cerradotrk   ##########################
Cerradotrk=rbind(Cerrado1trk,Cerrado2trk); Cerradotrk

###############################################################################################################################
                                   ###  COSTA RICA  ### 
### # 6) Costa Rica

#CRica <- read.delim(file="c:/RWorkDir/jaguardatapaper/CRica.txt")  # <-- if use read table it will be required to adjust variables again!!!
#age <- as.numeric(levels(CRica$age))[CRica$age]
#weight <- as.numeric(levels(CRica$weight))[CRica$weight]
#CRica$id <- as.factor(CRica$individual.local.identifier..ID.)  
#CRica$age <- age   ### converted to number
#CRica$weight <- weight   ### converted to number
#CRica <- read.delim(file="c:/RWorkDir/jaguardatapaper/CRica.txt")
#CRica$timestamp.posix <- as.character(CRica$timestamp.posix)
#CRica$timestamp.posix<- as.POSIXct(CRica$timestamp.posix, format ="%Y-%m-%d %H:%M:%S", tz = 'GMT') 
#CRica$date <- as.character(CRica$date)
#CRica$date <- as.POSIXct(CRica$date, format ="%Y-%m-%d %H:%M:%S", tz = 'GMT') 
head(CRica);str(CRica)

trk <- mk_track(CRica,.x=utm_x, .y=utm_y, .t=date, id=id,Event_ID=Event_ID,project_region=project_region, sex=sex, age=age, weight=weight, 
status=status, period=period,long_x=long_x, lat_y=lat_y, crs = CRS("+proj=utm +zone=16P +north +datum=WGS84 +units=m +no_defs"))   ###  WGS84
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
CRicatrk <-trk; CRicatrk

###############################################################################################################################
                                    ###   DRY CHACO  ###   

### # 7) Drych1 

#Drych1 <- read.delim(file="c:/RWorkDir/jaguardatapaper/Drych1.txt")  # <-- if use read table it will be required to adjust variables again!!!
#age <- as.numeric(levels(Drych1$age))[Drych1$age]
#weight <- as.numeric(levels(Drych1$weight))[Drych1$weight]
#Drych1$id <- as.factor(Drych1$individual.local.identifier..ID.)  
#Drych1$age <- age   ### converted to number
#Drych1$weight <- weight   ### converted to number
#Drych1 <- read.delim(file="c:/RWorkDir/jaguardatapaper/Drych1.txt")
#Drych1$timestamp.posix <- as.character(Drych1$timestamp.posix)
#Drych1$timestamp.posix<- as.POSIXct(Drych1$timestamp.posix, format ="%Y-%m-%d %H:%M:%S", tz = 'GMT') 
#Drych1$date <- as.character(Drych1$date)
#Drych1$date <- as.POSIXct(Drych1$date, format ="%Y-%m-%d %H:%M:%S", tz = 'GMT') 
head(Drych1);str(Drych1)

trk <- mk_track(Drych1,.x=utm_x, .y=utm_y, .t=date, id=id,Event_ID=Event_ID,project_region=project_region, sex=sex, age=age, weight=weight, 
status=status, period=period,long_x=long_x, lat_y=lat_y, crs = CRS("+proj=utm +zone=20K +south +datum=WGS84 +units=m +no_defs"))   ###  WGS84
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
Drych1trk <-trk
########################################

### # 8) Drych2   

#Drych2 <- read.delim(file="c:/RWorkDir/jaguardatapaper/Drych2.txt")  # <-- if use read table it will be required to adjust variables again!!!
#age <- as.numeric(levels(Drych2$age))[Drych2$age]
#weight <- as.numeric(levels(Drych2$weight))[Drych2$weight]
#Drych2$id <- as.factor(Drych2$individual.local.identifier..ID.)  
#Drych2$age <- age   ### converted to number
#Drych2$weight <- weight   ### converted to number
#Drych2 <- read.delim(file="c:/RWorkDir/jaguardatapaper/Drych2.txt")
#Drych2$timestamp.posix <- as.character(Drych2$timestamp.posix)
#Drych2$timestamp.posix<- as.POSIXct(Drych2$timestamp.posix, format ="%Y-%m-%d %H:%M:%S", tz = 'GMT') 
#Drych2$date <- as.character(Drych2$date)
#Drych2$date <- as.POSIXct(Drych2$date, format ="%Y-%m-%d %H:%M:%S", tz = 'GMT') 
head(Drych2);str(Drych2)

trk <- mk_track(Drych2,.x=utm_x, .y=utm_y, .t=date, id=id,Event_ID=Event_ID,project_region=project_region, sex=sex, age=age, weight=weight, 
status=status, period=period,long_x=long_x, lat_y=lat_y, crs = CRS("+proj=utm +zone=21K +south +datum=WGS84 +units=m +no_defs"))   ###  WGS84
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
Drych2trk <-trk

###  DDRY CHACO trk => Drychtrk   ##########################
Drychtrk=rbind(Drych1trk,Drych2trk); Drychtrk

###############################################################################################################################
                                 ###    HUMID CHACO  ###
### # 9) Humid Chaco       

#Hch <- read.delim(file="c:/RWorkDir/jaguardatapaper/Hch.txt")  # <-- if use read table it will be required to adjust variables again!!!
#age <- as.numeric(levels(Hch$age))[Hch$age]
#weight <- as.numeric(levels(Hch$weight))[Hch$weight]
#Hch$id <- as.factor(Hch$individual.local.identifier..ID.)  
#Hch$age <- age   ### converted to number
#Hch$weight <- weight   ### converted to number
#Hch <- read.delim(file="c:/RWorkDir/jaguardatapaper/Hch.txt")
#Hch$timestamp.posix <- as.character(Hch$timestamp.posix)
#Hch$timestamp.posix<- as.POSIXct(Hch$timestamp.posix, format ="%Y-%m-%d %H:%M:%S", tz = 'GMT') 
#Hch$date <- as.character(Hch$date)
#Hch$date <- as.POSIXct(Hch$date, format ="%Y-%m-%d %H:%M:%S", tz = 'GMT') 
head(Hch);str(Hch)

trk <- mk_track(Hch,.x=utm_x, .y=utm_y, .t=date, id=id,Event_ID=Event_ID, project_region=project_region,sex=sex, age=age, weight=weight, 
status=status, period=period,long_x=long_x, lat_y=lat_y, crs = CRS("+proj=utm +zone=21K +south +datum=WGS84 +units=m +no_defs"))   ###  WGS84
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
Hchtrk <-trk
###############################################################################################################################
                               #######   FOREST PARAGUAY ######
							   
### # 10) Forest Paraguay   

#FPy <- read.delim(file="c:/RWorkDir/jaguardatapaper/FPy.txt")  # <-- if use read table it will be required to adjust variables again!!!
#age <- as.numeric(levels(FPy$age))[FPy$age]
#weight <- as.numeric(levels(FPy$weight))[FPy$weight]
#FPy$id <- as.factor(FPy$individual.local.identifier..ID.)  
#FPy$age <- age   ### converted to number
#FPy$weight <- weight   ### converted to number
#FPy <- read.delim(file="c:/RWorkDir/jaguardatapaper/FPy.txt")
#FPy$timestamp.posix <- as.character(FPy$timestamp.posix)
#FPy$timestamp.posix<- as.POSIXct(FPy$timestamp.posix, format ="%Y-%m-%d %H:%M:%S", tz = 'GMT') 
#FPy$date <- as.character(FPy$date)
#FPy$date <- as.POSIXct(FPy$date, format ="%Y-%m-%d %H:%M:%S", tz = 'GMT') 
head(FPy);str(FPy)

trk <- mk_track(FPy,.x=utm_x, .y=utm_y, .t=date, id=id,Event_ID=Event_ID, project_region=project_region,sex=sex, age=age, weight=weight, 
status=status, period=period,long_x=long_x, lat_y=lat_y, crs = CRS("+proj=utm +zone=21J +south +datum=WGS84 +units=m +no_defs"))   ###  WGS84
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
FPytrk <-trk

###############################################################################################################################
                           ############# Iguazu  #################  
	   
### # 11      Iguazu1

#Iguazu1 <- read.delim(file="c:/RWorkDir/jaguardatapaper/Iguazu1.txt")  # <-- if use read table it will be required to adjust variables again!!!
#age <- as.numeric(levels(Iguazu1$age))[Iguazu1$age]
#weight <- as.numeric(levels(Iguazu1$weight))[Iguazu1$weight]
#Iguazu1$id <- as.factor(Iguazu1$individual.local.identifier..ID.)  
#Iguazu1$age <- age   ### converted to number
#Iguazu1$weight <- weight   ### converted to number
#Iguazu1 <- read.delim(file="c:/RWorkDir/jaguardatapaper/Iguazu1.txt")
#Iguazu1$timestamp.posix <- as.character(Iguazu1$timestamp.posix)
#Iguazu1$timestamp.posix<- as.POSIXct(Iguazu1$timestamp.posix, format ="%Y-%m-%d %H:%M:%S", tz = 'GMT') 
#Iguazu1$date <- as.character(Iguazu1$date)
#Iguazu1$date <- as.POSIXct(Iguazu1$date, format ="%Y-%m-%d %H:%M:%S", tz = 'GMT') 
head(Iguazu1);str(Iguazu1)

trk <- mk_track(Iguazu1,.x=utm_x, .y=utm_y, .t=date, id=id,Event_ID=Event_ID,project_region=project_region, sex=sex, age=age, weight=weight, 
status=status, period=period,long_x=long_x, lat_y=lat_y, crs = CRS("+proj=utm +zone=21J +south +datum=WGS84 +units=m +no_defs"))   ###  WGS84
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
Iguazu1trk <-trk


################################################# 
### # 12 Iguazu2

#Iguazu2 <- read.delim(file="c:/RWorkDir/jaguardatapaper/Iguazu2.txt")  # <-- if use read table it will be required to adjust variables again!!!
#age <- as.numeric(levels(Iguazu2$age))[Iguazu2$age]
#weight <- as.numeric(levels(Iguazu2$weight))[Iguazu2$weight]
#Iguazu2$id <- as.factor(Iguazu2$individual.local.identifier..ID.)  
#Iguazu2$age <- age   ### converted to number
#Iguazu2$weight <- weight   ### converted to number
#Iguazu2 <- read.delim(file="c:/RWorkDir/jaguardatapaper/Iguazu2.txt")
#Iguazu2$timestamp.posix <- as.character(Iguazu2$timestamp.posix)
#Iguazu2$timestamp.posix<- as.POSIXct(Iguazu2$timestamp.posix, format ="%Y-%m-%d %H:%M:%S", tz = 'GMT') 
#Iguazu2$date <- as.character(Iguazu2$date)
#Iguazu2$date <- as.POSIXct(Iguazu2$date, format ="%Y-%m-%d %H:%M:%S", tz = 'GMT') 
head(Iguazu2);str(Iguazu2)

trk <- mk_track(Iguazu2,.x=utm_x, .y=utm_y, .t=date, id=id,Event_ID=Event_ID, project_region=project_region,sex=sex, age=age, weight=weight, 
status=status, period=period,long_x=long_x, lat_y=lat_y, crs = CRS("+proj=utm +zone=22J +south +datum=WGS84 +units=m +no_defs"))   ###  WGS84
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
Iguazu2trk <-trk
########################################
###  Iguazutrk   ##########################
Iguazutrk=rbind(Iguazu1trk,Iguazu2trk); Iguazutrk

###################################################################################################################################
                                    ###   AMAZONIA  ###    
 ### # 13)  Amazonia Mamiraua (Brazil)    ####  Flooded  ########################
 
#Mamiraua <- read.delim(file="c:/RWorkDir/jaguardatapaper/Mamiraua.txt")  # <-- if use read table it will be required to adjust variables again!!!
#age <- as.numeric(levels(Mamiraua$age))[Mamiraua$age]
#weight <- as.numeric(levels(Mamiraua$weight))[Mamiraua$weight]
#Mamiraua$id <- as.factor(Mamiraua$individual.local.identifier..ID.)  
#Mamiraua$age <- age   ### converted to number
#Mamiraua$weight <- weight   ### converted to number
#Mamiraua <- read.delim(file="c:/RWorkDir/jaguardatapaper/Mamiraua.txt")
#Mamiraua$timestamp.posix <- as.character(Mamiraua$timestamp.posix)
#Mamiraua$timestamp.posix<- as.POSIXct(Mamiraua$timestamp.posix, format ="%Y-%m-%d %H:%M:%S", tz = 'GMT') 
#Mamiraua$date <- as.character(Mamiraua$date)
#Mamiraua$date <- as.POSIXct(Mamiraua$date, format ="%Y-%m-%d %H:%M:%S", tz = 'GMT') 
head(Mamiraua);str(Mamiraua)

trk <- mk_track(Mamiraua,.x=utm_x, .y=utm_y, .t=date, id=id,Event_ID=Event_ID,project_region=project_region, sex=sex, age=age, weight=weight, 
status=status, period=period,long_x=long_x, lat_y=lat_y, crs = CRS("+proj=utm +zone=20M +south +datum=WGS84 +units=m +no_defs"))   ###  WGS84
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
Mamirauatrk <-trk

################################################################################################################################# 
### # 14        Dry/ East Amazonia, IOP PA, translocated

#iopPA <- read.delim(file="c:/RWorkDir/jaguardatapaper/iopPA.txt")  # <-- if use read table it will be required to adjust variables again!!!
#age <- as.numeric(levels(iopPA$age))[iopPA$age]
#weight <- as.numeric(levels(iopPA$weight))[iopPA$weight]
#iopPA$id <- as.factor(iopPA$individual.local.identifier..ID.)  
#iopPA$age <- age   ### converted to number
#iopPA$weight <- weight   ### converted to number
#iopPA <- read.delim(file="c:/RWorkDir/jaguardatapaper/iopPA.txt")
#iopPA$timestamp.posix <- as.character(iopPA$timestamp.posix)
#iopPA$timestamp.posix<- as.POSIXct(iopPA$timestamp.posix, format ="%Y-%m-%d %H:%M:%S", tz = 'GMT') 
#iopPA$date <- as.character(iopPA$date)
#iopPA$date <- as.POSIXct(iopPA$date, format ="%Y-%m-%d %H:%M:%S", tz = 'GMT') 
head(iopPA);str(iopPA)

trk <- mk_track(iopPA,.x=utm_x, .y=utm_y, .t=date, id=id,Event_ID=Event_ID,project_region=project_region, sex=sex, age=age, weight=weight, 
status=status, period=period,long_x=long_x, lat_y=lat_y, crs = CRS("+proj=utm +zone=22M +south +datum=WGS84 +units=m +no_defs"))   ###  WGS84
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
iopPAtrk <-trk

###################################################################################################################################
                                  ### Greater Lacandona, Mexico ##
											 
 ### # 15 Lacandona
#Lacandona <- read.delim(file="c:/RWorkDir/jaguardatapaper/Lacandona.txt")  # <-- if use read table it will be required to adjust variables again!!!
#age <- as.numeric(levels(Lacandona$age))[Lacandona$age]
#weight <- as.numeric(levels(Lacandona$weight))[Lacandona$weight]
#Lacandona$id <- as.factor(Lacandona$individual.local.identifier..ID.)  
#Lacandona$age <- age   ### converted to number
#Lacandona$weight <- weight   ### converted to number
#Lacandona <- read.delim(file="c:/RWorkDir/jaguardatapaper/Lacandona.txt")
#Lacandona$timestamp.posix <- as.character(Lacandona$timestamp.posix)
#Lacandona$timestamp.posix<- as.POSIXct(Lacandona$timestamp.posix, format ="%Y-%m-%d %H:%M:%S", tz = 'GMT') 
#Lacandona$date <- as.character(Lacandona$date)
#Lacandona$date <- as.POSIXct(Lacandona$date, format ="%Y-%m-%d %H:%M:%S", tz = 'GMT') 
head(Lacandona);str(Lacandona)

trk <- mk_track(Lacandona,.x=utm_x, .y=utm_y, .t=date, id=id,Event_ID=Event_ID,project_region=project_region, sex=sex, age=age, weight=weight, 
status=status, period=period,long_x=long_x, lat_y=lat_y, crs = CRS("+proj=utm +zone=15Q +south +datum=WGS84 +units=m +no_defs"))   ###  WGS84
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
Lacandonatrk <-trk


######################################################################################################################## 
                                      ### # 16) MexEast, Mexico
 ### # 16) MexEast
#MexEast <- read.delim(file="c:/RWorkDir/jaguardatapaper/MexEast.txt")  # <-- if use read table it will be required to adjust variables again!!!
#age <- as.numeric(levels(MexEast$age))[MexEast$age]
#weight <- as.numeric(levels(MexEast$weight))[MexEast$weight]
#MexEast$id <- as.factor(MexEast$individual.local.identifier..ID.)  
#MexEast$age <- age   ### converted to number
#MexEast$weight <- weight   ### converted to number
#MexEast <- read.delim(file="c:/RWorkDir/jaguardatapaper/MexEast.txt")
#MexEast$timestamp.posix <- as.character(MexEast$timestamp.posix)
#MexEast$timestamp.posix<- as.POSIXct(MexEast$timestamp.posix, format ="%Y-%m-%d %H:%M:%S", tz = 'GMT') 
#MexEast$date <- as.character(MexEast$date)
#MexEast$date <- as.POSIXct(MexEast$date, format ="%Y-%m-%d %H:%M:%S", tz = 'GMT') 
head(MexEast);str(MexEast)

trk <- mk_track(MexEast,.x=utm_x, .y=utm_y, .t=date, id=id,Event_ID=Event_ID,project_region=project_region, sex=sex, age=age, weight=weight, 
status=status, period=period,long_x=long_x, lat_y=lat_y, crs = CRS("+proj=utm +zone=16Q +south +datum=WGS84 +units=m +no_defs"))   ###  WGS84
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
MexEasttrk <-trk


##################################################################################################################################### 
                              ### Sonora, Mexico ###

 ### # 17)  Sonora
#Sonora <- read.delim(file="c:/RWorkDir/jaguardatapaper/Sonora.txt")  # <-- if use read table it will be required to adjust variables again!!!
#age <- as.numeric(levels(Sonora$age))[Sonora$age]
#weight <- as.numeric(levels(Sonora$weight))[Sonora$weight]
#Sonora$id <- as.factor(Sonora$individual.local.identifier..ID.)  
#Sonora$age <- age   ### converted to number
#Sonora$weight <- weight   ### converted to number
#Sonora <- read.delim(file="c:/RWorkDir/jaguardatapaper/Sonora.txt")
#Sonora$timestamp.posix <- as.character(Sonora$timestamp.posix)
#Sonora$timestamp.posix<- as.POSIXct(Sonora$timestamp.posix, format ="%Y-%m-%d %H:%M:%S", tz = 'GMT') 
#Sonora$date <- as.character(Sonora$date)
#Sonora$date <- as.POSIXct(Sonora$date, format ="%Y-%m-%d %H:%M:%S", tz = 'GMT') 
head(Sonora);str(Sonora)

trk <- mk_track(Sonora,.x=utm_x, .y=utm_y, .t=date, id=id,Event_ID=Event_ID,project_region=project_region, sex=sex, age=age, weight=weight, 
status=status, period=period,long_x=long_x, lat_y=lat_y, crs = CRS("+proj=utm +zone=12R +south +datum=WGS84 +units=m +no_defs"))   ###  WGS84
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
Sonoratrk <-trk

###############################################################################################################################
                                  ######   PANTANAL  ######  
									
### # 18) PantanalTotal Brazil&Paraguay     7 Projects in total, all in the same UTM zone
#Pantanal <- read.delim(file="c:/RWorkDir/jaguardatapaper/Pantanal.txt")  # <-- if use read table it will be required to adjust variables again!!!
#age <- as.numeric(levels(Pantanal$age))[Pantanal$age]
#weight <- as.numeric(levels(Pantanal$weight))[Pantanal$weight]
#Pantanal$id <- as.factor(Pantanal$individual.local.identifier..ID.)  
#Pantanal$age <- age   ### converted to number
#Pantanal$weight <- weight   ### converted to number
#Pantanal <- read.delim(file="c:/RWorkDir/jaguardatapaper/Pantanal.txt")
#Pantanal$timestamp.posix <- as.character(Pantanal$timestamp.posix)
#Pantanal$timestamp.posix<- as.POSIXct(Pantanal$timestamp.posix, format ="%Y-%m-%d %H:%M:%S", tz = 'GMT') 
#Pantanal$date <- as.character(Pantanal$date)
#Pantanal$date <- as.POSIXct(Pantanal$date, format ="%Y-%m-%d %H:%M:%S", tz = 'GMT') 
head(Pantanal);str(Pantanal)

trk <- mk_track(Pantanal,.x=utm_x, .y=utm_y, .t=date, id=id,Event_ID=Event_ID,project_region=project_region, sex=sex, age=age, weight=weight, 
status=status, period=period,long_x=long_x, lat_y=lat_y, crs = CRS("+proj=utm +zone=21K +south +datum=WGS84 +units=m +no_defs"))   ###  WGS84
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
(Pantanaltrk <-trk)
 
########################################
 ### Objects for Pantanal project regions
 (Oncafaritrk=Pantanaltrk %>% filter(project_region == "Oncafari"))
 (Paraguaytrk=Pantanaltrk %>% filter(project_region == "Pantanal Paraguay"))
 (Panthera1trk=Pantanaltrk %>% filter(project_region == "Panthera1"))
 (Panthera2trk=Pantanaltrk %>% filter(project_region == "Panthera2"))
 (RioNegrotrk=Pantanaltrk %>% filter(project_region == "Rio Negro"))
 (SaoBentotrk=Pantanaltrk %>% filter(project_region == "Sao Bento"))
 (Taiamatrk=Pantanaltrk %>% filter(project_region == "Taiama"))
 
 (Ptrk=rbind(Oncafaritrk,Paraguaytrk,Panthera1trk,Panthera2trk,RioNegrotrk,SaoBentotrk,Taiamatrk))
 Pantanaltrk

################################################################################################################################
           ###################    ALL  JAGUARS  trk =>     jaguartrk           ########################
################################################################################################################################
		   
jaguartrk=rbind(AFWtrk,Caatingatrk,Cerrado1trk,Cerrado2trk,CRicatrk,Pantanaltrk,Drych1trk,Drych2trk,Hchtrk,FPytrk,
Iguazu1trk, Iguazu2trk, Mamirauatrk,iopPAtrk, Lacandonatrk, MexEasttrk, Sonoratrk)
jaguartrk
### All Project regions trk
jtrk=rbind(AFW1trk,AFW2trk,Caatingatrk,Cerrado1trk,Cerrado2trk,CRicatrk,Drychtrk, Hchtrk,FPytrk,Iguazutrk,Mamirauatrk,iopPAtrk,
			 Lacandonatrk, MexEasttrk, Sonoratrk,Oncafaritrk,Paraguaytrk,Panthera1trk,Panthera2trk,RioNegrotrk,SaoBentotrk,Taiamatrk) 
			 
### jaguartrk <- jtrk

#################################################################################################################################
 ### Dataframe for each individual (run only if need any specic individual)
 
 X1=subset(jaguar_df,id=='1')
 X2=subset(jaguar_df,id=='2')
 X3=subset(jaguar_df,id=='3')
 X4=subset(jaguar_df,id=='4')
 X5=subset(jaguar_df,id=='5')
 X6=subset(jaguar_df,id=='6')
 X7=subset(jaguar_df,id=='7')
 X8=subset(jaguar_df,id=='8')
 X9=subset(jaguar_df,id=='9')
X10=subset(jaguar_df,id=='10') 
X11=subset(jaguar_df,id=='11') 
X12=subset(jaguar_df,id=='12') 
X13=subset(jaguar_df,id=='13')
X14=subset(jaguar_df,id=='14') 
X15=subset(jaguar_df,id=='15') 
X16=subset(jaguar_df,id=='16') 
X17=subset(jaguar_df,id=='17') 
X18=subset(jaguar_df,id=='18') 
X19=subset(jaguar_df,id=='19') 
X20=subset(jaguar_df,id=='20')
X21=subset(jaguar_df,id=='21')
X22=subset(jaguar_df,id=='22')
X23=subset(jaguar_df,id=='23')
X24=subset(jaguar_df,id=='24')
X25=subset(jaguar_df,id=='25')
X26=subset(jaguar_df,id=='26')
X27=subset(jaguar_df,id=='27')
X28=subset(jaguar_df,id=='28')
X29=subset(jaguar_df,id=='29')
X30=subset(jaguar_df,id=='30')
X31=subset(jaguar_df,id=='31')
X32=subset(jaguar_df,id=='32')
X33=subset(jaguar_df,id=='33')
X34=subset(jaguar_df,id=='34')
X35=subset(jaguar_df,id=='35')
X36=subset(jaguar_df,id=='36')
X37=subset(jaguar_df,id=='37')
X38=subset(jaguar_df,id=='38')
X39=subset(jaguar_df,id=='39')
X40=subset(jaguar_df,id=='40')
X41=subset(jaguar_df,id=='41')
X42=subset(jaguar_df,id=='42')
X43=subset(jaguar_df,id=='43')
X44=subset(jaguar_df,id=='44')
X45=subset(jaguar_df,id=='45')
X46=subset(jaguar_df,id=='46')
X47=subset(jaguar_df,id=='47')
X48=subset(jaguar_df,id=='48')
X49=subset(jaguar_df,id=='49')
X50=subset(jaguar_df,id=='50')
X51=subset(jaguar_df,id=='51')
X52=subset(jaguar_df,id=='52')
X53=subset(jaguar_df,id=='53')
X54=subset(jaguar_df,id=='54')
X55=subset(jaguar_df,id=='55')
X56=subset(jaguar_df,id=='56')
X57=subset(jaguar_df,id=='57')
X58=subset(jaguar_df,id=='58')
X59=subset(jaguar_df,id=='59')
X60=subset(jaguar_df,id=='60')
X61=subset(jaguar_df,id=='61')
X62=subset(jaguar_df,id=='62')
X63=subset(jaguar_df,id=='63')
X64=subset(jaguar_df,id=='64')
X65=subset(jaguar_df,id=='65')
X66=subset(jaguar_df,id=='66')
X67=subset(jaguar_df,id=='67')
X68=subset(jaguar_df,id=='68')
X69=subset(jaguar_df,id=='69')
X70=subset(jaguar_df,id=='70')
X71=subset(jaguar_df,id=='71')
X72=subset(jaguar_df,id=='72')
X73=subset(jaguar_df,id=='73')
X74=subset(jaguar_df,id=='74')
X75=subset(jaguar_df,id=='75')
X76=subset(jaguar_df,id=='76')
X77=subset(jaguar_df,id=='77')
X78=subset(jaguar_df,id=='78')
X79=subset(jaguar_df,id=='79')
X80=subset(jaguar_df,id=='80')
X81=subset(jaguar_df,id=='81')
X82=subset(jaguar_df,id=='82')
X83=subset(jaguar_df,id=='83')
X84=subset(jaguar_df,id=='84')
X85=subset(jaguar_df,id=='85')
X86=subset(jaguar_df,id=='86')
X87=subset(jaguar_df,id=='87')
X88=subset(jaguar_df,id=='88')
X89=subset(jaguar_df,id=='89')
X90=subset(jaguar_df,id=='90')
X91=subset(jaguar_df,id=='91')
X92=subset(jaguar_df,id=='92')
X93=subset(jaguar_df,id=='93')
X94=subset(jaguar_df,id=='94')
X95=subset(jaguar_df,id=='95')
X96=subset(jaguar_df,id=='96')
X97=subset(jaguar_df,id=='97')
X98=subset(jaguar_df,id=='98')
X99=subset(jaguar_df,id=='99')
X100=subset(jaguar_df,id=='100')
X101=subset(jaguar_df,id=='101')
X102=subset(jaguar_df,id=='102')
X103=subset(jaguar_df,id=='103')
X104=subset(jaguar_df,id=='104')
X105=subset(jaguar_df,id=='105')
X106=subset(jaguar_df,id=='106')
X107=subset(jaguar_df,id=='107')
X108=subset(jaguar_df,id=='108')
X109=subset(jaguar_df,id=='109')
X110=subset(jaguar_df,id=='110')
X111=subset(jaguar_df,id=='111')
X112=subset(jaguar_df,id=='112')
X113=subset(jaguar_df,id=='113')
X114=subset(jaguar_df,id=='114')
X115=subset(jaguar_df,id=='115')
X116=subset(jaguar_df,id=='116')
X117=subset(jaguar_df,id=='117')

 