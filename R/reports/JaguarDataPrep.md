---
title: "Jaguar Data Preparation"
authors: "Alan E. de Barros,Bernardo Niebuhr, Vanesa Bejarano, Julia Oshima, Claudia Kanda, Milton Ribeiro, Ronaldo Morato, Paulo Prado"
date: ""
---


#  **Jaguar Data Preparation**

#### *Alan E. de Barros,Bernardo Niebuhr,Vanesa Bejarano,Julia Oshima,Claudia Kanda,Milton Ribeiro,Ronaldo Morato,Paulo Prado*
date: "March, 08 2019"


```r
# Adapted from Bernardo Niebuhr data preparation, Luca Borger lectures and John Fieberg's amt's movebank scripts
```






#### *Preamble*
For a fresh start, clean everything in working memory


```r
rm(list= ls())                                                 
```


Install and load packages


```r
if(!require(install.load)) install.packages('install.load'); library(install.load)
install.load::install_load("move", "adehabitatLT", "amt") # Movement packages
install.load::install_load("maptools", "raster", "rgdal","sp") # Spatial packages
install.load::install_load("colorspace","ggmap", "rgl", "lattice", "leaflet") # Visualization packages
install.load::install_load("RCurl", "dplyr", "readr", "lubridate", "tibble") # Aux packages
install.load::install_load("circular", "caTools") # Stats packages
install.load::install_load("knitr", "ezknitr") # To render documents 
```



### Source functions from GitHub local directory 


```r
# Do not required to set directory if Rscript have been opened from the GitHub local directory
source("DataPrepFunctions.R") 
```


#### Load the data and create a dataframe object:


```r
mov.data.org <- read.delim(file="../data/mov.data.org.txt")
mov.data.org <- dplyr::select(mov.data.org, -(individual.taxon.canonical.name:tag.local.identifier))
```


#### Add Individual info  (see older development files for details)


```r
info <- read.delim(file="../data/info.txt")
```


#### Merge movement with individual info/parameters


```r
merged<- merge(mov.data.org,info)
mov.data.org <- merged 
```


#### Organize data

Test
get.year(time.stamp = mov.data.org$timestamp[10000])

All individuals


```r
year <- as.numeric(sapply(mov.data.org$timestamp, get.year))
```

Add 1900/2000


```r
new.year <- as.character(ifelse(year > 50, year + 1900, year + 2000))
```

Test


```r
set.year(time.stamp = as.character(mov.data.org$timestamp[10000]), year = '2013')
```

```
## [1] "12/19/2013 21:00"
```

All individuals


```r
date.time <- as.character(mapply(set.year, as.character(mov.data.org$timestamp),
                                 new.year))
```

Date/Time as POSIXct object


```r
mov.data.org$timestamp.posix <- as.POSIXct(date.time, 
                                           format = "%m/%d/%Y %H:%M", tz = 'GMT')
mov.data.org$GMTtime <- mov.data.org$timestamp.posix
```



## Get local time
 


```r
# A column to represent the local timezone (already with the - signal) has been included to then multiply the timestamp and get the difference:
mov.data.org$local_time <- mov.data.org$timestamp.posix + mov.data.org$timezone*60*60
mov.data.org$timestamp.posix <- mov.data.org$local_time 
```

Now all the (timestamp.posix)'s calculations are based on local time


### Adjusting for Movement Packages
### adehabitatLT

Transforms in ltraj object


```r
coords <- data.frame(mov.data.org$location.long, mov.data.org$location.lat)
mov.traj <- as.ltraj(xy = coords, date=mov.data.org$timestamp.posix, 
                     id=mov.data.org$individual.local.identifier..ID., 
                     burst=mov.data.org$individual.local.identifier..ID., 
                     infolocs = mov.data.org[,-c(3:6, ncol(mov.data.org))])
mov.traj.df <- ld(mov.traj)
```



### move


```r
# Organize data as a move package format
move.data <- move(x = mov.traj.df$x, y = mov.traj.df$y, 
                  time = mov.traj.df$date, 
                  proj = CRS('+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'),
                  data = mov.traj.df, animal = mov.traj.df$id, sensor = 'GPS')
```

move.data  (moveStack)
Separate individual animals' trajectories 


```r
unstacked <- split(move.data)
```

```
## Error in split(as.data.frame(coordinates(x)), f): argument "f" is missing, with no default
```

```r
jaguar_df <- as(move.data, "data.frame")
```

Reclassifying variables


```r
age <- as.numeric(levels(jaguar_df$age))[jaguar_df$age]
```

```
## Warning: NAs introduced by coercion
```

```r
weight <- as.numeric(levels(jaguar_df$weight))[jaguar_df$weight]
```

```
## Warning: NAs introduced by coercion
```

```r
jaguar_df$id <- as.factor(jaguar_df$individual.local.identifier..ID.)  
jaguar_df$age <- age   
jaguar_df$weight <- weight  
```


Cleaning up columns which will be in excess due to repetition of analysis 


```r
jaguar_df$dx <- NULL
jaguar_df$dy <- NULL
jaguar_df$dist <- NULL
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
jaguar_df$Event_ID  <- NULL
#jaguar_df$timezone <- NULL
# jaguar_df$dt <- NULL
```




```r
jaguar_df$week <- as.numeric(strftime(as.POSIXlt(jaguar_df$timestamp.posix),format="%W"))
jaguar_df$day <- as.numeric(strftime(as.POSIXlt(jaguar_df$timestamp.posix),format="%j"))
jaguar_df$year <- as.numeric(strftime(as.POSIXlt(jaguar_df$timestamp.posix),format="%Y"))
jaguar_df$hour <- as.numeric(strftime(as.POSIXlt(jaguar_df$timestamp.posix),format="%H"))
jaguar_df$min <- as.numeric(strftime(as.POSIXlt(jaguar_df$timestamp.posix),format="%M"))
jaguar_df$time <- jaguar_df$hour + (jaguar_df$min)/60
```




```r
# day  7 to 16  and  night 19 to 4 and  sun riseset 5,6,17,18
#jaguar_df$month=as.numeric(substr(jaguar_df$date,6,7))
jaguar_df$period=ifelse(jaguar_df$hour==7,"day",
                        ifelse(jaguar_df$hour==8,"day",
                        ifelse(jaguar_df$hour==9,"day", 
                        ifelse(jaguar_df$hour==10,"day", 
                        ifelse(jaguar_df$hour==11,"day", 
                        ifelse(jaguar_df$hour==12,"day", 
                        ifelse(jaguar_df$hour==13,"day", 
                        ifelse(jaguar_df$hour==14,"day",
                        ifelse(jaguar_df$hour==15,"day",
                        ifelse(jaguar_df$hour==16,"day",
                        ifelse(jaguar_df$hour==19,"night",
                        ifelse(jaguar_df$hour==20,"night",
                        ifelse(jaguar_df$hour==21,"night",
                        ifelse(jaguar_df$hour==22,"night", 
                        ifelse(jaguar_df$hour==23,"night",
                        ifelse(jaguar_df$hour==0,"night",
                        ifelse(jaguar_df$hour==1,"night",
                        ifelse(jaguar_df$hour==2,"night",
                        ifelse(jaguar_df$hour==3,"night",
                        ifelse(jaguar_df$hour==4,"night","riseset"))))))))))))))))))))
```


# More Data checking and cleaning 

Delete observations where missing lat or long or a timestamp. (There are no missing observations but it is a good practice)


```r
ind<-complete.cases(jaguar_df[,c("y","x","date")])
jaguar_df<-jaguar_df[ind==TRUE,]
```




```r
#" Check order (data should already be ordered)
jaguar_ord <- jaguar_df[order(jaguar_df$id,jaguar_df$date),]
all.equal(jaguar_df,jaguar_ord)
```

```
## [1] TRUE
```


Check for duplicated observations (ones with same lat, long, timestamp,
 and individual identifier). No duplicate observations in this data set


```r
ind2<-jaguar_df %>% select("date","x","y","id") %>% duplicated
```

sum(ind2) # no duplicates


```r
jaguar_df<-jaguar_df[ind2!=TRUE,]
```


Clean suspectly close points!!!  Above 1200 secs or 20 min minimal interval 


```r
excludes <- filter(jaguar_df, dt < 1200)
removed<- anti_join(jaguar_df, excludes)
```

```
## Joining, by = c("x", "y", "date", "dt", "timestamp.posix", "GMTtime", "id", "country", "status", "model", "sex", "age", "weight", "schedule", "project_region", "project_bioveg", "timezone", "week", "day", "year", "hour", "min", "time", "period")
```

```r
jaguar_df <- removed 
```





## Add UTMs and Adjust Posix accordingly with timezone

#### Grouping project regions when they occur within the same UTM zone


```r
# Code with function => crs.convert
```



###    ATLANTIC FOREST WEST    (Project regions 1 and 2)

1) Atlantic Forest W1


```r
AFW1 <- crs.convert(data = subset(jaguar_df,project_region=='Atlantic Forest W1'),
                    crs.input = "+proj=longlat +datum=WGS84",
                    crs.output = "+proj=utm +zone=22K +south +datum=WGS84 +units=m +no_defs",
                    point.names = c("utm_x", "utm_y", "long_x", "lat_y"))
AFW1$date <- as.character(AFW1$date)
AFW1$date <- as.POSIXct(AFW1$date, format ="%Y-%m-%d %H:%M:%S", tz = "Etc/GMT+3") 
```


2) Atlantic Forest W2


```r
AFW2 <- crs.convert(data =subset(jaguar_df,project_region=='Atlantic Forest W2'),
                    crs.input = "+proj=longlat +datum=WGS84",
                    crs.output = "+proj=utm +zone=22K +south +datum=WGS84 +units=m +no_defs",
                    point.names = c("utm_x", "utm_y", "long_x", "lat_y"))
AFW2$date <- as.character(AFW2$date)
AFW2$date <- as.POSIXct(AFW2$date, format ="%Y-%m-%d %H:%M:%S", tz = "Etc/GMT+3") 
```


Atlantic Forest West


```r
AFW=rbind(AFW1,AFW2)
```



### CAATINGA   

3) Caatinga


```r
Caatinga <- crs.convert(data = subset(jaguar_df,project_region=='Caatinga'),
                    crs.input = "+proj=longlat +datum=WGS84",
                    crs.output = "+proj=utm +zone=23L +south +datum=WGS84 +units=m +no_defs",
                    point.names = c("utm_x", "utm_y", "long_x", "lat_y"))
Caatinga$date <- as.character(Caatinga$date)
Caatinga$date <- as.POSIXct(Caatinga$date, format ="%Y-%m-%d %H:%M:%S", tz = "Etc/GMT+3") 
```




### CERRADO   2 Projects (4 and 5)

4)  Cerrado1



```r
Cerrado1<- crs.convert(data = subset(jaguar_df,project_region=='Cerrado1'),
                       crs.input = "+proj=longlat +datum=WGS84",
                       crs.output = "+proj=utm +zone=22K +south +datum=WGS84 +units=m +no_defs",
                       point.names = c("utm_x", "utm_y", "long_x", "lat_y"))
Cerrado1$date <- as.character(Cerrado1$date)
Cerrado1$date <- as.POSIXct(Cerrado1$date, format ="%Y-%m-%d %H:%M:%S", tz = "Etc/GMT+3") 
```


5)  Cerrado2



```r
Cerrado2<- crs.convert(data = subset(jaguar_df,project_region=='Cerrado2'),
                       crs.input = "+proj=longlat +datum=WGS84",
                       crs.output = "+proj=utm +zone=22L +south +datum=WGS84 +units=m +no_defs",
                       point.names = c("utm_x", "utm_y", "long_x", "lat_y"))
Cerrado2$date <- as.character(Cerrado2$date)
Cerrado2$date <- as.POSIXct(Cerrado2$date, format ="%Y-%m-%d %H:%M:%S", tz = "Etc/GMT+3") 
Cerrado=rbind(Cerrado1,Cerrado2)
```



### COSTA RICA

6) Costa Rica




```r
CRica<- crs.convert(data = subset(jaguar_df,project_region=='Costa Rica'),
                       crs.input = "+proj=longlat +datum=WGS84",
                       crs.output = "+proj=utm +zone=16P +north +datum=WGS84 +units=m +no_defs",
                       point.names = c("utm_x", "utm_y", "long_x", "lat_y"))
CRica$date <- as.character(CRica$date)
CRica$date <- as.POSIXct(CRica$date, format ="%Y-%m-%d %H:%M:%S", tz = "Etc/GMT+6") 
```



###   DRY CHACO 
Transforming coordinate to UTM # ids: 16,70 => (most 20K), ids: 76,77 => (20 K); ids: 71,72,73 => (21 K) 



```r
Drych=subset(jaguar_df,project_region=='Dry chaco')
X16=subset(Drych,id=='16')
X70=subset(Drych,id=='70')
X71=subset(Drych,id=='71')
X72=subset(Drych,id=='72')
X73=subset(Drych,id=='73')
X76=subset(Drych,id=='76')
X77=subset(Drych,id=='77')
```

 Drych1 and Drych2


```r
Drych1=rbind(X16,X70,X76,X77)
Drych2=rbind(X71,X72,X73)
```


7) Drych1    


```r
# 76,77  Zone 20 K; 16 & 70 (most is in Zone 20 K too, but a little bit in 21K)
Drych1<- crs.convert(data = Drych1,
                    crs.input = "+proj=longlat +datum=WGS84",
                    crs.output = "+proj=utm +zone=20K +south +datum=WGS84 +units=m +no_defs",
                    point.names = c("utm_x", "utm_y", "long_x", "lat_y"))
Drych1$date <- as.character(Drych1$date)
Drych1$date <- as.POSIXct(Drych1$date, format ="%Y-%m-%d %H:%M:%S", tz = "Etc/GMT+4") 
```



8) Drych2  


```r
# 71,72,73 => (21 K)
Drych2<- crs.convert(data = Drych2,
                     crs.input = "+proj=longlat +datum=WGS84",
                     crs.output = "+proj=utm +zone=21K +south +datum=WGS84 +units=m +no_defs",
                     point.names = c("utm_x", "utm_y", "long_x", "lat_y"))
Drych2$date <- as.character(Drych2$date)
Drych2$date <- as.POSIXct(Drych2$date, format ="%Y-%m-%d %H:%M:%S", tz = "Etc/GMT+4") 
# Dry Chaco Total
Drych=rbind(Drych1,Drych2)
```



###  HUMID CHACO

9) Humid Chaco  
     


```r
Hch<- crs.convert(data = subset(jaguar_df,project_region=='Humid chaco'),
                    crs.input = "+proj=longlat +datum=WGS84",
                    crs.output = "+proj=utm +zone=21K +south +datum=WGS84 +units=m +no_defs",
                    point.names = c("utm_x", "utm_y", "long_x", "lat_y"))
Hch$date <- as.character(Hch$date)
Hch$date <- as.POSIXct(Hch$date, format ="%Y-%m-%d %H:%M:%S", tz = "Etc/GMT+4") 
```



###  FOREST PARAGUAY 

10) Forest Paraguay   



```r
FPy<- crs.convert(data = subset(jaguar_df,project_region=='Forest Paraguay'),
                  crs.input = "+proj=longlat +datum=WGS84",
                  crs.output = "+proj=utm +zone=21J +south +datum=WGS84 +units=m +no_defs",
                  point.names = c("utm_x", "utm_y", "long_x", "lat_y"))
FPy$date <- as.character(FPy$date)
FPy$date <- as.POSIXct(FPy$date, format ="%Y-%m-%d %H:%M:%S", tz = "Etc/GMT+4") 
```




### IGUAZU



```r
Iguazu=subset(jaguar_df,project_region=='Iguazu')
# Transforming coordinate to UTM   # 42,66,80,90 => (21 J);      83 => (22 J) 
X42=subset(Iguazu,id=='42')
X66=subset(Iguazu,id=='66')
X80=subset(Iguazu,id=='80')
X90=subset(Iguazu,id=='90')
X83=subset(Iguazu,id=='83')
```

 Iguazu1 and Iguazu2


```r
Iguazu1=rbind(X42,X66,X80,X90)
Iguazu2=(X83)
```


11) Iguazu1    


```r
# 42,66,80,90 => ( Zone 21 J)
Iguazu1<- crs.convert(data = Iguazu1,
                     crs.input = "+proj=longlat +datum=WGS84",
                     crs.output = "+proj=utm +zone=21J +south +datum=WGS84 +units=m +no_defs",
                     point.names = c("utm_x", "utm_y", "long_x", "lat_y"))
Iguazu1$date <- as.character(Iguazu1$date)
Iguazu1$date <- as.POSIXct(Iguazu1$date, format ="%Y-%m-%d %H:%M:%S", tz = "Etc/GMT+3") 
```



12) Iguazu2    


```r
#  83 => (22 J) 
Iguazu2<- crs.convert(data = Iguazu2,
                      crs.input = "+proj=longlat +datum=WGS84",
                      crs.output = "+proj=utm +zone=22J +south +datum=WGS84 +units=m +no_defs",
                      point.names = c("utm_x", "utm_y", "long_x", "lat_y"))
Iguazu2$date <- as.character(Iguazu2$date)
Iguazu2$date <- as.POSIXct(Iguazu2$date, format ="%Y-%m-%d %H:%M:%S", tz = "Etc/GMT+3") 
```

Iguazu


```r
Iguazu=rbind(Iguazu1,Iguazu2)
```



###  AMAZONIA

13)  Amazonia Mamiraua (Brazil) - Flooded Amazonia


```r
Mamiraua<- crs.convert(data = subset(jaguar_df,project_region=='Mamiraua'),
                  crs.input = "+proj=longlat +datum=WGS84",
                  crs.output = "+proj=utm +zone=20M +south +datum=WGS84 +units=m +no_defs",
                  point.names = c("utm_x", "utm_y", "long_x", "lat_y"))
Mamiraua$date <- as.character(Mamiraua$date)
Mamiraua$date <- as.POSIXct(Mamiraua$date, format ="%Y-%m-%d %H:%M:%S", tz = "Etc/GMT+4") 
```



Dry Amazonia, PA, translocated 
14) IOP Para Amazonia   


```r
# Translocated animal   
iopPA<- crs.convert(data = subset(jaguar_df,id=='24'),
                       crs.input = "+proj=longlat +datum=WGS84",
                       crs.output = "+proj=utm +zone=22M +south +datum=WGS84 +units=m +no_defs",
                       point.names = c("utm_x", "utm_y", "long_x", "lat_y"))
iopPA$date <- as.character(iopPA$date)
iopPA$date <- as.POSIXct(iopPA$date, format ="%Y-%m-%d %H:%M:%S", tz = "Etc/GMT+3") 
```



###   MEXICO

 Greater Lacandona  15) Mexico 


```r
Lacandona<- crs.convert(data = subset(jaguar_df,project_region=='Greater Lacandona'),
                    crs.input = "+proj=longlat +datum=WGS84",
                    crs.output = "+proj=utm +zone=15Q +south +datum=WGS84 +units=m +no_defs",
                    point.names = c("utm_x", "utm_y", "long_x", "lat_y"))
Lacandona$date <- as.character(Lacandona$date)
Lacandona$date <- as.POSIXct(Lacandona$date, format ="%Y-%m-%d %H:%M:%S", tz = "Etc/GMT+6") 
```

 

Mexico East  16) 


```r
MexEast<- crs.convert(data = subset(jaguar_df,project_region=='Mexico East'),
                        crs.input = "+proj=longlat +datum=WGS84",
                        crs.output = "+proj=utm +zone=16Q +south +datum=WGS84 +units=m +no_defs",
                        point.names = c("utm_x", "utm_y", "long_x", "lat_y"))
MexEast$date <- as.character(MexEast$date)
MexEast$date <- as.POSIXct(MexEast$date, format ="%Y-%m-%d %H:%M:%S", tz = "Etc/GMT+6") 
```

 

Mexico Sonora 17) 


```r
Sonora<- crs.convert(data = subset(jaguar_df,project_region=='Mexico Sonora'),
                      crs.input = "+proj=longlat +datum=WGS84",
                      crs.output = "+proj=utm +zone=12R +south +datum=WGS84 +units=m +no_defs",
                      point.names = c("utm_x", "utm_y", "long_x", "lat_y"))
Sonora$date <- as.character(Sonora$date)
Sonora$date <- as.POSIXct(Sonora$date, format ="%Y-%m-%d %H:%M:%S", tz = "Etc/GMT+7") 
```


 Mexico


```r
Mex=rbind(Lacandona,MexEast,Sonora)
```




###  PANTANAL 

18) PantanalTotal Brazil&Paraguay   -  7 Projects in total, all in the same UTM zone



```r
Pantanal<- crs.convert(data = subset(jaguar_df,project_bioveg=='Pantanal'),
                     crs.input = "+proj=longlat +datum=WGS84",
                     crs.output = "+proj=utm +zone=21K +south +datum=WGS84 +units=m +no_defs",
                     point.names = c("utm_x", "utm_y", "long_x", "lat_y"))
Pantanal$date <- as.character(Pantanal$date)
Pantanal$date <- as.POSIXct(Pantanal$date, format ="%Y-%m-%d %H:%M:%S", tz = "Etc/GMT+4") 
```



###  Jaguar Dataframe with UTMs 
head(AFW1);head(AFW2);head(Caatinga);head(Cerrado1);head(Cerrado2);head(CRica);head(Pantanal);head(Drych1);str(Drych)
head(Hch);head(FPy);head(Iguazu1);head(Iguazu2);head(Mamiraua);head(iopPA);head(Lacandona);head(MexEast);head(Sonora)


```r
jaguar_df=rbind(AFW1,AFW2,Caatinga,Cerrado1,Cerrado2,CRica,Pantanal,Drych1,Drych2,Hch,FPy,Iguazu1,Iguazu2,Mamiraua,iopPA,Lacandona,MexEast,Sonora)
```

Need re-order the data again


```r
jaguar_ord <- jaguar_df[order(jaguar_df$id,jaguar_df$date),]
jaguar_df <- jaguar_ord
```


And create a new Event_ID


```r
jaguar_df$Event_ID <- seq.int(nrow(jaguar_df))
```

Save as an additional object


```r
jaguar<-jaguar_df
```

In case we want save and read it as txt


```r
# write.table(jaguar_df,file="../data/jaguar_df.txt",row.names = F,quote=F,col.names=T,sep="\t")
# jaguar <- read.delim(file="../data/jaguar_df.txt")
```





## Creating a track in amt 


```r
# Commom to both RSF and SSF and based on UTM and project regions 
```


### ATLANTIC FOREST WEST  (Project regions 1 and 2)

Atlantic Forest West


```r
AFWtrk <- trk.convert(data = AFW,
                      .x=utm_x,
                      .y=utm_y,
                      .t=date,
                      id=id,
                      project_region=project_region,
                      sex=sex,
                      age=age,
                      weight=weight, 
                      status=status,
                      period=period,
                      long_x=long_x,
                      lat_y=lat_y,
                      crs = CRS("+proj=utm +zone=22K +south +datum=WGS84 +units=m +no_defs"))
```

```
## .t found, creating `track_xyt`.
```

Objects for AFW project regions


```r
 AFW1trk=AFWtrk %>% filter(project_region == "Atlantic Forest W1")
 AFW2trk=AFWtrk %>% filter(project_region == "Atlantic Forest W2")
```

###  CAATINGA 

3) Caatinga   


```r
Caatingatrk <- trk.convert(data = Caatinga,
                      .x=utm_x,
                      .y=utm_y,
                      .t=date,
                      id=id,
                      project_region=project_region,
                      sex=sex,
                      age=age,
                      weight=weight, 
                      status=status,
                      period=period,
                      long_x=long_x,
                      lat_y=lat_y,
                      crs = CRS("+proj=utm +zone=23L +south +datum=WGS84 +units=m +no_defs"))
```

```
## .t found, creating `track_xyt`.
```


###  CERRADO   (2 Projects: 4 and 5)

4) Cerrado1   


```r
Cerrado1trk <- trk.convert(data = Cerrado1,
                           .x=utm_x,
                           .y=utm_y,
                           .t=date,
                           id=id,
                           project_region=project_region,
                           sex=sex,
                           age=age,
                           weight=weight, 
                           status=status,
                           period=period,
                           long_x=long_x,
                           lat_y=lat_y,
                           crs = CRS("+proj=utm +zone=22K +south +datum=WGS84 +units=m +no_defs"))
```

```
## .t found, creating `track_xyt`.
```

5) Cerrado2   


```r
Cerrado2trk <- trk.convert(data = Cerrado2,
                           .x=utm_x,
                           .y=utm_y,
                           .t=date,
                           id=id,
                           project_region=project_region,
                           sex=sex,
                           age=age,
                           weight=weight, 
                           status=status,
                           period=period,
                           long_x=long_x,
                           lat_y=lat_y,
                           crs = CRS("+proj=utm +zone=22L +south +datum=WGS84 +units=m +no_defs"))
```

```
## .t found, creating `track_xyt`.
```

Cerradotrk


```r
Cerradotrk=rbind(Cerrado1trk,Cerrado2trk)
```



###  COSTA RICA  

6) Costa Rica


```r
CRicatrk <- trk.convert(data = CRica,
                          .x=utm_x,
                          .y=utm_y,
                          .t=date,
                          id=id,
                          project_region=project_region,
                          sex=sex,
                          age=age,
                          weight=weight, 
                          status=status,
                          period=period,
                          long_x=long_x,
                          lat_y=lat_y,
                          crs = CRS("+proj=utm +zone=16P +north +datum=WGS84 +units=m +no_defs"))
```

```
## .t found, creating `track_xyt`.
```

### DRY CHACO

7) Drych1


```r
Drych1trk <- trk.convert(data = Drych1,
                        .x=utm_x,
                        .y=utm_y,
                        .t=date,
                        id=id,
                        project_region=project_region,
                        sex=sex,
                        age=age,
                        weight=weight, 
                        status=status,
                        period=period,
                        long_x=long_x,
                        lat_y=lat_y,
                        crs = CRS("+proj=utm +zone=20K +south +datum=WGS84 +units=m +no_defs"))
```

```
## .t found, creating `track_xyt`.
```

8) Drych2


```r
Drych2trk <- trk.convert(data = Drych2,
                         .x=utm_x,
                         .y=utm_y,
                         .t=date,
                         id=id,
                         project_region=project_region,
                         sex=sex,
                         age=age,
                         weight=weight, 
                         status=status,
                         period=period,
                         long_x=long_x,
                         lat_y=lat_y,
                         crs = CRS("+proj=utm +zone=21K +south +datum=WGS84 +units=m +no_defs"))
```

```
## .t found, creating `track_xyt`.
```


 DRY CHACO trk => Drychtrk


```r
Drychtrk=rbind(Drych1trk,Drych2trk)
```



### HUMID CHACO  

9) Humid Chaco       


```r
Hchtrk <- trk.convert(data = Hch,
                         .x=utm_x,
                         .y=utm_y,
                         .t=date,
                         id=id,
                         project_region=project_region,
                         sex=sex,
                         age=age,
                         weight=weight, 
                         status=status,
                         period=period,
                         long_x=long_x,
                         lat_y=lat_y,
                         crs = CRS("+proj=utm +zone=21K +south +datum=WGS84 +units=m +no_defs"))
```

```
## .t found, creating `track_xyt`.
```



###   FOREST PARAGUAY 
10) Forest Paraguay   


```r
FPytrk <- trk.convert(data = FPy,
                      .x=utm_x,
                      .y=utm_y,
                      .t=date,
                      id=id,
                      project_region=project_region,
                      sex=sex,
                      age=age,
                      weight=weight, 
                      status=status,
                      period=period,
                      long_x=long_x,
                      lat_y=lat_y,
                      crs = CRS("+proj=utm +zone=21J +south +datum=WGS84 +units=m +no_defs"))
```

```
## .t found, creating `track_xyt`.
```


  
### IGUAZU
	   
11) Iguazu1


```r
Iguazu1trk <- trk.convert(data = Iguazu1,
                      .x=utm_x,
                      .y=utm_y,
                      .t=date,
                      id=id,
                      project_region=project_region,
                      sex=sex,
                      age=age,
                      weight=weight, 
                      status=status,
                      period=period,
                      long_x=long_x,
                      lat_y=lat_y,
                      crs = CRS("+proj=utm +zone=21J +south +datum=WGS84 +units=m +no_defs"))
```

```
## .t found, creating `track_xyt`.
```


12) Iguazu2


```r
Iguazu2trk <- trk.convert(data = Iguazu2,
                          .x=utm_x,
                          .y=utm_y,
                          .t=date,
                          id=id,
                          project_region=project_region,
                          sex=sex,
                          age=age,
                          weight=weight, 
                          status=status,
                          period=period,
                          long_x=long_x,
                          lat_y=lat_y,
                          crs = CRS("+proj=utm +zone=22J +south +datum=WGS84 +units=m +no_defs"))
```

```
## .t found, creating `track_xyt`.
```


 Iguazutrk   


```r
Iguazutrk=rbind(Iguazu1trk,Iguazu2trk)
```



###   AMAZONIA   
13) Flooded  Amazonia, Mamiraua (Brazil)   


```r
Mamirauatrk <- trk.convert(data = Mamiraua,
                          .x=utm_x,
                          .y=utm_y,
                          .t=date,
                          id=id,
                          project_region=project_region,
                          sex=sex,
                          age=age,
                          weight=weight, 
                          status=status,
                          period=period,
                          long_x=long_x,
                          lat_y=lat_y,
                          crs = CRS("+proj=utm +zone=20M +south +datum=WGS84 +units=m +no_defs"))
```

```
## .t found, creating `track_xyt`.
```



14   Dry/ East Amazonia, IOP PA, translocated


```r
iopPAtrk <- trk.convert(data = iopPA,
                           .x=utm_x,
                           .y=utm_y,
                           .t=date,
                           id=id,
                           project_region=project_region,
                           sex=sex,
                           age=age,
                           weight=weight, 
                           status=status,
                           period=period,
                           long_x=long_x,
                           lat_y=lat_y,
                           crs = CRS("+proj=utm +zone=22M +south +datum=WGS84 +units=m +no_defs"))
```

```
## .t found, creating `track_xyt`.
```



### Greater Lacandona, Mexico

15) Lacandona


```r
Lacandonatrk <- trk.convert(data = Lacandona,
                        .x=utm_x,
                        .y=utm_y,
                        .t=date,
                        id=id,
                        project_region=project_region,
                        sex=sex,
                        age=age,
                        weight=weight, 
                        status=status,
                        period=period,
                        long_x=long_x,
                        lat_y=lat_y,
                        crs = CRS("+proj=utm +zone=15Q +south +datum=WGS84 +units=m +no_defs"))
```

```
## .t found, creating `track_xyt`.
```



16) MexEast, Mexico


```r
MexEasttrk <- trk.convert(data = MexEast,
                            .x=utm_x,
                            .y=utm_y,
                            .t=date,
                            id=id,
                            project_region=project_region,
                            sex=sex,
                            age=age,
                            weight=weight, 
                            status=status,
                            period=period,
                            long_x=long_x,
                            lat_y=lat_y,
                            crs = CRS("+proj=utm +zone=16Q +south +datum=WGS84 +units=m +no_defs"))
```

```
## .t found, creating `track_xyt`.
```



### Sonora, Mexico

17)  Sonora


```r
Sonoratrk <- trk.convert(data = Sonora,
                          .x=utm_x,
                          .y=utm_y,
                          .t=date,
                          id=id,
                          project_region=project_region,
                          sex=sex,
                          age=age,
                          weight=weight, 
                          status=status,
                          period=period,
                          long_x=long_x,
                          lat_y=lat_y,
                          crs = CRS("+proj=utm +zone=12R +south +datum=WGS84 +units=m +no_defs"))
```

```
## .t found, creating `track_xyt`.
```



###   PANTANAL

18) PantanalTotal Brazil&Paraguay  -  7 Projects in total, all in the same UTM zone


```r
Pantanaltrk <- trk.convert(data = Pantanal,
                         .x=utm_x,
                         .y=utm_y,
                         .t=date,
                         id=id,
                         project_region=project_region,
                         sex=sex,
                         age=age,
                         weight=weight, 
                         status=status,
                         period=period,
                         long_x=long_x,
                         lat_y=lat_y,
                         crs = CRS("+proj=utm +zone=21K +south +datum=WGS84 +units=m +no_defs"))
```

```
## .t found, creating `track_xyt`.
```

 
Objects for Pantanal project regions


```r
 Oncafaritrk=Pantanaltrk %>% filter(project_region == "Oncafari")
 Paraguaytrk=Pantanaltrk %>% filter(project_region == "Pantanal Paraguay")
 Panthera1trk=Pantanaltrk %>% filter(project_region == "Panthera1")
 Panthera2trk=Pantanaltrk %>% filter(project_region == "Panthera2")
 RioNegrotrk=Pantanaltrk %>% filter(project_region == "Rio Negro")
 SaoBentotrk=Pantanaltrk %>% filter(project_region == "Sao Bento")
 Taiamatrk=Pantanaltrk %>% filter(project_region == "Taiama")
```



###    ALL  JAGUARS  trk =>     jaguartrk  
		   
All Project regions trk


```r
jaguartrk=rbind(AFW1trk,AFW2trk,Caatingatrk,Cerrado1trk,Cerrado2trk,CRicatrk,Drychtrk, Hchtrk,FPytrk,Iguazutrk,
                Mamirauatrk,iopPAtrk,Lacandonatrk, MexEasttrk, Sonoratrk,Oncafaritrk,Paraguaytrk,Panthera1trk,
                Panthera2trk,RioNegrotrk,SaoBentotrk,Taiamatrk) 
```



Dataframe for each individual (to run if we need any specic individual)


```r
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
```


Produce an output using ezknitr


```r
ezspin(file = "JaguarDataPrep.R", out_dir = "reports",
       params = list("DATASET_NAME" = "jaguar.dat"), 
       keep_html = TRUE, keep_rmd = TRUE)
```

```
## ezspin output in
## D:\Documents\GitHub\jaguar-codes\R\reports
```

```r
open_output_dir()
```

```
## Opening D:\Documents\GitHub\jaguar-codes\R\reports
```

Other info


```r
sessionInfo()	  
```

```
## R version 3.5.2 (2018-12-20)
## Platform: x86_64-w64-mingw32/x64 (64-bit)
## Running under: Windows >= 8 x64 (build 9200)
## 
## Matrix products: default
## 
## locale:
## [1] LC_COLLATE=English_United States.1252 
## [2] LC_CTYPE=English_United States.1252   
## [3] LC_MONETARY=English_United States.1252
## [4] LC_NUMERIC=C                          
## [5] LC_TIME=English_United States.1252    
## 
## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
##  [1] bindrcpp_0.2.2      ezknitr_0.6         knitr_1.21         
##  [4] caTools_1.17.1.1    circular_0.4-93     tibble_1.4.2       
##  [7] lubridate_1.7.4     readr_1.2.1         dplyr_0.7.8        
## [10] RCurl_1.95-4.11     bitops_1.0-6        leaflet_2.0.2      
## [13] lattice_0.20-38     rgl_0.99.16         ggmap_2.6.1        
## [16] ggplot2_3.1.0       colorspace_1.4-0    maptools_0.9-4     
## [19] amt_0.0.5.0         adehabitatLT_0.3.23 CircStats_0.2-6    
## [22] boot_1.3-20         MASS_7.3-51.1       adehabitatMA_0.3.12
## [25] ade4_1.7-13         move_3.1.0          rgdal_1.3-6        
## [28] raster_2.8-4        sp_1.3-1            geosphere_1.5-7    
## [31] install.load_1.2.1 
## 
## loaded via a namespace (and not attached):
##  [1] webshot_0.5.1           httr_1.4.0             
##  [3] tools_3.5.2             utf8_1.1.4             
##  [5] R6_2.3.0                lazyeval_0.2.1         
##  [7] manipulateWidget_0.10.0 withr_2.1.2            
##  [9] tidyselect_0.2.5        compiler_3.5.2         
## [11] cli_1.0.1               xml2_1.2.0             
## [13] scales_1.0.0            mvtnorm_1.0-8          
## [15] stringr_1.3.1           digest_0.6.18          
## [17] foreign_0.8-71          R.utils_2.7.0          
## [19] jpeg_0.1-8              pkgconfig_2.0.2        
## [21] htmltools_0.3.6         highr_0.7              
## [23] maps_3.3.0              htmlwidgets_1.3        
## [25] rlang_0.3.0.1           rstudioapi_0.8         
## [27] shiny_1.2.0             bindr_0.1.1            
## [29] jsonlite_1.5            crosstalk_1.0.0        
## [31] R.oo_1.22.0             magrittr_1.5           
## [33] Matrix_1.2-15           Rcpp_1.0.0             
## [35] munsell_0.5.0           fansi_0.4.0            
## [37] proto_1.0.0             R.methodsS3_1.7.1      
## [39] stringi_1.2.4           yaml_2.2.0             
## [41] plyr_1.8.4              grid_3.5.2             
## [43] parallel_3.5.2          promises_1.0.1         
## [45] crayon_1.3.4            miniUI_0.1.1.1         
## [47] splines_3.5.2           mapproj_1.2.6          
## [49] hms_0.4.2               pillar_1.3.0           
## [51] rjson_0.2.20            markdown_0.8           
## [53] reshape2_1.4.3          codetools_0.2-15       
## [55] glue_1.3.0              evaluate_0.12          
## [57] png_0.1-7               httpuv_1.4.5           
## [59] RgoogleMaps_1.4.3       gtable_0.2.0           
## [61] purrr_0.2.5             tidyr_0.8.2            
## [63] assertthat_0.2.0        xfun_0.4               
## [65] mime_0.6                xtable_1.8-3           
## [67] later_0.7.5             survival_2.43-3        
## [69] memoise_1.1.0
```

```r
proc.time()
```

```
##    user  system elapsed 
## 1055.92   37.62 5358.71
```

