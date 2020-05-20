
#  **Jaguar Data Preparation**

#### *Alan E. de Barros, Bernardo Niebuhr, Vanesa Bejarano, Julia Oshima,Claudia Kanda, Milton Ribeiro, Ronaldo Morato,Paulo Prado*
date: "March, 13 2019"
##### Scripts adapted from Bernardo Niebuhr data preparation, and Luca Borger and John Fieberg's lectures.



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
install.load::install_load("colorspace", "rgl", "lattice", "leaflet","ggplot2") # Visualization packages
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
mov.data.org$timestamp.posix <- as.POSIXct(date.time,format = "%m/%d/%Y %H:%M", tz = 'UTC')
mov.data.org$GMTtime <- mov.data.org$timestamp.posix
```



## Get local time
 
A column to represent the local timezone (already with the - signal) has been included to then multiply the timestamp and get the difference:


```r
mov.data.org$local_time <- mov.data.org$timestamp.posix + mov.data.org$timezone*60*60
mov.data.org$timestamp.posix <- mov.data.org$local_time 
```

Now all the (timestamp.posix)'s calculations are based on local time
#### *Note: UTC has also been assigned to the local time here. But the correct timezones of each region are assigned below*


### Adjusting for Movement Packages
#### Creating a track in amt (with lat long for now)


```r
trk1 <- mk_track(mov.data.org,.x=location.long, .y=location.lat, .t=timestamp.posix,
                id=individual.local.identifier..ID.,sex=sex,age=age, weight=weight, 
                status=status,model=model,project_region=project_region, project_bioveg=project_bioveg,
                country=country,GMTtime=GMTtime,local_time=local_time,crs = CRS("+init=epsg:4326"))
```

```
## .t found, creating `track_xyt`.
```

```r
trk1 <- trk1 %>% arrange(id)
jaguar_df <- as.data.frame(trk1)
```


Reclassifying variables


```r
jaguar_df$idf <- as.factor(jaguar_df$id)
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
jaguar_df$age <- age   
jaguar_df$weight <- weight  
jaguar_df<-rename(jaguar_df, x = x_)
jaguar_df<-rename(jaguar_df, y = y_)
jaguar_df<-rename(jaguar_df, date = t_)
```



### More Data checking and cleaning 

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
ind2<-duplicated(jaguar_df[,c("date","x","y","id")])
```

sum(ind2) # no duplicates


```r
jaguar_df<-jaguar_df[ind2!=TRUE,]
```


Clean suspectly close points!!!  Above 1200 secs or 20 min minimal interval 


```r
#  Lets add on the time difference to each obs to select and exclude them after
jaguar_df<-jaguar_df %>% group_by(id) %>% mutate(dt = date - lag(date, default = NA))
excludes <- filter(jaguar_df, dt < 1200)
removed<- anti_join(jaguar_df, excludes)
```

```
## Joining, by = c("x", "y", "date", "id", "sex", "age", "weight", "status", "model", "project_region", "project_bioveg", "country", "GMTtime", "local_time", "idf", "dt")
```

```r
jaguar_df <- removed
jaguar_df <- as.data.frame(jaguar_df)
```

Create a new Event_ID


```r
jaguar_ord <- jaguar_df[order(jaguar_df$id,jaguar_df$date),]
jaguar_df <- jaguar_ord
jaguar_df$Event_ID <- seq.int(nrow(jaguar_df))
```



#### **Need still add some individual cleaning of outliers (based on ctmm and large timelags)!!!**


### Add UTMs and POSIX timezone
#### Grouping project regions when they occur within the same UTM and time zone


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
J16=subset(Drych,id=='16')
J70=subset(Drych,id=='70')
J71=subset(Drych,id=='71')
J72=subset(Drych,id=='72')
J73=subset(Drych,id=='73')
J76=subset(Drych,id=='76')
J77=subset(Drych,id=='77')
```

 Drych1 and Drych2


```r
Drych1=rbind(J16,J70,J76,J77)
Drych2=rbind(J71,J72,J73)
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
J42=subset(Iguazu,id=='42')
J66=subset(Iguazu,id=='66')
J80=subset(Iguazu,id=='80')
J90=subset(Iguazu,id=='90')
J83=subset(Iguazu,id=='83')
```

 Iguazu1 and Iguazu2


```r
Iguazu1=rbind(J42,J66,J80,J90)
Iguazu2=(J83)
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

#### *Important Note: A single dataframe can NOT hold multiple POSIX tzs*
This will assign the POSIXct of the first region (AFW in this case). 
Hence it is required to store the regional POSIXct in multiple dataframes (or within a tibble for example)


```r
jaguar_df=rbind(AFW1,AFW2,Caatinga,Cerrado1,Cerrado2,CRica,Pantanal,
          Drych1,Drych2,Hch,FPy,Iguazu1,Iguazu2,Mamiraua,iopPA,Lacandona,MexEast,Sonora)
```

Need re-order the data again


```r
jaguar_ord <- jaguar_df[order(jaguar_df$id,jaguar_df$date),]
jaguar_df <- jaguar_ord
```

We assign UTC again just to store a full dataframe with UTMS but bellow we use regional POSIXct to create trks (amt package).


```r
jaguar_df$date <- jaguar_df$local_time
```

In case we want save and read it as txt


```r
# write.table(jaguar_df,file="../data/jaguar_df.txt",row.names = F,quote=F,col.names=T,sep="\t")
# jaguar <- read.delim(file="../data/jaguar_df.txt")
```





## Creating tracks in amt 
#### Based on UTM, project regions and timezones (Commom to both RSF and SSF)

### ATLANTIC FOREST WEST  (Project regions 1 and 2)

Atlantic Forest West


```r
AFWtrk <- trk.convert(data = AFW,
                      .x=utm_x,
                      .y=utm_y,
                      .t=date,
                      id=id,
                      dt=dt,
                      project_region=project_region,
                      sex=sex,
                      age=age,
                      weight=weight, 
                      status=status,
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
                      dt=dt,
                      project_region=project_region,
                      sex=sex,
                      age=age,
                      weight=weight, 
                      status=status,
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
                           dt=dt,
                           project_region=project_region,
                           sex=sex,
                           age=age,
                           weight=weight, 
                           status=status,
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
                           dt=dt,
                           project_region=project_region,
                           sex=sex,
                           age=age,
                           weight=weight, 
                           status=status,
                           long_x=long_x,
                           lat_y=lat_y,
                           crs = CRS("+proj=utm +zone=22L +south +datum=WGS84 +units=m +no_defs"))
```

```
## .t found, creating `track_xyt`.
```

```r
# Cerradotrk
# Cerradotrk=rbind(Cerrado1trk,Cerrado2trk) # This apparently would assign only the first UTM zone for both
```



###  COSTA RICA  

6) Costa Rica


```r
CRicatrk <- trk.convert(data = CRica,
                          .x=utm_x,
                          .y=utm_y,
                          .t=date,
                          id=id,
                          dt=dt,
                          project_region=project_region,
                          sex=sex,
                          age=age,
                          weight=weight, 
                          status=status,
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
                        dt=dt,
                        project_region=project_region,
                        sex=sex,
                        age=age,
                        weight=weight, 
                        status=status,
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
                         dt=dt,
                         project_region=project_region,
                         sex=sex,
                         age=age,
                         weight=weight, 
                         status=status,
                         long_x=long_x,
                         lat_y=lat_y,
                         crs = CRS("+proj=utm +zone=21K +south +datum=WGS84 +units=m +no_defs"))
```

```
## .t found, creating `track_xyt`.
```




```r
#  DRY CHACO trk => Drychtrk
# Drychtrk=rbind(Drych1trk,Drych2trk) # This apparently would assign only the first UTM zone for both
```



### HUMID CHACO  

9) Humid Chaco       


```r
Hchtrk <- trk.convert(data = Hch,
                         .x=utm_x,
                         .y=utm_y,
                         .t=date,
                         id=id,
                         dt=dt,
                         project_region=project_region,
                         sex=sex,
                         age=age,
                         weight=weight, 
                         status=status,
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
                      dt=dt,
                      project_region=project_region,
                      sex=sex,
                      age=age,
                      weight=weight, 
                      status=status,
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
                      dt=dt,
                      project_region=project_region,
                      sex=sex,
                      age=age,
                      weight=weight, 
                      status=status,
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
                          dt=dt,
                          project_region=project_region,
                          sex=sex,
                          age=age,
                          weight=weight, 
                          status=status,
                          long_x=long_x,
                          lat_y=lat_y,
                          crs = CRS("+proj=utm +zone=22J +south +datum=WGS84 +units=m +no_defs"))
```

```
## .t found, creating `track_xyt`.
```




```r
# Iguazutrk   
# Iguazutrk=rbind(Iguazu1trk,Iguazu2trk) # This apparently would assign only the first UTM zone for both
```



###   AMAZONIA   
13) Flooded  Amazonia, Mamiraua (Brazil)   


```r
Mamirauatrk <- trk.convert(data = Mamiraua,
                          .x=utm_x,
                          .y=utm_y,
                          .t=date,
                          id=id,
                          dt=dt,
                          project_region=project_region,
                          sex=sex,
                          age=age,
                          weight=weight, 
                          status=status,
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
                           dt=dt,
                           project_region=project_region,
                           sex=sex,
                           age=age,
                           weight=weight, 
                           status=status,
                           long_x=long_x,
                           lat_y=lat_y,
                           crs = CRS("+proj=utm +zone=22M +south +datum=WGS84 +units=m +no_defs"))
```

```
## .t found, creating `track_xyt`.
```



###   MEXICO

Greater Lacandona
15) Lacandona


```r
Lacandonatrk <- trk.convert(data = Lacandona,
                        .x=utm_x,
                        .y=utm_y,
                        .t=date,
                        id=id,
                        dt=dt,
                        project_region=project_region,
                        sex=sex,
                        age=age,
                        weight=weight, 
                        status=status,
                        long_x=long_x,
                        lat_y=lat_y,
                        crs = CRS("+proj=utm +zone=15Q +south +datum=WGS84 +units=m +no_defs"))
```

```
## .t found, creating `track_xyt`.
```



16) MexEast


```r
MexEasttrk <- trk.convert(data = MexEast,
                            .x=utm_x,
                            .y=utm_y,
                            .t=date,
                            id=id,
                            dt=dt,
                            project_region=project_region,
                            sex=sex,
                            age=age,
                            weight=weight, 
                            status=status,
                            long_x=long_x,
                            lat_y=lat_y,
                            crs = CRS("+proj=utm +zone=16Q +south +datum=WGS84 +units=m +no_defs"))
```

```
## .t found, creating `track_xyt`.
```



17)  Sonora


```r
Sonoratrk <- trk.convert(data = Sonora,
                          .x=utm_x,
                          .y=utm_y,
                          .t=date,
                          id=id,
                          dt=dt,
                          project_region=project_region,
                          sex=sex,
                          age=age,
                          weight=weight, 
                          status=status,
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
                         dt=dt,
                         project_region=project_region,
                         sex=sex,
                         age=age,
                         weight=weight, 
                         status=status,
                         long_x=long_x,
                         lat_y=lat_y,
                         crs = CRS("+proj=utm +zone=21K +south +datum=WGS84 +units=m +no_defs"))
```

```
## .t found, creating `track_xyt`.
```

 
Objects for Pantanal project regions (trk and dataframes)


```r
 Oncafaritrk=Pantanaltrk %>% filter(project_region == "Oncafari")
 Oncafari=Pantanal %>% filter(project_region == "Oncafari")
 Paraguaytrk=Pantanaltrk %>% filter(project_region == "Pantanal Paraguay")
 Paraguay=Pantanal %>% filter(project_region == "Pantanal Paraguay")
 Panthera1trk=Pantanaltrk %>% filter(project_region == "Panthera1")
 Panthera1=Pantanal %>% filter(project_region == "Panthera1")
 Panthera2trk=Pantanaltrk %>% filter(project_region == "Panthera2")
 Panthera2=Pantanal %>% filter(project_region == "Panthera2")
 RioNegrotrk=Pantanaltrk %>% filter(project_region == "Rio Negro")
 RioNegro=Pantanal %>% filter(project_region == "Rio Negro")
 SaoBentotrk=Pantanaltrk %>% filter(project_region == "Sao Bento")
 SaoBento=Pantanal %>% filter(project_region == "Sao Bento")
 Taiamatrk=Pantanaltrk %>% filter(project_region == "Taiama")
 Taiama=Pantanal %>% filter(project_region == "Taiama")
```




Dataframe and trk for each individual with adjusted timezones (to run if we need any specic individual)


```r
J1=subset(Hch,id=='1');J1$idloc <- seq.int(nrow(J1)) ####### Hch
J1trk<-Hchtrk %>% filter(id=="1")
J2=subset(FPy,id=='2');J2$idloc <- seq.int(nrow(J2)) ######### FPy
J2trk<-FPytrk %>% filter(id=="2")
J3=subset(Hch,id=='3');J3$idloc <- seq.int(nrow(J3)) ####### Hch
J3trk<-Hchtrk %>% filter(id=="3")
J4=subset(Hch,id=='4');J4$idloc <- seq.int(nrow(J4))   ####### Hch
J4trk<-Hchtrk %>% filter(id=="4")
J5=subset(Hch,id=='5');J5$idloc <- seq.int(nrow(J5))   ####### Hch
J5trk<-Hchtrk %>% filter(id=="5")
J6=subset(Hch,id=='6');J6$idloc <- seq.int(nrow(J6))   ####### Hch
J6trk<-Hchtrk %>% filter(id=="6")
J7=subset(Hch,id=='7');J7$idloc <- seq.int(nrow(J7))   ####### Hch
J7trk<-Hchtrk %>% filter(id=="7")
J8=subset(FPy,id=='8');J8$idloc <- seq.int(nrow(J8))   ######### FPy
J8trk<-FPytrk %>% filter(id=="8")
J9=subset(Hch,id=='9');J9$idloc <- seq.int(nrow(J9))   ####### Hch
J9trk<-Hchtrk %>% filter(id=="9")
J10=subset(Hch,id=='10');J10$idloc <- seq.int(nrow(J10)) ####### Hch
J10trk<-Hchtrk %>% filter(id=="10")
J11=subset(Hch,id=='11');J11$idloc <- seq.int(nrow(J11)) ####### Hch
J11trk<-Hchtrk %>% filter(id=="11")
J12=subset(Pantanal,id=='12');J12$idloc <- seq.int(nrow(J12)) ########## Pantanal
J12trk<-Pantanaltrk %>% filter(id=="12")
J13=subset(Pantanal,id=='13');J13$idloc <- seq.int(nrow(J13)) ########## Pantanal
J13trk<-Pantanaltrk %>% filter(id=="13")
J14=subset(Pantanal,id=='14');J14$idloc <- seq.int(nrow(J14)) ########## Pantanal
J14trk<-Pantanaltrk %>% filter(id=="14")
J15=subset(Pantanal,id=='15');J15$idloc <- seq.int(nrow(J15)) ########## Pantanal
J15trk<-Pantanaltrk %>% filter(id=="15")
J16=subset(Drych1,id=='16');J16$idloc <- seq.int(nrow(J16)) ######### Drych1
J16trk<-Drych1trk %>% filter(id=="16")
J17=subset(Cerrado1,id=='17');J17$idloc <- seq.int(nrow(J17)) ############### Cerrado1
J17trk<-Cerrado1trk %>% filter(id=="17")
J18=subset(Pantanal,id=='18');J18$idloc <- seq.int(nrow(J18)) ########## Pantanal
J18trk<-Pantanaltrk %>% filter(id=="18")
J19=subset(Pantanal,id=='19');J19$idloc <- seq.int(nrow(J19)) ########## Pantanal
J19trk<-Pantanaltrk %>% filter(id=="19")
J20=subset(Caatinga,id=='20');J20$idloc <- seq.int(nrow(J20)) ################# Caatinga 
J20trk<-Caatingatrk %>% filter(id=="20")
J21=subset(FPy,id=='21');J21$idloc <- seq.int(nrow(J21))  ######### FPy
J21trk<-FPytrk %>% filter(id=="21")
J22=subset(Pantanal,id=='22');J22$idloc <- seq.int(nrow(J22)) ########## Pantanal
J22trk<-Pantanaltrk %>% filter(id=="22")
J23=subset(Pantanal,id=='23');J23$idloc <- seq.int(nrow(J23)) ########## Pantanal
J23trk<-Pantanaltrk %>% filter(id=="23")
J24=subset(iopPA,id=='24');J24$idloc <- seq.int(nrow(J24)) ###### iopPA
J24trk<-iopPAtrk %>% filter(id=="24")
J25=subset(Pantanal,id=='25');J25$idloc <- seq.int(nrow(J25)) ########## Pantanal
J25trk<-Pantanaltrk %>% filter(id=="25")
J26=subset(CRica,id=='26');J26$idloc <- seq.int(nrow(J26))  ######### CRica
J26trk<-CRicatrk %>% filter(id=="26")
J27=subset(Pantanal,id=='27');J27$idloc <- seq.int(nrow(J27)) ########## Pantanal
J27trk<-Pantanaltrk %>% filter(id=="27")
J28=subset(Pantanal,id=='28');J28$idloc <- seq.int(nrow(J28)) ########## Pantanal
J28trk<-Pantanaltrk %>% filter(id=="28")
J29=subset(Pantanal,id=='29');J29$idloc <- seq.int(nrow(J29)) ########## Pantanal
J29trk<-Pantanaltrk %>% filter(id=="29")
J30=subset(Pantanal,id=='30');J30$idloc <- seq.int(nrow(J30)) ########## Pantanal
J30trk<-Pantanaltrk %>% filter(id=="30")
J31=subset(Pantanal,id=='31');J31$idloc <- seq.int(nrow(J31)) ########## Pantanal
J31trk<-Pantanaltrk %>% filter(id=="31")
J32=subset(Pantanal,id=='32');J32$idloc <- seq.int(nrow(J32)) ########## Pantanal
J32trk<-Pantanaltrk %>% filter(id=="32")
J33=subset(Pantanal,id=='33');J33$idloc <- seq.int(nrow(J33)) ########## Pantanal
J33trk<-Pantanaltrk %>% filter(id=="33")
J34=subset(AFW1,id=='34');J34$idloc <- seq.int(nrow(J34))   ############### AFW1
J34trk<-AFW1trk %>% filter(id=="34")
J35=subset(AFW1,id=='35');J35$idloc <- seq.int(nrow(J35))   ############### AFW1
J35trk<-AFW1trk %>% filter(id=="35")
J36=subset(AFW1,id=='36');J36$idloc <- seq.int(nrow(J36))   ############### AFW1
J36trk<-AFW1trk %>% filter(id=="36")
J37=subset(AFW1,id=='37');J37$idloc <- seq.int(nrow(J37))   ############### AFW1
J37trk<-AFW1trk %>% filter(id=="37")
J38=subset(AFW1,id=='38');J38$idloc <- seq.int(nrow(J38))   ############### AFW1
J38trk<-AFW1trk %>% filter(id=="38")
J39=subset(AFW2,id=='39');J39$idloc <- seq.int(nrow(J39))   ###############  AFW2
J39trk<-AFW2trk %>% filter(id=="39")
J40=subset(AFW2,id=='40');J40$idloc <- seq.int(nrow(J40))   ###############  AFW2
J40trk<-AFW2trk %>% filter(id=="40")
J41=subset(Pantanal,id=='41');J41$idloc <- seq.int(nrow(J41)) ########## Pantanal
J41trk<-Pantanaltrk %>% filter(id=="41")
J42=subset(Iguazu1,id=='42');J42$idloc <- seq.int(nrow(J42)) ######## Iguazu1
J42trk<-Iguazu1trk %>% filter(id=="42")
J43=subset(Sonora,id=='43');J43$idloc <- seq.int(nrow(J43)) ######## Sonora
J43trk<-Sonoratrk %>% filter(id=="43")
J44=subset(Lacandona,id=='44');J44$idloc <- seq.int(nrow(J44)) ######## Lacandona
J44trk<-Lacandonatrk %>% filter(id=="44")
J45=subset(Lacandona,id=='45');J45$idloc <- seq.int(nrow(J45)) ######## Lacandona
J45trk<-Lacandonatrk %>% filter(id=="45")
J46=subset(Lacandona,id=='46');J46$idloc <- seq.int(nrow(J46))  ######## Lacandona
J46trk<-Lacandonatrk %>% filter(id=="46")
J47=subset(Lacandona,id=='47');J47$idloc <- seq.int(nrow(J47))  ######## Lacandona
J47trk<-Lacandonatrk %>% filter(id=="47")
J48=subset(Lacandona,id=='48');J48$idloc <- seq.int(nrow(J48))  ######## Lacandona
J48trk<-Lacandonatrk %>% filter(id=="48")
J49=subset(MexEast,id=='49');J49$idloc <- seq.int(nrow(J49))  ######### MexEast
J49trk<-MexEasttrk %>% filter(id=="49")
J50=subset(Caatinga,id=='50');J50$idloc <- seq.int(nrow(J50))  ############### Caatinga
J50trk<-Caatingatrk %>% filter(id=="50")
J51=subset(Pantanal,id=='51');J51$idloc <- seq.int(nrow(J51)) ########## Pantanal
J51trk<-Pantanaltrk %>% filter(id=="51")
J52=subset(Pantanal,id=='52');J52$idloc <- seq.int(nrow(J52)) ########## Pantanal
J52trk<-Pantanaltrk %>% filter(id=="52")
J53=subset(Pantanal,id=='53');J53$idloc <- seq.int(nrow(J53)) ########## Pantanal
J53trk<-Pantanaltrk %>% filter(id=="53")
J54=subset(Pantanal,id=='54');J54$idloc <- seq.int(nrow(J54)) ########## Pantanal
J54trk<-Pantanaltrk %>% filter(id=="54")
J55=subset(Pantanal,id=='55');J55$idloc <- seq.int(nrow(J55)) ########## Pantanal
J55trk<-Pantanaltrk %>% filter(id=="55")
J56=subset(Pantanal,id=='56');J56$idloc <- seq.int(nrow(J56)) ########## Pantanal
J56trk<-Pantanaltrk %>% filter(id=="56")
J57=subset(Pantanal,id=='57');J57$idloc <- seq.int(nrow(J57)) ########## Pantanal
J57trk<-Pantanaltrk %>% filter(id=="57")
J58=subset(AFW1,id=='58');J58$idloc <- seq.int(nrow(J58))   ###############  AFW1
J58trk<-AFW1trk %>% filter(id=="58")
J59=subset(Pantanal,id=='59');J59$idloc <- seq.int(nrow(J59)) ########## Pantanal
J59trk<-Pantanaltrk %>% filter(id=="59")
J60=subset(Pantanal,id=='60');J60$idloc <- seq.int(nrow(J60)) ########## Pantanal
J60trk<-Pantanaltrk %>% filter(id=="60")
J61=subset(Pantanal,id=='61');J61$idloc <- seq.int(nrow(J61)) ########## Pantanal
J61trk<-Pantanaltrk %>% filter(id=="61")
J62=subset(AFW1,id=='62');J62$idloc <- seq.int(nrow(J62))   ###############  AFW1
J62trk<-AFW1trk %>% filter(id=="62")
J63=subset(AFW2,id=='63');J63$idloc <- seq.int(nrow(J63))   ###############  AFW2
J63trk<-AFW2trk %>% filter(id=="63")
J64=subset(Sonora,id=='64');J64$idloc <- seq.int(nrow(J64)) ######## Sonora
J64trk<-Sonoratrk %>% filter(id=="64")
J65=subset(Cerrado1,id=='65');J65$idloc <- seq.int(nrow(J65))  ########## Cerrado1
J65trk<-Cerrado1trk %>% filter(id=="65")
J66=subset(Iguazu1,id=='66');J66$idloc <- seq.int(nrow(J66))  ########## Iguazu1
J66trk<-Iguazu1trk %>% filter(id=="66")
J67=subset(Cerrado1,id=='67');J67$idloc <- seq.int(nrow(J67))  ########## Cerrado1
J67trk<-Cerrado1trk %>% filter(id=="67")
J68=subset(Pantanal,id=='68');J68$idloc <- seq.int(nrow(J68)) ########## Pantanal
J68trk<-Pantanaltrk %>% filter(id=="68")
J69=subset(Pantanal,id=='69');J69$idloc <- seq.int(nrow(J69)) ########## Pantanal
J69trk<-Pantanaltrk %>% filter(id=="69")
J70=subset(Drych1,id=='70');J70$idloc <- seq.int(nrow(J70))  ######## Drych1
J70trk<-Drych1trk %>% filter(id=="70")
J71=subset(Drych2,id=='71');J71$idloc <- seq.int(nrow(J71))  ####### Drych2
J71trk<-Drych2trk %>% filter(id=="71")
J72=subset(Drych2,id=='72');J72$idloc <- seq.int(nrow(J72))  ####### Drych2
J72trk<-Drych2trk %>% filter(id=="72")
J73=subset(Drych2,id=='73');J73$idloc <- seq.int(nrow(J73))  ####### Drych2
J73trk<-Drych2trk %>% filter(id=="73")
J74=subset(Pantanal,id=='74');J74$idloc <- seq.int(nrow(J74))   ########## Pantanal
J74trk<-Pantanaltrk %>% filter(id=="74")
J75=subset(Pantanal,id=='75');J75$idloc <- seq.int(nrow(J75))   ########## Pantanal
J75trk<-Pantanaltrk %>% filter(id=="75")
J76=subset(Drych1,id=='76');J76$idloc <- seq.int(nrow(J76))  ####### Drych1
J76trk<-Drych1trk %>% filter(id=="76")
J77=subset(Drych1,id=='77');J77$idloc <- seq.int(nrow(J77))  ####### Drych1
J77trk<-Drych1trk %>% filter(id=="77")
J78=subset(FPy,id=='78');J78$idloc <- seq.int(nrow(J78))   ####### FPy
J78trk<-FPytrk %>% filter(id=="78")
J79=subset(Pantanal,id=='79');J79$idloc <- seq.int(nrow(J79)) ########## Pantanal
J79trk<-Pantanaltrk %>% filter(id=="79")
J80=subset(Iguazu1,id=='80');J80$idloc <- seq.int(nrow(J80)) ######## Iguazu1
J80trk<-Iguazu1trk %>% filter(id=="80")
J81=subset(Pantanal,id=='81');J81$idloc <- seq.int(nrow(J81)) ########## Pantanal
J81trk<-Pantanaltrk %>% filter(id=="81")
J82=subset(Cerrado1,id=='82');J82$idloc <- seq.int(nrow(J82))  ########## Cerrado1
J82trk<-Cerrado1trk %>% filter(id=="82")
J83=subset(Iguazu2,id=='83');J83$idloc <- seq.int(nrow(J83))  ########## Iguazu2
J83trk<-Iguazu2trk %>% filter(id=="83")
J84=subset(Pantanal,id=='84');J84$idloc <- seq.int(nrow(J84))  ########## Pantanal
J84trk<-Pantanaltrk %>% filter(id=="84")
J85=subset(Cerrado1,id=='85');J85$idloc <- seq.int(nrow(J85))  ########## Cerrado1
J85trk<-Cerrado1trk %>% filter(id=="85")
J86=subset(Pantanal,id=='86');J86$idloc <- seq.int(nrow(J86)) ########## Pantanal
J86trk<-Pantanaltrk %>% filter(id=="86")
J87=subset(Pantanal,id=='87');J87$idloc <- seq.int(nrow(J87)) ########## Pantanal
J87trk<-Pantanaltrk %>% filter(id=="87")
J88=subset(Pantanal,id=='88');J88$idloc <- seq.int(nrow(J88)) ########## Pantanal
J88trk<-Pantanaltrk %>% filter(id=="88")
J89=subset(Cerrado2,id=='89');J89$idloc <- seq.int(nrow(J89)) ######### Cerrado2
J89trk<-Cerrado2trk %>% filter(id=="89")
J90=subset(Iguazu1,id=='90');J90$idloc <- seq.int(nrow(J90)) ######### Iguazu1
J90trk<-Iguazu1trk %>% filter(id=="90")
J91=subset(Pantanal,id=='91');J91$idloc <- seq.int(nrow(J91)) ########## Pantanal
J91trk<-Pantanaltrk %>% filter(id=="91")
J92=subset(Pantanal,id=='92');J92$idloc <- seq.int(nrow(J92)) ########## Pantanal
J92trk<-Pantanaltrk %>% filter(id=="92")
J93=subset(Mamiraua,id=='93');J93$idloc <- seq.int(nrow(J93)) ######### Mamiraua
J93trk<-Mamirauatrk %>% filter(id=="93")
J94=subset(Mamiraua,id=='94');J94$idloc <- seq.int(nrow(J94)) ######### Mamiraua
J94trk<-Mamirauatrk %>% filter(id=="94")
J95=subset(Mamiraua,id=='95');J95$idloc <- seq.int(nrow(J95)) ######### Mamiraua
J95trk<-Mamirauatrk %>% filter(id=="95")
J96=subset(Mamiraua,id=='96');J96$idloc <- seq.int(nrow(J96)) ######### Mamiraua
J96trk<-Mamirauatrk %>% filter(id=="96")
J97=subset(Mamiraua,id=='97');J97$idloc <- seq.int(nrow(J97)) ######### Mamiraua
J97trk<-Mamirauatrk %>% filter(id=="97")
J98=subset(Mamiraua,id=='98');J98$idloc <- seq.int(nrow(J98)) ######### Mamiraua
J98trk<-Mamirauatrk %>% filter(id=="98")
J99=subset(Mamiraua,id=='99');J99$idloc <- seq.int(nrow(J99)) ######### Mamiraua
J99trk<-Mamirauatrk %>% filter(id=="99")
J100=subset(Mamiraua,id=='100');J100$idloc <- seq.int(nrow(J100))######### Mamiraua
J100trk<-Mamirauatrk %>% filter(id=="100")
J101=subset(Pantanal,id=='101');J101$idloc <- seq.int(nrow(J101)) ########## Pantanal
J101trk<-Pantanaltrk %>% filter(id=="101")
J102=subset(Pantanal,id=='102');J102$idloc <- seq.int(nrow(J102)) ########## Pantanal
J102trk<-Pantanaltrk %>% filter(id=="102")
J103=subset(Pantanal,id=='103');J103$idloc <- seq.int(nrow(J103)) ########## Pantanal
J103trk<-Pantanaltrk %>% filter(id=="103")
J104=subset(Pantanal,id=='104');J104$idloc <- seq.int(nrow(J104)) ########## Pantanal
J104trk<-Pantanaltrk %>% filter(id=="104")
J105=subset(Pantanal,id=='105');J105$idloc <- seq.int(nrow(J105)) ########## Pantanal
J105trk<-Pantanaltrk %>% filter(id=="105")
J106=subset(Pantanal,id=='106');J106$idloc <- seq.int(nrow(J106)) ########## Pantanal
J106trk<-Pantanaltrk %>% filter(id=="106")
J107=subset(Pantanal,id=='107');J107$idloc <- seq.int(nrow(J107))########## Pantanal
J107trk<-Pantanaltrk %>% filter(id=="107")
J108=subset(Pantanal,id=='108');J108$idloc <- seq.int(nrow(J108))########## Pantanal
J108trk<-Pantanaltrk %>% filter(id=="108")
J109=subset(Pantanal,id=='109');J109$idloc <- seq.int(nrow(J109))########## Pantanal
J109trk<-Pantanaltrk %>% filter(id=="109")
J110=subset(Pantanal,id=='110');J110$idloc <- seq.int(nrow(J110))########## Pantanal
J110trk<-Pantanaltrk %>% filter(id=="110")
J111=subset(Pantanal,id=='111');J111$idloc <- seq.int(nrow(J111))########## Pantanal
J111trk<-Pantanaltrk %>% filter(id=="111")
J112=subset(Pantanal,id=='112');J112$idloc <- seq.int(nrow(J112))########## Pantanal
J112trk<-Pantanaltrk %>% filter(id=="112")
J113=subset(Pantanal,id=='113');J113$idloc <- seq.int(nrow(J113))########## Pantanal
J113trk<-Pantanaltrk %>% filter(id=="113")
J114=subset(Pantanal,id=='114');J114$idloc <- seq.int(nrow(J114))########## Pantanal
J114trk<-Pantanaltrk %>% filter(id=="114")
J115=subset(Pantanal,id=='115');J115$idloc <- seq.int(nrow(J115))########## Pantanal
J115trk<-Pantanaltrk %>% filter(id=="115")
J116=subset(Pantanal,id=='116');J116$idloc <- seq.int(nrow(J116))########## Pantanal
J116trk<-Pantanaltrk %>% filter(id=="116")
J117=subset(Pantanal,id=='117');J117$idloc <- seq.int(nrow(J117))########## Pantanal
J117trk<-Pantanaltrk %>% filter(id=="117")
```





```r
#    ALL  JAGUARS  =>  jaguartrk  #This apparently would assign only the first timezone appearing for all
# All Project regions trk
# jaguartrk=bind_rows(AFW1trk,AFW2trk,Caatingatrk,Cerrado1trk,Cerrado2trk,CRicatrk,Drychtrk, Hchtrk,FPytrk,
# Iguazutrk,Mamirauatrk,iopPAtrk,Lacandonatrk, MexEasttrk, Sonoratrk,Oncafaritrk,Paraguaytrk,Panthera1trk,
# Panthera2trk,RioNegrotrk,SaoBentotrk,Taiamatrk) 
# to get it ordered
# jaguartrk<-arrange(jaguartrk,id,t_)
```

   ALL  JAGUARS  => jaguarlist 


```r
# To put individuals together while keeping respective UTMs and Timezones, just in case we need it
jaguarlist = list(J1trk,J2trk,J3trk,J4trk,J5trk,J6trk,J7trk,J8trk,J9trk,J10trk,
                  J11trk,J12trk,J13trk,J14trk,J15trk,J16trk,J17trk,J18trk,J19trk,J20trk,
                  J21trk,J22trk,J23trk,J24trk,J25trk,J26trk,J27trk,J28trk,J29trk,J30trk,
                  J31trk,J32trk,J33trk,J34trk,J35trk,J36trk,J37trk,J38trk,J39trk,J40trk,
                  J41trk,J42trk,J43trk,J44trk,J45trk,J46trk,J47trk,J48trk,J49trk,J50trk,
                  J51trk,J52trk,J53trk,J54trk,J55trk,J56trk,J57trk,J58trk,J59trk,J60trk,
                  J61trk,J62trk,J63trk,J64trk,J65trk,J66trk,J67trk,J68trk,J69trk,J70trk,
                  J71trk,J72trk,J73trk,J74trk,J75trk,J76trk,J77trk,J78trk,J79trk,J80trk,
                  J81trk,J82trk,J83trk,J84trk,J85trk,J86trk,J87trk,J88trk,J89trk,J90trk,
                  J91trk,J92trk,J93trk,J94trk,J95trk,J96trk,J97trk,J98trk,J99trk,J100trk,
              J101trk,J102trk,J103trk,J104trk,J105trk,J106trk,J107trk,J108trk,J109trk,J110trk,
              J111trk,J112trk,J113trk,J114trk,J115trk,J116trk,J117trk)
```

Other info


```r
sessionInfo()	  
```

```
## R version 3.5.3 (2019-03-11)
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
## [10] RCurl_1.95-4.12     bitops_1.0-6        ggplot2_3.1.0      
## [13] leaflet_2.0.2       lattice_0.20-38     rgl_0.100.19       
## [16] colorspace_1.4-0    maptools_0.9-5      amt_0.0.5.0        
## [19] adehabitatLT_0.3.23 CircStats_0.2-6     boot_1.3-20        
## [22] MASS_7.3-51.1       adehabitatMA_0.3.12 ade4_1.7-13        
## [25] move_3.1.0          rgdal_1.4-3         raster_2.8-19      
## [28] sp_1.3-1            geosphere_1.5-7     install.load_1.2.1 
## 
## loaded via a namespace (and not attached):
##  [1] nlme_3.1-137            webshot_0.5.1          
##  [3] RColorBrewer_1.1-2      httr_1.4.0             
##  [5] tools_3.5.3             utf8_1.1.4             
##  [7] R6_2.3.0                rgeos_0.4-2            
##  [9] lazyeval_0.2.1          mgcv_1.8-27            
## [11] manipulateWidget_0.10.0 npsurv_0.4-0           
## [13] withr_2.1.2             tidyselect_0.2.5       
## [15] compiler_3.5.3          cli_1.0.1              
## [17] xml2_1.2.0              labeling_0.3           
## [19] scales_1.0.0            mvtnorm_1.0-8          
## [21] stringr_1.4.0           digest_0.6.18          
## [23] foreign_0.8-71          R.utils_2.7.0          
## [25] pkgconfig_2.0.2         htmltools_0.3.6        
## [27] highr_0.7               htmlwidgets_1.3        
## [29] rlang_0.3.0.1           rstudioapi_0.8         
## [31] shiny_1.2.0             bindr_0.1.1            
## [33] jsonlite_1.5            crosstalk_1.0.0        
## [35] R.oo_1.22.0             magrittr_1.5           
## [37] Matrix_1.2-15           Rcpp_1.0.0             
## [39] munsell_0.5.0           fansi_0.4.0            
## [41] R.methodsS3_1.7.1       stringi_1.2.4          
## [43] yaml_2.2.0              plyr_1.8.4             
## [45] grid_3.5.3              parallel_3.5.3         
## [47] promises_1.0.1          crayon_1.3.4           
## [49] miniUI_0.1.1.1          splines_3.5.3          
## [51] hms_0.4.2               pillar_1.3.0           
## [53] codetools_0.2-16        glue_1.3.0             
## [55] evaluate_0.12           lsei_1.2-0             
## [57] httpuv_1.4.5            gtable_0.2.0           
## [59] purrr_0.2.5             tidyr_0.8.2            
## [61] assertthat_0.2.0        xfun_0.4               
## [63] mime_0.6                xtable_1.8-3           
## [65] later_0.7.5             survival_2.43-3        
## [67] memoise_1.1.0           fitdistrplus_1.0-11
```

```r
proc.time()
```

```
##    user  system elapsed 
##  175.18  132.51 8191.76
```
