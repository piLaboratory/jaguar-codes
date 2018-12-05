############ Adapted from MoveBank ("John Fieberg")& Morato et al.2018 / Bernardo scripts
rm(list= ls())                                                  ### For a fresh start
setwd("C:/RWorkDir/jaguardatapaper")                            ### Set directory
###########################################################################################################         																					
###Enter jaguar data from Morato et al. 2018
## Load packages
if(!require(install.load)) install.packages('install.load'); library(install.load)
install.load::install_load("ggmap","maptools",'move',"circular","RCurl","dplyr","readr","caTools","adehabitatLT","ggsn","magick","rgl","tidyverse","htmltools","rglwidget","knitr","lubridate","raster","amt","tibble","leaflet","ezknitr")
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
info <- read.delim(file="c:/RWorkDir/jaguardatapaper/info.txt")

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
#hist(mov.traj.df$dist)

## move
# Organize data as a move package format
move.data <- move(x = mov.traj.df$x, y = mov.traj.df$y, 
                  time = mov.traj.df$date, 
                  proj = CRS('+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'),
                  data = mov.traj.df, animal = mov.traj.df$id, sensor = 'GPS')

move.data  ### moveStack
summary(move.data)
###Separate individual animals' trajectories 
#unstacked <- split(move.data)
#head(unstacked)
jaguar_df <- as(move.data, "data.frame")
id <- as.integer(levels(jaguar_df$id))[jaguar_df$id]
jaguar_df$id <- id   ### does the id need to be an integer?  Or factor works fine?
head(jaguar_df)
str(jaguar_df)
#########################################################################################################							   							   									  
#### Preamble  ##########################################################################################
library(knitr)
options(width=150) #in merge file
options(width=165,digits.secs = 3)  #show figures in vignette test
opts_chunk$set(fig.width=12,fig.height=4.5, error=TRUE,cache = F) #show figures in the vignette test
#' Record time for running all code
ptm<-proc.time()
#' Set the seed for the random number generator, so it will be possible
#' to reproduce the random points
set.seed(10299)

#' ### Data cleaning ##########################################
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

##############################################################################################################
#' ## Creating a track in amt (commom to both RSF and SSF)  
#' ###########################################################################################################
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
trk <- mk_track(jaguar_df, .x=x, .y=y, .t=date, id=id, crs = CRS("+init=epsg:4326"))
trk <- trk %>% arrange(id)
trk
#trk <- mk_track(jaguar_df, .x=x, .y=y, .t=date, id =id, 
                #crs = CRS('+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'))	# same as

# Now it is easy to calculate day/night with either movement track
trk <- trk %>% time_of_day()  ### Data says it is GMT, so it would need to be corrected to local??????

#' Now, we can transform back to geographic coordinates (Perhaps there is no need since it looks the same)
#trk <- transform_coords(trk, CRS('+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'))

#' If we create trk.temp as above, we can make sure that we got back the same coordinates 
#' by comparing to the original track
#mean(trk$x_ - trk.temp$x_) # look good, there is virtually no difference
#mean(trk$y_ - trk.temp$y_) # same for y

#' Save the class here (and apply it later after adding columns to the 
#' object)
trk.class<-class(trk)

#' ## Movement Characteristics
#' 
#' - dir_abs will calculate absolute angles for steps
#' - dir_rel will calculate turning angles (relative angles)
#' - step_lengths will calculate distances between points
#' - nsd = Net Squared Displacement (distance from first point)
#' 
#' Arguments direction_abs:
#' 
#' - full_circle will calculate between 0 and 360 (rather than -180 and 180)
#' - zero gives the direction = 0 for absolute angle
#' 
#' Note:  we have to calculate these characteristics separately for each 
#' individual (to avoid calculating a distance between say the last observation
#' of the first individual and the first observation of the second individual).
#' 
#' 
#' To do this, we could loop through individuals, calculate these
#' characteristics for each individual, then rbind the data 
#' back together.  Or, use nested data frames and the map function
#' in the purrr library to do this with very little code. 
#' 
#' To see how nesting works, we can create a nested object by individual
nesttrk<-trk%>%nest(-id)
nesttrk

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

###################################################################################################
###########            Applied it to each Area/ Project       ############################################
### An Example
#  1) IOP Cerrado (including Para)
animals <- as(move.data, 'data.frame')
iop=subset(animals,study.name=='IOP')
iopPA=subset(iop,id=='24')
iop17=subset(iop,id=='17')
iop65=subset(iop,id=='65')
iop67=subset(iop,id=='67')
iop82=subset(iop,id=='82')
iop85=subset(iop,id=='85')

#IOP Cerrado  (without the individual from Para)
iopCerrado=rbind(iop17, iop65, iop67, iop82, iop85)
head(iopCerrado)
str(iopCerrado)

### iop includes Cerrado and Pauapebas, iopCerrado includes only Cerrado animals
trk <- mk_track(iopCerrado, .x=x, .y=y, .t=date, id=id, crs = CRS("+init=epsg:4326")) 
trk <- trk %>% arrange(id)
trk
trk <- trk %>% time_of_day()  ### I believe it is GMT so it would need to be corrected
trk.class<-class(trk)
nesttrk<-trk%>%nest(-id)
nesttrk

trk<-trk %>% nest(-id) %>% 
  mutate(dir_abs = map(data, direction_abs,full_circle=TRUE, zero="N"), 
         dir_rel = map(data, direction_rel), 
         sl = map(data, step_lengths),
         nsd_=map(data, nsd))%>%unnest()
		 
trk<-trk%>% mutate(week=week(t_),
    month = month(t_, label=TRUE), 
    year=year(t_),
    hour = hour(t_))
	
class(trk)
class(trk)<-trk.class

#' Lets take a look at what we created
trk
### trk is corresponding to iopCerrado in the examples
################################################################################################################
#' ## Some plots of movement characteristics

#' ### Absolute angles (for each movement) relative to North 
#' We could use a rose diagram (below) to depict the distribution of angles. 
#+fig.height=12, fig.width=12
ggplot(trk, aes(x = dir_abs, y=..density..)) + geom_histogram(breaks = seq(0,360, by=20))+
  coord_polar(start = 0) + theme_minimal() + 
  scale_fill_brewer() + ylab("Density") + ggtitle("Angles Direct") + 
  scale_x_continuous("", limits = c(0, 360), breaks = seq(0, 360, by=20), 
                     labels = seq(0, 360, by=20))+
  facet_wrap(~id)

#' ### Turning angles 
#' 
#' Note: a 0 indicates the animal continued to move in a straight line, a 180 
#' indicates the animal turned around (but note, resting + measurement error often can
#' make it look like the animal turned around).
#+fig.height=12, fig.width=12
ggplot(trk, aes(x = dir_rel, y=..density..)) + geom_histogram(breaks = seq(-180,180, by=20))+
  coord_polar(start = 0) + theme_minimal() + 
  scale_fill_brewer() + ylab("Density") + ggtitle("Angles Direct") + 
  scale_x_continuous("", limits = c(-180, 180), breaks = seq(-180, 180, by=20), 
                     labels = seq(-180, 180, by=20))+
  facet_wrap(~id)

#' ### Turning angles as histograms
#+fig.height=12, fig.width=12
ggplot(trk, aes(x = dir_rel)) +  geom_histogram(breaks = seq(-180,180, by=20))+
  theme_minimal() + 
  scale_fill_brewer() + ylab("Count") + ggtitle("Angles Relative") + 
  scale_x_continuous("", limits = c(-180, 180), breaks = seq(-180, 180, by=20),
                     labels = seq(-180, 180, by=20))+facet_wrap(~id, scales="free")

#' ### Net-squared displacement over time for each individual
#+fig.height=12, fig.width=12
ggplot(trk, aes(x = t_, y=nsd_)) + geom_point()+
  facet_wrap(~id, scales="free")


#' ## Explore movement characteristics by (day/night, hour, month)
#' 
#' ### step length distribution by day/night
#' 
#+fig.height=12, fig.width=12, warning=FALSE, message=FALSE
ggplot(trk, aes(x = tod_, y = log(sl))) + 
  geom_boxplot()+geom_smooth()+facet_wrap(~id)

#' ### Space use (MCP or KDE) by week, month, and year
#' 
#' Note:  this code will only work for your critters if you
#' have multiple observations for each combination of
#' (month, year).  If you don't have many observations, you could
#' try:  nest(-id, -year) and unnest(id,year)
mcps.week<-trk %>% nest(-id,-year,  -month, -week) %>%  
  mutate(mcparea = map(data, ~hr_mcp(., levels = c(0.95)) %>% hr_area)) %>% 
  select(id, year, month, week, mcparea) %>% unnest()

#+fig.height=12, fig.width=12, warning=FALSE, message=FALSE
ggplot(mcps.week, aes(x = week, y = area, colour=as.factor(year))) + geom_point()+
  geom_smooth()+ facet_wrap(~id, scales="free")

#' Same for KDE   (for some reason not working...)
#kde.week<-trk %>% nest(-id,-year, -month, -week) %>%  
  #mutate(kdearea = map(data, ~hr_kde(., levels=c(0.95)) %>% hr_area)) %>%
  #select(id, year, month, week,  kdearea) %>% unnest()

#+fig.height=12, fig.width=12, warning=FALSE, message=FALSE
#ggplot(kde.week, aes(x = week, y = kdearea, colour=as.factor(year))) + geom_point()+
 # geom_smooth()+ facet_wrap(~id, scales="free")


#' ## RSF prep
#' 
#' Generate random points within MCPs for a single individual using amt functions.
#' Notes:
#' 
#' - It is common to generate points randomly, but other options are possible. 
#' - In particular, it can beneficial to generate a systematically placed sample
#' - Samples can also be generated using the *spsample* function in the sp library or
#' using a GIS (note: amt uses the spsample function within its function random_points)
#' - Other home range polygons could be used (e.g., kernel density, local convex hull
#' etc.)
#' 
#' #### Random points: illustrate for 1 individual
trk %>% filter(id=="17") %>%
  random_points(.,factor = 20) %>% plot

#' Illustrate systematic points (to do this, we need to create the mcp first)
trk%>%filter(id=="17") %>% 
  random_points(., factor = 20, type="regular") %>% 
  plot() 

#' Now, lets generate points for all individuals. We can do this
#' efficiently by making use of pipes (%>%),nested data frames, and
#' then by adding a new column -- a list-column -- to trks
avail.pts <- trk %>% nest(-id) %>% 
  mutate(rnd_pts = map(data, ~ random_points(., factor = 20, type="regular"))) %>% 
  select(id, rnd_pts) %>%  # you dont want to have the original point twice, hence drop data
  unnest()

#' Or, we could do this using a loop (commented out, below)
#avail.pts<-NULL
#uid<-unique(trk$id) # individual identifiers
#luid<-length(uid) # number of unique individuals
#for(i in 1:luid){
# random_points will generate random points within mcp
# Add on the individual id and combine all data
#  temp<-cbind(id=uid[i],trk%>%filter(id==uid[i])%>%random_points)
#  avail.pts<-rbind(avail.pts, temp)
#}
#avail.pts<-as_tibble(avail.pts)
#avail.pts

#' ## Write out data for further annotating
#' 
#' Need to rename variables so everything is in the format Movebank requires for annotation of generic time-location 
#' records (see https://www.movebank.org/node/6608#envdata_generic_request). This means, we need the following variables:
#' 
#' - location-lat (perhaps with addition of Easting/Northing in UTMs)
#' - location-long (perhaps with addition of Easting/Northing in UTMs)
#' - timestamp (in Movebank format)
#' 
#' Need to project to lat/long, while also keeping lat/long. Then rename
#' variables and write out the data sets.
#avail <- SpatialPointsDataFrame(avail.pts[,c("x_","y_")], avail.pts, 
                              #  proj4string=CRS("+proj=utm +zone=18N +datum=WGS84"))  
#avail.df <- data.frame(spTransform(avail, CRS("+proj=longlat +datum=WGS84")))[,1:6]
#names(avail.df)<-c("idr", "case_", "utm.easting", "utm.northing", "location-long", "location-lat")


avail.df <-data.frame(avail.pts)
avail.df$long<-avail.df$x_
avail.df$lat<-avail.df$y_


#' Check to make sure everything looks right ### ('location-lat', 'location-long', utm.easting, utm.northing)
test<-subset(avail.df, case_==TRUE)
test %>% select('long','lat') %>% 
  summarise_all(mean)

 #For jaguar_df (entire dataset) or iopCerrado (current example)

#jaguar_df %>% summarize(meanloc.long=mean(x), meanloc.lat=mean(y)) 
iopCerrado %>% summarize(meanloc.long=mean(x), meanloc.lat=mean(y))

#' Add a timestamp to annotate these data with environmental covariates in Movebank using Env-DATA (https://www.movebank.org/node/6607).
#' Here we just use the first timestamp, however meaningful timestamps are needed if annotating variables that vary in time.
#avail.df$timestamp<-jaguar_df$date[1]
avail.df$timestamp<-iopCerrado$date[1]

#' These points then need to be annotated prior to fitting rsfs. Let's 
#' write out 2 files:
#' 
#' - FisherRSF2018.csv will contain all points and identifying information. 
#' - FisherRSFannotate.csv will contain only the columns used to create the annotation.
#' 
#' The latter file will take up less space, making it easier to annotate (and also possible to upload to github)
#write.csv(avail.df, file="data/FisherRSF2018.csv", row.names = FALSE)
write.csv(avail.df, file="c:/RWorkDir/jaguardatapaper/JaguarRSF_iopCerrado.csv",row.names = FALSE)
available<-avail.df %>% select("timestamp", "long", "lat")
write.csv(avail.df, file="c:/RWorkDir/jaguardatapaper/RSFannotateiopCerrado.csv", row.names = FALSE)


#############################################################################################################
   ################  TO BE ADJUSTED !!!!!!  #################################################################
   #################      MERGE WITH GIS #################################################################### 
   
## RSF data

#' Read in original data (used and available points) and merge on environmental
#' data.
#rsfdattemp<-read.csv("data/FisherRSF2018.csv")
#annotated<-read.csv("data/FisherRSFannotate.csv-6951661009794024140.csv")
rsfdattemp<-read.csv("c:/RWorkDir/AnimoveRscripts/FisherRSF2018.csv")
annotated<-read.csv("c:/RWorkDir/AnimoveRscripts/FisherRSFannotate.csv-6951661009794024140.csv")

#martesStack<- getMovebankData(study="Martes pennanti LaPoint New York", login=login)
#martesStack
#Martes<-read.csv("C:/RWorkDir/AnimoveRscripts/Martes.csv")

#' Now, merge these by "timestamp", "location-long", "location-lat")
rsfdat<-merge(rsfdattemp, annotated)

#' Write out data for use in FisherRSF2018.R
#write.csv(rsfdat, "data/FisherRSF2018-EnvDATA-results.csv", row.names = FALSE)
write.csv(rsfdat, "C:/RWorkDir/AnimoveRscripts/FisherRSF2018-EnvDATA-results.csv", row.names = FALSE)
   
##############################################################################################################
#############################################################################################################
   ################  TO BE ADJUSTED !!!!!!  #################################################################
   #################      ANALYSIS       #################################################################### 
   
   
#' Read in annotated data. 
#rsfdat<-read.csv("data/FisherRSF2018-EnvDATA-results.csv")
rsfdat<-read.csv("C:/RWorkDir/AnimoveRscripts/FisherRSF2018-EnvDATA-results.csv")


#' Simplify some variable names
names(rsfdat)[c(4,8:10)]<-c("id","LandClass", "Elevation", "PopDens")
rsfdat$case_<-as.numeric(rsfdat$case_)

#' Create landcover classes (as suggested by Scott Lapoint :)
rsfdat$LandClass<-as.character(rsfdat$LandClass)
rsfdat<-rsfdat %>% mutate(landC = fct_collapse(LandClass,
      agri = c("11", "14", "30"),
      forest =c("30","40","50","60", "70","80", "90","100"),
      shrub= c("110", "130", "150"),
      grass = c("120", "140"),
      wet= c("160"),
      other = c("170", "180", "190", "200", "210", "220")))

#' Center and scale variables
rsfdat<-rsfdat %>% mutate(elev=as.numeric(scale(Elevation)), 
                          popD=as.numeric(scale(PopDens)))

#' ## Explore the data:

#' Look Distribution of variables for used and available
#+fig.width=8, fig.height=8
binomial_smooth <- function(...) {
  geom_smooth(method = "glm", method.args = list(family = "binomial"), ...)
}
ggplot(rsfdat,aes(x=Elevation, y=case_))+
  stat_smooth(method="glm", method.args = list(family = "binomial"))+
  binomial_smooth(formula = y ~ splines::ns(x, 5), colour="red")+
  facet_wrap(~id, scales="free")

ggplot(rsfdat,aes(x=PopDens, y=case_))+
  stat_smooth(method="glm", method.args = list(family = "binomial"))+
  binomial_smooth(formula = y ~ splines::ns(x, 5), colour="red")+
  facet_wrap(~id, scales="free")

ggplot(rsfdat, aes(x=landC, y=..prop..,group=case_, colour=case_))+
  geom_bar(position="dodge", aes(fill=case_))+facet_wrap(~id, scales="free")

#' ## RSF fitting

#' Weight available data 
rsfdat$w<-ifelse(rsfdat$case_==1, 1, 5000)

#' We can fit an RSF model to a single animal using logistic regression
summary(glm(case_ ~ elev+popD+landC, data = subset(rsfdat, id=="M2"), weight=w,family = binomial))

#' Note, this individual did not experience all landcover classes
rsfdat %>% filter(id=="M2") %>% with(table(case_, landC))  
rsfdat$used<-as.factor(rsfdat$case_)
rsfdat$used<-fct_recode(rsfdat$used, "avail"="0", "used"="1")

#+fig.width=6, fig.height=4
ggplot(subset(rsfdat, id=="M2"),  aes(x=landC,group=used))+
  geom_bar(position=position_dodge(), aes(y=..prop.., fill = used), stat="count") +
  scale_fill_brewer(palette="Paired")+
  geom_text(aes( label = scales::percent(..prop..),
                 y= ..prop.. ), stat= "count", vjust = -.3, position=position_dodge(0.9)) +
  labs(y = "Proportion", fill="used", x="Landcover") 

#' Now, fit an RSF model to data from each animal.  Since not all animals experience
#' all habitat types, lets just explore forest versus non-forest
rsfdat$forest<-ifelse(rsfdat$landC=="forest", 1, 0)

fit_rsf <- function(data){
  mod <- glm(case_ ~ elev+popD+forest, data = data, weight=w,family = binomial)
  return(mod)
}
rsffits <- rsfdat %>%  nest(-id) %>% 
  dplyr::mutate(mod = purrr::map(data, fit_rsf)) 

#' This stores a list of model fits
rsffits 

#' Look at first model
rsffits$mod[[1]]
 
#' Now, use tidy to extract information about the model fits into a nested
#' data frame
rsffits <- rsffits %>%
   dplyr::mutate(tidy = purrr::map(mod, broom::tidy),
                 n = purrr::map(data, nrow) %>% simplify())
rsffits 
rsffits$tidy
 
#' Now, create data frame w/ the coefficients, etc
 rsf_coefs <- rsffits %>%
 tidyr::unnest(tidy) %>%
   dplyr::select(-(std.error:p.value))
 
 rsf_coefs %>% tidyr::spread(term, estimate)
  
#' Plot coefficients
#+fig.width=12, fig.heigh=4
 rsf_coefs %>% filter(term!="(Intercept)") %>%
   ggplot(., aes(x=1, y=estimate)) + 
      geom_dotplot(binaxis="y", stackdir="center")+geom_hline(yintercept=0)+
      facet_wrap(~term, scales="free")
 
#' Write out coefficients for MultipleAnimals.R
# save(rsf_coefs, file="data/rsfcoefs.Rdata")
save(rsf_coefs, file="C:/RWorkDir/AnimoveRscripts/rsfcoefs.Rdata")

#' ## Document Footer	

#' Session Information:	
#' 	
 sessionInfo()	  
    