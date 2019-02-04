#' --- setwd("C:/RWorkDir/AnimoveRscripts")    
#' title: "TestVignette"
#' author: "John Fieberg"
#' date: ""
#' ---
#' 
#' ## Purpose
#' 
#' - test whether you have everything set up for the workshop
#' - retrieve and explore animal movement data from Movebank
#' - illustrate how to work with functions in the animal movement tools package (amt)
#' - write out available points for further annotating in EnvDATA
#' 
#' 
#' #### Preamble
#' 
#' Load libraries
#+warning=FALSE, message=FALSE
library(knitr)
library(lubridate)
library(maptools)
library(raster)
library(move)
library(amt) 
library(ggmap)
library(tibble)
library(leaflet)
library(dplyr)
options(width=165,digits.secs = 3)
opts_chunk$set(fig.width=12,fig.height=4.5, error=TRUE,cache = F)

#' Record time for running all code
ptm<-proc.time()

#' Set the seed for the random number generator, so it will be possible
#' to reproduce the random points
set.seed(10299)

#' Create a login object for a user account at movebank.org
loginStored <- movebankLogin(username="MovebankWorkshop", password="genericuserforthisexercise")

#' Get overview information about a Movebank study. Be sure to check the citation and license terms if not using your own data.

getMovebankStudy(study="Martes pennanti LaPoint New York", login=loginStored) # see study-level info

#' You may receive a message that you do not have access to data in a study. If you are not a 
#' data manager for a study, the owner may require that you read and accept the license terms 
#' once before you may download it. Currently license terms cannot be accepted directly from 
#' the move package. To view and accept them, go to movebank.org, search for the study and begin a 
#' download, and view and accept the license terms. After that you will be able to load from R.

#' Load data from a study in Movebank and create a MoveStack object. For more details and options see https://cran.r-project.org/web/packages/move/index.html.
fisher.move <- getMovebankData(study="Martes pennanti LaPoint New York", login=loginStored)
head(fisher.move)

#' Note: If there are duplicate animal-timestamp records in Movebank, you will get a warning. 
#' You can exclude duplicate records on import using removeDuplicatedTimestamps=T. If you are 
#' a data manager for a study in Movebank you can also filter them out directly in the study 
#' so they are excluded by default in downloads (see https://www.movebank.org/node/27252).

#' Create a data frame from the MoveStack object
fisher.dat <- as(fisher.move, "data.frame")

#' ### Data cleaning

#' Delete observations where missing lat or long or a timestamp.  There are no missing
#' observations in this data set, but it is still good practice to check.
ind<-complete.cases(fisher.dat[,c("location_lat", "location_long", "timestamp")])
fisher.dat<-fisher.dat[ind==TRUE,]

#' Check for duplicated observations (ones with same lat, long, timestamp,
#'  and individual identifier). There are no duplicate
#' observations in this data set, but it is still good practice to check.
ind2<-fisher.dat %>% select(timestamp, location_long, location_lat, local_identifier) %>%
  duplicated
sum(ind2) # no duplicates
fisher.dat<-fisher.dat[ind2!=TRUE,]

#' Make timestamp a date/time variable
fisher.dat$timestamp<-as.POSIXct(fisher.dat$timestamp, format="%Y-%m-%d %H:%M:%OS", tz="UTC")

#' Look at functions in the move package.
plot(fisher.move)
show(fisher.move)
summary(fisher.move)

#' ## Plots of the data
#' 
#' Note: there are lots of ways to plot data using R.  I have included code that 
#' illustrates 2 simple plotting methods for a single individual:
#' 
#' - using ggmap (with google maps)
#' - using leaflet
#' 
#' Both can become cumbersome and slow with too many observations.
#' Lets look at the data from F2. Note: When reading data from Movebank always be 
#' sure to refer to the animal-id and individual-local-identifer
#' and not the tag-id or tag-local-identifier in order to exclude pre- and 
#' post-deployment records and correctly separate out individual animals when the 
#' study includes redeployments.

fisherF2<-fisher.dat %>% filter(local_identifier=="F2")

#' We can plot the data using ggmap and ggplot;
#' ggmap extends ggplot2 by adding maps from e.g. google maps as background

#' Get the map that covers the study area. Zoom out a bit from what calc_zoom
#' suggests.  For your own data, you may have to change the zoom.

z<-(calc_zoom(location_long, location_lat, fisherF2))
map <- get_map(location = c(lon = mean(fisherF2$location_long), 
                            lat = mean(fisherF2$location_lat)), zoom = 12,
               maptype = "hybrid", source = "google")

ggmap(map) + 
  geom_point(data=fisherF2, aes(x=location_long, y=location_lat), size=2.5)

# plot the track on the map

#' Now, using leaflet
leaflet(fisherF2)%>%addTiles()%>%
  addCircles(fisherF2$location_long, fisherF2$location_lat)

#' ### Using ggplot without a background
#' 
#' Use separate axes for each individual (add scales="free" to facet_wrap)
#+fig.height=12, fig.width=12
ggplot(fisher.dat, aes(x=location_long, y=location_lat))+geom_point()+
  facet_wrap(~local_identifier, scales="free")

#' Now, all on 1 plot
#+fig.height=6, fig.width=12
ggplot(fisher.dat, aes(x=location_long, y=location_lat, color=as.factor(local_identifier)))+
  geom_point() 

#' ## Creating a track in amt
#' 
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
trk <- mk_track(fisher.dat, .x=location_long, .y=location_lat, .t=timestamp, id = local_identifier, 
                crs = CRS("+init=epsg:4326"))

# Now it is easy to calculate day/night with either movement track
trk <- trk %>% time_of_day()

#' Now, we can transform back to geographic coordinates
trk <- transform_coords(trk, CRS("+init=epsg:32618"))

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
temp<-trk %>% filter(id=="M1") %>% direction_rel
head(temp)

#' Or, we can add a columns to each nested column of data using purrr::map
trk<-trk %>% nest(-id) %>% 
  mutate(dir_abs = map(data, direction_abs,full_circle=TRUE, zero="N"), 
         dir_rel = map(data, direction_rel), 
         sl = map(data, step_lengths),
         nsd_=map(data, nsd))%>%unnest()

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

#' Same for KDE
kde.week<-trk %>% nest(-id,-year,  -month, -week) %>%  
  mutate(kdearea = map(data, ~hr_kde(., levels=c(0.95)) %>% hr_area)) %>%
  select(id, year, month, week,  kdearea) %>% unnest()

#+fig.height=12, fig.width=12, warning=FALSE, message=FALSE
ggplot(kde.week, aes(x = week, y = kdearea, colour=as.factor(year))) + geom_point()+
  geom_smooth()+ facet_wrap(~id, scales="free")


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
trk %>% filter(id=="F1") %>%
  random_points(.,factor = 20) %>% plot

#' Illustrate systematic points (to do this, we need to create the mcp first)
trk%>%filter(id=="F1") %>% 
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
avail <- SpatialPointsDataFrame(avail.pts[,c("x_","y_")], avail.pts, 
                                proj4string=CRS("+proj=utm +zone=18N +datum=WGS84"))  
avail.df <- data.frame(spTransform(avail, CRS("+proj=longlat +datum=WGS84")))[,1:6]
names(avail.df)<-c("idr", "case_", "utm.easting", "utm.northing", "location-long", "location-lat")

#' Check to make sure everything looks right
test<-subset(avail.df, case_==TRUE)
test %>% select('location-lat', 'location-long', utm.easting, utm.northing) %>% 
  summarise_all(mean)

fisher.dat %>% summarize(meanloc.lat=mean(location_lat), 
                         meanloc.long=mean(location_long))

#' Add a timestamp to annotate these data with environmental covariates in Movebank using Env-DATA (https://www.movebank.org/node/6607).
#' Here we just use the first timestamp, however meaningful timestamps are needed if annotating variables that vary in time.
avail.df$timestamp<-fisher.dat$timestamp[1]

#' These points then need to be annotated prior to fitting rsfs. Let's 
#' write out 2 files:
#' 
#' - FisherRSF2018.csv will contain all points and identifying information. 
#' - FisherRSFannotate.csv will contain only the columns used to create the annotation.
#' 
#' The latter file will take up less space, making it easier to annotate (and also possible to upload to github)
#write.csv(avail.df, file="data/FisherRSF2018.csv", row.names = FALSE)
write.csv(avail.df, file="c:/RWorkDir/AnimoveRscripts/FisherRSF2018.csv",row.names = FALSE)

avail.df<-avail.df %>% select("timestamp", "location-long", "location-lat")
write.csv(avail.df, file="c:/RWorkDir/AnimoveRscripts/FisherRSFannotate.csv", row.names = FALSE)

#' ## SSF prep
#' 
#' SSFs assume that data have been collected at regular time intervals.
#' We can use the track_resample function to regularize the trajectory so that
#' all points are located within some tolerence of each other in time. To figure
#' out a meaningful tolerance range, we should calculate time differences between
#' locations & look at as a function of individual.
(timestats<-trk %>% nest(-id) %>% mutate(sr = map(data, summarize_sampling_rate)) %>%
    select(id, sr) %>% unnest)

#' Time intervals range from every 2 to 15 minutes on average, depending
#' on the individual.  Lets add on the time difference to each obs.
trk<-trk %>% group_by(id) %>% mutate(dt_ = t_ - lag(t_, default = NA))

#' Let's illustrate track regularization with ID = F2. Let's
#' go from every 2 minutes to every 10.
tempF2<-trk %>% filter(id=="F2") %>% track_resample(rate=minutes(10), tolerance=minutes(2))
tempF2 %>% select(id, x_, y_, t_, burst_)

#' Now loop over individuals and do the following:
#' 
#' - Regularize trajectories using an appropriate time window (see e.g., below) 
#' - calculate new dt values
#' - Create bursts using individual-specific time intervals
#' - Generate random steps within each burst
#' 
#' The random steps are generated using the following approach:
#' 
#' 1. Fit a gamma distribution to step lenghts
#' 2. Fit a von mises distribution to turn angles
#' 3. Use these distribution to draw new turns and step lengths, form new simulated steps
#' and generate random x,y values.
#' 

#+warning=FALSE
ssfdat<-NULL
temptrk<-with(trk, track(x=x_, y=y_, t=t_, id=id))
uid<-unique(trk$id) # individual identifiers
luid<-length(uid) # number of unique individuals
for(i in 1:luid){
  # Subset individuals & regularize track
  temp<-temptrk%>% filter(id==uid[i]) %>% 
    track_resample(rate=minutes(round(timestats$median[i])), 
                   tolerance=minutes(max(10,round(timestats$median[i]/5))))
  
  # Get rid of any bursts without at least 2 points
  temp<-filter_min_n_burst(temp, 2)
  
  # burst steps
  stepstemp<-steps_by_burst(temp)
  
  # create random steps using fitted gamma and von mises distributions and append
  rnd_stps <- stepstemp %>%  random_steps(n = 15)
  
  # append id
  rnd_stps<-rnd_stps%>%mutate(id=uid[i])
  
  # append new data to data from other individuals
  ssfdat<-rbind(rnd_stps, ssfdat)
}
ssfdat<-as_tibble(ssfdat)
ssfdat

###Just to check the medians of what we did in comparison with the totals
ssfdatemp<-ssfdat %>% filter(case_==TRUE) 
ssfdatdf <- data.frame(ssfdatemp)
# A summary applied to ungrouped tbl returns a single row group by id,burst_
bursts <-ssfdatdf %>%
group_by(id,burst_) %>%
summarise(n_distinct = n_distinct(burst_), n = n())%>% data.frame
bursts  

#SSFs sampled data (x1, y1, t1)
ssfdattrk <- mk_track(ssfdatemp, .x=x1_, .y=y1_, .t=t1_, id=id, crs = CRS("+init=epsg:4326")) 
(timestatsssf<-ssfdattrk %>% nest(-id) %>% mutate(sr = map(data, summarize_sampling_rate, time_unit = "hour")) %>%
    select(id, sr) %>% unnest)
### SSFs sampled data (x2, y2, t2)
ssfdattrk <- mk_track(ssfdatemp, .x=x2_, .y=y2_, .t=t2_, id=id, crs = CRS("+init=epsg:4326")) 
(timestatsssf<-ssfdattrk %>% nest(-id) %>% mutate(sr = map(data, summarize_sampling_rate, time_unit = "hour")) %>%
    select(id, sr) %>% unnest)
###summarize_sampling_rate(ssfdattrk %>% filter(id=="F1")) ### if we want to look one specific animal

##Summary of totals from the original trk(as we did above just to compare)
(timestats<-trk %>% nest(-id) %>% mutate(sr = map(data, summarize_sampling_rate)) %>%
    select(id, sr) %>% unnest)





#' Now, lets plot the data for random and matched points
#' 
#+fig.height=12, fig.width=12, warning=FALSE
ggplot(ssfdat, aes(x2_, y2_, color=case_))+geom_point()+facet_wrap(~id, scales="free")


#' Relabel as utms 
ssfdat$utm.easting<-ssfdat$x2_
ssfdat$utm.northing<-ssfdat$y2_

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
#' variables and write out the data sets. With the SSFs, we have the extra complication of
#' having a time and location at both the start and end of the step.  
#' 
#' For the time being, we will assume we want to annotate variables at the end of the step
#' but use the starting point of the step as the timestamp.
#' 
#' You could also calculate the midpoint of the timestep like this:
#' data$timestamp.midpoint <- begintime + (endtime-begintime)/2

ssfdat2 <- SpatialPointsDataFrame(ssfdat[,c("x2_","y2_")], ssfdat, 
                                  proj4string=CRS("+proj=utm +zone=18N +datum=WGS84"))  
ssf.df <- data.frame(spTransform(ssfdat2, CRS("+proj=longlat +datum=WGS84"))) 
names(ssf.df)[c(13,16,17)] <-c("individual.local.identifier", "location-long", "location-lat")
ssf.df$timestamp<-ssf.df$t1_
ssf.df %>% select('location-lat', utm.easting, x1_, x2_, y1_, y2_, 'location-long', utm.northing) %>% head


#' These points then need to be annotated prior to fitting ssfs. Let's 
#' write out 2 files:
#' 
#' - FisherSSF2018.csv will contain all points and identifying information. 
#' - FisherSSFannotate.csv will contain only the columns used to create the annotation.
#' 
#' The latter file will take up less space, making it easier to annotate (and also possible to upload to github)
#write.csv(ssf.df, file="data/AllStepsFisher2018.csv", row.names=FALSE)
#ssf.df<-ssf.df %>% select("timestamp", "location-long", "location-lat")
#write.csv(ssf.df, file="data/FisherSSFannotate.csv", row.names = FALSE)

write.csv(ssf.df, file="c:/RWorkDir/AnimoveRscripts/AllStepsFisher2018.csv", row.names=FALSE)
ssf.df<-ssf.df %>% select("timestamp", "location-long", "location-lat")
write.csv(ssf.df, file="c:/RWorkDir/AnimoveRscripts/FisherSSFannotate.csv", row.names = FALSE)


#' ## Using nested data frames
#' 
#' This works if all animals are sampled at a constant sampling rate.
#' Again, you need to create a new column (using `mutate`) where you save the random steps to
#' Not all animals have a 15 min sampling rate, so we might drop the first
#' animal that has a 15 minute sampling rate. 
#trk %>% nest(-id) %>% mutate(sr = map(.$data, summarize_sampling_rate)) %>% 
#  select(id, sr) %>% unnest()

# Then, we could avoid a loop with:
#ssfdat <- trk %>% filter(id != "M1") %>% nest(-id) %>% 
#  mutate(ssf = map(data, function(d) {
#    d %>%
#      track_resample(rate = minutes(10), tolerance = minutes(2)) %>% 
#      filter_min_n_burst(min_n = 3) %>% 
#      steps_by_burst() %>% random_steps()
#  })) %>% select(id, ssf) %>% unnest()

#' ## Document Footer	
#' 	
#' Session Information:	
#' 	
sessionInfo()	  

proc.time()-ptm