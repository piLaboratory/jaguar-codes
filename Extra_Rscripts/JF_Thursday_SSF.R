#' # Fitting SSFs to individual animals
#'
#' Goals: 
#' 
#' - Illustrate data development for fitting SSF models
#' - Fit models to individual animals and then generate population-level summaries
#' 
#' ### Preamble
#'
#' First, delete everything in working memory
rm(list=ls(all=TRUE))
#' 
#' Load libraries
#+ warning=FALSE, message=FALSE 
library(ezknitr)
library(knitr)
library(lubridate)
library(raster)
library(move)
library(spatstat)
library(maptools)
library(amt) 

#' Load amt last, to avoid conflicts (with functions having the same
#' name in multiple packages)
#+ warning=FALSE, message=FALSE 
library(amt)  

#' Other options
options(width=150)
opts_chunk$set(fig.width=12,fig.height=4.5, dev="png")


#' Load processed data
load("Datasets/ProcessedData/Hyena.Rdata")
load("Datasets/ProcessedData/GeogrLayers.Rdata")

#' Turn into normal data frame
SHd<-as.data.frame(SH) 

#' Sort by id and timestamp
SHd<-SHd[order(SHd$IND, SHd$TIMESTAMP),]

#' ## Data development for distance to dens
#' 
#' Before we can calculate the distance between points and den sites, 
#' we need to project the den site locations into the correct coordinate system.  
#' We will create a data frame w/ these coordinates & then merge this information
#' back onto the original data set.
ind<-is.na(SHd$DEN_LAT)!=TRUE # observations with den locations

#' Create a SpatialPointsDataFrame w/ the locations with den sites, give it the 
#' current projection
SHDen<-SHd%>%filter(ind==TRUE) %>% select("DEN_LONG","DEN_LAT") 
SHDen<-SpatialPointsDataFrame(coords = SHDen[c("DEN_LONG","DEN_LAT")],
                              data=SHDen,
                              proj=CRS("+proj=longlat +datum=WGS84"))
                  
#' Transform the coordinates system 
#' -> FIRST we define it with proj4string, THEN we transform it
SHDen <- spTransform(SHDen, CRS("+proj=utm +zone=34 +south +ellps=WGS84"))
DenDat<-data.frame(SHDen@coords)
names(DenDat)<-c("UTM.x", "UTM.y")
head(DenDat)

#' Merge on and keep only data including den sites
SHd<-SHd[is.na(SHd$DEN_LAT)!=TRUE,]
SHd<-cbind(SHd, DenDat)

#' ## SSF data prep
#'  
#' First, create a track.   
trk<-with(SHd,track(x=LONGITUDE , y=LATITUDE , t=TIMESTAMP, id=IND,
                    DEN.x=UTM.x, DEN.y=UTM.y, den=DENNING))  

#' Add on week, month, hour, year variables to help explore sampling patterns.
trk<-trk%>% 
  mutate(
    week=week(t_),
    month = month(t_, label=TRUE), 
    year=year(t_),
    hour = hour(t_)
  )

#' ## SSF prep
#' 
#' #### Data prep: we need data collected regularly in time
#'  
#' SSFs assume that data have been collected at regular time intervals.  Lets see 
#' how close are data fit this assumption.
trk<-trk %>% group_by(id) %>% mutate(dt = t_ - lag(t_, default = NA))

#' Plot duration between points by hour.  Note: sampling frequency was reduced mid-day
#+ fig.height=12, fig.width=12
ggplot(trk, aes(as.factor(hour), dt))+geom_boxplot()+ facet_wrap(~id)+ylim(0,24)
table(trk$hour)

#' Lets subset the data so we only consider hours 18-6
trk2<- filter(trk, hour%in%c(0,2,4,6,8,18,20,22))

#' ### Regularize the trajectory further
#' 
#' We can use the track_resample function to regularize the trajectory so that
#' all points are located within some tolerence of each other in time. To figure
#' out a meaningful tolerance range, we could calculate time differences between
#' locations & look at as a function of individual.
trk
timestats<- trk2 %>% nest(-id) %>% 
  mutate(sr = map(data, summarize_sampling_rate)) %>% 
  dplyr::select(sr, id) %>% unnest()
timestats

#' Or, alternatively...
timestats<-trk2%>%group_by(id)%>%summarize(meandt=mean(dt, na.rm=T), mediandt=median(dt, na.rm=T), 
                                          mindt=min(dt, na.rm=T), IQR1=IQR(dt, na.rm=T)) 
timestats

#' Note, shortest time difference is 2 hours.  Lets regularize the data
#' so that we get observations every 2 hours. 

#' We can use track_resample to regularize the trajectories.  Lets explore how this
#' works, just using data for Sonny. 
temp2<-trk2%>% filter(id=="Sonny") %>% 
    track_resample(rate=hours(2), tolerance=minutes(30)) 
print(temp2, width=Inf)

#' Now, can get rid of any observations without at least 2 consecutive observations
temp3<-temp2%>% filter_min_n_burst(2)

#' Now loop over individuals and do the following:
#' 
#' - Regularize trajectories using a two hour time window 
#' - calculate new dt values
#' - Generate random steps within each burst (defined by a set of consecutive 
#' locations all 2 hours apart)
#' 
#' The random steps are generated using the following approach:
#' 
#' 1. Fit a gamma distribution to step lengths
#' 2. Fit a von mises distribution to turn angles
#' 3. Use these distribution to draw new turns and step lengths, form new simulated steps
#' and generate random x,y values.
#' 
#+warning=FALSE
ssfdat<-NULL # set up object to hold data
uid<-unique(trk$id) # unique individuals
luid<-length(uid) # number of unique individuals
temptrk<-with(trk, track(x=x_, y=y_, t=t_, id=id, den=den))
for(i in 1:luid){
    # Subset individuals & regularize track
    temp<-temptrk%>% filter(id==uid[i]) %>% 
      track_resample(rate=hours(2), 
                     tolerance=minutes(30))
    
    # Get rid of any bursts without at least 2 points
    temp<-filter_min_n_burst(temp, 2)
    
    # Create bursts and steps
    stepstemp<-steps_by_burst(temp)
    
    # create random steps using fitted gamma and von mises distributions and append
    rnd_stps <- stepstemp %>%  random_steps(n = 15)
    
    # append id
    rnd_stps<-rnd_stps%>%mutate(id=uid[i])
    
    # append new data to data from other individuals
    ssfdat<-rbind(rnd_stps, ssfdat)
}
ssfdat.t<-as_tibble(ssfdat)
ssfdat.t

#' Merge back on den locations whether or not the individual was denning
trktemp<-trk%>%select("id", "t_", "DEN.x", "DEN.y", "den") 
ssfdat<-base::merge(ssfdat.t, trktemp, by.x=c("id", "t1_"), by.y=c("id", "t_"))
print(sfsdat, width=Inf) 

#' Calculate distance from points to den
ssfdat$dendist<-with(ssfdat, sqrt((x2_-DEN.x)^2+(y2_-DEN.y)^2))
ssfdat.t<-ssfdat

#' Plot the data
#+fig.height=12, fig.width=12, warning=FALSE
ggplot(ssfdat.t, aes(x2_, y2_, color=case_))+geom_point(alpha=0.1)+facet_wrap(~id, scales="free")

#' ### Extract covariate data at the points
#' 
#' Now, extract covariates at the end of the used and avialable steps
ssfdat<-with(ssfdat.t, track(x=x2_, y=y2_, id=id, case_=case_, sl=sl_,ta=ta_, 
                             step_id=step_id_,den=den,
                             dendist=dendist))
ssfdat<-ssfdat %>% extract_covariates(env)
print(sfsdat, width=Inf)

#' Now create "distance to" covariates (also takes time!).  To do this
#' we will make use of the nncross function in the spatstat library.
ssf.ppp<-as.ppp(as_sp(ssfdat)) 
ssfdat$distroad<- nncross(ssf.ppp, as.psp(roads))$dist
ssfdat$distcamps<-  nncross(ssf.ppp, as.ppp(camps))$dist
ssfdat$distrivers<- nncross(ssf.ppp, as.psp(rivers))$dist
ssfdat$distcattlep<-nncross(ssf.ppp, as.ppp(cattleP))$dist

#' Scale and center all covariates by subtracting their mean and dividing by
#' their standard deviation.  Regression coefficients will give the relative 
#' risk of use per 1 sd change in each covariate.
ssfdat<-mutate(ssfdat, dendist.sc=scale(dendist),
               distroad.sc=scale(distroad),
               distcamps.sc=scale(distcamps),
               distrivers.sc=scale(distrivers), 
               distcattlep.sc=scale(distcattlep),
               ndvi.sc=scale(NDVI))

#' Have to be careful with categorical variables, make sure 
#' all animals are found in all categories.  Note: veg categories 1 and 5 are
#' problematic for some individuals (e.g., Gin, Fly, Bumble, Tori).
with(ssfdat, table(id, vegetation, case_))

#' ## Explore the data:

#' Look Distribution of variables for used (case_ ==TRUE) and available (case_==FALSE)
#' 
#' ### Vegetation
#+fig.width=8, fig.height=8
ggplot(ssfdat, aes(as.factor(vegetation), y=..prop..,group=case_, colour=case_))+
  geom_bar(position="dodge", aes(fill=case_))+facet_wrap(~id, scales="free")+
  xlab("Vegatation")

#' ### NDVI
#+fig.width=8, fig.height=8
ggplot(ssfdat, aes(ndvi.sc, colour=case_))+geom_density()+facet_wrap(~id, scales="free")+
  xlab("NDVI")

#' ### Distance to nearest Road
#+fig.width=8, fig.height=8
ggplot(ssfdat, aes(distroad.sc, colour=case_))+geom_density()+facet_wrap(~id, scales="free")+
  xlab("Distance to Road")

#' ### Distance to nearest Camp
#+fig.width=8, fig.height=8
ggplot(ssfdat, aes(distcamps.sc, colour=case_))+geom_density()+facet_wrap(~id, scales="free")+
  xlab("Distance to Camp")

#' ### Distance to nearest River
#+fig.width=8, fig.height=8
ggplot(ssfdat, aes(distrivers.sc, colour=case_))+geom_density()+facet_wrap(~id, scales="free")+
  xlab("Distance to River")

#' ### Distance to nearest Cattle Post
#+fig.width=8, fig.height=8
ggplot(ssfdat, aes(distcattlep.sc, colour=case_))+geom_density()+facet_wrap(~id, scales="free")+
  xlab("Distance to Cattle Posts")

#' ### Distance to nearest den
#+fig.width=8, fig.height=8
ggplot(ssfdat, aes(dendist.sc, colour=case_))+geom_density()+facet_wrap(~id, scales="free")+
  xlab("Distance to nearest den")


#' Fit an SSF to a single animal
summary(clogit(case_ ~ as.factor(vegetation)+dendist.sc+ndvi.sc+distroad.sc+distcamps.sc+
                 distrivers.sc+distcattlep.sc+sl+log(sl)+strata(step_id), data = subset(ssfdat, id=="Apollo")))

#' Fit an SSF model to data from each animal
fit_ssf <- function(data){
  mod <- clogit(case_ ~ ndvi.sc+distroad.sc+distcamps.sc+dendist.sc+ dendist.sc+
                  distrivers.sc+distcattlep.sc+sl+log(sl)+strata(step_id), data = data)
  return(mod)
}
ssffits <-ssfdat %>%  nest(-id) %>% 
  dplyr::mutate(mod = purrr::map(data, fit_ssf)) 

#' Look at first model
ssffits$mod[[1]]

#' Now, use tidy to extract information about the model fits into a nested
#' data frame
ssffits <- ssffits %>%
  dplyr::mutate(tidy = purrr::map(mod, broom::tidy),
                n = purrr::map(data, nrow) %>% simplify())
ssffits 
ssffits$tidy

#' Now, create data frame w/ the coefficients, etc
ssf_coefs2 <- ssffits %>%
  tidyr::unnest(tidy) %>%
  dplyr::select(-(std.error:p.value))  
ssf_coefs2


#' ## Summary of Individual-level Models
#' 
#' Summarize mean and SE of the coefficients (treating them as data)
pop.coef<-ssf_coefs2 %>% group_by(term) %>% 
  summarize(mean.beta=mean(estimate, na.rm=TRUE),
            se.mean.beta=sd(estimate, na.rm=TRUE)/sqrt(n()))

#' Look to see if response to den depends on proportion of time denning
denintense<-SHd %>% group_by(IND) %>% summarize(intense=mean(I(DENNING=="Y")))
                            
dencoef <-  ssf_coefs2 %>% filter(term=="dendist.sc")
dencoef<-base::merge(dencoef, denintense, by.x="id", by.y="IND")

#' Plot the relationship between the estimated coefficients for distance
#' to den versus proportion of time denning.
#' 
#+fig.width=5, fig.height=4
ggplot(dencoef, aes(intense, estimate)) + geom_point()+geom_smooth(method="lm") +
  xlab("Proportion of time denning")

#' Not a lot to see here, but perhaps individuals that are dennning most of the time
#' (intensity = 1) tend to have more negative selection coefficients (they more strongly)
#' select for locations near the den.
#' 
#' Plot coefficients
#+ fig.height=16, fig.width=16
ggplot(ssf_coefs2, aes(x=term, y=estimate, col=id)) +
  geom_dotplot(binaxis="y", stackdir="center", size=3, position=position_dodge(0.2))+
  geom_hline(yintercept=0) +
  theme(text = element_text(size = 25))+ 
  theme(axis.text.x = element_text(angle = 90))+
    scale_x_discrete(breaks= unique(pop.coef$term),
                   labels=c("Dist Den", "Dist Camp",
                            "Dist CattleP","Dist Rivers","Dist Roads","logsl", "NDVI", "sl"))


#' ## Population-level Summaries
#' 
#' Plot of mean coefficient in the population
#+ fig.height=12, fig.width=16
ggplot(pop.coef, aes(x=term, y=mean.beta)) +
  geom_errorbar(aes(ymin=mean.beta-1.96*se.mean.beta, ymax=mean.beta+1.96*se.mean.beta), width=0.1)+
  geom_point()+xlab("Predictor")+ylab(expression(paste("Mean  ", hat(beta))))+
  geom_hline(yintercept=0) + theme(text = element_text(size = 25))+
  theme(axis.text.x = element_text(angle = 90))+
  scale_x_discrete(breaks= unique(pop.coef$term),
                   labels=c("Dist Den", "Dist Camp",
                            "Dist CattleP","Dist Rivers","Dist Roads","logsl", "NDVI", "sl"))


#' ## Document Footer	
#' 	
#' 	
#' Session Information:	
#' 	
devtools::session_info()  
