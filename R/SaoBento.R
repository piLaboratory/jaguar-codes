#########################################################################################################################################               
                ###   JAGUAR DATASET - PANTANAL  (Sao  Bento) ### 
###########################################################################################################
### Preliminary tests for RSF & SSF 
###########################################################################################################
 
    ### RUN JAGUAR DATA PREPARATION FIRST !!!   NEED THE trk files!!!	

   SaoBentotrk -> trk	
################################################################################################################
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
trk %>% filter(id=="105") %>%
  random_points(.,factor = 20) %>% plot

#' Illustrate systematic points (to do this, we need to create the mcp first)
trk%>%filter(id=="105") %>% 
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
avail.df <-data.frame(avail.pts)
avail.df$x<-avail.df$x_
avail.df$y<-avail.df$y_


#' Check to make sure everything looks right 
test<-subset(avail.df, case_==TRUE)
test %>% select('x','y') %>% 
  summarise_all(mean)

SaoBento=subset(jaguar,project_region=='Sao Bento')
write.table(SaoBento,file="c:/RWorkDir/jaguardatapaper/SaoBento.txt",row.names = F,quote=F,col.names=T,sep="\t")

SaoBento %>% summarize(meanloc.x=mean(utm_x), meanloc.y=mean(utm_y))

#' Add a timestamp to annotate these data with environmental covariates in Movebank using Env-DATA (https://www.movebank.org/node/6607).
#' Here we just use the first timestamp, however meaningful timestamps are needed if annotating variables that vary in time.
#avail.df$timestamp<-jaguar_df$date[1]
avail.df$timestamp.posix<-SaoBento$date[1]

#' These points then need to be annotated prior to fitting rsfs. Let's 
#' write out 2 files:
#' 
#' - FisherRSF2018.csv will contain all points and identifying information. 
#' - FisherRSFannotate.csv will contain only the columns used to create the annotation.
#' 
#' The latter file will take up less space, making it easier to annotate (and also possible to upload to github)
#write.csv(avail.df, file="data/FisherRSF2018.csv", row.names = FALSE)
#write.csv(avail.df, file="c:/RWorkDir/jaguardatapaper/JaguarRSF_SaoBento.csv",row.names = FALSE)
#available<-avail.df %>% select("timestamp", "x", "y")
write.csv(avail.df, file="c:/RWorkDir/jaguardatapaper/RSFannotateSaoBento.csv", row.names = FALSE)
RSFannotateSaoBento<-read.csv("c:/RWorkDir/jaguardatapaper/RSFannotateSaoBento.csv")


#############################################################################################################
  #################      GIS LAYERS ###  RSF  ####################################################### 
#############################################################################################################
#rm(list= ls())                           
#unloadNamespace("amt"); unloadNamespace("tidyverse"); unloadNamespace("modelr"); unloadNamespace("broom"); unloadNamespace("tidyr")

## Load packages a few more packages                      
if(!require(install.load)) install.packages('install.load'); library(install.load)
install.load::install_load("raster","RCurl","rts","rvest","rgl","lubridate","lattice","rgdal","sp","stringr","methods",
"maptools", "vegan","sp","spatialEco","installr")



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




 #############################################################################################################
  #################      EXTRACTING RASTER VALUES ###  RSF  ####################################################### 
#############################################################################################################
################################# csv file ####################################################
 
RSFannotateSaoBento<-read.csv("c:/RWorkDir/jaguardatapaper/RSFannotateSaoBento.csv")
head(RSFannotateSaoBento)

### Coverting to points
RSFpts <- SpatialPointsDataFrame(coords = RSFannotateSaoBento[,c("x_","y_")], data = RSFannotateSaoBento,
                                            proj4string = CRS("+proj=utm +zone=21K +south +datum=WGS84 +units=m +no_defs"))
#plot(pts)
crs(RSFpts) #vendo o sistema de projecao do shape de ptos
compareCRS(tree_cover,RSFpts)
#plot(tree_cover)
#plot(RSFpts,add=T)

#### Extracting values
RSFpts$tree_cover<-extract_tree_cover<-extract(tree_cover,RSFpts)
head(RSFpts)
RSFpts$land_cover<-extract_land_cover<-extract(land_cover,RSFpts)
head(RSFpts)
RSFpts$livestock<-extract_livestock<-extract(livestock,RSFpts)
head(RSFpts)
#RSFpts$human_foot<-extract_human_foot<-extract(human_foot,RSFpts)

#add others

##### Exporting to txt file
RSFsig<-data.frame (RSFpts)
write.table(RSFsig, "RSFsig" , sep="\t", append=F, col.names = T, row.names = F, quote = F)
       
## RSF data
install.packages("tidyverse"); library("tidyverse")
#' Read in original data (used and available points) and merge on environmental data.
#call SaoBento object
#rsfdattemp <-SaoBento; head(rsfdattemp)
# or a previously created file
rsfdattemp<- read.delim(file="c:/RWorkDir/jaguardatapaper/SaoBento.txt")

#rsfdattemp <- Pantanal
annotated<-read.csv("c:/RWorkDir/jaguardatapaper/RSFsig.csv")
#rsfdat<-annotated  <-RSFsig


#' Now, merge these by "timestamp.posix", "x", "y")
#rsfdat<-merge(rsfdattemp, annotated)    <----------------------------------------------- PROBLEMS with MERGE !!!


#' Simplify some variable names
names(rsfdat)[c(8:10)]<-c("tree_cover", "land_cover","livestock") #,"human_foot", "d_drain","neo", "pop","elevation","d_water","d_roads","d_forest")
rsfdat$case_<-as.numeric(rsfdat$case_)

#' Create landcover classes (as suggested by Scott Lapoint :)
rsfdat$land_cover<-as.character(rsfdat$land_cover)

rsfdat<-rsfdat %>% mutate(landC = fct_collapse(land_cover,
      agri = c("11", "14", "30"),
      forest =c("30","40","50","60", "70","80", "90","100"),
      shrub= c("110", "130", "150"),
      grass = c("120", "140"),
      wet= c("160"),
      other = c("170", "180", "190", "200", "210", "220")))

#' Center and scale variables
rsfdat<-rsfdat %>% mutate(livestock=as.numeric(scale(livestock)), tree_coverSC=as.numeric(scale(tree_cover)))

 #d_forestSC=as.numeric(scale(d_forest)),d_roadsSC=as.numeric(scale(d_roads)),d_waterSC=as.numeric(scale(d_water)),d_drainSC=as.numeric(scale(d_drain)),popSC=as.numeric(scale(pop)))

						  
#' ## Explore the data:

#' Look Distribution of variables for used and available
#+fig.width=8, fig.height=8
binomial_smooth <- function(...) {
  geom_smooth(method = "glm", method.args = list(family = "binomial"), ...)
}

#ggplot(rsfdat,aes(x=human_foot, y=case_))+stat_smooth(method="glm", method.args = list(family = "binomial"))+binomial_smooth(formula = y ~ splines::ns(x, 5), colour="red")+facet_wrap(~id, scales="free")

ggplot(rsfdat,aes(x=livestock, y=case_))+
  stat_smooth(method="glm", method.args = list(family = "binomial"))+
  binomial_smooth(formula = y ~ splines::ns(x, 5), colour="red")+
  facet_wrap(~id, scales="free")
  
#ggplot(rsfdat,aes(x=elevation, y=case_))+stat_smooth(method="glm", method.args = list(family = "binomial"))+binomial_smooth(formula = y ~ splines::ns(x, 5), colour="red")+facet_wrap(~id, scales="free")

#ggplot(rsfdat,aes(x=pop, y=case_))+stat_smooth(method="glm", method.args = list(family = "binomial"))+binomial_smooth(formula = y ~ splines::ns(x, 5), colour="red")+facet_wrap(~id, scales="free")
  
ggplot(rsfdat,aes(x=tree_cover, y=case_))+
  stat_smooth(method="glm", method.args = list(family = "binomial"))+
  binomial_smooth(formula = y ~ splines::ns(x, 5), colour="red")+
  facet_wrap(~id, scales="free")
 
  
#ggplot(rsfdat,aes(x=d_forest, y=case_))+ stat_smooth(method="glm", method.args = list(family = "binomial"))+ binomial_smooth(formula = y ~ splines::ns(x, 5), colour="red")+ facet_wrap(~id, scales="free")
  
#ggplot(rsfdat,aes(x=d_roads, y=case_))+stat_smooth(method="glm", method.args = list(family = "binomial"))+binomial_smooth(formula = y ~ splines::ns(x, 5), colour="red")+facet_wrap(~id, scales="free")
  
#ggplot(rsfdat,aes(x=d_water, y=case_))+stat_smooth(method="glm", method.args = list(family = "binomial"))+binomial_smooth(formula = y ~ splines::ns(x, 5), colour="red")+facet_wrap(~id, scales="free")
  
#ggplot(rsfdat,aes(x=d_drain, y=case_))+stat_smooth(method="glm", method.args = list(family = "binomial"))+binomial_smooth(formula = y ~ splines::ns(x, 5), colour="red")+facet_wrap(~id, scales="free")




  

ggplot(rsfdat, aes(x=landC, y=..prop..,group=case_, colour=case_))+
  geom_bar(position="dodge", aes(fill=case_))+facet_wrap(~id, scales="free")




#' ## RSF fitting

#' Weight available data 
rsfdat$w<-ifelse(rsfdat$case_==1, 1, 5000)

#' We can fit an RSF model to a single animal using logistic regression
summary(glm(case_ ~tree_coverSC+livestock +landC, data = subset(rsfdat, id=="115"), weight=w,family = binomial))

 #elevSC+popSC+d_forestSC+d_roadsSC+d_waterSC+d_drainSC


#' landcover classes use
rsfdat %>% filter(id=="115") %>% with(table(case_, landC))  
rsfdat$used<-as.factor(rsfdat$case_)
rsfdat$used<-fct_recode(rsfdat$used, "avail"="0", "used"="1")

#+fig.width=6, fig.height=4
ggplot(subset(rsfdat, id=="115"),  aes(x=landC,group=used))+
  geom_bar(position=position_dodge(), aes(y=..prop.., fill = used), stat="count") +
  scale_fill_brewer(palette="Paired")+
  geom_text(aes( label = scales::percent(..prop..),
                 y= ..prop.. ), stat= "count", vjust = -.3, position=position_dodge(0.9)) +
  labs(y = "Proportion", fill="used", x="Landcover") 

#' Now, fit an RSF model to data from each animal.  Since not all animals experience
#' all habitat types, lets just explore forest versus non-forest
rsfdat$forest<-ifelse(rsfdat$landC=="forest", 1, 0)

fit_rsf <- function(data){
  mod <- glm(case_ ~ tree_cover+livestock+forest, data = data, weight=w,family = binomial)
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

#########################################################################################################################################
################################################################################################################
#' ## SSF prep
#' 
SaoBentotrk -> trk	
####' SSFs assume that data have been collected at regular time intervals.
#' We can use the track_resample function to regularize the trajectory so that
#' all points are located within some tolerence of each other in time. To figure
#' out a meaningful tolerance range, we should calculate time differences between
#' locations & look at as a function of individual.

(timestats<-trk %>% nest(-id,-sex,-age,-weight,-status) %>% mutate(sr = map(data, summarize_sampling_rate, time_unit = "hour")) %>%
    select(id,sex,age,weight,status,sr) %>% unnest)
	
## Lets add on the time difference to each obs.
trk<-trk %>% group_by(id) %>% mutate(dt_ = t_ - lag(t_, default = NA))
trk
#' Let's illustrate track regularization with ID = 115. Let's
#' keeping is as 3 hours with tolerance of 3 hours
temp115<-trk %>% filter(id=="115") %>% track_resample(rate=hours(3), tolerance=hours(0))
temp115 %>% select(id, x_, y_, t_, burst_)
#Just checking it
summarize_sampling_rate(temp115)

#' #' Regularization of ID = 110 which had the highest median intervals
# Changing from 6 hours to 3 hours, with tolerance of 3 hours
temp110<-trk %>% filter(id=="110") %>% track_resample(rate=hours(3), tolerance=hours(0))
temp110 %>% select(id, x_, y_, t_, burst_)
#Just checking it
summarize_sampling_rate(temp110)


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

ssfdat<-NULL
temptrk<-with(trk, track(x=x_, y=y_, t=t_, id=id)) ###,Event_ID=Event_ID,project_region=project_region, sex=sex, age=age, weight=weight,status=status, period=period))
uid<-unique(trk$id) # individual identifiers
luid<-length(uid) # number of unique individuals
for(i in 1:luid){
  # Subset individuals & regularize track 
  temp<-temptrk%>% filter(id==uid[i]) %>% 
  track_resample(rate=hours(round(timestats$median[i])), 
                   tolerance=hours(max(3,round(timestats$median[i]/6))))
				   
  #temp<-temptrk%>% filter(id==uid[i]) %>% 
    #track_resample(rate=hours(8), tolerance=hours(4)) 
### This type of approach would discard many points and increase the median time differences
	
  # Get rid of any bursts without at least 3 points
  temp<-filter_min_n_burst(temp, min_n =3) 
  
  # burst steps
  stepstemp<-steps_by_burst(temp)
  
  # create random steps using fitted gamma and von mises distributions and append
  rnd_stps <- stepstemp %>% random_steps(n = 15)
  
###Here we may include the adjust_param functions from amt package
###(adjust_shape, adjust_scale, adjust_scale, adjust_kappa) to create a iSSF
  
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
head(bursts)  

  	
#SSFs sampled data (x1, y1, t1)
ssfdattrk <- mk_track(ssfdatemp, .x=x1_, .y=y1_, .t=t1_, id=id, crs = CRS("+proj=utm +zone=21K +south +datum=WGS84 +units=m +no_defs")) 
##Summary of totals from the original trk(as we did above just to compare)
(timestatsssf<-ssfdattrk %>% nest(-id) %>% mutate(sr = map(data, summarize_sampling_rate, time_unit = "hour")) %>%
    select(id, sr) %>% unnest)
	
### SSFs sampled data (x2, y2, t2)
ssfdattrk <- mk_track(ssfdatemp, .x=x2_, .y=y2_, .t=t2_, id=id, crs = CRS("+proj=utm +zone=21K +south +datum=WGS84 +units=m +no_defs")) 
##Summary of totals from the original trk(as we did above just to compare)
(timestatsssf<-ssfdattrk %>% nest(-id) %>% mutate(sr = map(data, summarize_sampling_rate, time_unit = "hour")) %>%
    select(id, sr) %>% unnest)
## And comparing that to the original file
SaoBentotrk -> trk	
(timestats<-trk %>% nest(-id,-sex,-age,-weight,-status) %>% mutate(sr = map(data, summarize_sampling_rate, time_unit = "hour")) %>%
    select(id,sex,age,weight,status,sr) %>% unnest)


#' Now, lets plot the data for random and matched points
#' 
#+fig.height=12, fig.width=12, warning=FALSE
ggplot(ssfdat, aes(x2_, y2_, color=case_))+geom_point()+facet_wrap(~id,scales="free")

##########################################################################################################
#' ## Write out data for further annotating
#' #######################################################################################################

#'  With the SSFs, we have the extra complication of
#' having a time and location at both the start and end of the step.  
#' 
#' For the time being, we will assume we want to annotate variables at the end of the step
#' but use the starting point of the step as the timestamp.

ssf.df <- data.frame(ssfdat)
head(ssf.df)
ssf.df$timestamp <-ssf.df$t1_
ssf.df$x <-ssf.df$x2_
ssf.df$y <-ssf.df$y2_
head(ssf.df)
ssfALL<-ssf.df
head(ssfALL)
#ssfALL$sl_ <- NULL; ssfALL$ta_<- NULL
ssfALL<-ssfALL%>% mutate(obs = 1:n()) %>% select(obs, everything())
head(ssfALL)

stepSB<-ssfALL %>% select("obs","id","case_", "x", "y","timestamp","step_id_")
head(stepSB)

#AllSteps <- ssf.df  ### include name of the area/ project
#' These points then need to be annotated prior to fitting ssfs. 2 files All and to annotate

write.csv(ssfALL, file="C:/RWorkDir/jaguardatapaper/ssfALL.csv", row.names = FALSE) ###annotated_projectname
write.csv(stepSB, file="C:/RWorkDir/jaguardatapaper/stepSB.csv", row.names = FALSE) ###annotated_projectname

#stepSB<- read.csv(file="c:/RWorkDir/jaguardatapaper/stepSB.csv")
#ssfALL<- read.csv(file="c:/RWorkDir/jaguardatapaper/ssfALL.csv")



#############################################################################################################
  #################      GIS LAYERS ###  SSF  ####################################################### 
#############################################################################################################
#rm(list= ls())                                                
## Commented because is not generic
#setwd("C:/RWorkDir/jaguardatapaper")                           ### Set directory

require(raster)
require (spatialEco)
require(sp)
require(rgdal)

### Check if raster have been uploaded and projected to UTM
land_cover; crs(land_cover) 
plot(land_cover)

tree_cover; crs(tree_cover)
plot(tree_cover)

livestock; crs(livestock)
plot(livestock)

###  add others

### Import table points
stepSB<- read.csv(file="c:/RWorkDir/jaguardatapaper/stepSB.csv")
head(stepSB)

### Coverting to points
SSFpts <- SpatialPointsDataFrame(coords = stepSB[,c("x","y")], data = stepSB,
                                            proj4string = CRS("+proj=utm +zone=21K +south +datum=WGS84 +units=m +no_defs"))

crs(SSFpts) #vendo o sistema de projecao do shape de ptos
compareCRS(tree_cover,SSFpts)
#plot(tree_cover)
#plot(RSFpts,add=T)

#### Extracting values  ###  IF USE EXTRACTION DID NOT WORK IT MAY BE REQUIRED EITHER TO USE rm(list= ls()) or restart R.
#  and then run only a few required packages such as raster, sp, rgdal, spatialEco...  RUN amt package ONLY AFTER THE EXTRACTIONS,
# BECAUSE CONFLICYS!!!
### Another possible solution may be the use of extract_covariates command from amt package 


SSFpts$tree_cover<-extract_tree_cover<-extract(tree_cover,SSFpts)
head(SSFpts)
SSFpts$land_cover<-extract_land_cover<-extract(land_cover,SSFpts)
head(SSFpts)
SSFpts$livestock<-extract_livestock<-extract(livestock,SSFpts)
head(SSFpts)

# ADD OTHERS
#SSFpts$human_foot<-extract_human_foot<-extract(human_foot,SSFpts)
#shp_albers$land_use<-extract_land_use<-extract(land_use,shp_albers)
#head(shp_albers)
#shp_albers$tree_cover<-extract_tree_cover<-extract(tree_cover,shp_albers)
#head(shp_albers)
#shp_albers$land_cover<-extract_land_cover<-extract(land_cover,shp_albers)
#head(shp_albers)
#shp_albers$human_foot<-extract_human_foot<-extract(human_foot,shp_albers)
#head(shp_albers)
#shp_albers$livestock<-extract_livestock<-extract(livestock,shp_albers)
#head(shp_albers)
#shp_albers$pop<-extract_pop<-extract(pop,shp_albers)
#head(shp_albers)
#shp_albers$elevat<-extract_elevat<-extract(elevat,shp_albers)
#head(shp_albers)
#shp_albers$water30m<-extract_water30m<-extract(water30m,shp_albers)
#head(shp_albers)
#shp_albers$distroad<-extract_distroad<-extract(distroad,shp_albers)
#head(shp_albers)
#shp_albers$distdrain<-extract_distdrain<-extract(distdrain,shp_albers)
#head(shp_albers)
#shp_albers$distwater<-extract_distwater<-extract(distwater,shp_albers)
#head(shp_albers)
#shp_albers$disttreecover<-extract_disttreecover<-extract(disttreecover,shp_albers)
#head(shp_albers)

ssf.df$timestamp<-ssf.df$t1_
ssf.df <- data.frame(ssfdat)
head(ssf.df)
AllSteps <- ssf.df  ### include name of the area/ project  (e.g. iopCerradoAllSteps)
AllSteps ->ssfALL

write.csv(ssf.df, file="C:/RWorkDir/jaguardatapaper/AllSteps.csv", row.names=FALSE) ###AllSteps_projectname
#ssf.df<-ssf.df %>% select("timestamp", "long", "lat")
#write.csv(ssf.df, file="C:/RWorkDir/jaguardatapaper/annotated_iopCerrado.csv", row.names = FALSE) ###annotated_projectname



##### Exporting to txt file
ssfSB<-data.frame (SSFpts)
write.table(ssfSB,file="c:/RWorkDir/jaguardatapaper/ssfSB.txt",row.names = F,quote=F,col.names=T,sep="\t")

head(ssfSB)
#add_row_numbers(ssfSB, name = "row_number", zero_based = FALSE)

################################################################################################################
                                   ###  MERGE ###     
################################################################################################################
## SSF data

#' Read in original data (used and available points) and merge on environmental data.
ssfSB<- read.delim(file="c:/RWorkDir/jaguardatapaper/ssfSB.txt") ###annotated_projectname 
ssfALL<- read.csv(file="c:/RWorkDir/jaguardatapaper/ssfALL.csv") ###AllSteps_projectname

ssfSB$id <- as.character(ssfSB$id)
ssfALL$id <- as.character(ssfALL$id)
ssfSB$obs <- as.integer(ssfSB$obs)
ssfALL$obs <- as.integer(ssfALL$obs)
ssfSB$timestamp <-as.character(ssfSB$timestamp)
ssfSB$timestamp <-as.POSIXct(ssfSB$timestamp, format ="%Y-%m-%d %H:%M:%S", tz = 'GMT') 
ssfALL$timestamp <-as.character(ssfSB$timestamp)
ssfALL$timestamp <-as.POSIXct(ssfALL$timestamp, format ="%Y-%m-%d %H:%M:%S", tz = 'GMT') 



# Now, merge these by "id","timestamp", "x", "y")
ssfSBGIS<-merge (ssfSB,ssfALL)  ### NEED check this!!! #, by = c("id","timestamp"), all.x = TRUE)        #c("id","timestamp","long","lat","obs"))   
head(ssfSBGIS)
str(ssfSBGIS)

ssfSBGIS$id.y <- NULL; ssfSBGIS$timestamp.y<- NULL; ssfSBGIS$long.y<- NULL; 
ssfSBGIS$lat.y<- NULL; ssfSBGIS$step_id_.y<- NULL
ssfSBGIS$case_.y <- NULL

ssfSBGIS <- ssfSBGIS %>% rename(case=case_.x)
ssfSBGIS <- ssfSBGIS %>% rename(long=long.x)
ssfSBGIS <- ssfSBGIS %>% rename(lat=lat.x)
ssfSBGIS <- ssfSBGIS %>% rename(timestamp=timestamp.x)
ssfSBGIS <- ssfSBGIS %>% rename(step_id=step_id_.x)
ssfSBGIS <- ssfSBGIS %>% rename(id=id.x)

 
#' Write out data for use in SSF analises
write.csv(ssfSBGIS, file="C:/RWorkDir/jaguardatapaper/ssfSBGIS.csv", row.names=FALSE)
ssfSBGIS<- read.csv(file="c:/RWorkDir/jaguardatapaper/ssfSBGIS.csv")

#' Now, merge these by "timestamp", "location-long", "location-lat")
#ssfdat<-merge(ssfdattemp,annotated)
#' Write out data for use in SSF analises
#write.csv(ssfdat, file="C:/RWorkDir/jaguardatapaper/iopCerrado_AllSteps_GIS.csv", row.names=FALSE)


###############################################################################################################
 ######                          SSF script  AND Analysis                                            #########
###############################################################################################################
### use libraries as above

#' Read in annotated available data for SSF modeling
#ssfdat<-read.csv("C:/RWorkDir/jaguardatapaper/ssfSBGIS.csv")  ### NEED CHECK MERGE FIRST!
ssfdat<-read.csv("C:/RWorkDir/jaguardatapaper/ssfSB.csv")


#' Convert time variables
ssfdat$t1_<-as.POSIXct(ssfdat$t1_, format ="%Y-%m-%d %H:%M:%S", tz="UTC")
ssfdat$t2_<-as.POSIXct(ssfdat$t2_, format ="%Y-%m-%d %H:%M:%S", tz="UTC")
ssfdat$timestamp<-as.POSIXct(ssfdat$timestamp, format ="%Y-%m-%d %H:%M:%S", tz="UTC")
head(ssfdat)
str(ssfdat)
   
# Organize variables
#numeric
ssfdat$case_<-as.numeric(ssfdat$case)
#ssfdat$hansen<-as.numeric(ssfdat$hansen)
ssfdat$tree_cover<-as.numeric(ssfdat$tree_cover)
#ssfdat$water30m<-as.numeric(ssfdat$water30m)
#ssfdat$human_foot<-as.numeric(ssfdat$human_foot)
ssfdat$livestock<-as.numeric(ssfdat$livestock)
#ssfdat$pop<-as.numeric(ssfdat$pop)
#ssfdat$elevat<-as.numeric(ssfdat$elevat)
#ssfdat$distroad<-as.numeric(ssfdat$distroad)
#ssfdat$distdrain<-as.numeric(ssfdat$distdrain)
#ssfdat$distwater<-as.numeric(ssfdat$distwater)
#ssfdat$disttreecover<-as.numeric(ssfdat$disttreecover)

head(ssfdat)


#Character
ssfdat$id<-as.character(ssfdat$id)

#' landcover classes 
ssfdat$land_cover <-as.factor(ssfdat$land_cover)


ssfdat<-ssfdat %>% mutate(land_C = fct_collapse(land_use,
                                               crop = c("10", "20"),
											     mosaic_crop = c("30"),
                                               forest =c("60","70","80"),
											     mosaic_treeshrub =c("100"),
												 mosaic_herb =c("110"),
												 sparse_tsh =c("150"),
											   flood_tree =c("160","170"),
											   flood_shrubherb = c("180"),
                                               shrub= c("120"),
                                               grass = c("130"),
                                               water= c("210"),
											   urban = c("190"),
                                               bare = c("200")))
											   

#Value Label: 0 No Data 
#      10 Cropland, rainfed 
 #     20 Cropland, irrigated or post-flooding 
  #    30 Mosaic cropland (>50%) / natural vegetation (tree, shrub, herbaceous cover) (50%) / cropland (15%) 
   #   60 Tree cover, broadleaved, deciduous, closed to open (>15%)
    #  70 Tree cover, needleleaved, evergreen, closed to open (>15%) 
     # 80 Tree cover, needleleaved, deciduous, closed to open (>15%) 90 Tree cover, mixed leaf type (broadleaved and needleleaved)
 #    100 Mosaic tree and shrub (>50%) / herbaceous cover (<50%) 
  #   110 Mosaic herbaceous cover (>50%) / tree and shrub (<50%)
   #   120 Shrubland 
    #  130 Grassland 
     #                                      140 Lichens and mosses 
      #150 Sparse vegetation (tree, shrub, herbaceous cover) (<15%) 
      #160 Tree cover, flooded, fresh or brakish water
      #170 Tree cover, flooded, saline water 
      #180 Shrub or herbaceous cover, flooded, fresh/saline/brakish water 
    #190 Urban areas 
    #200 Bare areas 
    #210 Water bodies 
    #                                 220 Permanent snow and ice
											   
											   

											   

#' Center and scale variables
ssfdat<-ssfdat %>% mutate(f_coverSC=as.numeric(scale(tree_cover)),livestockSC=as.numeric(scale(livestock))

#elevSC=as.numeric(scale(elevation)), f_coverSC=as.numeric(scale(tree_cover)),
                          human_footSC=as.numeric(scale(human_foot)), waterSC=as.numeric(scale(water30m)), 
                          d_forestSC=as.numeric(scale(disttreecover)),d_roadsSC=as.numeric(scale(distroad)),
						  d_waterSC=as.numeric(scale(distwater)),d_drainSC=as.numeric(scale(distdrain)),
                          popSC=as.numeric(scale(pop))) 
						  


#' ## Explore the data:

#' Look Distribution of variables for used and available
#+fig.width=8, fig.height=8
binomial_smooth <- function(...) {
  geom_smooth(method = "glm", method.args = list(family = "binomial"), ...)
}
#ggplot(ssfdat,aes(x=elevat, y=case_))+
  stat_smooth(method="glm", method.args = list(family = "binomial"))+
  binomial_smooth(formula = y ~ splines::ns(x, 5), colour="red")+
  facet_wrap(~id, scales="free")

#ggplot(ssfdat,aes(x=pop, y=case_))+
  stat_smooth(method="glm", method.args = list(family = "binomial"))+
  binomial_smooth(formula = y ~ splines::ns(x, 5), colour="red")+
  facet_wrap(~id, scales="free")

ggplot(ssfdat, aes(x=land_cover, y=..prop..,group=case_, colour=case_))+
  geom_bar(position="dodge", aes(fill=case_))+facet_wrap(~id, scales="free")


#' Now, fit an SSF model to data from each animal.  Since not all animals experience
#' all habitat types, lets just explore forest versus non-forest
#ssfdat$forest<-ifelse(ssfdat$land_cover=="forest", 1, 0)

#' Fit an SSF to a single animal
summary(fit_issf(case_ ~ elevSC+popSC+f_coverSC+human_footSC+waterSC+d_forestSC+d_roadsSC+d_waterSC+d_drainSC+strata(step_id), 
                 data = subset(ssfdat, id=="105")))
summary(fit_issf(case_ ~ elevSC+popSC+f_coverSC+human_footSC+waterSC+d_forestSC+d_roadsSC+d_waterSC+d_drainSC+strata(step_id), 
                 data = subset(ssfdat, id=="111")))
summary(fit_issf(case_ ~ elevSC+popSC+f_coverSC+human_footSC+waterSC+d_forestSC+d_roadsSC+d_waterSC+d_drainSC+strata(step_id), 
                 data = subset(ssfdat, id=="113")))
summary(fit_issf(case_ ~ elevSC+popSC+f_coverSC+human_footSC+waterSC+d_forestSC+d_roadsSC+d_waterSC+d_drainSC+strata(step_id), 
                 data = subset(ssfdat, id=="114")))
summary(fit_issf(case_ ~ elevSC+popSC+f_coverSC+human_footSC+waterSC+d_forestSC+d_roadsSC+d_waterSC+d_drainSC+strata(step_id), 
                 data = subset(ssfdat, id=="115")))
        

#' Note, animal 1587 was always in forest, so we can't include forest 
#' for this individual.

#' Fit an SSF model to data from each animal
fitted_ssf <- function(data){
  n0<-sum(data$forest==0) 
  if(n0==0){ #  this will be id = M5
    mod <- fit_issf(case_ ~ Elevation+PopDens+sl_+log(sl_)+strata(step_id_), data=data)
  }
  if(n0>0){
  mod <- fit_issf(case_ ~ Elevation+PopDens+forest+sl_+log(sl_)+strata(step_id_), data=data)
  }
  return(mod)
}
ssffits <-ssfdat %>%  nest(-id) %>% 
  dplyr::mutate(mod = purrr::map(data, fitted_ssf)) 

#' Look at first model
ssffits$mod[[1]]

#' Now, use tidy to extract information about the model fits
ssffits <- ssffits %>%
  dplyr::mutate(tidy = purrr::map(mod, ~broom::tidy(.x$model)),
                n = purrr::map(data, nrow) %>% simplify())

ssffits$tidy

#' Now, create data frame w/ the coefficients, etc
ssf_coefs <- ssffits %>%
  tidyr::unnest(tidy) %>%
  dplyr::select(-(std.error:conf.high)) 
  
ssf_coefs %>% tidyr::spread(term, estimate)

#' Plot coefficients
#+ fig.width=12, fig.height= 8
ssf_coefs %>% 
  ggplot(., aes(x=1, y=estimate)) + 
  geom_dotplot(binaxis="y", stackdir="center")+geom_hline(yintercept=0)+
  facet_wrap(~term, scales="free")


#' ## Document Footer	
#' 	
#' Document spun with:  ezspin("FisherSSF.R",  fig_dir = "figures", keep_md=FALSE)  	
#' 	
#' Session Information:	
#' 	
sessionInfo()	  


################# Function to extrat points #################

### Covertendo a tabela para arquivo shape de pontos
#' @param points tabela com os pontos gerados ...
#' @param paths.to.raster vector with the path names for the rasters files of each environmental variable
#' @param crs.points CRS object, with the crs of the point data
#' @details this function ...
env.points=function(points, paths.to.raster, 
                    crs.points=CRS("+proj=longlat +datum=WGS84")){
  pts <- SpatialPointsDataFrame(coords = points[,c("long","lat")],
                                data = points,
                                proj4string = crs.points )
  output=matrix(nrow=nrow(points), ncol=length(paths.to.raster))
  for(i in 1:length(paths.to.raster)){
    env.raster = raster(paths.to.raster[i])
    shp.pts<- spTransform(pts, crs(env.raster))
    output[,i] = extract(env.raster,shp.pts)
  }
data.frame(cbind(shp.pts,output))  
}

teste=env.points(tabela, 
                 paths.to.raster = c("C:/SIG/b115_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif", "C:/SIG/b115_Elevation_md_SRTM_1km_neotropic_albers_tif_exp.tif"))

pts <- SpatialPointsDataFrame(coords = tabela[,c("long","lat")], data = tabela,
                                            proj4string = CRS("+proj=longlat +datum=WGS84"))
#plot(pts)
crs(pts) #vendo o sistema de projecao do shape de ptos

###Colocando o shape de points na mesma projecao do raster

shp_albers<- spTransform(pts, crs(hansen))

compareCRS(hansen,shp_albers)

###ploting
plot(shp_albers, add=TRUE)

### exportando o shape com a projecao do raster
writeOGR(shp_albers, dsn = '.',"points_albers", driver="ESRI Shapefile") ## exportar para arcgis so pontos selecionados

#### Extraindo o valor do raster em cada pto
shp_albers$hansen<-extract_hansen<-extract(hansen,shp_albers)
head(shp_albers)
shp_albers$elevation<-extract_elevation<-extract(elevation,shp_albers)
head(shp_albers)
shp_albers$land_use<-extract_land_use<-extract(land_use,shp_albers)
head(shp_albers)

##### Exportando o dbf do shape como uma tabela txt
data_final<-data.frame (shp_albers)
write.table(data_final, "Results_table_atributos_all_points.txt" , sep="\t", append=F, col.names = T, row.names = F, quote = F)




