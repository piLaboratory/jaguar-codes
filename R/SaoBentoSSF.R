#'
#' #  **JAGUAR DATASET - PANTANAL  (Sao  Bento)**
#' 
#' #### *Alan E. de Barros, Bernardo Niebuhr, Vanesa Bejarano, Julia Oshima,Claudia Kanda, Milton Ribeiro, Ronaldo Morato,Paulo Prado*
#' date: "April, 04 2019"
#' ##### Scripts adapted from Johannes Signer and John John Fieberg's lectures.
#' 
#' #### Run JaguarDataPrep first !!! 
source("JaguarDataPrep.R")
#'
#' ### Preliminary tests for SSF using data from Sao Bento, Pantanal 
options(width=165,digits.secs = 3)
opts_chunk$set(fig.width=12,fig.height=4.5, error=TRUE,cache = F)

#' Record time for running all code
ptm<-proc.time()

#' Set the seed for the random number generator, so it will be possible
#' to reproduce the random points
set.seed(10299)
options(max.print=1000000)

SaoBentotrk -> trk	
trk # %>% print(trk, n = Inf)
summary(trk)
colSums(ifelse(table(trk$id,trk$sex)>0,1,0)) #table(droplevels(trk$id),trk$sex)  # sex
colSums(ifelse(table(trk$id,trk$age)>0,1,0)) #table(droplevels(trk$id),trk$age)  # age
colSums(ifelse(table(trk$id,trk$weight)>0,1,0)) #table(droplevels(trk$id),trk$weight) # weight
colSums(ifelse(table(trk$id,trk$status)>0,1,0)) #table(droplevels(trk$id),trk$status) # residence status


#' ### Visualize all data with leaflet andd plots
pal <- colorFactor(palette = 'Paired',domain = SaoBento$id)
leaflet(SaoBento)%>%addTiles()%>%addCircles(SaoBento$x,SaoBento$y,color=~pal(id))%>% addProviderTiles(providers$Esri.WorldImagery)
   
#' ##### 3D plot
open3d()
# To get a bigger window than the default
par3d(windowRect = c(100, 100, 612, 612))
Sys.sleep(0.1) # Allow sluggish window managers to catch up
with(SaoBento, plot3d(x,y,date, type="l", col=as.factor(SaoBento$id)))
(stcube<-with(SaoBento, plot3d(x,y,date, type="l",col=as.numeric(cut(SaoBento$weight,11)), alpha=0.4)))
# Or with points
open3d()
# To get a bigger window than the default
par3d(windowRect = c(100, 100, 612, 612))
Sys.sleep(0.1) # Allow sluggish window managers to catch up

with(SaoBento, plot3d(x,y,date, type="p", col=as.factor(SaoBento$id)))
(stcube<-with(SaoBento, plot3d(x,y,date, type="p",col=as.numeric(cut(SaoBento$weight,11)), alpha=0.4)))

#' ## Some plots of movement characteristics
#' 
#' We can select id and steps, unnest the new data_frame and create a plot of the step-length distributions.
trk %>% select(id, sl) %>% unnest %>% 
  ggplot(aes(sl, fill = factor(id))) + geom_density(alpha = 0.4)

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
ggplot(trk, aes(x = tod_, y = log(sl)))+ geom_boxplot()+geom_smooth()+facet_wrap(~id)


#' ## GIS LAYERS in UTM  **These NEED BE ADJUSTED TO EACH INDIVIDUAL!!!**
#  For now we use bases from ID 115 because its 70 Km radius emcompasses all the other individuals
#' 
#' Distances **These NEED BE RE-CALCULATED!!!**
dist2drainage <-"D:/GISUTM/J115/UTM_b115_dist2drainage_exp.tif"
dist2water<-"D:/GISUTM/J115/UTM_b115_dist2waterbodies_exp.tif"
dist2forest <-"D:/GISUTM/J115/UTM_b115_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist.tif"
dist2roads <-"D:/GISUTM/J115/UTM_b115_dist2roads_exp.tif"
#' Anthropic
human_footprint<-"D:/GISUTM/J115/UTM_b115_human_footprint_2009_1km_tif_exp.tif"
Livestock<-"D:/GISUTM/J115/UTM_b115_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp.tif"
popdensity<-"D:/GISUTM/J115/UTM_b115_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp.tif"
#' Landscape classificatiom and physical
Landcover <-"D:/GISUTM/J115/UTM_b115_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp.tif"
treecover <-"D:/GISUTM/J115/UTM_b115_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp.tif"
topo<-"D:/GISUTM/J115/UTM_b115_Neotropic_Earthenv_dem90m_tif_exp.tif"
water<-"D:/GISUTM/J115/UTM_b115_water_frequency_2010_30m_tif_exp.tif"
#' Check projection and images
dist2drainage=raster(dist2drainage);crs(dist2drainage); plot(dist2drainage)
dist2water=raster(dist2water); crs(dist2water); plot(dist2water)
dist2forest=raster(dist2forest); crs(dist2forest); plot(dist2forest)
dist2roads=raster(dist2roads); crs(dist2roads); plot(dist2roads)
human_footprint=raster(human_footprint); crs(human_footprint); plot(human_footprint)
Livestock=raster(Livestock); crs(Livestock); plot(Livestock)
popdensity=raster(popdensity); crs(popdensity); plot(popdensity)
Landcover=raster(Landcover); crs(Landcover); plot(Landcover)
treecover=raster(treecover); crs(treecover); plot(treecover)
topo=raster(topo); crs(topo); plot(topo)
water=raster(water); crs(water); plot(water)
#'

#' ## SSF prep
#' 
SaoBentotrk -> trk	
####' SSFs assume that data have been collected at regular time intervals.
#' We can use the track_resample function to regularize the trajectory so that
#' all points are located within some tolerence of each other in time. To figure
#' out a meaningful tolerance range, we should calculate time differences between
#' locations & look at as a function of individual.
(timestats<-trk %>% nest(-id,-sex,-age,-weight,-status,-project_region) %>% mutate(sr = map(data, summarize_sampling_rate, time_unit = "hour")) %>%
    select(id,sex,age,weight,status,project_region,sr) %>% unnest)
	
## Lets add on the time difference to each obs.  Alerady DONE!!! 
#trk<-trk %>% group_by(id) %>% mutate(dt_ = t_ - lag(t_, default = NA))
#trk

#' Jaguar 105
#' Let's illustrate track regularization with ID = 105. Let's
#' keeping is as 3 hours with tolerance of 1 hours
temp105<-trk %>% filter(id=="105") %>% track_resample(rate=hours(3), tolerance=hours(1))
temp105 %>% select(id, x_, y_, t_, burst_)

temp105trk<-with(trk, track(x=x_, y=y_, t=t_, id=id)) ###,Event_ID=Event_ID,project_region=project_region, sex=sex, age=age, weight=weight,status=status, period

#Just checking it
summarize_sampling_rate(temp105)

#' Before fitting a step selection, the data well need to prepared. First, we change from a point representation to a step representation, using the function steps_by_burst, which in contrast to the steps function accounts for bursts.
ssf105 <- temp105 %>% steps_by_burst()

#' Next, we generate random steps with the function random_steps. This function fits by default a Gamma distribution to the step lengths and a von Mises distribution to the turn angles, and then pairs each observed step with n random steps.
ssf105 <- ssf105 %>% random_steps(n = 15)

#' As a last step, we have to extract the covariates at the end point of each step. We can do this with extract_covariates.
ssf105 <- ssf105 %>% extract_covariates(dist2drainage); ssf105 <- rename(ssf105, dist2drainage = UTM_b115_dist2drainage_exp)
ssf105 <- ssf105 %>% extract_covariates(dist2water); ssf105 <- rename(ssf105, dist2water = UTM_b115_dist2waterbodies_exp) 
ssf105 <- ssf105 %>% extract_covariates(dist2forest);ssf105 <- rename(ssf105, dist2forest =UTM_b115_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist)  
ssf105 <- ssf105 %>% extract_covariates(dist2roads); ssf105 <- rename(ssf105, dist2roads = UTM_b115_dist2roads_exp) 
ssf105 <- ssf105 %>% extract_covariates(human_footprint); ssf105 <- rename(ssf105, human_footprint = UTM_b115_human_footprint_2009_1km_tif_exp) 
ssf105 <- ssf105 %>% extract_covariates(Livestock); ssf105 <- rename(ssf105, Livestock = UTM_b115_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp) 
ssf105 <- ssf105 %>% extract_covariates(popdensity); ssf105 <- rename(ssf105, popdensity =UTM_b115_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp) 
ssf105 <- ssf105 %>% extract_covariates(Landcover); ssf105 <- rename(ssf105, Landcover = UTM_b115_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp) 
ssf105 <- ssf105 %>% extract_covariates(treecover); ssf105 <- rename(ssf105, treecover = UTM_b115_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp) 
ssf105 <- ssf105 %>% extract_covariates(topo); ssf105 <- rename(ssf105, topo = UTM_b115_Neotropic_Earthenv_dem90m_tif_exp)  
ssf105 <- ssf105 %>% extract_covariates(water); ssf105 <- rename(ssf105, water = UTM_b115_water_frequency_2010_30m_tif_exp)  

print(ssf105, width=Inf)
#'

#' Jaguar 106
#' Let's illustrate track regularization with ID = 106. Let's
#' keeping is as 3 hours with tolerance of 1 hours
temp106<-trk %>% filter(id=="106") %>% track_resample(rate=hours(3), tolerance=hours(1))
temp106 %>% select(id, x_, y_, t_, burst_)

#Just checking it
summarize_sampling_rate(temp106)
#' Before fitting a step selection, the data well need to prepared. First, we change from a point representation to a step representation, using the function steps_by_burst, which in contrast to the steps function accounts for bursts.
ssf106 <- temp106 %>% steps_by_burst()

#' Next, we generate random steps with the function random_steps. This function fits by default a Gamma distribution to the step lengths and a von Mises distribution to the turn angles, and then pairs each observed step with n random steps.
ssf106 <- ssf106 %>% random_steps(n = 15)

#' As a last step, we have to extract the covariates at the end point of each step. We can do this with extract_covariates.
ssf106 <- ssf106 %>% extract_covariates(dist2drainage); ssf106 <- rename(ssf106, dist2drainage = UTM_b115_dist2drainage_exp)
ssf106 <- ssf106 %>% extract_covariates(dist2water); ssf106 <- rename(ssf106, dist2water = UTM_b115_dist2waterbodies_exp) 
ssf106 <- ssf106 %>% extract_covariates(dist2forest);ssf106 <- rename(ssf106, dist2forest =UTM_b115_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist)  
ssf106 <- ssf106 %>% extract_covariates(dist2roads); ssf106 <- rename(ssf106, dist2roads = UTM_b115_dist2roads_exp) 
ssf106 <- ssf106 %>% extract_covariates(human_footprint); ssf106 <- rename(ssf106, human_footprint = UTM_b115_human_footprint_2009_1km_tif_exp) 
ssf106 <- ssf106 %>% extract_covariates(Livestock); ssf106 <- rename(ssf106, Livestock = UTM_b115_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp) 
ssf106 <- ssf106 %>% extract_covariates(popdensity); ssf106 <- rename(ssf106, popdensity =UTM_b115_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp) 
ssf106 <- ssf106 %>% extract_covariates(Landcover); ssf106 <- rename(ssf106, Landcover = UTM_b115_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp) 
ssf106 <- ssf106 %>% extract_covariates(treecover); ssf106 <- rename(ssf106, treecover = UTM_b115_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp) 
ssf106 <- ssf106 %>% extract_covariates(topo); ssf106 <- rename(ssf106, topo = UTM_b115_Neotropic_Earthenv_dem90m_tif_exp)  
ssf106 <- ssf106 %>% extract_covariates(water); ssf106 <- rename(ssf106, water = UTM_b115_water_frequency_2010_30m_tif_exp)  

print(ssf106, width=Inf)
#'


#' Jaguar 107
#' Let's illustrate track regularization with ID = 107. Let's
#' keeping is as 3 hours with tolerance of 1 hours
temp107<-trk %>% filter(id=="107") %>% track_resample(rate=hours(3), tolerance=hours(1))
temp107 %>% select(id, x_, y_, t_, burst_)

#Just checking it
summarize_sampling_rate(temp107)
#' Before fitting a step selection, the data well need to prepared. First, we change from a point representation to a step representation, using the function steps_by_burst, which in contrast to the steps function accounts for bursts.
ssf107 <- temp107 %>% steps_by_burst()

#' Next, we generate random steps with the function random_steps. This function fits by default a Gamma distribution to the step lengths and a von Mises distribution to the turn angles, and then pairs each observed step with n random steps.
ssf107 <- ssf107 %>% random_steps(n = 15)

#' As a last step, we have to extract the covariates at the end point of each step. We can do this with extract_covariates.
ssf107 <- ssf107 %>% extract_covariates(dist2drainage); ssf107 <- rename(ssf107, dist2drainage = UTM_b115_dist2drainage_exp)
ssf107 <- ssf107 %>% extract_covariates(dist2water); ssf107 <- rename(ssf107, dist2water = UTM_b115_dist2waterbodies_exp) 
ssf107 <- ssf107 %>% extract_covariates(dist2forest);ssf107 <- rename(ssf107, dist2forest =UTM_b115_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist)  
ssf107 <- ssf107 %>% extract_covariates(dist2roads); ssf107 <- rename(ssf107, dist2roads = UTM_b115_dist2roads_exp) 
ssf107 <- ssf107 %>% extract_covariates(human_footprint); ssf107 <- rename(ssf107, human_footprint = UTM_b115_human_footprint_2009_1km_tif_exp) 
ssf107 <- ssf107 %>% extract_covariates(Livestock); ssf107 <- rename(ssf107, Livestock = UTM_b115_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp) 
ssf107 <- ssf107 %>% extract_covariates(popdensity); ssf107 <- rename(ssf107, popdensity =UTM_b115_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp) 
ssf107 <- ssf107 %>% extract_covariates(Landcover); ssf107 <- rename(ssf107, Landcover = UTM_b115_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp) 
ssf107 <- ssf107 %>% extract_covariates(treecover); ssf107 <- rename(ssf107, treecover = UTM_b115_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp) 
ssf107 <- ssf107 %>% extract_covariates(topo); ssf107 <- rename(ssf107, topo = UTM_b115_Neotropic_Earthenv_dem90m_tif_exp)  
ssf107 <- ssf107 %>% extract_covariates(water); ssf107 <- rename(ssf107, water = UTM_b115_water_frequency_2010_30m_tif_exp)  

print(ssf107, width=Inf)
#'

#' Jaguar 108
#' Let's illustrate track regularization with ID = 108. Let's
#' keeping is as 3 hours with tolerance of 1 hours
temp108<-trk %>% filter(id=="108") %>% track_resample(rate=hours(3), tolerance=hours(1))
temp108 %>% select(id, x_, y_, t_, burst_)

#Just checking it
summarize_sampling_rate(temp108)
#' Before fitting a step selection, the data well need to prepared. First, we change from a point representation to a step representation, using the function steps_by_burst, which in contrast to the steps function accounts for bursts.
ssf108 <- temp108 %>% steps_by_burst()

#' Next, we generate random steps with the function random_steps. This function fits by default a Gamma distribution to the step lengths and a von Mises distribution to the turn angles, and then pairs each observed step with n random steps.
ssf108 <- ssf108 %>% random_steps(n = 15)

#' As a last step, we have to extract the covariates at the end point of each step. We can do this with extract_covariates.
ssf108 <- ssf108 %>% extract_covariates(dist2drainage); ssf108 <- rename(ssf108, dist2drainage = UTM_b115_dist2drainage_exp)
ssf108 <- ssf108 %>% extract_covariates(dist2water); ssf108 <- rename(ssf108, dist2water = UTM_b115_dist2waterbodies_exp) 
ssf108 <- ssf108 %>% extract_covariates(dist2forest);ssf108 <- rename(ssf108, dist2forest =UTM_b115_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist)  
ssf108 <- ssf108 %>% extract_covariates(dist2roads); ssf108 <- rename(ssf108, dist2roads = UTM_b115_dist2roads_exp) 
ssf108 <- ssf108 %>% extract_covariates(human_footprint); ssf108 <- rename(ssf108, human_footprint = UTM_b115_human_footprint_2009_1km_tif_exp) 
ssf108 <- ssf108 %>% extract_covariates(Livestock); ssf108 <- rename(ssf108, Livestock = UTM_b115_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp) 
ssf108 <- ssf108 %>% extract_covariates(popdensity); ssf108 <- rename(ssf108, popdensity =UTM_b115_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp) 
ssf108 <- ssf108 %>% extract_covariates(Landcover); ssf108 <- rename(ssf108, Landcover = UTM_b115_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp) 
ssf108 <- ssf108 %>% extract_covariates(treecover); ssf108 <- rename(ssf108, treecover = UTM_b115_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp) 
ssf108 <- ssf108 %>% extract_covariates(topo); ssf108 <- rename(ssf108, topo = UTM_b115_Neotropic_Earthenv_dem90m_tif_exp)  
ssf108 <- ssf108 %>% extract_covariates(water); ssf108 <- rename(ssf108, water = UTM_b115_water_frequency_2010_30m_tif_exp)  

print(ssf108, width=Inf)
#'


#' Jaguar 109
#' Let's illustrate track regularization with ID = 109. Let's
#' keeping is as 3 hours with tolerance of 1 hours
temp109<-trk %>% filter(id=="109") %>% track_resample(rate=hours(3), tolerance=hours(1))
temp109 %>% select(id, x_, y_, t_, burst_)

#Just checking it
summarize_sampling_rate(temp109)
#' Before fitting a step selection, the data well need to prepared. First, we change from a point representation to a step representation, using the function steps_by_burst, which in contrast to the steps function accounts for bursts.
ssf109 <- temp109 %>% steps_by_burst()

#' Next, we generate random steps with the function random_steps. This function fits by default a Gamma distribution to the step lengths and a von Mises distribution to the turn angles, and then pairs each observed step with n random steps.
ssf109 <- ssf109 %>% random_steps(n = 15)

#' As a last step, we have to extract the covariates at the end point of each step. We can do this with extract_covariates.
ssf109 <- ssf109 %>% extract_covariates(dist2drainage); ssf109 <- rename(ssf109, dist2drainage = UTM_b115_dist2drainage_exp)
ssf109 <- ssf109 %>% extract_covariates(dist2water); ssf109 <- rename(ssf109, dist2water = UTM_b115_dist2waterbodies_exp) 
ssf109 <- ssf109 %>% extract_covariates(dist2forest);ssf109 <- rename(ssf109, dist2forest =UTM_b115_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist)  
ssf109 <- ssf109 %>% extract_covariates(dist2roads); ssf109 <- rename(ssf109, dist2roads = UTM_b115_dist2roads_exp) 
ssf109 <- ssf109 %>% extract_covariates(human_footprint); ssf109 <- rename(ssf109, human_footprint = UTM_b115_human_footprint_2009_1km_tif_exp) 
ssf109 <- ssf109 %>% extract_covariates(Livestock); ssf109 <- rename(ssf109, Livestock = UTM_b115_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp) 
ssf109 <- ssf109 %>% extract_covariates(popdensity); ssf109 <- rename(ssf109, popdensity =UTM_b115_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp) 
ssf109 <- ssf109 %>% extract_covariates(Landcover); ssf109 <- rename(ssf109, Landcover = UTM_b115_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp) 
ssf109 <- ssf109 %>% extract_covariates(treecover); ssf109 <- rename(ssf109, treecover = UTM_b115_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp) 
ssf109 <- ssf109 %>% extract_covariates(topo); ssf109 <- rename(ssf109, topo = UTM_b115_Neotropic_Earthenv_dem90m_tif_exp)  
ssf109 <- ssf109 %>% extract_covariates(water); ssf109 <- rename(ssf109, water = UTM_b115_water_frequency_2010_30m_tif_exp)  

print(ssf109, width=Inf)
#'

#' Jaguar 110
#' Let's illustrate track regularization with ID = 110. Let's
#' keeping is as 3 hours with tolerance of 1 hours
temp110<-trk %>% filter(id=="110") %>% track_resample(rate=hours(3), tolerance=hours(1))
temp110 %>% select(id, x_, y_, t_, burst_)

#Just checking it
summarize_sampling_rate(temp110)
#' Before fitting a step selection, the data well need to prepared. First, we change from a point representation to a step representation, using the function steps_by_burst, which in contrast to the steps function accounts for bursts.
ssf110 <- temp110 %>% steps_by_burst()

#' Next, we generate random steps with the function random_steps. This function fits by default a Gamma distribution to the step lengths and a von Mises distribution to the turn angles, and then pairs each observed step with n random steps.
ssf110 <- ssf110 %>% random_steps(n = 15)

#' As a last step, we have to extract the covariates at the end point of each step. We can do this with extract_covariates.
ssf110 <- ssf110 %>% extract_covariates(dist2drainage); ssf110 <- rename(ssf110, dist2drainage = UTM_b115_dist2drainage_exp)
ssf110 <- ssf110 %>% extract_covariates(dist2water); ssf110 <- rename(ssf110, dist2water = UTM_b115_dist2waterbodies_exp) 
ssf110 <- ssf110 %>% extract_covariates(dist2forest);ssf110 <- rename(ssf110, dist2forest =UTM_b115_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist)  
ssf110 <- ssf110 %>% extract_covariates(dist2roads); ssf110 <- rename(ssf110, dist2roads = UTM_b115_dist2roads_exp) 
ssf110 <- ssf110 %>% extract_covariates(human_footprint); ssf110 <- rename(ssf110, human_footprint = UTM_b115_human_footprint_2009_1km_tif_exp) 
ssf110 <- ssf110 %>% extract_covariates(Livestock); ssf110 <- rename(ssf110, Livestock = UTM_b115_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp) 
ssf110 <- ssf110 %>% extract_covariates(popdensity); ssf110 <- rename(ssf110, popdensity =UTM_b115_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp) 
ssf110 <- ssf110 %>% extract_covariates(Landcover); ssf110 <- rename(ssf110, Landcover = UTM_b115_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp) 
ssf110 <- ssf110 %>% extract_covariates(treecover); ssf110 <- rename(ssf110, treecover = UTM_b115_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp) 
ssf110 <- ssf110 %>% extract_covariates(topo); ssf110 <- rename(ssf110, topo = UTM_b115_Neotropic_Earthenv_dem90m_tif_exp)  
ssf110 <- ssf110 %>% extract_covariates(water); ssf110 <- rename(ssf110, water = UTM_b115_water_frequency_2010_30m_tif_exp)  

print(ssf110, width=Inf)
#'


#' Jaguar 111
#' Let's illustrate track regularization with ID = 111. Let's
#' keeping it as 3 hours with tolerance of 1 hours
temp111<-trk %>% filter(id=="111") %>% track_resample(rate=hours(3), tolerance=hours(1))
temp111 %>% select(id, x_, y_, t_, burst_)

#Just checking it
summarize_sampling_rate(temp111)
#' Before fitting a step selection, the data well need to prepared. First, we change from a point representation to a step representation, using the function steps_by_burst, which in contrast to the steps function accounts for bursts.
ssf111 <- temp111 %>% steps_by_burst()

#' Next, we generate random steps with the function random_steps. This function fits by default a Gamma distribution to the step lengths and a von Mises distribution to the turn angles, and then pairs each observed step with n random steps.
ssf111 <- ssf111 %>% random_steps(n = 15)

#' As a last step, we have to extract the covariates at the end point of each step. We can do this with extract_covariates.
ssf111 <- ssf111 %>% extract_covariates(dist2drainage); ssf111 <- rename(ssf111, dist2drainage = UTM_b115_dist2drainage_exp)
ssf111 <- ssf111 %>% extract_covariates(dist2water); ssf111 <- rename(ssf111, dist2water = UTM_b115_dist2waterbodies_exp) 
ssf111 <- ssf111 %>% extract_covariates(dist2forest);ssf111 <- rename(ssf111, dist2forest =UTM_b115_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist)  
ssf111 <- ssf111 %>% extract_covariates(dist2roads); ssf111 <- rename(ssf111, dist2roads = UTM_b115_dist2roads_exp) 
ssf111 <- ssf111 %>% extract_covariates(human_footprint); ssf111 <- rename(ssf111, human_footprint = UTM_b115_human_footprint_2009_1km_tif_exp) 
ssf111 <- ssf111 %>% extract_covariates(Livestock); ssf111 <- rename(ssf111, Livestock = UTM_b115_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp) 
ssf111 <- ssf111 %>% extract_covariates(popdensity); ssf111 <- rename(ssf111, popdensity =UTM_b115_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp) 
ssf111 <- ssf111 %>% extract_covariates(Landcover); ssf111 <- rename(ssf111, Landcover = UTM_b115_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp) 
ssf111 <- ssf111 %>% extract_covariates(treecover); ssf111 <- rename(ssf111, treecover = UTM_b115_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp) 
ssf111 <- ssf111 %>% extract_covariates(topo); ssf111 <- rename(ssf111, topo = UTM_b115_Neotropic_Earthenv_dem90m_tif_exp)  
ssf111 <- ssf111 %>% extract_covariates(water); ssf111 <- rename(ssf111, water = UTM_b115_water_frequency_2010_30m_tif_exp)  

print(ssf111, width=Inf)
#'

#' Jaguar 112
#' Let's illustrate track regularization with ID = 112. Let's
#' keeping is as 3 hours with tolerance of 1 hours
temp112<-trk %>% filter(id=="112") %>% track_resample(rate=hours(3), tolerance=hours(1))
temp112 %>% select(id, x_, y_, t_, burst_)

#Just checking it
summarize_sampling_rate(temp112)
#' Before fitting a step selection, the data well need to prepared. First, we change from a point representation to a step representation, using the function steps_by_burst, which in contrast to the steps function accounts for bursts.
ssf112 <- temp112 %>% steps_by_burst()

#' Next, we generate random steps with the function random_steps. This function fits by default a Gamma distribution to the step lengths and a von Mises distribution to the turn angles, and then pairs each observed step with n random steps.
ssf112 <- ssf112 %>% random_steps(n = 15)

#' As a last step, we have to extract the covariates at the end point of each step. We can do this with extract_covariates.
ssf112 <- ssf112 %>% extract_covariates(dist2drainage); ssf112 <- rename(ssf112, dist2drainage = UTM_b115_dist2drainage_exp)
ssf112 <- ssf112 %>% extract_covariates(dist2water); ssf112 <- rename(ssf112, dist2water = UTM_b115_dist2waterbodies_exp) 
ssf112 <- ssf112 %>% extract_covariates(dist2forest);ssf112 <- rename(ssf112, dist2forest =UTM_b115_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist)  
ssf112 <- ssf112 %>% extract_covariates(dist2roads); ssf112 <- rename(ssf112, dist2roads = UTM_b115_dist2roads_exp) 
ssf112 <- ssf112 %>% extract_covariates(human_footprint); ssf112 <- rename(ssf112, human_footprint = UTM_b115_human_footprint_2009_1km_tif_exp) 
ssf112 <- ssf112 %>% extract_covariates(Livestock); ssf112 <- rename(ssf112, Livestock = UTM_b115_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp) 
ssf112 <- ssf112 %>% extract_covariates(popdensity); ssf112 <- rename(ssf112, popdensity =UTM_b115_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp) 
ssf112 <- ssf112 %>% extract_covariates(Landcover); ssf112 <- rename(ssf112, Landcover = UTM_b115_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp) 
ssf112 <- ssf112 %>% extract_covariates(treecover); ssf112 <- rename(ssf112, treecover = UTM_b115_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp) 
ssf112 <- ssf112 %>% extract_covariates(topo); ssf112 <- rename(ssf112, topo = UTM_b115_Neotropic_Earthenv_dem90m_tif_exp)  
ssf112 <- ssf112 %>% extract_covariates(water); ssf112 <- rename(ssf112, water = UTM_b115_water_frequency_2010_30m_tif_exp)  

print(ssf112, width=Inf)
#'


#' Jaguar 113
#' Let's illustrate track regularization with ID = 113. Let's
#' keeping is as 3 hours with tolerance of 1 hours
temp113<-trk %>% filter(id=="113") %>% track_resample(rate=hours(3), tolerance=hours(1))
temp113 %>% select(id, x_, y_, t_, burst_)

#Just checking it
summarize_sampling_rate(temp113)
#' Before fitting a step selection, the data well need to prepared. First, we change from a point representation to a step representation, using the function steps_by_burst, which in contrast to the steps function accounts for bursts.
ssf113 <- temp113 %>% steps_by_burst()

#' Next, we generate random steps with the function random_steps. This function fits by default a Gamma distribution to the step lengths and a von Mises distribution to the turn angles, and then pairs each observed step with n random steps.
ssf113 <- ssf113 %>% random_steps(n = 15)

#' As a last step, we have to extract the covariates at the end point of each step. We can do this with extract_covariates.
ssf113 <- ssf113 %>% extract_covariates(dist2drainage); ssf113 <- rename(ssf113, dist2drainage = UTM_b115_dist2drainage_exp)
ssf113 <- ssf113 %>% extract_covariates(dist2water); ssf113 <- rename(ssf113, dist2water = UTM_b115_dist2waterbodies_exp) 
ssf113 <- ssf113 %>% extract_covariates(dist2forest);ssf113 <- rename(ssf113, dist2forest =UTM_b115_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist)  
ssf113 <- ssf113 %>% extract_covariates(dist2roads); ssf113 <- rename(ssf113, dist2roads = UTM_b115_dist2roads_exp) 
ssf113 <- ssf113 %>% extract_covariates(human_footprint); ssf113 <- rename(ssf113, human_footprint = UTM_b115_human_footprint_2009_1km_tif_exp) 
ssf113 <- ssf113 %>% extract_covariates(Livestock); ssf113 <- rename(ssf113, Livestock = UTM_b115_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp) 
ssf113 <- ssf113 %>% extract_covariates(popdensity); ssf113 <- rename(ssf113, popdensity =UTM_b115_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp) 
ssf113 <- ssf113 %>% extract_covariates(Landcover); ssf113 <- rename(ssf113, Landcover = UTM_b115_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp) 
ssf113 <- ssf113 %>% extract_covariates(treecover); ssf113 <- rename(ssf113, treecover = UTM_b115_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp) 
ssf113 <- ssf113 %>% extract_covariates(topo); ssf113 <- rename(ssf113, topo = UTM_b115_Neotropic_Earthenv_dem90m_tif_exp)  
ssf113 <- ssf113 %>% extract_covariates(water); ssf113 <- rename(ssf113, water = UTM_b115_water_frequency_2010_30m_tif_exp)  

print(ssf113, width=Inf)
#'

#' Jaguar 114
#' Let's illustrate track regularization with ID = 114. Let's
#' keeping is as 3 hours with tolerance of 1 hours
temp114<-trk %>% filter(id=="114") %>% track_resample(rate=hours(3), tolerance=hours(1))
temp114 %>% select(id, x_, y_, t_, burst_)

#Just checking it
summarize_sampling_rate(temp114)
#' Before fitting a step selection, the data well need to prepared. First, we change from a point representation to a step representation, using the function steps_by_burst, which in contrast to the steps function accounts for bursts.
ssf114 <- temp114 %>% steps_by_burst()

#' Next, we generate random steps with the function random_steps. This function fits by default a Gamma distribution to the step lengths and a von Mises distribution to the turn angles, and then pairs each observed step with n random steps.
ssf114 <- ssf114 %>% random_steps(n = 15)

#' As a last step, we have to extract the covariates at the end point of each step. We can do this with extract_covariates.
ssf114 <- ssf114 %>% extract_covariates(dist2drainage); ssf114 <- rename(ssf114, dist2drainage = UTM_b115_dist2drainage_exp)
ssf114 <- ssf114 %>% extract_covariates(dist2water); ssf114 <- rename(ssf114, dist2water = UTM_b115_dist2waterbodies_exp) 
ssf114 <- ssf114 %>% extract_covariates(dist2forest);ssf114 <- rename(ssf114, dist2forest =UTM_b115_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist)  
ssf114 <- ssf114 %>% extract_covariates(dist2roads); ssf114 <- rename(ssf114, dist2roads = UTM_b115_dist2roads_exp) 
ssf114 <- ssf114 %>% extract_covariates(human_footprint); ssf114 <- rename(ssf114, human_footprint = UTM_b115_human_footprint_2009_1km_tif_exp) 
ssf114 <- ssf114 %>% extract_covariates(Livestock); ssf114 <- rename(ssf114, Livestock = UTM_b115_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp) 
ssf114 <- ssf114 %>% extract_covariates(popdensity); ssf114 <- rename(ssf114, popdensity =UTM_b115_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp) 
ssf114 <- ssf114 %>% extract_covariates(Landcover); ssf114 <- rename(ssf114, Landcover = UTM_b115_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp) 
ssf114 <- ssf114 %>% extract_covariates(treecover); ssf114 <- rename(ssf114, treecover = UTM_b115_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp) 
ssf114 <- ssf114 %>% extract_covariates(topo); ssf114 <- rename(ssf114, topo = UTM_b115_Neotropic_Earthenv_dem90m_tif_exp)  
ssf114 <- ssf114 %>% extract_covariates(water); ssf114 <- rename(ssf114, water = UTM_b115_water_frequency_2010_30m_tif_exp)  

print(ssf114, width=Inf)
#'



#' Jaguar 115
#' Let's illustrate track regularization with ID = 115. Let's
#' keeping it as 3 hours with tolerance of 1 hours
temp115<-trk %>% filter(id=="115") %>% track_resample(rate=hours(3), tolerance=hours(1))
temp115 %>% select(id, x_, y_, t_, burst_)

#Just checking it
summarize_sampling_rate(temp115)
#' Before fitting a step selection, the data well need to prepared. First, we change from a point representation to a step representation, using the function steps_by_burst, which in contrast to the steps function accounts for bursts.
ssf115 <- temp115 %>% steps_by_burst()

#' Next, we generate random steps with the function random_steps. This function fits by default a Gamma distribution to the step lengths and a von Mises distribution to the turn angles, and then pairs each observed step with n random steps.
ssf115 <- ssf115 %>% random_steps(n = 15)

#' As a last step, we have to extract the covariates at the end point of each step. We can do this with extract_covariates.
ssf115 <- ssf115 %>% extract_covariates(dist2drainage); ssf115 <- rename(ssf115, dist2drainage = UTM_b115_dist2drainage_exp)
ssf115 <- ssf115 %>% extract_covariates(dist2water); ssf115 <- rename(ssf115, dist2water = UTM_b115_dist2waterbodies_exp) 
ssf115 <- ssf115 %>% extract_covariates(dist2forest);ssf115 <- rename(ssf115, dist2forest =UTM_b115_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist)  
ssf115 <- ssf115 %>% extract_covariates(dist2roads); ssf115 <- rename(ssf115, dist2roads = UTM_b115_dist2roads_exp) 
ssf115 <- ssf115 %>% extract_covariates(human_footprint); ssf115 <- rename(ssf115, human_footprint = UTM_b115_human_footprint_2009_1km_tif_exp) 
ssf115 <- ssf115 %>% extract_covariates(Livestock); ssf115 <- rename(ssf115, Livestock = UTM_b115_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp) 
ssf115 <- ssf115 %>% extract_covariates(popdensity); ssf115 <- rename(ssf115, popdensity =UTM_b115_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp) 
ssf115 <- ssf115 %>% extract_covariates(Landcover); ssf115 <- rename(ssf115, Landcover = UTM_b115_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp) 
ssf115 <- ssf115 %>% extract_covariates(treecover); ssf115 <- rename(ssf115, treecover = UTM_b115_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp) 
ssf115 <- ssf115 %>% extract_covariates(topo); ssf115 <- rename(ssf115, topo = UTM_b115_Neotropic_Earthenv_dem90m_tif_exp)  
ssf115 <- ssf115 %>% extract_covariates(water); ssf115 <- rename(ssf115, water = UTM_b115_water_frequency_2010_30m_tif_exp)  

print(ssf115, width=Inf)

# Organize variables
#numeric
ssf105$case_<-as.numeric(ssf105$case_)
#' landcover classes 
ssf105$landcover <-as.factor(ssf105$Landcover)
ssf105<-ssf105 %>% mutate(landcover = fct_collapse(landcover,
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
print(ssf105, width=Inf)

#' Center and scale variables
ssf105<-ssf105 %>% mutate(tcoverSC=as.numeric(scale(treecover)),livestockSC=as.numeric(scale(Livestock)),
                          topoSC=as.numeric(scale(topo)),
                          humanfootSC=as.numeric(scale(human_footprint)), waterSC=as.numeric(scale(water)), 
                          dist2forestSC=as.numeric(scale(dist2forest)),dist2roadSC=as.numeric(scale(dist2roads)),
                          dist2waterSC=as.numeric(scale(dist2water)),distdrainSC=as.numeric(scale(dist2drainage)),
                          popSC=as.numeric(scale(popdensity))) 
#' ## Explore the data:

#' Look Distribution of variables for used and available
#+fig.width=8, fig.height=8
binomial_smooth <- function(...) {
  geom_smooth(method = "glm", method.args = list(family = "binomial"), ...)
}
ggplot(ssf105,aes(x=elevat, y=case_))+
stat_smooth(method="glm", method.args = list(family = "binomial"))+
  binomial_smooth(formula = y ~ splines::ns(x, 5), colour="red")+
  facet_wrap(~id, scales="free")



#' ALTERNATIVELY

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
                   tolerance=hours(max(3,round(timestats$median[i]/2))))
				   
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
ssfdat   


xyplot(id~date, data = J105, groups = id,xlab = "Sampling period",ylab = "ID Number")
ssf105<-ssfdat %>% filter(id=="105")
xyplot(id~t1_, data = as.data.frame(ssf105), groups = id,xlab = "Sampling period",ylab = "ID Number")

xyplot(id~date, data = J106, groups = id,xlab = "Sampling period",ylab = "ID Number")
ssf106<-ssfdat %>% filter(id=="106")
xyplot(id~t1_, data = as.data.frame(ssf106), groups = id,xlab = "Sampling period",ylab = "ID Number")

xyplot(id~date, data = J107, groups = id,xlab = "Sampling period",ylab = "ID Number")
ssf107<-ssfdat %>% filter(id=="107")
xyplot(id~t1_, data = as.data.frame(ssf107), groups = id,xlab = "Sampling period",ylab = "ID Number")

xyplot(id~date, data = J108, groups = id,xlab = "Sampling period",ylab = "ID Number")
ssf108<-ssfdat %>% filter(id=="108")
xyplot(id~t1_, data = as.data.frame(ssf108), groups = id,xlab = "Sampling period",ylab = "ID Number")

xyplot(id~date, data = J109, groups = id,xlab = "Sampling period",ylab = "ID Number")
ssf109<-ssfdat %>% filter(id=="109")
xyplot(id~t1_, data = as.data.frame(ssf109), groups = id,xlab = "Sampling period",ylab = "ID Number")

xyplot(id~date, data = J110, groups = id,xlab = "Sampling period",ylab = "ID Number")
ssf110<-ssfdat %>% filter(id=="110")
xyplot(id~t1_, data = as.data.frame(ssf110), groups = id,xlab = "Sampling period",ylab = "ID Number")

xyplot(id~date, data = J111, groups = id,xlab = "Sampling period",ylab = "ID Number")
ssf111<-ssfdat %>% filter(id=="111")
xyplot(id~t1_, data = as.data.frame(ssf111), groups = id,xlab = "Sampling period",ylab = "ID Number")

xyplot(id~date, data = J112, groups = id,xlab = "Sampling period",ylab = "ID Number")
ssf112<-ssfdat %>% filter(id=="112")
xyplot(id~t1_, data = as.data.frame(ssf112), groups = id,xlab = "Sampling period",ylab = "ID Number")

xyplot(id~date, data = J113, groups = id,xlab = "Sampling period",ylab = "ID Number")
ssf113<-ssfdat %>% filter(id=="113")
xyplot(id~t1_, data = as.data.frame(ssf113), groups = id,xlab = "Sampling period",ylab = "ID Number")

xyplot(id~date, data = J114, groups = id,xlab = "Sampling period",ylab = "ID Number")
ssf114<-ssfdat %>% filter(id=="114")
xyplot(id~t1_, data = as.data.frame(ssf114), groups = id,xlab = "Sampling period",ylab = "ID Number")

xyplot(id~date, data = J115, groups = id,xlab = "Sampling period",ylab = "ID Number")
ssf115<-ssfdat %>% filter(id=="115")
xyplot(id~t1_, data = as.data.frame(ssf115), groups = id,xlab = "Sampling period",ylab = "ID Number")


#' As a last step, we have to extract the covariates at the end point of each step. We can do this with extract_covariates.
ssfdat <- ssfdat %>% extract_covariates(dist2drainage); ssfdat <- rename(ssfdat, dist2drainage = UTM_b115_dist2drainage_exp)
ssfdat <- ssfdat %>% extract_covariates(dist2water); ssfdat <- rename(ssfdat, dist2water = UTM_b115_dist2waterbodies_exp) 
ssfdat <- ssfdat %>% extract_covariates(dist2forest);ssfdat <- rename(ssfdat, dist2forest =UTM_b115_dist2forestedges_Neotropic_Hansen_forest1_0_95percenttreecover_2000_30_tif_exp_edge_dist)  
ssfdat <- ssfdat %>% extract_covariates(dist2roads); ssfdat <- rename(ssfdat, dist2roads = UTM_b115_dist2roads_exp) 
ssfdat <- ssfdat %>% extract_covariates(human_footprint); ssfdat <- rename(ssfdat, human_footprint = UTM_b115_human_footprint_2009_1km_tif_exp) 
ssfdat <- ssfdat %>% extract_covariates(Livestock); ssfdat <- rename(ssfdat, Livestock = UTM_b115_Livestock_Cattle_CC2006_AD_1km_neotropic_albers_tif_exp) 
ssfdat <- ssfdat %>% extract_covariates(popdensity); ssfdat <- rename(ssfdat, popdensity =UTM_b115_Population_density_gpw_v4_rev10_2015_1km_neotropic_albers_tif_exp) 
ssfdat <- ssfdat %>% extract_covariates(Landcover); ssfdat <- rename(ssfdat, Landcover = UTM_b115_Landcover_ESACCI_2015_300m_neotropic_albers_tif_exp) 
ssfdat <- ssfdat %>% extract_covariates(treecover); ssfdat <- rename(ssfdat, treecover = UTM_b115_Neotropic_Hansen_percenttreecover_2000_30m_tif_exp) 
ssfdat <- ssfdat %>% extract_covariates(topo); ssfdat <- rename(ssfdat, topo = UTM_b115_Neotropic_Earthenv_dem90m_tif_exp)  
ssfdat <- ssfdat %>% extract_covariates(water); ssfdat <- rename(ssfdat, water = UTM_b115_water_frequency_2010_30m_tif_exp)  

print(ssfdat, width=Inf)








###Just to check the medians of what we did in comparison with the totals
ssfdatemp<-ssfdat %>% filter(case_==TRUE) 
ssfdatdf <- data.frame(ssfdatemp)

# A summary applied to ungrouped tbl returns a single row group by id,burst_
bursts <-ssfdatdf %>%
group_by(id,burst_) %>%
summarise(n_distinct = n_distinct(burst_), n = n())%>% data.frame
head(bursts)  

  	
#SSFs sampled data (x1, y1, t1)
#ssfdattrk <- mk_track(ssfdatemp, .x=x1_, .y=y1_, .t=t1_, id=id, crs = CRS("+proj=utm +zone=21K +south +datum=WGS84 +units=m +no_defs")) 
##Summary of totals from the original trk(as we did above just to compare)
#(timestatsssf<-ssfdattrk %>% nest(-id) %>% mutate(sr = map(data, summarize_sampling_rate, time_unit = "hour")) %>% select(id, sr) %>% unnest)
### SSFs sampled data (x2, y2, t2)
ssfdattrk <- mk_track(ssfdatemp, .x=x2_, .y=y2_, .t=t2_, id=id, crs = CRS("+proj=utm +zone=21K +south +datum=WGS84 +units=m +no_defs")) 
##Summary of totals from the original trk(as we did above just to compare)
(timestatsssf<-ssfdattrk %>% nest(-id) %>% mutate(sr = map(data, summarize_sampling_rate, time_unit = "hour")) %>% select(id, sr) %>% unnest)
## And comparing that to the original file
#SaoBentotrk -> trk	
#(timestats<-trk %>% nest(-id,-sex,-age,-weight,-status) %>% mutate(sr = map(data, summarize_sampling_rate, time_unit = "hour")) %>%select(id,sex,age,weight,status,sr) %>% unnest)

#' Now, lets plot the data for random and matched points
#' 
#+fig.height=12, fig.width=12, warning=FALSE
ggplot(ssfdat, aes(x2_, y2_, color=case_))+geom_point()+facet_wrap(~id,scales="free")

#' ## Write out data for further annotating

#'  With the SSFs, we have the extra complication of
#' having a time and location at both the start and end of the step.  
#' 
#' For the time being, we will assume we want to annotate variables at the end of the step
#' but use the starting point of the step as the timestamp.


#ssfALL<-ssfdat %>% select(t1_,x2_)


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
stepSB<-as_tibble(stepSB)
#print(ssfdat, width=Inf)



#AllSteps <- ssf.df  ### include name of the area/ project
#' These points then need to be annotated prior to fitting ssfs. 2 files All and to annotate
#write.csv(ssfALL, file="C:/RWorkDir/jaguardatapaper/ssfALL.csv", row.names = FALSE) ###annotated_projectname
#write.csv(stepSB, file="C:/RWorkDir/jaguardatapaper/stepSB.csv", row.names = FALSE) ###annotated_projectname

#'







#" Simplify some variable names and make case a numeric variable
ssfdat$case_<-as.numeric(ssfdat$case)
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
ssfdat<-ssfdat %>% mutate(f_coverSC=as.numeric(scale(tree_cover)),livestockSC=as.numeric(scale(livestock)),
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




