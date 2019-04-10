#' #  **Jaguar Map Move**
#' 
#' #### *Alan E. de Barros, Bernardo Niebuhr, Vanesa Bejarano, Julia Oshima,Claudia Kanda, Milton Ribeiro, Ronaldo Morato,Paulo Prado*
#' date: "March, 13 2019"
#' #### Run JaguarDataPrep first !!! 
#'
#' ### Visualize all data with leaflet
pal <- colorFactor(palette = 'Dark2',domain = jaguar_df$id)
leaflet(jaguar_df)%>%addTiles()%>%addCircles(jaguar_df$x,jaguar_df$y,color=~pal(id))%>% addProviderTiles(providers$Esri.WorldImagery)
xyplot(id~date, data = jaguar_df, groups = project_region, auto.key=list(columns = 3))
options(max.print=1000000)
#^ Summary
summary(jaguar_df)
#' ##### All individuals basic counts
colSums(ifelse(table(jaguar_df$id,jaguar_df$sex)>0,1,0)) #table(droplevels(jaguar_df$id),jaguar_df$sex)  # sex
colSums(ifelse(table(jaguar_df$id,jaguar_df$age)>0,1,0)) #table(droplevels(jaguar_df$id),jaguar_df$age)  # age
colSums(ifelse(table(jaguar_df$id,jaguar_df$weight)>0,1,0)) #table(droplevels(jaguar_df$id),jaguar_df$weight) # weight
colSums(ifelse(table(jaguar_df$id,jaguar_df$status)>0,1,0)) #table(droplevels(jaguar_df$id),jaguar_df$status) # residence status
colSums(ifelse(table(jaguar_df$id,jaguar_df$country)>0,1,0)) #table(droplevels(jaguar_df$id),jaguar_df$country) # country
colSums(ifelse(table(jaguar_df$id,jaguar_df$project_region)>0,1,0)) #table(droplevels(jaguar_df$id),jaguar_df$project_region) 
colSums(ifelse(table(jaguar_df$id,jaguar_df$project_bioveg)>0,1,0)) #table(droplevels(jaguar_df$id),jaguar_df$project_bioveg) 

#'
#' #### 1) Atlantic Forest W1
pal <- colorFactor(palette = 'Dark2',domain = AFW1$id)
#pal <- colorFactor(palette = 'Paired',domain = AFW1$id)
leaflet(AFW1)%>%addTiles()%>%addCircles(AFW1$x, AFW1$y, color =~pal(id))%>% addProviderTiles(providers$Esri.WorldImagery)

ft=ftable(AFW1$id); f<-as.data.frame.table(ft);f$Var1 <- NULL;f$Var2 <- NULL;subset(f,Freq>0)
xyplot(id~date, data = AFW1, groups = id)

#' ##### 3D plot
open3d()
# To get a bigger window than the default
par3d(windowRect = c(100, 100, 612, 612))
Sys.sleep(0.1) # Allow sluggish window managers to catch up
with(AFW1, plot3d(x,y,date, type="l", col=as.factor(AFW1$idf)))
(stcube<-with(AFW1, plot3d(x,y,date, type="l",col=as.numeric(cut(AFW1$id,7)), alpha=0.4)))
#' ##### 3D plot
open3d()
# To get a bigger window than the default
par3d(windowRect = c(100, 100, 612, 612))
Sys.sleep(0.1) # Allow sluggish window managers to catch up
with(AFW1, plot3d(x,y,date, type="p", col=as.numeric(AFW1$weight)))
(stcube<-with(AFW1, plot3d(x,y,date, type="p",col=as.numeric(cut(AFW1$weight,7)), alpha=0.4)))


#' ##### Individuals basic counts
summary(AFW1)
table(droplevels(AFW1$id),AFW1$sex); colSums(ifelse(table(AFW1$id,AFW1$sex)>0,1,0))  # sex
table(droplevels(AFW1$id),AFW1$age); colSums(ifelse(table(AFW1$id,AFW1$age)>0,1,0)) # age
table(droplevels(AFW1$id),AFW1$weight); colSums(ifelse(table(AFW1$id,AFW1$weight)>0,1,0))  # weight
table(droplevels(AFW1$id),AFW1$status); colSums(ifelse(table(AFW1$id,AFW1$status)>0,1,0)) # residence status

#' ##### Working with amt
AFW1trk -> trk; get_crs(AFW1trk) # double checking CRS
AFW1trk %>% dplyr::select(id,sl) %>% unnest %>% 
  ggplot(aes(sl, fill = factor(id))) + geom_density(alpha = 0.4)

AFW1trk %>% dplyr::select(id,dt) %>% unnest %>% 
  ggplot(aes(dt, fill = factor(id))) + geom_density(alpha = 0.4)


trk <- trk %>%dplyr::select(x_,y_,t_,id,tod_, everything())
trk <- trk%>%dplyr::select(-sl,-project_region,-dir_abs,-dir_rel,-nsd_,-long_x,-lat_y, 
                   -week, -month, -year,-hour, everything())
trk # %>% print(n = Inf)

trk <- trk%>%select(-sl,-project_region,-dir_abs,-dir_rel,-nsd_,-long_x,-lat_y, 
                    -week, -month, -year,-hour, everything())


#' Check the current sampling rate
(timestats<-trk %>% nest(-id) %>% mutate(sr = map(data, summarize_sampling_rate, time_unit = "hour")) %>%
    select(id, sr) %>% unnest)
#  Lets add on the time difference to each obs.
trk<-trk %>% group_by(id) %>% mutate(dt_ = t_ - lag(t_, default = NA))
trk <- trk %>%select(x_,y_,t_,id,tod_, everything())
trk



#' ##### J34
leaflet(J34)%>%addTiles()%>%addCircles(J34$x, J34$y,color = "red")%>% addProviderTiles(providers$Esri.WorldImagery)

leaflet(J34)%>%addTiles()%>% addPolylines(J34$x, J34$y, smoothFactor = 0.01,stroke=0.01,weight = 2,color = "red")%>% addProviderTiles(providers$Esri.WorldImagery)


xyplot(id~date, data = J34, groups = id,xlab = "Sampling period",ylab = "ID Number")
J34trk %>% dplyr::select(id,dt) %>% unnest %>% 
  ggplot(aes(dt, fill = factor(id))) + geom_density(alpha = 0.4)

head(J34)
nrow(J34);range(J34$date)
max(J34$dt, na.rm = TRUE)/60/60/24;mean(J34$dt, na.rm = TRUE)/60/60/24;median(J34$dt, na.rm = TRUE)/60/60/24 # days
min(J34$dt, na.rm = TRUE)/60/60 #hours

J34trk 

#' ## Display space-time cube. 
open3d()
# To get a bigger window than the default
par3d(windowRect = c(100, 100, 612, 612))
Sys.sleep(0.1) # Allow sluggish window managers to catch up
with(J34, plot3d(x,y,date, type="l", col=as.numeric(J34$date)))
(stcube<-with(J34, plot3d(x,y,date, type="l",col=as.numeric(cut(J34$date,20)), alpha=0.4)))

#Histograms and plots
hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 34],
     main = 'Individual 34', xlab = "Time between fixes (h)")
plot(density.default(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 34],
                     na.rm = T), main = 'Individual 34')

breaks <- c(0, 0.5, 1.5, 2.5, 4.5, 8.5, 16.5, 32.5, 64.5, 128.5, 256.5, 1000)
breaks <- c(0, 1, 2*(1:256))
hh <- hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 34 &
                                     mov.data.diff$time.diff > 0 & 
                                     mov.data.diff$time.diff < 512],
           main = 'Individual 34',xlab = 'Time between fixes (h)',
           plot = T, breaks = breaks)







leaflet(J35 )%>%addTiles()%>%addCircles(J35$x, J35$y, color = "red" )%>% addProviderTiles(providers$Esri.WorldImagery)
leaflet(J36 )%>%addTiles()%>%addCircles(J36$x, J36$y,color = "red" )%>% addProviderTiles(providers$Esri.WorldImagery)
leaflet(J37 )%>%addTiles()%>%addCircles(J37 $x, J37 $y,color = "red")%>% addProviderTiles(providers$Esri.WorldImagery)
leaflet(J38 )%>%addTiles()%>%addCircles(J38 $x, J38 $y,color = "red")%>% addProviderTiles(providers$Esri.WorldImagery)
leaflet(J58 )%>%addTiles()%>%addCircles(J58 $x, J58 $y,color = "blue")%>% addProviderTiles(providers$Esri.WorldImagery)
leaflet(J62 )%>%addTiles()%>%addCircles(J62 $x, J62 $y,color = "blue")%>% addProviderTiles(providers$Esri.WorldImagery)

#'
#' #### 2) Atlantic Forest W2
ft=ftable(AFW2$id); f<-as.data.frame.table(ft);f$Var1 <- NULL;f$Var2 <- NULL;subset(f,Freq>0)
leaflet(AFW2)%>%addTiles()%>%addCircles(AFW2$x, AFW2$y, color =~pal(id))%>% addProviderTiles(providers$Esri.WorldImagery)

AFW2trk %>% dplyr::select(id,sl) %>% unnest %>% 
  ggplot(aes(sl, fill = factor(id))) + geom_density(alpha = 0.4)

AFW2trk %>% dplyr::select(id,dt) %>% unnest %>% 
  ggplot(aes(dt, fill = factor(id))) + geom_density(alpha = 0.4)


#' ##### AFW2 individuals

#' Sex  # table(AFW2$id,AFW2$sex);colSums(ifelse(table(AFW2$id,AFW2$sex)>0,1,0))
Male=subset(AFW2,sex=='Male')  # Males
f<-as.data.frame.table(ftable(Male$id));f$Var1 <- NULL;f$Var2 <- NULL;subset(f,Freq>0) 
# colSums(ifelse(table(Male$id,Male$age)>0,1,0));table(Male$id,Male$age)
Female=subset(AFW2,sex=='Female')  # Females
f<-as.data.frame.table(ftable(Female$id));f$Var1 <- NULL;f$Var2 <- NULL;subset(f,Freq>0) 
# colSums(ifelse(table(Female$id,Female$age)>0,1,0));table(Female$id,Female$age)


new_df <- AFW2trk %>%select(period,tod_, everything())
new_df%>% print(n = Inf)


leaflet(J39 )%>%addTiles()%>%addCircles(J39 $x, J39 $y)%>% addProviderTiles(providers$Esri.WorldImagery)
leaflet(J40 )%>%addTiles()%>%addCircles(J40 $x, J40 $y)%>% addProviderTiles(providers$Esri.WorldImagery)
leaflet(J63 )%>%addTiles()%>%addCircles(J63 $x, J63 $y)%>% addProviderTiles(providers$Esri.WorldImagery)

#'
#' #### 3) Caatinga
leaflet(Caatinga)%>%addTiles()%>%addCircles(Caatinga$x, Caatinga$y)%>% addProviderTiles(providers$Esri.WorldImagery)


Caatingatrk %>% dplyr::select(id,sl) %>% unnest %>% 
  ggplot(aes(sl, fill = factor(id))) + geom_density(alpha = 0.4)

Caatingatrk %>% dplyr::select(id,dt) %>% unnest %>% 
  ggplot(aes(dt, fill = factor(id))) + geom_density(alpha = 0.4)

#'Caatinga individuals
ft=ftable(Caatinga$id); f<-as.data.frame.table(ft);f$Var1 <- NULL;f$Var2 <- NULL;subset(f,Freq>0)



leaflet(J20 )%>%addTiles()%>%addCircles(J20 $x, J20 $y)%>% addProviderTiles(providers$Esri.WorldImagery)
leaflet(J50 )%>%addTiles()%>%addCircles(J50 $x, J50 $y)%>% addProviderTiles(providers$Esri.WorldImagery)

#'
#' #### 4) Cerrado1
leaflet(Cerrado1)%>%addTiles()%>%addCircles(Cerrado1$x, Cerrado1$y)%>% addProviderTiles(providers$Esri.WorldImagery)

Cerrado1trk %>% dplyr::select(id,sl) %>% unnest %>% 
  ggplot(aes(sl, fill = factor(id))) + geom_density(alpha = 0.4)

Cerrado1trk %>% dplyr::select(id,dt) %>% unnest %>% 
  ggplot(aes(dt, fill = factor(id))) + geom_density(alpha = 0.4)


#'Cerrado1 individuals
ft=ftable(Cerrado1$id); f<-as.data.frame.table(ft);f$Var1 <- NULL;f$Var2 <- NULL;subset(f,Freq>0)
leaflet(J17 )%>%addTiles()%>%addCircles(J17 $x, J17 $y)%>% addProviderTiles(providers$Esri.WorldImagery)
leaflet(J65 )%>%addTiles()%>%addCircles(J65 $x, J65 $y)%>% addProviderTiles(providers$Esri.WorldImagery)
leaflet(J67 )%>%addTiles()%>%addCircles(J67 $x, J67 $y)%>% addProviderTiles(providers$Esri.WorldImagery)
leaflet(J82 )%>%addTiles()%>%addCircles(J82 $x, J82 $y)%>% addProviderTiles(providers$Esri.WorldImagery)
leaflet(J85 )%>%addTiles()%>%addCircles(J85 $x, J85 $y)%>% addProviderTiles(providers$Esri.WorldImagery)

#'
#' #### 5) Cerrado2
leaflet(Cerrado2)%>%addTiles()%>%addCircles(Cerrado2$x, Cerrado2$y)%>% addProviderTiles(providers$Esri.WorldImagery)

Cerrado2trk %>% dplyr::select(id,sl) %>% unnest %>% 
  ggplot(aes(sl, fill = factor(id))) + geom_density(alpha = 0.4)

Cerrado2trk %>% dplyr::select(id,dt) %>% unnest %>% 
  ggplot(aes(dt, fill = factor(id))) + geom_density(alpha = 0.4)


#'Cerrado2 individuals
ft=ftable(Cerrado2$id); f<-as.data.frame.table(ft);f$Var1 <- NULL;f$Var2 <- NULL;subset(f,Freq>0)
leaflet(J89 )%>%addTiles()%>%addCircles(J89 $x, J89 $y)%>% addProviderTiles(providers$Esri.WorldImagery)

#'
#' #### 6) Costa Rica
leaflet(CRica)%>%addTiles()%>%addCircles(CRica$x, CRica$y)%>% addProviderTiles(providers$Esri.WorldImagery)

CRicatrk %>% dplyr::select(id,sl) %>% unnest %>% 
  ggplot(aes(sl, fill = factor(id))) + geom_density(alpha = 0.4)

CRicatrk %>% dplyr::select(id,dt) %>% unnest %>% 
  ggplot(aes(dt, fill = factor(id))) + geom_density(alpha = 0.4)


#'Costa Rica individuals
ft=ftable(CRica$id); f<-as.data.frame.table(ft);f$Var1 <- NULL;f$Var2 <- NULL;subset(f,Freq>0)
leaflet(J26 )%>%addTiles()%>%addCircles(J26 $x, J26 $y)%>% addProviderTiles(providers$Esri.WorldImagery)

#'
#' #### 7) Drych1 
leaflet(Drych1)%>%addTiles()%>%addCircles(Drych1$x, Drych1$y)%>% addProviderTiles(providers$Esri.WorldImagery)

Drych1trk %>% dplyr::select(id,sl) %>% unnest %>% 
  ggplot(aes(sl, fill = factor(id))) + geom_density(alpha = 0.4)

Drych1trk %>% dplyr::select(id,dt) %>% unnest %>% 
  ggplot(aes(dt, fill = factor(id))) + geom_density(alpha = 0.4)


#' Drych1 individuals
ft=ftable(Drych1$id); f<-as.data.frame.table(ft);f$Var1 <- NULL;f$Var2 <- NULL;subset(f,Freq>0)
leaflet(J16 )%>%addTiles()%>%addCircles(J16 $x, J16 $y)%>% addProviderTiles(providers$Esri.WorldImagery)
leaflet(J70 )%>%addTiles()%>%addCircles(J70 $x, J70 $y)%>% addProviderTiles(providers$Esri.WorldImagery)
leaflet(J76 )%>%addTiles()%>%addCircles(J76 $x, J76 $y)%>% addProviderTiles(providers$Esri.WorldImagery)
leaflet(J77 )%>%addTiles()%>%addCircles(J77 $x, J77 $y)%>% addProviderTiles(providers$Esri.WorldImagery)

#' 
#' #### 8) Drych2
leaflet(Drych2)%>%addTiles()%>%addCircles(Drych2$x, Drych2$y)%>% addProviderTiles(providers$Esri.WorldImagery)

Drych2trk %>% dplyr::select(id,sl) %>% unnest %>% 
  ggplot(aes(sl, fill = factor(id))) + geom_density(alpha = 0.4)

Drych2trk %>% dplyr::select(id,dt) %>% unnest %>% 
  ggplot(aes(dt, fill = factor(id))) + geom_density(alpha = 0.4)

#'
#' #### 9) Humid Chaco 
leaflet(Hch)%>%addTiles()%>%addCircles(Hch$x, Hch$y)%>% addProviderTiles(providers$Esri.WorldImagery)

Hchtrk %>% dplyr::select(id,sl) %>% unnest %>% 
  ggplot(aes(sl, fill = factor(id))) + geom_density(alpha = 0.4)

Hchtrk %>% dplyr::select(id,dt) %>% unnest %>% 
  ggplot(aes(dt, fill = factor(id))) + geom_density(alpha = 0.4)

#' 
J1

#' #### 10) Forest Paraguay
leaflet(FPy)%>%addTiles()%>%addCircles(FPy$x, FPy$y)%>% addProviderTiles(providers$Esri.WorldImagery)
#' 
#' #### 11) Iguazu1
leaflet(Iguazu1)%>%addTiles()%>%addCircles(Iguazu1$x, Iguazu1$y)%>% addProviderTiles(providers$Esri.WorldImagery)
#' 
#'#### 12) Iguazu2
leaflet(Iguazu2)%>%addTiles()%>%addCircles(Iguazu2$x, Iguazu2$y)%>% addProviderTiles(providers$Esri.WorldImagery)
#' 
#' #### 13) Amazonia Mamiraua (Brazil) - Flooded Amazonia
leaflet(Mamiraua)%>%addTiles()%>%addCircles(Mamiraua$x, Mamiraua$y)%>% addProviderTiles(providers$Esri.WorldImagery)
#' 
#' #### 14) IOP Para Amazonia   
# Translocated animal   
leaflet(iopPA)%>%addTiles()%>%addCircles(iopPA$x, iopPA$y)%>% addProviderTiles(providers$Esri.WorldImagery)
#' 
#' #### 15) Greater Lacandona, Mexico 
leaflet(Lacandona)%>%addTiles()%>%addCircles(Lacandona$x, Lacandona$y)%>% addProviderTiles(providers$Esri.WorldImagery)
#' 
#' #### 16) Mexico East
leaflet(MexEast)%>%addTiles()%>%addCircles(MexEast$x, MexEast$y)%>% addProviderTiles(providers$Esri.WorldImagery)
#' 
#' #### 17) Mexico Sonora
leaflet(Sonora)%>%addTiles()%>%addCircles(Sonora$x, Sonora$y)%>% addProviderTiles(providers$Esri.WorldImagery)
#' 
#' #### 18) PantanalTotal Brazil&Paraguay   -  7 Projects in total
leaflet(Pantanal)%>%addTiles()%>%addCircles(Pantanal$x, Pantanal$y)%>% addProviderTiles(providers$Esri.WorldImagery)

ft=ftable(Pantanal$id); f<-as.data.frame.table(ft);f$Var1 <- NULL;f$Var2 <- NULL;subset(f,Freq>0)
x11()
xyplot(id~date, data = Pantanal, groups = id)
xyplot(id~date, data = Pantanal, groups = id, auto.key=list(columns = 3))



#' 
leaflet(Oncafari)%>%addTiles()%>%addCircles(Oncafari$x, Oncafari$y)%>% addProviderTiles(providers$Esri.WorldImagery)
leaflet(Paraguay)%>%addTiles()%>%addCircles(Paraguay$x, Paraguay$y)%>% addProviderTiles(providers$Esri.WorldImagery)
leaflet(Panthera1)%>%addTiles()%>%addCircles(Panthera1$x, Panthera1$y)%>% addProviderTiles(providers$Esri.WorldImagery)
leaflet(Panthera2)%>%addTiles()%>%addCircles(Panthera2$x, Panthera2$y)%>% addProviderTiles(providers$Esri.WorldImagery)
leaflet(RioNegro)%>%addTiles()%>%addCircles(RioNegro$x, RioNegro$y)%>% addProviderTiles(providers$Esri.WorldImagery)
leaflet(SaoBento)%>%addTiles()%>%addCircles(SaoBento$x, SaoBento$y)%>% addProviderTiles(providers$Esri.WorldImagery)
leaflet(Taiama)%>%addTiles()%>%addCircles(Taiama$x, Taiama$y)%>% addProviderTiles(providers$Esri.WorldImagery)




#' 
#' 




#' 
#' 


#' ### Organize data again accordingly with packages for visualization/visual analysis
#' 
#' ####  amt
#'  already set in data prep (trk files for region or individuals e.g. Pantanaltrk or )


#' 
#' #### Adjust again as move object
move.data <- move(x = jaguar_df$x, y = jaguar_df$y, 
                  time = jaguar_df$date, 
                  proj = CRS('+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'),
                  data = jaguar_df, animal = jaguar_df$id, sensor = 'GPS')

move.data  ### moveStack
###Separate individual animals' trajectories 
unstacked <- split(move.data)
head(unstacked)


#-----------------------------
# Distribution of fix rates
mov.data.diff <- jaguar_df

# Calculation of time between fixes
time.between <- function(individual, dat) {
  # Select individual
  ind <- dat[dat$id == individual,]
  # Calculate difference in time, in hours, and return
  c(as.numeric(diff.POSIXt(ind$timestamp.posix), units = 'hours'), NA)
}
# All individuals
inds <- unique(mov.data.diff$id)
mov.data.diff$time.diff <- unlist(sapply(inds, FUN = time.between, dat = mov.data.diff))

# Histogram
for(i in inds) {
  print(paste(i, paste(range(mov.data.diff$time.diff[mov.data.diff$id == i], na.rm = T), collapse = ', '), sep = ': '))
}

hist(mov.data.diff$time.diff[mov.data.diff$time.diff > 0 & mov.data.diff$time.diff < 48],
     main = 'All individuals pooled', xlab = "Time between fixes (h)")

hist(mov.data.diff$time.diff[mov.data.diff$time.diff > 0 & mov.data.diff$time.diff > 720],
     main = 'All individuals pooled', xlab = "Time between fixes (h)")

table(mov.data.diff$time.diff[mov.data.diff$time.diff  > 720])
z=(mov.data.diff$time.diff[mov.data.diff$time.diff  > 720])

# Checking the range
(summary_range=tapply(Pantanal$date,list(Pantanal$id),range))
# Checking the range in terms of days
(dr=with(Pantanal, tapply(date, list(id), function(date) max(date)-min(date))))


#-----------------------------


	## Plot individual animals' trajectories (Maximum Zoom to include all points)
#X1	
unstacked <- split(move.data)
X1=unstacked$X1
X1_df <- as(X1, "data.frame")
head(X1)
nrow(X1_df)
range(X1_df$date)
#' Now, using leaflet
leaflet(X1)%>%addTiles()%>%addCircles(X1$x, X1$y)

leaflet(X1)%>%addTiles()%>%addCircles(X1$x, X1$y)%>% addProviderTiles(providers$Esri.NatGeoWorldMap)

leaflet(X1)%>%addTiles()%>%addCircles(X1$x, X1$y)%>% addProviderTiles(providers$Esri.WorldImagery)



m <- get_map(bbox(extent(X1)*1.1), source="google", zoom=12, maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=X1_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
#aes(colour = as.numeric(date)
#geom_point
#' ## Display space-time cube. 
open3d()
# To get a bigger window than the default
par3d(windowRect = c(100, 100, 612, 612))
Sys.sleep(0.1) # Allow sluggish window managers to catch up
with(X1_df, plot3d(x,y,date, type="l", col=as.numeric(X1_df$Event_ID)))
(stcube<-with(X1_df, plot3d(x,y,date, type="l",col=as.numeric(cut(X1_df$dist,3)), alpha=0.4)))
#Timelag
head(timeLag(X1, units="mins"))
head(timeLag(X1, units="hours"))
median(timeLag(X1, units="hours"))
max(timeLag(X1, units="hours"))
max(timeLag(X1, units="days"))
plot(X1_df$individual.local.identifier..ID.~X1_df $date,xlab = "Sampling period",ylab = "Individuals ID Number")
#write.table(X1,file="c:/RWorkDir/jaguardatapaper/ID/X1.txt",row.names = F,quote=F,col.names=T,sep="\t")
#Histograms and plots
hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 1],
     main = 'Individual 1', xlab = "Time between fixes (h)")
plot(density.default(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 1],
     na.rm = T), main = 'Individual 1')
	 
breaks <- c(0, 0.5, 1.5, 2.5, 4.5, 8.5, 16.5, 32.5, 64.5, 128.5, 256.5, 1000)
breaks <- c(0, 1, 2*(1:256))
hh <- hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 1 &
                                     mov.data.diff$time.diff > 0 & 
                                     mov.data.diff$time.diff < 512],
									 main = 'Individual 1',xlab = 'Time between fixes (h)',
           plot = T, breaks = breaks)

#X2 
X2=unstacked$X2
X2_df <- as(X2, "data.frame")
head(X2_df)
nrow(X2_df)
range(X2_df$date)
m <- get_map(bbox(extent(X2)*1.1), source="google", zoom=13, maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=X2_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
#m <- get_map(bbox(extent(X2)*1.1), source="google", zoom=12, maptype="satellite") #define the box for map based on the trajectory coordinates
#ggmap(m)+geom_path(data=X2_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
#' ## Display space-time cube. 
open3d()
# To get a bigger window than the default
par3d(windowRect = c(100, 100, 612, 612))
Sys.sleep(0.1) # Allow sluggish window managers to catch up
with(X2_df, plot3d(x,y,date, type="l", col=as.numeric(X2_df$Event_ID)))
(stcube<-with(X2_df, plot3d(x,y,date, type="l",col=as.numeric(cut(X2_df$dist,3)), alpha=0.4)))
#Timelag
head(timeLag(X2, units="mins"))
head(timeLag(X2, units="hours"))
median(timeLag(X2, units="hours"))
max(timeLag(X2, units="hours"))
max(timeLag(X2, units="days"))
plot(X2_df$individual.local.identifier..ID.~X2_df$date,xlab = "Sampling period",ylab = "Individuals ID Number")
#write.table(X2,file="c:/RWorkDir/jaguardatapaper/ID/X2.txt",row.names = F,quote=F,col.names=T,sep="\t")
#Histograms and plots
hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 2],
     main = 'Individual 2', xlab = "Time between fixes (h)")
plot(density.default(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 2],
     na.rm = T), main = 'Individual 2')
	 
breaks <- c(0, 0.5, 1.5, 2.5, 4.5, 8.5, 16.5, 32.5, 64.5, 128.5, 256.5, 1000)
breaks <- c(0, 1, 2*(1:256))
hh <- hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 2 &
                                     mov.data.diff$time.diff > 0 & 
                                     mov.data.diff$time.diff < 512],
									 main = 'Individual 2',xlab = 'Time between fixes (h)',
           plot = T, breaks = breaks)

#X3
X3=unstacked$X3
X3_df <- as(X3, "data.frame")
head(X3_df)
nrow(X3_df)
range(X3_df$date)
m <- get_map(bbox(extent(X3)*1.1), source="google", zoom=10, maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=X3_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
#m <- get_map(bbox(extent(X3)*1.1), source="google", zoom=11, maptype="satellite") #define the box for map based on the trajectory coordinates
#ggmap(m)+geom_path(data=X3_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
#Miss 1 row
#' ## Display space-time cube. 
open3d()
# To get a bigger window than the default
par3d(windowRect = c(100, 100, 612, 612))
Sys.sleep(0.1) # Allow sluggish window managers to catch up
with(X3_df, plot3d(x,y,date, type="l", col=as.numeric(X3_df$Event_ID)))
(stcube<-with(X3_df, plot3d(x,y,date, type="l",col=as.numeric(cut(X3_df$dist,3)), alpha=0.4)))
#Timelag
head(timeLag(X3, units="mins"))
head(timeLag(X3, units="hours"))
median(timeLag(X3, units="hours"))
max(timeLag(X3, units="hours"))
max(timeLag(X3, units="days"))
plot(X3_df$individual.local.identifier..ID.~X3_df$date,xlab = "Sampling period",ylab = "Individuals ID Number")
#write.table(X3,file="c:/RWorkDir/jaguardatapaper/ID/X3.txt",row.names = F,quote=F,col.names=T,sep="\t")
#Histograms and plots
hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 3],
     main = 'Individual 3', xlab = "Time between fixes (h)")
plot(density.default(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 3],
     na.rm = T), main = 'Individual 3')
	 
breaks <- c(0, 0.5, 1.5, 2.5, 4.5, 8.5, 16.5, 32.5, 64.5, 128.5, 256.5, 1000)
breaks <- c(0, 1, 2*(1:256))
hh <- hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 3 &
                                     mov.data.diff$time.diff > 0 & 
                                     mov.data.diff$time.diff < 512],
									 main = 'Individual 3',xlab = 'Time between fixes (h)',
           plot = T, breaks = breaks)

#X4
X4=unstacked$X4
X4_df <- as(X4, "data.frame")
head(X4_df)
nrow(X4_df)
range(X4_df$date)
m <- get_map(bbox(extent(X4)*1.1), source="google", zoom=12, maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=X4_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
#' ## Display space-time cube. 
open3d()
# To get a bigger window than the default
par3d(windowRect = c(100, 100, 612, 612))
Sys.sleep(0.1) # Allow sluggish window managers to catch up
with(X4_df, plot3d(x,y,date, type="l", col=as.numeric(X4_df$Event_ID)))
(stcube<-with(X4_df, plot3d(x,y,date, type="l",col=as.numeric(cut(X4_df$dist,3)), alpha=0.4)))
#Timelag
head(timeLag(X4, units="mins"))
head(timeLag(X4, units="hours"))
median(timeLag(X4, units="hours"))
max(timeLag(X4, units="hours"))
max(timeLag(X4, units="days"))
plot(X4_df$individual.local.identifier..ID.~X4_df$date,xlab = "Sampling period",ylab = "Individuals ID Number")
#write.table(X4,file="c:/RWorkDir/jaguardatapaper/ID/X4.txt",row.names = F,quote=F,col.names=T,sep="\t")
#Histograms and plots
hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 4],
     main = 'Individual 4', xlab = "Time between fixes (h)")
plot(density.default(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 4],
     na.rm = T), main = 'Individual 4')
	 
breaks <- c(0, 0.5, 1.5, 2.5, 4.5, 8.5, 16.5, 32.5, 64.5, 128.5, 256.5, 1000)
breaks <- c(0, 1, 2*(1:256))
hh <- hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 4 &
                                     mov.data.diff$time.diff > 0 & 
                                     mov.data.diff$time.diff < 512],
									 main = 'Individual 4',xlab = 'Time between fixes (h)',
           plot = T, breaks = breaks)


#X5
X5=unstacked$X5
X5_df <- as(X5, "data.frame")
head(X5_df)
nrow(X5_df)
range(X5_df$date)
m <- get_map(bbox(extent(X5)*1.1), source="google", zoom=13, maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=X5_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
m <- get_map(bbox(extent(X5)*1.1), source="google", zoom=12, maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=X5_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
#' ## Display space-time cube. 
open3d()
# To get a bigger window than the default
par3d(windowRect = c(100, 100, 612, 612))
Sys.sleep(0.1) # Allow sluggish window managers to catch up
with(X5_df, plot3d(x,y,date, type="l", col=as.numeric(X5_df$Event_ID)))
(stcube<-with(X5_df, plot3d(x,y,date, type="l",col=as.numeric(cut(X5_df$dist,3)), alpha=0.4)))
#Timelag
head(timeLag(X5, units="mins"))
head(timeLag(X5, units="hours"))
median(timeLag(X5, units="hours"))
max(timeLag(X5, units="hours"))
max(timeLag(X5, units="days"))
plot(X5_df$individual.local.identifier..ID.~X5_df$date,xlab = "Sampling period",ylab = "Individuals ID Number")
#write.table(X5,file="c:/RWorkDir/jaguardatapaper/ID/X5.txt",row.names = F,quote=F,col.names=T,sep="\t")
#Histograms and plots
hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 5],
     main = 'Individual 5', xlab = "Time between fixes (h)")
plot(density.default(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 5],
     na.rm = T), main = 'Individual 5')
	 
breaks <- c(0, 0.5, 1.5, 2.5, 4.5, 8.5, 16.5, 32.5, 64.5, 128.5, 256.5, 1000)
breaks <- c(0, 1, 2*(1:256))
hh <- hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 5 &
                                     mov.data.diff$time.diff > 0 & 
                                     mov.data.diff$time.diff < 512],
									 main = 'Individual 5',xlab = 'Time between fixes (h)',
           plot = T, breaks = breaks)

#X6
X6=unstacked$X6
X6_df <- as(X6, "data.frame")
head(X6_df)
nrow(X6_df)
range(X6_df$date)
m <- get_map(bbox(extent(X6)*1.1), source="google", zoom=11, maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=X6_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
m <- get_map(bbox(extent(X6)*1.1), source="google", zoom=12, maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=X6_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
#' ## Display space-time cube. 
open3d()
# To get a bigger window than the default
par3d(windowRect = c(100, 100, 612, 612))
Sys.sleep(0.1) # Allow sluggish window managers to catch up
with(X6_df, plot3d(x,y,date, type="l", col=as.numeric(X6_df$Event_ID)))
(stcube<-with(X6_df, plot3d(x,y,date, type="l",col=as.numeric(cut(X6_df$dist,3)), alpha=0.4)))
#Timelag
head(timeLag(X6, units="mins"))
head(timeLag(X6, units="hours"))
median(timeLag(X6, units="hours"))
max(timeLag(X6, units="hours"))
max(timeLag(X6, units="days"))
plot(X6_df$individual.local.identifier..ID.~X6_df$date,xlab = "Sampling period",ylab = "Individuals ID Number")
#write.table(X6,file="c:/RWorkDir/jaguardatapaper/ID/X6.txt",row.names = F,quote=F,col.names=T,sep="\t")
#Histograms and plots
hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 6],
     main = 'Individual 6', xlab = "Time between fixes (h)")
plot(density.default(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 6],
     na.rm = T), main = 'Individual 6')
	 
breaks <- c(0, 0.5, 1.5, 2.5, 4.5, 8.5, 16.5, 32.5, 64.5, 128.5, 256.5, 1000)
breaks <- c(0, 1, 2*(1:256))
hh <- hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 6 &
                                     mov.data.diff$time.diff > 0 & 
                                     mov.data.diff$time.diff < 512],
									 main = 'Individual 6',xlab = 'Time between fixes (h)',
           plot = T, breaks = breaks)

#X7
X7=unstacked$X7
X7_df <- as(X7, "data.frame")
head(X7_df)
nrow(X7_df)
range(X7_df$date)
m <- get_map(bbox(extent(X7)*1.1), source="google", zoom=12, maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=X7_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
m <- get_map(bbox(extent(X7)*1.1), source="google", zoom=11, maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=X7_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
#' ## Display space-time cube. 
open3d()
# To get a bigger window than the default
par3d(windowRect = c(100, 100, 612, 612))
Sys.sleep(0.1) # Allow sluggish window managers to catch up
with(X7_df, plot3d(x,y,date, type="l", col=as.numeric(X7_df$Event_ID)))
(stcube<-with(X7_df, plot3d(x,y,date, type="l",col=as.numeric(cut(X7_df$dist,3)), alpha=0.4)))
#Timelag
head(timeLag(X7, units="mins"))
head(timeLag(X7, units="hours"))
median(timeLag(X7, units="hours"))
max(timeLag(X7, units="hours"))
max(timeLag(X7, units="days"))
plot(X7_df$individual.local.identifier..ID.~X7_df$date,xlab = "Sampling period",ylab = "Individuals ID Number")
#write.table(X7,file="c:/RWorkDir/jaguardatapaper/ID/X7.txt",row.names = F,quote=F,col.names=T,sep="\t")
#Histograms and plots
hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 7],
     main = 'Individual 7', xlab = "Time between fixes (h)")
plot(density.default(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 7],
     na.rm = T), main = 'Individual 7')
	 
breaks <- c(0, 0.5, 1.5, 2.5, 4.5, 8.5, 16.5, 32.5, 64.5, 128.5, 256.5, 1000)
breaks <- c(0, 1, 2*(1:256))
hh <- hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 7 &
                                     mov.data.diff$time.diff > 0 & 
                                     mov.data.diff$time.diff < 512],
									 main = 'Individual 7',xlab = 'Time between fixes (h)',
           plot = T, breaks = breaks)

#X8
X8=unstacked$X8
X8_df <- as(X8, "data.frame")
head(X8_df)
nrow(X8_df)
range(X8_df$date)
m <- get_map(bbox(extent(X8)*1.1), source="google", zoom=11, maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=X8_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
#m <- get_map(bbox(extent(X8)*1.1), source="google", zoom=12, maptype="satellite") #define the box for map based on the trajectory coordinates
#ggmap(m)+geom_path(data=X8_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
#6 row warning
#' ## Display space-time cube. 
open3d()
# To get a bigger window than the default
par3d(windowRect = c(100, 100, 612, 612))
Sys.sleep(0.1) # Allow sluggish window managers to catch up
with(X8_df, plot3d(x,y,date, type="l", col=as.numeric(X8_df$Event_ID)))
(stcube<-with(X8_df, plot3d(x,y,date, type="l",col=as.numeric(cut(X8_df$dist,3)), alpha=0.4)))
#Timelag
head(timeLag(X8, units="mins"))
head(timeLag(X8, units="hours"))
median(timeLag(X8, units="hours"))
max(timeLag(X8, units="hours"))
max(timeLag(X8, units="days"))
plot(X8_df$individual.local.identifier..ID.~X8_df$date,xlab = "Sampling period",ylab = "Individuals ID Number")
#write.table(X8,file="c:/RWorkDir/jaguardatapaper/ID/X8.txt",row.names = F,quote=F,col.names=T,sep="\t")
#Histograms and plots
hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 8],
     main = 'Individual 8', xlab = "Time between fixes (h)")
plot(density.default(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 8],
     na.rm = T), main = 'Individual 8')
	 
breaks <- c(0, 0.5, 1.5, 2.5, 4.5, 8.5, 16.5, 32.5, 64.5, 128.5, 256.5, 1000)
breaks <- c(0, 1, 2*(1:256))
hh <- hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 8 &
                                     mov.data.diff$time.diff > 0 & 
                                     mov.data.diff$time.diff < 512],
									 main = 'Individual 8',xlab = 'Time between fixes (h)',
           plot = T, breaks = breaks)

#X9
X9=unstacked$X9
X9_df <- as(X9, "data.frame")
head(X9_df)
nrow(X9_df)
range(X9_df$date)
m <- get_map(bbox(extent(X9)*1.1), source="google", zoom=13, maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=X9_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
m <- get_map(bbox(extent(X9)*1.1), source="google", zoom=12, maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=X9_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
#' ## Display space-time cube. 
open3d()
# To get a bigger window than the default
par3d(windowRect = c(100, 100, 612, 612))
Sys.sleep(0.1) # Allow sluggish window managers to catch up
with(X9_df, plot3d(x,y,date, type="l", col=as.numeric(X9_df$Event_ID)))
(stcube<-with(X9_df, plot3d(x,y,date, type="l",col=as.numeric(cut(X9_df$dist,3)), alpha=0.4)))
#Timelag
head(timeLag(X9, units="mins"))
head(timeLag(X9, units="hours"))
median(timeLag(X9, units="hours"))
max(timeLag(X9, units="hours"))
max(timeLag(X9, units="days"))
plot(X9_df$individual.local.identifier..ID.~X9_df$date,xlab = "Sampling period",ylab = "Individuals ID Number")
#write.table(X9,file="c:/RWorkDir/jaguardatapaper/ID/X9.txt",row.names = F,quote=F,col.names=T,sep="\t")
#Histograms and plots
hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 9],
     main = 'Individual 9', xlab = "Time between fixes (h)")
plot(density.default(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 9],
     na.rm = T), main = 'Individual 9')
	 
breaks <- c(0, 0.5, 1.5, 2.5, 4.5, 8.5, 16.5, 32.5, 64.5, 128.5, 256.5, 1000)
breaks <- c(0, 1, 2*(1:256))
hh <- hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 9 &
                                     mov.data.diff$time.diff > 0 & 
                                     mov.data.diff$time.diff < 512],
									 main = 'Individual 9',xlab = 'Time between fixes (h)',
           plot = T, breaks = breaks)


#X10	
unstacked <- split(move.data)
X10=unstacked$X10
X10_df <- as(X10, "data.frame")
head(X10_df)
nrow(X10_df)
range(X10_df$date)
m <- get_map(bbox(extent(X10)*1.1), source="google", zoom=12, maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=X10_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
#' ## Display space-time cube. 
open3d()
# To get a bigger window than the default
par3d(windowRect = c(100, 100, 612, 612))
Sys.sleep(0.1) # Allow sluggish window managers to catch up
with(X10_df, plot3d(x,y,date, type="l", col=as.numeric(X10_df$Event_ID)))
(stcube<-with(X10_df, plot3d(x,y,date, type="l",col=as.numeric(cut(X10_df$dist,3)), alpha=0.4)))
#Timelag
head(timeLag(X10, units="mins"))
head(timeLag(X10, units="hours"))
median(timeLag(X10, units="hours"))
max(timeLag(X10, units="hours"))
max(timeLag(X10, units="days"))
plot(X10_df$individual.local.identifier..ID.~X10_df $date,xlab = "Sampling period",ylab = "Individuals ID Number")
#write.table(X10,file="c:/RWorkDir/jaguardatapaper/ID/X10.txt",row.names = F,quote=F,col.names=T,sep="\t")
#Histograms and plots
hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 10],
     main = 'Individual 10', xlab = "Time between fixes (h)")
plot(density.default(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 10],
     na.rm = T), main = 'Individual 10')
	 
breaks <- c(0, 0.5, 1.5, 2.5, 4.5, 8.5, 16.5, 32.5, 64.5, 128.5, 256.5, 1000)
breaks <- c(0, 1, 2*(1:256))
hh <- hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 10 &
                                     mov.data.diff$time.diff > 0 & 
                                     mov.data.diff$time.diff < 512],
									 main = 'Individual 10',xlab = 'Time between fixes (h)',
           plot = T, breaks = breaks)


#X11	
unstacked <- split(move.data)
X11=unstacked$X11
X11_df <- as(X11, "data.frame")
head(X11_df)
nrow(X11_df)
range(X11_df$date)
m <- get_map(bbox(extent(X11)*1.1), source="google", zoom=12, maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=X11_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
#' ## Display space-time cube. 
open3d()
# To get a bigger window than the default
par3d(windowRect = c(100, 100, 612, 612))
Sys.sleep(0.1) # Allow sluggish window managers to catch up
with(X11_df, plot3d(x,y,date, type="l", col=as.numeric(X11_df$Event_ID)))
(stcube<-with(X11_df, plot3d(x,y,date, type="l",col=as.numeric(cut(X11_df$dist,3)), alpha=0.4)))
#Timelag
head(timeLag(X11, units="mins"))
head(timeLag(X11, units="hours"))
median(timeLag(X11, units="hours"))
max(timeLag(X11, units="hours"))
max(timeLag(X11, units="days"))
plot(X11_df$individual.local.identifier..ID.~X11_df $date,xlab = "Sampling period",ylab = "Individuals ID Number")
#write.table(X11,file="c:/RWorkDir/jaguardatapaper/ID/X11.txt",row.names = F,quote=F,col.names=T,sep="\t")
#Histograms and plots
hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 11],
     main = 'Individual 11', xlab = "Time between fixes (h)")
plot(density.default(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 11],
     na.rm = T), main = 'Individual 11')
	 
breaks <- c(0, 0.5, 1.5, 2.5, 4.5, 8.5, 16.5, 32.5, 64.5, 128.5, 256.5, 1000)
breaks <- c(0, 1, 2*(1:256))
hh <- hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 11 &
                                     mov.data.diff$time.diff > 0 & 
                                     mov.data.diff$time.diff < 512],
									 main = 'Individual 11',xlab = 'Time between fixes (h)',
           plot = T, breaks = breaks)


#X12	
unstacked <- split(move.data)
X12=unstacked$X12
X12_df <- as(X12, "data.frame")
head(X12_df)
nrow(X12_df)
range(X12_df$date)
m <- get_map(bbox(extent(X12)*1.1), source="google", zoom=11, maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=X12_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
#m <- get_map(bbox(extent(X12)*1.1), source="google", zoom=12, maptype="satellite") #define the box for map based on the trajectory coordinates
#ggmap(m)+geom_path(data=X12_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
#' ## Display space-time cube. 
open3d()
# To get a bigger window than the default
par3d(windowRect = c(100, 100, 612, 612))
Sys.sleep(0.1) # Allow sluggish window managers to catch up
with(X12_df, plot3d(x,y,date, type="l", col=as.numeric(X12_df$Event_ID)))
(stcube<-with(X12_df, plot3d(x,y,date, type="l",col=as.numeric(cut(X12_df$dist,3)), alpha=0.4)))
#Timelag
head(timeLag(X12, units="mins"))
head(timeLag(X12, units="hours"))
median(timeLag(X12, units="hours"))
max(timeLag(X12, units="hours"))
max(timeLag(X12, units="days"))
plot(X12_df$individual.local.identifier..ID.~X12_df $date,xlab = "Sampling period",ylab = "Individuals ID Number")
#write.table(X12,file="c:/RWorkDir/jaguardatapaper/ID/X12.txt",row.names = F,quote=F,col.names=T,sep="\t")
#Histograms and plots
hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 12],
     main = 'Individual 12', xlab = "Time between fixes (h)")
plot(density.default(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 12],
     na.rm = T), main = 'Individual 12')
	 
breaks <- c(0, 0.5, 1.5, 2.5, 4.5, 8.5, 16.5, 32.5, 64.5, 128.5, 256.5, 1000)
breaks <- c(0, 1, 2*(1:256))
hh <- hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 12 &
                                     mov.data.diff$time.diff > 0 & 
                                     mov.data.diff$time.diff < 512],
									 main = 'Individual 12',xlab = 'Time between fixes (h)',
           plot = T, breaks = breaks)

#X13	
unstacked <- split(move.data)
X13=unstacked$X13
X13_df <- as(X13, "data.frame")
head(X13_df)
nrow(X13_df)
range(X13_df$date)
m <- get_map(bbox(extent(X13)*1.1), source="google", zoom=12, maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=X13_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
#m <- get_map(bbox(extent(X13)*1.1), source="google", zoom=13, maptype="satellite") #define the box for map based on the trajectory coordinates
#ggmap(m)+geom_path(data=X13_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
#Removed 15 rows containing missing values (geom_path).
#' ## Display space-time cube. 
open3d()
# To get a bigger window than the default
par3d(windowRect = c(100, 100, 612, 612))
Sys.sleep(0.1) # Allow sluggish window managers to catch up
with(X13_df, plot3d(x,y,date, type="l", col=as.numeric(X13_df$Event_ID)))
(stcube<-with(X13_df, plot3d(x,y,date, type="l",col=as.numeric(cut(X13_df$dist,3)), alpha=0.4)))
#Timelag
head(timeLag(X13, units="mins"))
head(timeLag(X13, units="hours"))
median(timeLag(X13, units="hours"))
max(timeLag(X13, units="hours"))
max(timeLag(X13, units="days"))
plot(X13_df$individual.local.identifier..ID.~X13_df $date,xlab = "Sampling period",ylab = "Individuals ID Number")
#write.table(X13,file="c:/RWorkDir/jaguardatapaper/ID/X13.txt",row.names = F,quote=F,col.names=T,sep="\t")
#Histograms and plots
hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 13],
     main = 'Individual 13', xlab = "Time between fixes (h)")
plot(density.default(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 13],
     na.rm = T), main = 'Individual 13')
	 
breaks <- c(0, 0.5, 1.5, 2.5, 4.5, 8.5, 16.5, 32.5, 64.5, 128.5, 256.5, 1000)
breaks <- c(0, 1, 2*(1:256))
hh <- hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 13 &
                                     mov.data.diff$time.diff > 0 & 
                                     mov.data.diff$time.diff < 512],
									 main = 'Individual 13',xlab = 'Time between fixes (h)',
           plot = T, breaks = breaks)

#X14
unstacked <- split(move.data)
X14=unstacked$X14
X14_df <- as(X14, "data.frame")
head(X14_df)
nrow(X14_df)
range(X14_df$date)
m <- get_map(bbox(extent(X14)*1.1), source="google", zoom=12, maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=X14_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
#' ## Display space-time cube. 
open3d()
# To get a bigger window than the default
par3d(windowRect = c(100, 100, 612, 612))
Sys.sleep(0.1) # Allow sluggish window managers to catch up
with(X14_df, plot3d(x,y,date, type="l", col=as.numeric(X14_df$Event_ID)))
(stcube<-with(X14_df, plot3d(x,y,date, type="l",col=as.numeric(cut(X14_df$dist,3)), alpha=0.4)))
#Timelag
head(timeLag(X14, units="mins"))
head(timeLag(X14, units="hours"))
median(timeLag(X14, units="hours"))
max(timeLag(X14, units="hours"))
max(timeLag(X14, units="days"))
plot(X14_df$individual.local.identifier..ID.~X14_df $date,xlab = "Sampling period",ylab = "Individuals ID Number")
#write.table(X14,file="c:/RWorkDir/jaguardatapaper/ID/X14.txt",row.names = F,quote=F,col.names=T,sep="\t")
#Histograms and plots
hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 14],
     main = 'Individual 14', xlab = "Time between fixes (h)")
plot(density.default(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 14],
     na.rm = T), main = 'Individual 14')
	 
breaks <- c(0, 0.5, 1.5, 2.5, 4.5, 8.5, 16.5, 32.5, 64.5, 128.5, 256.5, 1000)
breaks <- c(0, 1, 2*(1:256))
hh <- hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 14 &
                                     mov.data.diff$time.diff > 0 & 
                                     mov.data.diff$time.diff < 512],
									 main = 'Individual 14',xlab = 'Time between fixes (h)',
           plot = T, breaks = breaks)

#X15
unstacked <- split(move.data)
X15=unstacked$X15
X15_df <- as(X15, "data.frame")
head(X15_df)
nrow(X15_df)
range(X15_df$date)
m <- get_map(bbox(extent(X15)*1.1), source="google", zoom=11, maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=X15_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
m <- get_map(bbox(extent(X15)*1.1), source="google", zoom=12, maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=X15_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
#' ## Display space-time cube. 
open3d()
# To get a bigger window than the default
par3d(windowRect = c(100, 100, 612, 612))
Sys.sleep(0.1) # Allow sluggish window managers to catch up
with(X1_df, plot3d(x,y,date, type="l", col=as.numeric(X1_df$Event_ID)))
(stcube<-with(X15_df, plot3d(x,y,date, type="l",col=as.numeric(cut(X15_df$dist,3)), alpha=0.4)))
#Timelag
head(timeLag(X15, units="mins"))
head(timeLag(X15, units="hours"))
median(timeLag(X15, units="hours"))
max(timeLag(X15, units="hours"))
max(timeLag(X15, units="days"))
plot(X15_df$individual.local.identifier..ID.~X15_df $date,xlab = "Sampling period",ylab = "Individuals ID Number")
#write.table(X15,file="c:/RWorkDir/jaguardatapaper/ID/X15.txt",row.names = F,quote=F,col.names=T,sep="\t")
#Histograms and plots
hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 15],
     main = 'Individual 15', xlab = "Time between fixes (h)")
plot(density.default(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 15],
     na.rm = T), main = 'Individual 15')
	 
breaks <- c(0, 0.5, 1.5, 2.5, 4.5, 8.5, 16.5, 32.5, 64.5, 128.5, 256.5, 1000)
breaks <- c(0, 1, 2*(1:256))
hh <- hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 15 &
                                     mov.data.diff$time.diff > 0 & 
                                     mov.data.diff$time.diff < 512],
									 main = 'Individual 15',xlab = 'Time between fixes (h)',
           plot = T, breaks = breaks)

#X16
unstacked <- split(move.data)
X16=unstacked$X16
X16_df <- as(X16, "data.frame")
head(X16_df)
nrow(X16_df)
range(X16_df$date)
m <- get_map(bbox(extent(X16)*1.1), source="google", zoom=11, maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=X16_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
m <- get_map(bbox(extent(X16)*1.1), source="google", zoom=12, maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=X16_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
#' ## Display space-time cube. 
open3d()
# To get a bigger window than the default
par3d(windowRect = c(100, 100, 612, 612))
Sys.sleep(0.1) # Allow sluggish window managers to catch up
with(X16_df, plot3d(x,y,date, type="l", col=as.numeric(X16_df$Event_ID)))
(stcube<-with(X16_df, plot3d(x,y,date, type="l",col=as.numeric(cut(X16_df$dist,3)), alpha=0.4)))
#Timelag
head(timeLag(X16, units="mins"))
head(timeLag(X16, units="hours"))
median(timeLag(X16, units="hours"))
max(timeLag(X16, units="hours"))
max(timeLag(X16, units="days"))
plot(X16_df$individual.local.identifier..ID.~X16_df $date,xlab = "Sampling period",ylab = "Individuals ID Number")
#write.table(X16,file="c:/RWorkDir/jaguardatapaper/ID/X16.txt",row.names = F,quote=F,col.names=T,sep="\t")
#Histograms and plots
hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 16],
     main = 'Individual 16', xlab = "Time between fixes (h)")
plot(density.default(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 16],
     na.rm = T), main = 'Individual 16')
	 
breaks <- c(0, 0.5, 1.5, 2.5, 4.5, 8.5, 16.5, 32.5, 64.5, 128.5, 256.5, 1000)
breaks <- c(0, 1, 2*(1:256))
hh <- hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 16 &
                                     mov.data.diff$time.diff > 0 & 
                                     mov.data.diff$time.diff < 512],
									 main = 'Individual 16',xlab = 'Time between fixes (h)',
           plot = T, breaks = breaks)

#X17
unstacked <- split(move.data)
X17=unstacked$X17
X17_df <- as(X17, "data.frame")
head(X17_df)
nrow(X17_df)
range(X17_df$date)
m <- get_map(bbox(extent(X17)*1.1), source="google", zoom=11, maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=X17_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
#' ## Display space-time cube. 
open3d()
# To get a bigger window than the default
par3d(windowRect = c(100, 100, 612, 612))
Sys.sleep(0.1) # Allow sluggish window managers to catch up
with(X17_df, plot3d(x,y,date, type="l", col=as.numeric(X17_df$Event_ID)))
(stcube<-with(X17_df, plot3d(x,y,date, type="l",col=as.numeric(cut(X17_df$dist,3)), alpha=0.4)))
#Timelag
head(timeLag(X17, units="mins"))
head(timeLag(X17, units="hours"))
median(timeLag(X17, units="hours"))
max(timeLag(X17, units="hours"))
max(timeLag(X17, units="days"))
plot(X17_df$individual.local.identifier..ID.~X17_df $date,xlab = "Sampling period",ylab = "Individuals ID Number")
#write.table(X17,file="c:/RWorkDir/jaguardatapaper/ID/X17.txt",row.names = F,quote=F,col.names=T,sep="\t")
#Histograms and plots
hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 17],
     main = 'Individual 17', xlab = "Time between fixes (h)")
plot(density.default(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 17],
     na.rm = T), main = 'Individual 17')
	 
breaks <- c(0, 0.5, 1.5, 2.5, 4.5, 8.5, 16.5, 32.5, 64.5, 128.5, 256.5, 1000)
breaks <- c(0, 1, 2*(1:256))
hh <- hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 17 &
                                     mov.data.diff$time.diff > 0 & 
                                     mov.data.diff$time.diff < 512],
									 main = 'Individual 17',xlab = 'Time between fixes (h)',
           plot = T, breaks = breaks)

#X18
unstacked <- split(move.data)
X18=unstacked$X18
X18_df <- as(X18, "data.frame")
head(X18_df)
nrow(X18_df)
range(X18_df$date)
m <- get_map(bbox(extent(X18)*1.1), source="google", zoom=12, maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=X18_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
#' ## Display space-time cube. 
open3d()
# To get a bigger window than the default
par3d(windowRect = c(100, 100, 612, 612))
Sys.sleep(0.1) # Allow sluggish window managers to catch up
with(X18_df, plot3d(x,y,date, type="l", col=as.numeric(X18_df$Event_ID)))
(stcube<-with(X18_df, plot3d(x,y,date, type="l",col=as.numeric(cut(X18_df$dist,3)), alpha=0.4)))
#Timelag
head(timeLag(X18, units="mins"))
head(timeLag(X18, units="hours"))
median(timeLag(X18, units="hours"))
max(timeLag(X18, units="hours"))
max(timeLag(X18, units="days"))
plot(X18_df$individual.local.identifier..ID.~X18_df $date,xlab = "Sampling period",ylab = "Individuals ID Number")
#write.table(X18,file="c:/RWorkDir/jaguardatapaper/ID/X18.txt",row.names = F,quote=F,col.names=T,sep="\t")
#Histograms and plots
hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 18],
     main = 'Individual 18', xlab = "Time between fixes (h)")
plot(density.default(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 18],
     na.rm = T), main = 'Individual 18')
	 
breaks <- c(0, 0.5, 1.5, 2.5, 4.5, 8.5, 16.5, 32.5, 64.5, 128.5, 256.5, 1000)
breaks <- c(0, 1, 2*(1:256))
hh <- hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 18 &
                                     mov.data.diff$time.diff > 0 & 
                                     mov.data.diff$time.diff < 512],
									 main = 'Individual 18',xlab = 'Time between fixes (h)',
           plot = T, breaks = breaks)

#X19
unstacked <- split(move.data)
X19=unstacked$X19
X19_df <- as(X19, "data.frame")
head(X19_df)
nrow(X19_df)
range(X19_df$date)
m <- get_map(bbox(extent(X19)*1.1), source="google", zoom=13, maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=X19_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
#' ## Display space-time cube. 
open3d()
# To get a bigger window than the default
par3d(windowRect = c(100, 100, 612, 612))
Sys.sleep(0.1) # Allow sluggish window managers to catch up
with(X19_df, plot3d(x,y,date, type="l", col=as.numeric(X19_df$Event_ID)))
(stcube<-with(X19_df, plot3d(x,y,date, type="l",col=as.numeric(cut(X19_df$dist,3)), alpha=0.4)))
#Timelag
head(timeLag(X19, units="mins"))
head(timeLag(X19, units="hours"))
median(timeLag(X19, units="hours"))
max(timeLag(X19, units="hours"))
max(timeLag(X19, units="days"))
plot(X19_df$individual.local.identifier..ID.~X19_df $date,xlab = "Sampling period",ylab = "Individuals ID Number")
#write.table(X19,file="c:/RWorkDir/jaguardatapaper/ID/X19.txt",row.names = F,quote=F,col.names=T,sep="\t")
#Histograms and plots
hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 19],
     main = 'Individual 19', xlab = "Time between fixes (h)")
plot(density.default(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 19],
     na.rm = T), main = 'Individual 19')
	 
breaks <- c(0, 0.5, 1.5, 2.5, 4.5, 8.5, 16.5, 32.5, 64.5, 128.5, 256.5, 1000)
breaks <- c(0, 1, 2*(1:256))
hh <- hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 19 &
                                     mov.data.diff$time.diff > 0 & 
                                     mov.data.diff$time.diff < 512],
									 main = 'Individual 19',xlab = 'Time between fixes (h)',
           plot = T, breaks = breaks)

#X20
unstacked <- split(move.data)
X20=unstacked$X20
X20_df <- as(X20, "data.frame")
head(X20_df)
nrow(X20_df)
range(X20_df$date)
m <- get_map(bbox(extent(X20)*1.1), source="google", zoom=10, maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=X20_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
m <- get_map(bbox(extent(X20)*1.1), source="google", zoom=11, maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=X20_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
#' ## Display space-time cube. 
open3d()
# To get a bigger window than the default
par3d(windowRect = c(100, 100, 612, 612))
Sys.sleep(0.1) # Allow sluggish window managers to catch up
with(X20_df, plot3d(x,y,date, type="l", col=as.numeric(X20_df$Event_ID)))
(stcube<-with(X20_df, plot3d(x,y,date, type="l",col=as.numeric(cut(X20_df$dist,3)), alpha=0.4)))
#Timelag
head(timeLag(X20, units="mins"))
head(timeLag(X20, units="hours"))
median(timeLag(X20, units="hours"))
max(timeLag(X20, units="hours"))
max(timeLag(X20, units="days"))
plot(X20_df$individual.local.identifier..ID.~X20_df $date,xlab = "Sampling period",ylab = "Individuals ID Number")
#write.table(X20,file="c:/RWorkDir/jaguardatapaper/ID/X20.txt",row.names = F,quote=F,col.names=T,sep="\t")
#Histograms and plots
hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 20],
     main = 'Individual 20', xlab = "Time between fixes (h)")
plot(density.default(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 20],
     na.rm = T), main = 'Individual 20')
	 
breaks <- c(0, 0.5, 1.5, 2.5, 4.5, 8.5, 16.5, 32.5, 64.5, 128.5, 256.5, 1000)
breaks <- c(0, 1, 2*(1:256))
hh <- hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 20 &
                                     mov.data.diff$time.diff > 0 & 
                                     mov.data.diff$time.diff < 512],
									 main = 'Individual 20',xlab = 'Time between fixes (h)',
           plot = T, breaks = breaks)

#X21
unstacked <- split(move.data)
X21=unstacked$X21
X21_df <- as(X21, "data.frame")
head(X21_df)
nrow(X21_df)
range(X21_df$date)
m <- get_map(bbox(extent(X21)*1.1), source="google", zoom=12, maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=X21_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
m <- get_map(bbox(extent(X21)*1.1), source="google", zoom=13, maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=X21_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
#' ## Display space-time cube. 
open3d()
# To get a bigger window than the default
par3d(windowRect = c(100, 100, 612, 612))
Sys.sleep(0.1) # Allow sluggish window managers to catch up
with(X21_df, plot3d(x,y,date, type="l", col=as.numeric(X21_df$Event_ID)))
(stcube<-with(X21_df, plot3d(x,y,date, type="l",col=as.numeric(cut(X21_df$dist,3)), alpha=0.4)))
#Timelag
head(timeLag(X21, units="mins"))
head(timeLag(X21, units="hours"))
median(timeLag(X21, units="hours"))
max(timeLag(X21, units="hours"))
max(timeLag(X21, units="days"))
plot(X21_df$individual.local.identifier..ID.~X21_df $date,xlab = "Sampling period",ylab = "Individuals ID Number")
#write.table(X21,file="c:/RWorkDir/jaguardatapaper/ID/X21.txt",row.names = F,quote=F,col.names=T,sep="\t")
#Histograms and plots
hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 21],
     main = 'Individual 21', xlab = "Time between fixes (h)")
plot(density.default(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 21],
     na.rm = T), main = 'Individual 21')
	 
breaks <- c(0, 0.5, 1.5, 2.5, 4.5, 8.5, 16.5, 32.5, 64.5, 128.5, 256.5, 1000)
breaks <- c(0, 1, 2*(1:256))
hh <- hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 21 &
                                     mov.data.diff$time.diff > 0 & 
                                     mov.data.diff$time.diff < 512],
									 main = 'Individual 21',xlab = 'Time between fixes (h)',
           plot = T, breaks = breaks)


#X22
unstacked <- split(move.data)
X22=unstacked$X22
X22_df <- as(X22, "data.frame")
head(X22_df)
nrow(X22_df)
range(X22_df$date)
m <- get_map(bbox(extent(X22)*1.1), source="google", zoom=12, maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=X22_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
m <- get_map(bbox(extent(X22)*1.1), source="google", zoom=13, maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=X22_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
#' ## Display space-time cube. 
open3d()
# To get a bigger window than the default
par3d(windowRect = c(100, 100, 612, 612))
Sys.sleep(0.1) # Allow sluggish window managers to catch up
with(X22_df, plot3d(x,y,date, type="l", col=as.numeric(X22_df$Event_ID)))
(stcube<-with(X22_df, plot3d(x,y,date, type="l",col=as.numeric(cut(X22_df$dist,3)), alpha=0.4)))
#Timelag
head(timeLag(X22, units="mins"))
head(timeLag(X22, units="hours"))
median(timeLag(X22, units="hours"))
max(timeLag(X22, units="hours"))
max(timeLag(X22, units="days"))
plot(X22_df$individual.local.identifier..ID.~X22_df $date,xlab = "Sampling period",ylab = "Individuals ID Number")
#write.table(X22,file="c:/RWorkDir/jaguardatapaper/ID/X22.txt",row.names = F,quote=F,col.names=T,sep="\t")
#Histograms and plots
hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 22],
     main = 'Individual 22', xlab = "Time between fixes (h)")
plot(density.default(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 22],
     na.rm = T), main = 'Individual 22')
	 
breaks <- c(0, 0.5, 1.5, 2.5, 4.5, 8.5, 16.5, 32.5, 64.5, 128.5, 256.5, 1000)
breaks <- c(0, 1, 2*(1:256))
hh <- hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 22 &
                                     mov.data.diff$time.diff > 0 & 
                                     mov.data.diff$time.diff < 512],
									 main = 'Individual 22',xlab = 'Time between fixes (h)',
           plot = T, breaks = breaks)

#X23
unstacked <- split(move.data)
X23=unstacked$X23
X23_df <- as(X23, "data.frame")
head(X23_df)
nrow(X23_df)
range(X23_df$date)
m <- get_map(bbox(extent(X23)*1.1), source="google", zoom=13, maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=X23_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
#' ## Display space-time cube. 
open3d()
# To get a bigger window than the default
par3d(windowRect = c(100, 100, 612, 612))
Sys.sleep(0.1) # Allow sluggish window managers to catch up
with(X23_df, plot3d(x,y,date, type="l", col=as.numeric(X23_df$Event_ID)))
(stcube<-with(X23_df, plot3d(x,y,date, type="l",col=as.numeric(cut(X23_df$dist,3)), alpha=0.4)))
#Timelag
head(timeLag(X23, units="mins"))
head(timeLag(X23, units="hours"))
median(timeLag(X23, units="hours"))
max(timeLag(X23, units="hours"))
max(timeLag(X23, units="days"))
plot(X23_df$individual.local.identifier..ID.~X23_df $date,xlab = "Sampling period",ylab = "Individuals ID Number")
#write.table(X23,file="c:/RWorkDir/jaguardatapaper/ID/X23.txt",row.names = F,quote=F,col.names=T,sep="\t")
#Histograms and plots
hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 23],
     main = 'Individual 23', xlab = "Time between fixes (h)")
plot(density.default(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 23],
     na.rm = T), main = 'Individual 23')
	 
breaks <- c(0, 0.5, 1.5, 2.5, 4.5, 8.5, 16.5, 32.5, 64.5, 128.5, 256.5, 1000)
breaks <- c(0, 1, 2*(1:256))
hh <- hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 23 &
                                     mov.data.diff$time.diff > 0 & 
                                     mov.data.diff$time.diff < 512],
									 main = 'Individual 23',xlab = 'Time between fixes (h)',
           plot = T, breaks = breaks)

#X24
unstacked <- split(move.data)
X24=unstacked$X24
X24_df <- as(X24, "data.frame")
head(X24_df)
nrow(X24_df)
range(X24_df$date)
m <- get_map(bbox(extent(X24)*1.1), source="google", zoom=11, maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=X24_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
#' ## Display space-time cube. 
open3d()
# To get a bigger window than the default
par3d(windowRect = c(100, 100, 612, 612))
Sys.sleep(0.1) # Allow sluggish window managers to catch up
with(X24_df, plot3d(x,y,date, type="l", col=as.numeric(X24_df$Event_ID)))
(stcube<-with(X24_df, plot3d(x,y,date, type="l",col=as.numeric(cut(X24_df$dist,3)), alpha=0.4)))
#Timelag
head(timeLag(X24, units="mins"))
head(timeLag(X24, units="hours"))
median(timeLag(X24, units="hours"))
max(timeLag(X24, units="hours"))
max(timeLag(X24, units="days"))
plot(X24_df$individual.local.identifier..ID.~X24_df $date,xlab = "Sampling period",ylab = "Individuals ID Number")
#write.table(X24,file="c:/RWorkDir/jaguardatapaper/ID/X24.txt",row.names = F,quote=F,col.names=T,sep="\t")
#Histograms and plots
hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 24],
     main = 'Individual 24', xlab = "Time between fixes (h)")
plot(density.default(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 24],
     na.rm = T), main = 'Individual 24')
	 
breaks <- c(0, 0.5, 1.5, 2.5, 4.5, 8.5, 16.5, 32.5, 64.5, 128.5, 256.5, 1000)
breaks <- c(0, 1, 2*(1:256))
hh <- hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 24 &
                                     mov.data.diff$time.diff > 0 & 
                                     mov.data.diff$time.diff < 512],
									 main = 'Individual 24',xlab = 'Time between fixes (h)',
           plot = T, breaks = breaks)

#X25
unstacked <- split(move.data)
X25=unstacked$X25
X25_df <- as(X25, "data.frame")
head(X25_df)
nrow(X25_df)
range(X25_df$date)
m <- get_map(bbox(extent(X25)*1.1), source="google", zoom=12, maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=X25_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
#m <- get_map(bbox(extent(X25)*1.1), source="google", zoom=13, maptype="satellite") #define the box for map based on the trajectory coordinates
#ggmap(m)+geom_path(data=X25_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
#' ## Display space-time cube. 
open3d()
# To get a bigger window than the default
par3d(windowRect = c(100, 100, 612, 612))
Sys.sleep(0.1) # Allow sluggish window managers to catch up
with(X25_df, plot3d(x,y,date, type="l", col=as.numeric(X25_df$Event_ID)))
(stcube<-with(X25_df, plot3d(x,y,date, type="l",col=as.numeric(cut(X25_df$dist,3)), alpha=0.4)))
#Timelag
head(timeLag(X25, units="mins"))
head(timeLag(X25, units="hours"))
median(timeLag(X25, units="hours"))
max(timeLag(X25, units="hours"))
max(timeLag(X25, units="days"))
plot(X25_df$individual.local.identifier..ID.~X25_df $date,xlab = "Sampling period",ylab = "Individuals ID Number")
#write.table(X25,file="c:/RWorkDir/jaguardatapaper/ID/X25.txt",row.names = F,quote=F,col.names=T,sep="\t")
#Histograms and plots
hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 25],
     main = 'Individual 25', xlab = "Time between fixes (h)")
plot(density.default(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 25],
     na.rm = T), main = 'Individual 25')
	 
breaks <- c(0, 0.5, 1.5, 2.5, 4.5, 8.5, 16.5, 32.5, 64.5, 128.5, 256.5, 1000)
breaks <- c(0, 1, 2*(1:256))
hh <- hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 25 &
                                     mov.data.diff$time.diff > 0 & 
                                     mov.data.diff$time.diff < 512],
									 main = 'Individual 25',xlab = 'Time between fixes (h)',
           plot = T, breaks = breaks)

#X26  ###  Costa Rica 
unstacked <- split(move.data)
X26=unstacked$X26
X26_df <- as(X26, "data.frame")
head(X26_df)
nrow(X26_df)
range(X26_df$date)
m <- get_map(bbox(extent(X26)*1.1), source="google", zoom=12, maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=X26_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
#m <- get_map(bbox(extent(X26)*1.1), source="google", zoom=13, maptype="satellite") #define the box for map based on the trajectory coordinates
#ggmap(m)+geom_path(data=X26_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
#' ## Display space-time cube. 
open3d()
# To get a bigger window than the default
par3d(windowRect = c(100, 100, 612, 612))
Sys.sleep(0.1) # Allow sluggish window managers to catch up
with(X26_df, plot3d(x,y,date, type="l", col=as.numeric(X26_df$Event_ID)))
(stcube<-with(X26_df, plot3d(x,y,date, type="l",col=as.numeric(cut(X26_df$dist,3)), alpha=0.4)))
#Timelag
head(timeLag(X26, units="mins"))
head(timeLag(X26, units="hours"))
median(timeLag(X26, units="hours"))
max(timeLag(X26, units="hours"))
max(timeLag(X26, units="days"))
plot(X26_df$individual.local.identifier..ID.~X26_df $date,xlab = "Sampling period",ylab = "Individuals ID Number")
#write.table(X26,file="c:/RWorkDir/jaguardatapaper/ID/X26.txt",row.names = F,quote=F,col.names=T,sep="\t")
#Histograms and plots
hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 26],
     main = 'Individual 26', xlab = "Time between fixes (h)")
plot(density.default(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 26],
     na.rm = T), main = 'Individual 26')
	 
breaks <- c(0, 0.5, 1.5, 2.5, 4.5, 8.5, 16.5, 32.5, 64.5, 128.5, 256.5, 1000)
breaks <- c(0, 1, 2*(1:256))
hh <- hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 26 &
                                     mov.data.diff$time.diff > 0 & 
                                     mov.data.diff$time.diff < 512],
									 main = 'Individual 26',xlab = 'Time between fixes (h)',
           plot = T, breaks = breaks)

#X27
unstacked <- split(move.data)
X27=unstacked$X27
X27_df <- as(X27, "data.frame")
head(X27_df)
nrow(X27_df)
range(X27_df$date)
m <- get_map(bbox(extent(X27)*1.1), source="google", zoom=13, maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=X27_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
#' ## Display space-time cube. 
open3d()
# To get a bigger window than the default
par3d(windowRect = c(100, 100, 612, 612))
Sys.sleep(0.1) # Allow sluggish window managers to catch up
with(X27_df, plot3d(x,y,date, type="l", col=as.numeric(X27_df$Event_ID)))
(stcube<-with(X27_df, plot3d(x,y,date, type="l",col=as.numeric(cut(X27_df$dist,3)), alpha=0.4)))
#Timelag
head(timeLag(X27, units="mins"))
head(timeLag(X27, units="hours"))
median(timeLag(X27, units="hours"))
max(timeLag(X27, units="hours"))
max(timeLag(X27, units="days"))
plot(X27_df$individual.local.identifier..ID.~X27_df $date,xlab = "Sampling period",ylab = "Individuals ID Number")
#write.table(X27,file="c:/RWorkDir/jaguardatapaper/ID/X27.txt",row.names = F,quote=F,col.names=T,sep="\t")
#Histograms and plots
hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 27],
     main = 'Individual 27', xlab = "Time between fixes (h)")
plot(density.default(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 27],
     na.rm = T), main = 'Individual 27')
	 
breaks <- c(0, 0.5, 1.5, 2.5, 4.5, 8.5, 16.5, 32.5, 64.5, 128.5, 256.5, 1000)
breaks <- c(0, 1, 2*(1:256))
hh <- hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 27 &
                                     mov.data.diff$time.diff > 0 & 
                                     mov.data.diff$time.diff < 512],
									 main = 'Individual 27',xlab = 'Time between fixes (h)',
           plot = T, breaks = breaks)

#X28
unstacked <- split(move.data)
X28=unstacked$X28
X28_df <- as(X28, "data.frame")
head(X28_df)
nrow(X28_df)
range(X28_df$date)
m <- get_map(bbox(extent(X28)*1.1), source="google", zoom=13, maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=X28_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
#' ## Display space-time cube. 
open3d()
# To get a bigger window than the default
par3d(windowRect = c(100, 100, 612, 612))
Sys.sleep(0.1) # Allow sluggish window managers to catch up
with(X28_df, plot3d(x,y,date, type="l", col=as.numeric(X28_df$Event_ID)))
(stcube<-with(X28_df, plot3d(x,y,date, type="l",col=as.numeric(cut(X28_df$dist,3)), alpha=0.4)))
#Timelag
head(timeLag(X28, units="mins"))
head(timeLag(X28, units="hours"))
median(timeLag(X28, units="hours"))
max(timeLag(X28, units="hours"))
max(timeLag(X28, units="days"))
plot(X28_df$individual.local.identifier..ID.~X28_df $date,xlab = "Sampling period",ylab = "Individuals ID Number")
#write.table(X28,file="c:/RWorkDir/jaguardatapaper/ID/X28.txt",row.names = F,quote=F,col.names=T,sep="\t")
#Histograms and plots
hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 28],
     main = 'Individual 28', xlab = "Time between fixes (h)")
plot(density.default(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 28],
     na.rm = T), main = 'Individual 28')
	 
breaks <- c(0, 0.5, 1.5, 2.5, 4.5, 8.5, 16.5, 32.5, 64.5, 128.5, 256.5, 1000)
breaks <- c(0, 1, 2*(1:256))
hh <- hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 28 &
                                     mov.data.diff$time.diff > 0 & 
                                     mov.data.diff$time.diff < 512],
									 main = 'Individual 28',xlab = 'Time between fixes (h)',
           plot = T, breaks = breaks)

#X29
unstacked <- split(move.data)
X29=unstacked$X29
X29_df <- as(X29, "data.frame")
head(X29_df)
nrow(X29_df)
range(X29_df$date)
m <- get_map(bbox(extent(X29)*1.1), source="google", zoom=14, maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=X29_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
m <- get_map(bbox(extent(X29)*1.1), source="google", zoom=13, maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=X29_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
#' ## Display space-time cube. 
open3d()
# To get a bigger window than the default
par3d(windowRect = c(100, 100, 612, 612))
Sys.sleep(0.1) # Allow sluggish window managers to catch up
with(X29_df, plot3d(x,y,date, type="l", col=as.numeric(X29_df$Event_ID)))
(stcube<-with(X29_df, plot3d(x,y,date, type="l",col=as.numeric(cut(X29_df$dist,3)), alpha=0.4)))
#Timelag
head(timeLag(X29, units="mins"))
head(timeLag(X29, units="hours"))
median(timeLag(X29, units="hours"))
max(timeLag(X29, units="hours"))
max(timeLag(X29, units="days"))
plot(X29_df$individual.local.identifier..ID.~X29_df $date,xlab = "Sampling period",ylab = "Individuals ID Number")
#write.table(X29,file="c:/RWorkDir/jaguardatapaper/ID/X29.txt",row.names = F,quote=F,col.names=T,sep="\t")
#Histograms and plots
hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 29],
     main = 'Individual 29', xlab = "Time between fixes (h)")
plot(density.default(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 29],
     na.rm = T), main = 'Individual 29')
	 
breaks <- c(0, 0.5, 1.5, 2.5, 4.5, 8.5, 16.5, 32.5, 64.5, 128.5, 256.5, 1000)
breaks <- c(0, 1, 2*(1:256))
hh <- hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 29 &
                                     mov.data.diff$time.diff > 0 & 
                                     mov.data.diff$time.diff < 512],
									 main = 'Individual 29',xlab = 'Time between fixes (h)',
           plot = T, breaks = breaks)

#X30
unstacked <- split(move.data)
X30=unstacked$X30
X30_df <- as(X30, "data.frame")
head(X30_df)
nrow(X30_df)
range(X30_df$date)
m <- get_map(bbox(extent(X30)*1.1), source="google", zoom=13, maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=X30_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
m <- get_map(bbox(extent(X30)*1.1), source="google", zoom=12, maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=X30_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
#' ## Display space-time cube. 
open3d()
# To get a bigger window than the default
par3d(windowRect = c(100, 100, 612, 612))
Sys.sleep(0.1) # Allow sluggish window managers to catch up
with(X30_df, plot3d(x,y,date, type="l", col=as.numeric(X30_df$Event_ID)))
(stcube<-with(X30_df, plot3d(x,y,date, type="l",col=as.numeric(cut(X30_df$dist,3)), alpha=0.4)))
#Timelag
head(timeLag(X30, units="mins"))
head(timeLag(X30, units="hours"))
median(timeLag(X30, units="hours"))
max(timeLag(X30, units="hours"))
max(timeLag(X30, units="days"))
plot(X30_df$individual.local.identifier..ID.~X30_df $date,xlab = "Sampling period",ylab = "Individuals ID Number")
##write.table(X30,file="c:/RWorkDir/jaguardatapaper/ID/X30.txt",row.names = F,quote=F,col.names=T,sep="\t")
#Histograms and plots
hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 30],
     main = 'Individual 30', xlab = "Time between fixes (h)")
plot(density.default(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 30],
     na.rm = T), main = 'Individual 30')
	 
breaks <- c(0, 0.5, 1.5, 2.5, 4.5, 8.5, 16.5, 32.5, 64.5, 128.5, 256.5, 1000)
breaks <- c(0, 1, 2*(1:256))
hh <- hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 30 &
                                     mov.data.diff$time.diff > 0 & 
                                     mov.data.diff$time.diff < 512],
									 main = 'Individual 30',xlab = 'Time between fixes (h)',
           plot = T, breaks = breaks)

#X31
unstacked <- split(move.data)
X31=unstacked$X31
X31_df <- as(X31, "data.frame")
head(X31_df)
nrow(X31_df)
range(X31_df$date)
m <- get_map(bbox(extent(X31)*1.1), source="google", zoom=13, maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=X31_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
#' ## Display space-time cube. 
open3d()
# To get a bigger window than the default
par3d(windowRect = c(100, 100, 612, 612))
Sys.sleep(0.1) # Allow sluggish window managers to catch up
with(X31_df, plot3d(x,y,date, type="l", col=as.numeric(X31_df$Event_ID)))
(stcube<-with(X31_df, plot3d(x,y,date, type="l",col=as.numeric(cut(X31_df$dist,3)), alpha=0.4)))
#Timelag
head(timeLag(X31, units="mins"))
head(timeLag(X31, units="hours"))
median(timeLag(X31, units="hours"))
max(timeLag(X31, units="hours"))
max(timeLag(X31, units="days"))
plot(X31_df$individual.local.identifier..ID.~X31_df $date,xlab = "Sampling period",ylab = "Individuals ID Number")
#write.table(X31,file="c:/RWorkDir/jaguardatapaper/ID/X31.txt",row.names = F,quote=F,col.names=T,sep="\t")
#Histograms and plots
hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 31],
     main = 'Individual 31', xlab = "Time between fixes (h)")
plot(density.default(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 31],
     na.rm = T), main = 'Individual 31')
	 
breaks <- c(0, 0.5, 1.5, 2.5, 4.5, 8.5, 16.5, 32.5, 64.5, 128.5, 256.5, 1000)
breaks <- c(0, 1, 2*(1:256))
hh <- hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 31 &
                                     mov.data.diff$time.diff > 0 & 
                                     mov.data.diff$time.diff < 512],
									 main = 'Individual 31',xlab = 'Time between fixes (h)',
           plot = T, breaks = breaks)

#X32
unstacked <- split(move.data)
X32=unstacked$X32
X32_df <- as(X32, "data.frame")
head(X32_df)
nrow(X32_df)
range(X32_df$date)
m <- get_map(bbox(extent(X32)*1.1), source="google", zoom=12, maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=X32_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
#' ## Display space-time cube. 
open3d()
# To get a bigger window than the default
par3d(windowRect = c(100, 100, 612, 612))
Sys.sleep(0.1) # Allow sluggish window managers to catch up
with(X32_df, plot3d(x,y,date, type="l", col=as.numeric(X32_df$Event_ID)))
(stcube<-with(X32_df, plot3d(x,y,date, type="l",col=as.numeric(cut(X32_df$dist,3)), alpha=0.4)))
#Timelag
head(timeLag(X32, units="mins"))
head(timeLag(X32, units="hours"))
median(timeLag(X32, units="hours"))
max(timeLag(X32, units="hours"))
max(timeLag(X32, units="days"))
plot(X32_df$individual.local.identifier..ID.~X32_df $date,xlab = "Sampling period",ylab = "Individuals ID Number")
#write.table(X32,file="c:/RWorkDir/jaguardatapaper/ID/X32.txt",row.names = F,quote=F,col.names=T,sep="\t")
#Histograms and plots
hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 32],
     main = 'Individual 32', xlab = "Time between fixes (h)")
plot(density.default(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 32],
     na.rm = T), main = 'Individual 32')
	 
breaks <- c(0, 0.5, 1.5, 2.5, 4.5, 8.5, 16.5, 32.5, 64.5, 128.5, 256.5, 1000)
breaks <- c(0, 1, 2*(1:256))
hh <- hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 32 &
                                     mov.data.diff$time.diff > 0 & 
                                     mov.data.diff$time.diff < 512],
									 main = 'Individual 32',xlab = 'Time between fixes (h)',
           plot = T, breaks = breaks)

#X33
unstacked <- split(move.data)
X33=unstacked$X33
X33_df <- as(X33, "data.frame")
head(X33_df)
nrow(X33_df)
range(X33_df$date)
m <- get_map(bbox(extent(X33)*1.1), source="google", zoom=12, maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=X33_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
m <- get_map(bbox(extent(X33)*1.1), source="google", zoom=13, maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=X33_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
#' ## Display space-time cube. 
open3d()
# To get a bigger window than the default
par3d(windowRect = c(100, 100, 612, 612))
Sys.sleep(0.1) # Allow sluggish window managers to catch up
with(X33_df, plot3d(x,y,date, type="l", col=as.numeric(X33_df$Event_ID)))
(stcube<-with(X33_df, plot3d(x,y,date, type="l",col=as.numeric(cut(X33_df$dist,3)), alpha=0.4)))
#Timelag
head(timeLag(X33, units="mins"))
head(timeLag(X33, units="hours"))
median(timeLag(X33, units="hours"))
max(timeLag(X33, units="hours"))
max(timeLag(X33, units="days"))
plot(X33_df$individual.local.identifier..ID.~X33_df $date,xlab = "Sampling period",ylab = "Individuals ID Number")
#write.table(X33,file="c:/RWorkDir/jaguardatapaper/ID/X33.txt",row.names = F,quote=F,col.names=T,sep="\t")
#Histograms and plots
hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 33],
     main = 'Individual 33', xlab = "Time between fixes (h)")
plot(density.default(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 33],
     na.rm = T), main = 'Individual 33')
	 
breaks <- c(0, 0.5, 1.5, 2.5, 4.5, 8.5, 16.5, 32.5, 64.5, 128.5, 256.5, 1000)
breaks <- c(0, 1, 2*(1:256))
hh <- hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 33 &
                                     mov.data.diff$time.diff > 0 & 
                                     mov.data.diff$time.diff < 512],
									 main = 'Individual 33',xlab = 'Time between fixes (h)',
           plot = T, breaks = breaks)

#X34
X34=unstacked$X34
X34_df <- as(X34, "data.frame")
head(X34_df)
nrow(X34_df)
range(X34_df$date)
range(X34_df$date)
m <- get_map(bbox(extent(X34)*1.1), source="google", zoom=12, maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=X34_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
#' ## Display space-time cube. 
open3d()
# To get a bigger window than the default
par3d(windowRect = c(100, 100, 612, 612))
Sys.sleep(0.1) # Allow sluggish window managers to catch up
with(X34_df, plot3d(x,y,date, type="l", col=as.numeric(X34_df$Event_ID)))
(stcube<-with(X34_df, plot3d(x,y,date, type="l",col=as.numeric(cut(X34_df$dist,3)), alpha=0.4)))
#Timelag
head(timeLag(X34, units="mins"))
head(timeLag(X34, units="hours"))
median(timeLag(X34, units="hours"))
max(timeLag(X34, units="hours"))
max(timeLag(X34, units="days"))
plot(X34_df$individual.local.identifier..ID.~X34_df $date,xlab = "Sampling period",ylab = "Individuals ID Number")
#write.table(X34,file="c:/RWorkDir/jaguardatapaper/ID/X34.txt",row.names = F,quote=F,col.names=T,sep="\t")
#Histograms and plots
hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 34],
     main = 'Individual 34', xlab = "Time between fixes (h)")
plot(density.default(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 34],
     na.rm = T), main = 'Individual 34')
	 
breaks <- c(0, 0.5, 1.5, 2.5, 4.5, 8.5, 16.5, 32.5, 64.5, 128.5, 256.5, 1000)
breaks <- c(0, 1, 2*(1:256))
hh <- hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 34 &
                                     mov.data.diff$time.diff > 0 & 
                                     mov.data.diff$time.diff < 512],
									 main = 'Individual 34',xlab = 'Time between fixes (h)',
           plot = T, breaks = breaks)

#X35
X35=unstacked$X35
X35_df <- as(X35, "data.frame")
head(X35_df)
nrow(X35_df)
range(X35_df$date)
m <- get_map(bbox(extent(X35)*1.1), source="google", zoom=12, maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=X35_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
m <- get_map(bbox(extent(X35)*1.1), source="google", zoom=11, maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=X35_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
#' ## Display space-time cube. 
open3d()
# To get a bigger window than the default
par3d(windowRect = c(100, 100, 612, 612))
Sys.sleep(0.1) # Allow sluggish window managers to catch up
with(X35_df, plot3d(x,y,date, type="l", col=as.numeric(X35_df$Event_ID)))
(stcube<-with(X35_df, plot3d(x,y,date, type="l",col=as.numeric(cut(X35_df$dist,3)), alpha=0.4)))
#Timelag
head(timeLag(X35, units="mins"))
head(timeLag(X35, units="hours"))
median(timeLag(X35, units="hours"))
max(timeLag(X35, units="hours"))
max(timeLag(X35, units="days"))
plot(X35_df$individual.local.identifier..ID.~X35_df $date,xlab = "Sampling period",ylab = "Individuals ID Number")
#write.table(X35,file="c:/RWorkDir/jaguardatapaper/ID/X35.txt",row.names = F,quote=F,col.names=T,sep="\t")
#Histograms and plots
hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 35],
     main = 'Individual 35', xlab = "Time between fixes (h)")
plot(density.default(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 35],
     na.rm = T), main = 'Individual 35')
	 
breaks <- c(0, 0.5, 1.5, 2.5, 4.5, 8.5, 16.5, 32.5, 64.5, 128.5, 256.5, 1000)
breaks <- c(0, 1, 2*(1:256))
hh <- hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 35 &
                                     mov.data.diff$time.diff > 0 & 
                                     mov.data.diff$time.diff < 512],
									 main = 'Individual 35',xlab = 'Time between fixes (h)',
           plot = T, breaks = breaks)

#X36
unstacked <- split(move.data)
X36=unstacked$X36
X36_df <- as(X36, "data.frame")
head(X36_df)
nrow(X36_df)
range(X36_df$date)
m <- get_map(bbox(extent(X36)*1.1), source="google", zoom=12, maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=X36_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
m <- get_map(bbox(extent(X36)*1.1), source="google", zoom=11, maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=X36_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
#' ## Display space-time cube. 
open3d()
# To get a bigger window than the default
par3d(windowRect = c(100, 100, 612, 612))
Sys.sleep(0.1) # Allow sluggish window managers to catch up
with(X36_df, plot3d(x,y,date, type="l", col=as.numeric(X36_df$Event_ID)))
(stcube<-with(X36_df, plot3d(x,y,date, type="l",col=as.numeric(cut(X36_df$dist,3)), alpha=0.4)))
#Timelag
head(timeLag(X36, units="mins"))
head(timeLag(X36, units="hours"))
median(timeLag(X36, units="hours"))
max(timeLag(X36, units="hours"))
max(timeLag(X36, units="days"))
plot(X36_df$individual.local.identifier..ID.~X36_df $date,xlab = "Sampling period",ylab = "Individuals ID Number")
#write.table(X36,file="c:/RWorkDir/jaguardatapaper/ID/X36.txt",row.names = F,quote=F,col.names=T,sep="\t")
#Histograms and plots
hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 36],
     main = 'Individual 36', xlab = "Time between fixes (h)")
plot(density.default(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 36],
     na.rm = T), main = 'Individual 36')
	 
breaks <- c(0, 0.5, 1.5, 2.5, 4.5, 8.5, 16.5, 32.5, 64.5, 128.5, 256.5, 1000)
breaks <- c(0, 1, 2*(1:256))
hh <- hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 36 &
                                     mov.data.diff$time.diff > 0 & 
                                     mov.data.diff$time.diff < 512],
									 main = 'Individual 36',xlab = 'Time between fixes (h)',
           plot = T, breaks = breaks)

#X37
unstacked <- split(move.data)
X37=unstacked$X37
X37_df <- as(X37, "data.frame")
head(X37_df)
nrow(X37_df)
range(X37_df$date)
m <- get_map(bbox(extent(X37)*1.1), source="google", zoom=12, maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=X37_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
#' ## Display space-time cube. 
open3d()
# To get a bigger window than the default
par3d(windowRect = c(100, 100, 612, 612))
Sys.sleep(0.1) # Allow sluggish window managers to catch up
with(X37_df, plot3d(x,y,date, type="l", col=as.numeric(X37_df$Event_ID)))
(stcube<-with(X37_df, plot3d(x,y,date, type="l",col=as.numeric(cut(X37_df$dist,3)), alpha=0.4)))
#Timelag
head(timeLag(X37, units="mins"))
head(timeLag(X37, units="hours"))
median(timeLag(X37, units="hours"))
max(timeLag(X37, units="hours"))
max(timeLag(X37, units="days"))
plot(X37_df$individual.local.identifier..ID.~X37_df $date,xlab = "Sampling period",ylab = "Individuals ID Number")
#write.table(X37,file="c:/RWorkDir/jaguardatapaper/ID/X37.txt",row.names = F,quote=F,col.names=T,sep="\t")
#Histograms and plots
hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 37],
     main = 'Individual 37', xlab = "Time between fixes (h)")
plot(density.default(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 37],
     na.rm = T), main = 'Individual 37')
	 
breaks <- c(0, 0.5, 1.5, 2.5, 4.5, 8.5, 16.5, 32.5, 64.5, 128.5, 256.5, 1000)
breaks <- c(0, 1, 2*(1:256))
hh <- hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 37 &
                                     mov.data.diff$time.diff > 0 & 
                                     mov.data.diff$time.diff < 512],
									 main = 'Individual 37',xlab = 'Time between fixes (h)',
           plot = T, breaks = breaks)

#X38
unstacked <- split(move.data)
X38=unstacked$X38
X38_df <- as(X38, "data.frame")
head(X38_df)
nrow(X38_df)
range(X38_df$date)
m <- get_map(bbox(extent(X38)*1.1), source="google", zoom=12, maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=X38_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
#' ## Display space-time cube. 
open3d()
# To get a bigger window than the default
par3d(windowRect = c(100, 100, 612, 612))
Sys.sleep(0.1) # Allow sluggish window managers to catch up
with(X38_df, plot3d(x,y,date, type="l", col=as.numeric(X38_df$Event_ID)))
(stcube<-with(X38_df, plot3d(x,y,date, type="l",col=as.numeric(cut(X38_df$dist,3)), alpha=0.4)))
#Timelag
head(timeLag(X38, units="mins"))
head(timeLag(X38, units="hours"))
median(timeLag(X38, units="hours"))
max(timeLag(X38, units="hours"))
max(timeLag(X38, units="days"))
plot(X38_df$individual.local.identifier..ID.~X38_df $date,xlab = "Sampling period",ylab = "Individuals ID Number")
#write.table(X38,file="c:/RWorkDir/jaguardatapaper/ID/X38.txt",row.names = F,quote=F,col.names=T,sep="\t")
#Histograms and plots
hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 38],
     main = 'Individual 38', xlab = "Time between fixes (h)")
plot(density.default(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 38],
     na.rm = T), main = 'Individual 38')
	 
breaks <- c(0, 0.5, 1.5, 2.5, 4.5, 8.5, 16.5, 32.5, 64.5, 128.5, 256.5, 1000)
breaks <- c(0, 1, 2*(1:256))
hh <- hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 38 &
                                     mov.data.diff$time.diff > 0 & 
                                     mov.data.diff$time.diff < 512],
									 main = 'Individual 38',xlab = 'Time between fixes (h)',
           plot = T, breaks = breaks)

#X39
unstacked <- split(move.data)
X39=unstacked$X39
X39_df <- as(X39, "data.frame")
head(X39_df)
nrow(X39_df)
range(X39_df$date)
m <- get_map(bbox(extent(X39)*1.1), source="google", zoom=12, maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=X39_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
m <- get_map(bbox(extent(X39)*1.1), source="google", zoom=11, maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=X39_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
#' ## Display space-time cube. 
open3d()
# To get a bigger window than the default
par3d(windowRect = c(100, 100, 612, 612))
Sys.sleep(0.1) # Allow sluggish window managers to catch up
with(X39_df, plot3d(x,y,date, type="l", col=as.numeric(X39_df$Event_ID)))
(stcube<-with(X39_df, plot3d(x,y,date, type="l",col=as.numeric(cut(X39_df$dist,3)), alpha=0.4)))
#Timelag
head(timeLag(X39, units="mins"))
head(timeLag(X39, units="hours"))
median(timeLag(X39, units="hours"))
max(timeLag(X39, units="hours"))
max(timeLag(X39, units="days"))
plot(X39_df$individual.local.identifier..ID.~X39_df $date,xlab = "Sampling period",ylab = "Individuals ID Number")
#write.table(X39,file="c:/RWorkDir/jaguardatapaper/ID/X39.txt",row.names = F,quote=F,col.names=T,sep="\t")
#Histograms and plots
hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 39],
     main = 'Individual 39', xlab = "Time between fixes (h)")
plot(density.default(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 39],
     na.rm = T), main = 'Individual 39')
	 
breaks <- c(0, 0.5, 1.5, 2.5, 4.5, 8.5, 16.5, 32.5, 64.5, 128.5, 256.5, 1000)
breaks <- c(0, 1, 2*(1:256))
hh <- hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 39 &
                                     mov.data.diff$time.diff > 0 & 
                                     mov.data.diff$time.diff < 512],
									 main = 'Individual 39',xlab = 'Time between fixes (h)',
           plot = T, breaks = breaks)

#X40
unstacked <- split(move.data)
X40=unstacked$X40
X40_df <- as(X40, "data.frame")
head(X40_df)
nrow(X40_df)
range(X40_df$date)
m <- get_map(bbox(extent(X40)*1.1), source="google", zoom=11, maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=X40_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
#' ## Display space-time cube. 
open3d()
# To get a bigger window than the default
par3d(windowRect = c(100, 100, 612, 612))
Sys.sleep(0.1) # Allow sluggish window managers to catch up
with(X40_df, plot3d(x,y,date, type="l", col=as.numeric(X40_df$Event_ID)))
(stcube<-with(X40_df, plot3d(x,y,date, type="l",col=as.numeric(cut(X40_df$dist,3)), alpha=0.4)))
#Timelag
head(timeLag(X40, units="mins"))
head(timeLag(X40, units="hours"))
median(timeLag(X40, units="hours"))
max(timeLag(X40, units="hours"))
max(timeLag(X40, units="days"))
plot(X40_df$individual.local.identifier..ID.~X40_df $date,xlab = "Sampling period",ylab = "Individuals ID Number")
#write.table(X40,file="c:/RWorkDir/jaguardatapaper/ID/X40.txt",row.names = F,quote=F,col.names=T,sep="\t")
#Histograms and plots
hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 40],
     main = 'Individual 40', xlab = "Time between fixes (h)")
plot(density.default(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 40],
     na.rm = T), main = 'Individual 40')
	 
breaks <- c(0, 0.5, 1.5, 2.5, 4.5, 8.5, 16.5, 32.5, 64.5, 128.5, 256.5, 1000)
breaks <- c(0, 1, 2*(1:256))
hh <- hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 40 &
                                     mov.data.diff$time.diff > 0 & 
                                     mov.data.diff$time.diff < 512],
									 main = 'Individual 40',xlab = 'Time between fixes (h)',
           plot = T, breaks = breaks)
#X41
unstacked <- split(move.data)
X41=unstacked$X41
X41_df <- as(X41, "data.frame")
head(X41_df)
nrow(X41_df)
range(X41_df$date)
m <- get_map(bbox(extent(X41)*1.1), source="google", zoom=13, maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=X41_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
m <- get_map(bbox(extent(X41)*1.1), source="google", zoom=14, maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=X41_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
#' ## Display space-time cube. 
open3d()
# To get a bigger window than the default
par3d(windowRect = c(100, 100, 612, 612))
Sys.sleep(0.1) # Allow sluggish window managers to catch up
with(X41_df, plot3d(x,y,date, type="l", col=as.numeric(X41_df$Event_ID)))
(stcube<-with(X41_df, plot3d(x,y,date, type="l",col=as.numeric(cut(X41_df$dist,3)), alpha=0.4)))
#Timelag
head(timeLag(X41, units="mins"))
head(timeLag(X41, units="hours"))
median(timeLag(X41, units="hours"))
max(timeLag(X41, units="hours"))
max(timeLag(X41, units="days"))
plot(X41_df$individual.local.identifier..ID.~X41_df $date,xlab = "Sampling period",ylab = "Individuals ID Number")
#write.table(X41,file="c:/RWorkDir/jaguardatapaper/ID/X41.txt",row.names = F,quote=F,col.names=T,sep="\t")
#Histograms and plots
hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 41],
     main = 'Individual 41', xlab = "Time between fixes (h)")
plot(density.default(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 41],
     na.rm = T), main = 'Individual 41')
	 
breaks <- c(0, 0.5, 1.5, 2.5, 4.5, 8.5, 16.5, 32.5, 64.5, 128.5, 256.5, 1000)
breaks <- c(0, 1, 2*(1:256))
hh <- hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 41 &
                                     mov.data.diff$time.diff > 0 & 
                                     mov.data.diff$time.diff < 512],
									 main = 'Individual 41',xlab = 'Time between fixes (h)',
           plot = T, breaks = breaks)


#X42
unstacked <- split(move.data)
X42=unstacked$X42
X42_df <- as(X42, "data.frame")
head(X42_df)
nrow(X42_df)
range(X42_df$date)
m <- get_map(bbox(extent(X42)*1.1), source="google", zoom=11, maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=X42_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
#' ## Display space-time cube. 
open3d()
# To get a bigger window than the default
par3d(windowRect = c(100, 100, 612, 612))
Sys.sleep(0.1) # Allow sluggish window managers to catch up
with(X42_df, plot3d(x,y,date, type="l", col=as.numeric(X42_df$Event_ID)))
(stcube<-with(X42_df, plot3d(x,y,date, type="l",col=as.numeric(cut(X42_df$dist,3)), alpha=0.4)))
#Timelag
head(timeLag(X42, units="mins"))
head(timeLag(X42, units="hours"))
median(timeLag(X42, units="hours"))
max(timeLag(X42, units="hours"))
max(timeLag(X42, units="days"))
plot(X42_df$individual.local.identifier..ID.~X42_df $date,xlab = "Sampling period",ylab = "Individuals ID Number")
#write.table(X42,file="c:/RWorkDir/jaguardatapaper/ID/X42.txt",row.names = F,quote=F,col.names=T,sep="\t")
#Histograms and plots
hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 42],
     main = 'Individual 42', xlab = "Time between fixes (h)")
plot(density.default(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 42],
     na.rm = T), main = 'Individual 42')
	 
breaks <- c(0, 0.5, 1.5, 2.5, 4.5, 8.5, 16.5, 32.5, 64.5, 128.5, 256.5, 1000)
breaks <- c(0, 1, 2*(1:256))
hh <- hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 42 &
                                     mov.data.diff$time.diff > 0 & 
                                     mov.data.diff$time.diff < 512],
									 main = 'Individual 42',xlab = 'Time between fixes (h)',
           plot = T, breaks = breaks)

#X43
X43=unstacked$X43
X43_df <- as(X43, "data.frame")
head(X43_df)
nrow(X43_df)
range(X43_df$date)
m <- get_map(bbox(extent(X43)*1.1), source="google", zoom=11, maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=X43_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
#' ## Display space-time cube. 
open3d()
# To get a bigger window than the default
par3d(windowRect = c(100, 100, 612, 612))
Sys.sleep(0.1) # Allow sluggish window managers to catch up
with(X43_df, plot3d(x,y,date, type="l", col=as.numeric(X43_df$Event_ID)))
(stcube<-with(X43_df, plot3d(x,y,date, type="l",col=as.numeric(cut(X43_df$dist,3)), alpha=0.4)))
#Timelag
head(timeLag(X43, units="mins"))
head(timeLag(X43, units="hours"))
median(timeLag(X43, units="hours"))
max(timeLag(X43, units="hours"))
max(timeLag(X43, units="days"))
plot(X43_df$individual.local.identifier..ID.~X43_df $date,xlab = "Sampling period",ylab = "Individuals ID Number")
#write.table(X43,file="c:/RWorkDir/jaguardatapaper/ID/X43.txt",row.names = F,quote=F,col.names=T,sep="\t")
#Histograms and plots
hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 43],
     main = 'Individual 43', xlab = "Time between fixes (h)")
plot(density.default(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 43],
     na.rm = T), main = 'Individual 43')
	 
breaks <- c(0, 0.5, 1.5, 2.5, 4.5, 8.5, 16.5, 32.5, 64.5, 128.5, 256.5, 1000)
breaks <- c(0, 1, 2*(1:256))
hh <- hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 43 &
                                     mov.data.diff$time.diff > 0 & 
                                     mov.data.diff$time.diff < 512],
									 main = 'Individual 43',xlab = 'Time between fixes (h)',
           plot = T, breaks = breaks)


#X44
X44=unstacked$X44
X44_df <- as(X44, "data.frame")
head(X44_df)
nrow(X44_df)
range(X44_df$date)
m <- get_map(bbox(extent(X44)*1.1), source="google", zoom=11, maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=X44_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
#' ## Display space-time cube. 
open3d()
# To get a bigger window than the default
par3d(windowRect = c(100, 100, 612, 612))
Sys.sleep(0.1) # Allow sluggish window managers to catch up
with(X44_df, plot3d(x,y,date, type="l", col=as.numeric(X44_df$Event_ID)))
(stcube<-with(X44_df, plot3d(x,y,date, type="l",col=as.numeric(cut(X44_df$dist,3)), alpha=0.4)))
#Timelag
head(timeLag(X44, units="mins"))
head(timeLag(X44, units="hours"))
median(timeLag(X44, units="hours"))
max(timeLag(X44, units="hours"))
max(timeLag(X44, units="days"))
plot(X44_df$individual.local.identifier..ID.~X44_df $date,xlab = "Sampling period",ylab = "Individuals ID Number")
#write.table(X44,file="c:/RWorkDir/jaguardatapaper/ID/X44.txt",row.names = F,quote=F,col.names=T,sep="\t")
#Histograms and plots
hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 44],
     main = 'Individual 44', xlab = "Time between fixes (h)")
plot(density.default(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 44],
     na.rm = T), main = 'Individual 44')
	 
breaks <- c(0, 0.5, 1.5, 2.5, 4.5, 8.5, 16.5, 32.5, 64.5, 128.5, 256.5, 1000)
breaks <- c(0, 1, 2*(1:256))
hh <- hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 44 &
                                     mov.data.diff$time.diff > 0 & 
                                     mov.data.diff$time.diff < 512],
									 main = 'Individual 44',xlab = 'Time between fixes (h)',
           plot = T, breaks = breaks)

#X45
X45=unstacked$X45
X45_df <- as(X45, "data.frame")
head(X45_df)
nrow(X45_df)
range(X45_df$date)
m <- get_map(bbox(extent(X45)*1.1), source="google", zoom=14, maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=X45_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
#' ## Display space-time cube. 
open3d()
# To get a bigger window than the default
par3d(windowRect = c(100, 100, 612, 612))
Sys.sleep(0.1) # Allow sluggish window managers to catch up
with(X45_df, plot3d(x,y,date, type="l", col=as.numeric(X45_df$Event_ID)))
(stcube<-with(X45_df, plot3d(x,y,date, type="l",col=as.numeric(cut(X45_df$dist,3)), alpha=0.4)))
#Timelag
head(timeLag(X45, units="mins"))
head(timeLag(X45, units="hours"))
median(timeLag(X45, units="hours"))
max(timeLag(X45, units="hours"))
max(timeLag(X45, units="days"))
plot(X45_df$individual.local.identifier..ID.~X45_df $date,xlab = "Sampling period",ylab = "Individuals ID Number")
#write.table(X45,file="c:/RWorkDir/jaguardatapaper/ID/X45.txt",row.names = F,quote=F,col.names=T,sep="\t")
#Histograms and plots
hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 45],
     main = 'Individual 45', xlab = "Time between fixes (h)")
plot(density.default(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 45],
     na.rm = T), main = 'Individual 45')
	 
breaks <- c(0, 0.5, 1.5, 2.5, 4.5, 8.5, 16.5, 32.5, 64.5, 128.5, 256.5, 1000)
breaks <- c(0, 1, 2*(1:256))
hh <- hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 45 &
                                     mov.data.diff$time.diff > 0 & 
                                     mov.data.diff$time.diff < 512],
									 main = 'Individual 45',xlab = 'Time between fixes (h)',
           plot = T, breaks = breaks)

#X46
X46=unstacked$X46
X46_df <- as(X46, "data.frame")
head(X46_df)
nrow(X46_df)
range(X46_df$date)
m <- get_map(bbox(extent(X46)*1.1), source="google", zoom=12, maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=X46_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
m <- get_map(bbox(extent(X46)*1.1), source="google", zoom=11, maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=X46_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
#' ## Display space-time cube. 
open3d()
# To get a bigger window than the default
par3d(windowRect = c(100, 100, 612, 612))
Sys.sleep(0.1) # Allow sluggish window managers to catch up
with(X46_df, plot3d(x,y,date, type="l", col=as.numeric(X46_df$Event_ID)))
(stcube<-with(X46_df, plot3d(x,y,date, type="l",col=as.numeric(cut(X46_df$dist,3)), alpha=0.4)))
#Timelag
head(timeLag(X46, units="mins"))
head(timeLag(X46, units="hours"))
median(timeLag(X46, units="hours"))
max(timeLag(X46, units="hours"))
max(timeLag(X46, units="days"))
plot(X46_df$individual.local.identifier..ID.~X46_df $date,xlab = "Sampling period",ylab = "Individuals ID Number")
#write.table(X46,file="c:/RWorkDir/jaguardatapaper/ID/X46.txt",row.names = F,quote=F,col.names=T,sep="\t")
#Histograms and plots
hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 46],
     main = 'Individual 46', xlab = "Time between fixes (h)")
plot(density.default(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 46],
     na.rm = T), main = 'Individual 46')
	 
breaks <- c(0, 0.5, 1.5, 2.5, 4.5, 8.5, 16.5, 32.5, 64.5, 128.5, 256.5, 1000)
breaks <- c(0, 1, 2*(1:256))
hh <- hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 46 &
                                     mov.data.diff$time.diff > 0 & 
                                     mov.data.diff$time.diff < 512],
									 main = 'Individual 46',xlab = 'Time between fixes (h)',
           plot = T, breaks = breaks)

#X47
X47=unstacked$X47
X47_df <- as(X47, "data.frame")
head(X47_df)
nrow(X47_df)
range(X47_df$date)
m <- get_map(bbox(extent(X47)*1.1), source="google", zoom=12, maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=X47_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
m <- get_map(bbox(extent(X47)*1.1), source="google", zoom=11, maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=X47_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
#' ## Display space-time cube. 
open3d()
# To get a bigger window than the default
par3d(windowRect = c(100, 100, 612, 612))
Sys.sleep(0.1) # Allow sluggish window managers to catch up
with(X47_df, plot3d(x,y,date, type="l", col=as.numeric(X47_df$Event_ID)))
(stcube<-with(X47_df, plot3d(x,y,date, type="l",col=as.numeric(cut(X47_df$dist,3)), alpha=0.4)))
#Timelag
head(timeLag(X47, units="mins"))
head(timeLag(X47, units="hours"))
median(timeLag(X47, units="hours"))
max(timeLag(X47, units="hours"))
max(timeLag(X47, units="days"))
plot(X47_df$individual.local.identifier..ID.~X47_df $date,xlab = "Sampling period",ylab = "Individuals ID Number")
#write.table(X47,file="c:/RWorkDir/jaguardatapaper/ID/X47.txt",row.names = F,quote=F,col.names=T,sep="\t")
#Histograms and plots
hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 47],
     main = 'Individual 47', xlab = "Time between fixes (h)")
plot(density.default(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 47],
     na.rm = T), main = 'Individual 47')
	 
breaks <- c(0, 0.5, 1.5, 2.5, 4.5, 8.5, 16.5, 32.5, 64.5, 128.5, 256.5, 1000)
breaks <- c(0, 1, 2*(1:256))
hh <- hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 47 &
                                     mov.data.diff$time.diff > 0 & 
                                     mov.data.diff$time.diff < 512],
									 main = 'Individual 47',xlab = 'Time between fixes (h)',
           plot = T, breaks = breaks)


#X48
X48=unstacked$X48
X48_df <- as(X48, "data.frame")
head(X48_df)
nrow(X48_df)
range(X48_df$date)
m <- get_map(bbox(extent(X48)*1.1), source="google", zoom=11, maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=X48_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
#' ## Display space-time cube. 
open3d()
# To get a bigger window than the default
par3d(windowRect = c(100, 100, 612, 612))
Sys.sleep(0.1) # Allow sluggish window managers to catch up
with(X48_df, plot3d(x,y,date, type="l", col=as.numeric(X48_df$Event_ID)))
(stcube<-with(X48_df, plot3d(x,y,date, type="l",col=as.numeric(cut(X48_df$dist,3)), alpha=0.4)))
#Timelag
head(timeLag(X48, units="mins"))
head(timeLag(X48, units="hours"))
median(timeLag(X48, units="hours"))
max(timeLag(X48, units="hours"))
max(timeLag(X48, units="days"))
plot(X48_df$individual.local.identifier..ID.~X48_df $date,xlab = "Sampling period",ylab = "Individuals ID Number")
#write.table(X48,file="c:/RWorkDir/jaguardatapaper/ID/X48.txt",row.names = F,quote=F,col.names=T,sep="\t")
#Histograms and plots
hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 48],
     main = 'Individual 48', xlab = "Time between fixes (h)")
plot(density.default(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 48],
     na.rm = T), main = 'Individual 48')
	 
breaks <- c(0, 0.5, 1.5, 2.5, 4.5, 8.5, 16.5, 32.5, 64.5, 128.5, 256.5, 1000)
breaks <- c(0, 1, 2*(1:256))
hh <- hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 48 &
                                     mov.data.diff$time.diff > 0 & 
                                     mov.data.diff$time.diff < 512],
									 main = 'Individual 48',xlab = 'Time between fixes (h)',
           plot = T, breaks = breaks)


#X49
X49=unstacked$X49
X49_df <- as(X49, "data.frame")
head(X49_df)
nrow(X49_df)
range(X49_df$date)
m <- get_map(bbox(extent(X49)*1.1), source="google", zoom=11, maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=X49_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
#' ## Display space-time cube. 
open3d()
# To get a bigger window than the default
par3d(windowRect = c(100, 100, 612, 612))
Sys.sleep(0.1) # Allow sluggish window managers to catch up
with(X49_df, plot3d(x,y,date, type="l", col=as.numeric(X49_df$Event_ID)))
(stcube<-with(X49_df, plot3d(x,y,date, type="l",col=as.numeric(cut(X49_df$dist,3)), alpha=0.4)))
#Timelag
head(timeLag(X49, units="mins"))
head(timeLag(X49, units="hours"))
median(timeLag(X49, units="hours"))
max(timeLag(X49, units="hours"))
max(timeLag(X49, units="days"))
plot(X49_df$individual.local.identifier..ID.~X49_df $date,xlab = "Sampling period",ylab = "Individuals ID Number")
#write.table(X49,file="c:/RWorkDir/jaguardatapaper/ID/X49.txt",row.names = F,quote=F,col.names=T,sep="\t")
#Histograms and plots
hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 49],
     main = 'Individual 49', xlab = "Time between fixes (h)")
plot(density.default(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 49],
     na.rm = T), main = 'Individual 49')
	 
breaks <- c(0, 0.5, 1.5, 2.5, 4.5, 8.5, 16.5, 32.5, 64.5, 128.5, 256.5, 1000)
breaks <- c(0, 1, 2*(1:256))
hh <- hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 49 &
                                     mov.data.diff$time.diff > 0 & 
                                     mov.data.diff$time.diff < 512],
									 main = 'Individual 49',xlab = 'Time between fixes (h)',
           plot = T, breaks = breaks)


#X50
X50=unstacked$X50
X50_df <- as(X50, "data.frame")
head(X50_df)
nrow(X50_df)
range(X50_df$date)
m <- get_map(bbox(extent(X50)*1.1), source="google", zoom=12, maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=X50_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
#' ## Display space-time cube. 
open3d()
# To get a bigger window than the default
par3d(windowRect = c(100, 100, 612, 612))
Sys.sleep(0.1) # Allow sluggish window managers to catch up
with(X50_df, plot3d(x,y,date, type="l", col=as.numeric(X50_df$Event_ID)))
(stcube<-with(X50_df, plot3d(x,y,date, type="l",col=as.numeric(cut(X50_df$dist,3)), alpha=0.4)))
#Timelag
head(timeLag(X50, units="mins"))
head(timeLag(X50, units="hours"))
median(timeLag(X50, units="hours"))
max(timeLag(X50, units="hours"))
max(timeLag(X50, units="days"))
plot(X50_df$individual.local.identifier..ID.~X50_df $date,xlab = "Sampling period",ylab = "Individuals ID Number")
#write.table(X50,file="c:/RWorkDir/jaguardatapaper/ID/X50.txt",row.names = F,quote=F,col.names=T,sep="\t")
#Histograms and plots
hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 50],
     main = 'Individual 50', xlab = "Time between fixes (h)")
plot(density.default(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 50],
     na.rm = T), main = 'Individual 50')
	 
breaks <- c(0, 0.5, 1.5, 2.5, 4.5, 8.5, 16.5, 32.5, 64.5, 128.5, 256.5, 1000)
breaks <- c(0, 1, 2*(1:256))
hh <- hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 50 &
                                     mov.data.diff$time.diff > 0 & 
                                     mov.data.diff$time.diff < 512],
									 main = 'Individual 50',xlab = 'Time between fixes (h)',
           plot = T, breaks = breaks)

#X51
X51=unstacked$X51
X51_df <- as(X51, "data.frame")
head(X51_df)
nrow(X51_df)
range(X51_df$date)
m <- get_map(bbox(extent(X51)*1.1), source="google", zoom=11, maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=X51_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
#m <- get_map(bbox(extent(X51)*1.1), source="google", zoom=12, maptype="satellite") #define the box for map based on the trajectory coordinates
#ggmap(m)+geom_path(data=X51_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
#' ## Display space-time cube. 
open3d()
# To get a bigger window than the default
par3d(windowRect = c(100, 100, 612, 612))
Sys.sleep(0.1) # Allow sluggish window managers to catch up
with(X51_df, plot3d(x,y,date, type="l", col=as.numeric(X51_df$Event_ID)))
(stcube<-with(X51_df, plot3d(x,y,date, type="l",col=as.numeric(cut(X51_df$dist,3)), alpha=0.4)))
#Timelag
head(timeLag(X51, units="mins"))
head(timeLag(X51, units="hours"))
median(timeLag(X51, units="hours"))
max(timeLag(X51, units="hours"))
max(timeLag(X51, units="days"))
plot(X51_df$individual.local.identifier..ID.~X51_df $date,xlab = "Sampling period",ylab = "Individuals ID Number")
#write.table(X51,file="c:/RWorkDir/jaguardatapaper/ID/X51.txt",row.names = F,quote=F,col.names=T,sep="\t")
#Histograms and plots
hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 51],
     main = 'Individual 51', xlab = "Time between fixes (h)")
plot(density.default(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 51],
     na.rm = T), main = 'Individual 51')
	 
breaks <- c(0, 0.5, 1.5, 2.5, 4.5, 8.5, 16.5, 32.5, 64.5, 128.5, 256.5, 1000)
breaks <- c(0, 1, 2*(1:256))
hh <- hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 51 &
                                     mov.data.diff$time.diff > 0 & 
                                     mov.data.diff$time.diff < 512],
									 main = 'Individual 51',xlab = 'Time between fixes (h)',
           plot = T, breaks = breaks)

#X52
X52=unstacked$X52
X52_df <- as(X52, "data.frame")
head(X52_df)
nrow(X52_df)
range(X52_df$date)
m <- get_map(bbox(extent(X52)*1.1), source="google", zoom=14, maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=X52_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
#' ## Display space-time cube. 
open3d()
# To get a bigger window than the default
par3d(windowRect = c(100, 100, 612, 612))
Sys.sleep(0.1) # Allow sluggish window managers to catch up
with(X52_df, plot3d(x,y,date, type="l", col=as.numeric(X52_df$Event_ID)))
(stcube<-with(X52_df, plot3d(x,y,date, type="l",col=as.numeric(cut(X52_df$dist,3)), alpha=0.4)))
#Timelag
head(timeLag(X52, units="mins"))
head(timeLag(X52, units="hours"))
median(timeLag(X52, units="hours"))
max(timeLag(X52, units="hours"))
max(timeLag(X52, units="days"))
plot(X52_df$individual.local.identifier..ID.~X52_df $date,xlab = "Sampling period",ylab = "Individuals ID Number")
#write.table(X52,file="c:/RWorkDir/jaguardatapaper/ID/X52.txt",row.names = F,quote=F,col.names=T,sep="\t")
#Histograms and plots
hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 52],
     main = 'Individual 52', xlab = "Time between fixes (h)")
plot(density.default(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 52],
     na.rm = T), main = 'Individual 52')
	 
breaks <- c(0, 0.5, 1.5, 2.5, 4.5, 8.5, 16.5, 32.5, 64.5, 128.5, 256.5, 1000)
breaks <- c(0, 1, 2*(1:256))
hh <- hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 52 &
                                     mov.data.diff$time.diff > 0 & 
                                     mov.data.diff$time.diff < 512],
									 main = 'Individual 52',xlab = 'Time between fixes (h)',
           plot = T, breaks = breaks)

#X53
X53=unstacked$X53
X53_df <- as(X53, "data.frame")
head(X53_df)
nrow(X53_df)
range(X53_df$date)
m <- get_map(bbox(extent(X53)*1.1), source="google", zoom=12, maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=X53_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
#' ## Display space-time cube. 
open3d()
# To get a bigger window than the default
par3d(windowRect = c(100, 100, 612, 612))
Sys.sleep(0.1) # Allow sluggish window managers to catch up
with(X53_df, plot3d(x,y,date, type="l", col=as.numeric(X53_df$Event_ID)))
(stcube<-with(X53_df, plot3d(x,y,date, type="l",col=as.numeric(cut(X53_df$dist,3)), alpha=0.4)))
#Timelag
head(timeLag(X53, units="mins"))
head(timeLag(X53, units="hours"))
median(timeLag(X53, units="hours"))
max(timeLag(X53, units="hours"))
max(timeLag(X53, units="days"))
plot(X53_df$individual.local.identifier..ID.~X53_df $date,xlab = "Sampling period",ylab = "Individuals ID Number")
#write.table(X53,file="c:/RWorkDir/jaguardatapaper/ID/X53.txt",row.names = F,quote=F,col.names=T,sep="\t")
#Histograms and plots
hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 53],
     main = 'Individual 53', xlab = "Time between fixes (h)")
plot(density.default(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 53],
     na.rm = T), main = 'Individual 53')
	 
breaks <- c(0, 0.5, 1.5, 2.5, 4.5, 8.5, 16.5, 32.5, 64.5, 128.5, 256.5, 1000)
breaks <- c(0, 1, 2*(1:256))
hh <- hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 53 &
                                     mov.data.diff$time.diff > 0 & 
                                     mov.data.diff$time.diff < 512],
									 main = 'Individual 53',xlab = 'Time between fixes (h)',
           plot = T, breaks = breaks)

#X54
X54=unstacked$X54
X54_df <- as(X54, "data.frame")
head(X54_df)
nrow(X54_df)
range(X54_df$date)
m <- get_map(bbox(extent(X54)*1.1), source="google", zoom=13, maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=X54_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
m <- get_map(bbox(extent(X54)*1.1), source="google", zoom=12, maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=X54_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
#' ## Display space-time cube. 
open3d()
# To get a bigger window than the default
par3d(windowRect = c(100, 100, 612, 612))
Sys.sleep(0.1) # Allow sluggish window managers to catch up
with(X54_df, plot3d(x,y,date, type="l", col=as.numeric(X54_df$Event_ID)))
(stcube<-with(X54_df, plot3d(x,y,date, type="l",col=as.numeric(cut(X54_df$dist,3)), alpha=0.4)))
#Timelag
head(timeLag(X54, units="mins"))
head(timeLag(X54, units="hours"))
median(timeLag(X54, units="hours"))
max(timeLag(X54, units="hours"))
max(timeLag(X54, units="days"))
plot(X54_df$individual.local.identifier..ID.~X54_df $date,xlab = "Sampling period",ylab = "Individuals ID Number")
#write.table(X54,file="c:/RWorkDir/jaguardatapaper/ID/X54.txt",row.names = F,quote=F,col.names=T,sep="\t")
#Histograms and plots
hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 54],
     main = 'Individual 54', xlab = "Time between fixes (h)")
plot(density.default(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 54],
     na.rm = T), main = 'Individual 54')
	 
breaks <- c(0, 0.5, 1.5, 2.5, 4.5, 8.5, 16.5, 32.5, 64.5, 128.5, 256.5, 1000)
breaks <- c(0, 1, 2*(1:256))
hh <- hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 54 &
                                     mov.data.diff$time.diff > 0 & 
                                     mov.data.diff$time.diff < 512],
									 main = 'Individual 54',xlab = 'Time between fixes (h)',
           plot = T, breaks = breaks)

#X55
X55=unstacked$X55
X55_df <- as(X55, "data.frame")
head(X55_df)
nrow(X55_df)
range(X55_df$date)
m <- get_map(bbox(extent(X55)*1.1), source="google", zoom=12, maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=X55_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
#' ## Display space-time cube. 
open3d()
# To get a bigger window than the default
par3d(windowRect = c(100, 100, 612, 612))
Sys.sleep(0.1) # Allow sluggish window managers to catch up
with(X55_df, plot3d(x,y,date, type="l", col=as.numeric(X55_df$Event_ID)))
(stcube<-with(X55_df, plot3d(x,y,date, type="l",col=as.numeric(cut(X55_df$dist,3)), alpha=0.4)))
#Timelag
head(timeLag(X55, units="mins"))
head(timeLag(X55, units="hours"))
median(timeLag(X55, units="hours"))
max(timeLag(X55, units="hours"))
max(timeLag(X55, units="days"))
plot(X55_df$individual.local.identifier..ID.~X55_df $date,xlab = "Sampling period",ylab = "Individuals ID Number")
#write.table(X55,file="c:/RWorkDir/jaguardatapaper/ID/X55.txt",row.names = F,quote=F,col.names=T,sep="\t")
#Histograms and plots
hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 55],
     main = 'Individual 55', xlab = "Time between fixes (h)")
plot(density.default(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 55],
     na.rm = T), main = 'Individual 55')
	 
breaks <- c(0, 0.5, 1.5, 2.5, 4.5, 8.5, 16.5, 32.5, 64.5, 128.5, 256.5, 1000)
breaks <- c(0, 1, 2*(1:256))
hh <- hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 55 &
                                     mov.data.diff$time.diff > 0 & 
                                     mov.data.diff$time.diff < 512],
									 main = 'Individual 55',xlab = 'Time between fixes (h)',
           plot = T, breaks = breaks)

#X56
X56=unstacked$X56
X56_df <- as(X56, "data.frame")
head(X56_df)
nrow(X56_df)
range(X56_df$date)
m <- get_map(bbox(extent(X56)*1.1), source="google", zoom=13, maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=X56_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
m <- get_map(bbox(extent(X56)*1.1), source="google", zoom=12, maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=X56_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
#' ## Display space-time cube. 
open3d()
# To get a bigger window than the default
par3d(windowRect = c(100, 100, 612, 612))
Sys.sleep(0.1) # Allow sluggish window managers to catch up
with(X56_df, plot3d(x,y,date, type="l", col=as.numeric(X56_df$Event_ID)))
(stcube<-with(X56_df, plot3d(x,y,date, type="l",col=as.numeric(cut(X56_df$dist,3)), alpha=0.4)))
#Timelag
head(timeLag(X56, units="mins"))
head(timeLag(X56, units="hours"))
median(timeLag(X56, units="hours"))
max(timeLag(X56, units="hours"))
max(timeLag(X56, units="days"))
plot(X56_df$individual.local.identifier..ID.~X56_df $date,xlab = "Sampling period",ylab = "Individuals ID Number")
#write.table(X56,file="c:/RWorkDir/jaguardatapaper/ID/X56.txt",row.names = F,quote=F,col.names=T,sep="\t")
#Histograms and plots
hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 56],
     main = 'Individual 56', xlab = "Time between fixes (h)")
plot(density.default(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 56],
     na.rm = T), main = 'Individual 56')
	 
breaks <- c(0, 0.5, 1.5, 2.5, 4.5, 8.5, 16.5, 32.5, 64.5, 128.5, 256.5, 1000)
breaks <- c(0, 1, 2*(1:256))
hh <- hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 56 &
                                     mov.data.diff$time.diff > 0 & 
                                     mov.data.diff$time.diff < 512],
									 main = 'Individual 56',xlab = 'Time between fixes (h)',
           plot = T, breaks = breaks)

#X57
X57=unstacked$X57
X57_df <- as(X57, "data.frame")
head(X57_df)
nrow(X57_df)
range(X57_df$date)
m <- get_map(bbox(extent(X57)*1.1), source="google", zoom=14, maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=X57_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
#' ## Display space-time cube. 
open3d()
# To get a bigger window than the default
par3d(windowRect = c(100, 100, 612, 612))
Sys.sleep(0.1) # Allow sluggish window managers to catch up
with(X57_df, plot3d(x,y,date, type="l", col=as.numeric(X57_df$Event_ID)))
(stcube<-with(X57_df, plot3d(x,y,date, type="l",col=as.numeric(cut(X57_df$dist,3)), alpha=0.4)))
#Timelag
head(timeLag(X57, units="mins"))
head(timeLag(X57, units="hours"))
median(timeLag(X57, units="hours"))
max(timeLag(X57, units="hours"))
max(timeLag(X57, units="days"))
plot(X57_df$individual.local.identifier..ID.~X57_df $date,xlab = "Sampling period",ylab = "Individuals ID Number")
#write.table(X57,file="c:/RWorkDir/jaguardatapaper/ID/X57.txt",row.names = F,quote=F,col.names=T,sep="\t")
#Histograms and plots
hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 57],
     main = 'Individual 57', xlab = "Time between fixes (h)")
plot(density.default(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 57],
     na.rm = T), main = 'Individual 57')
	 
breaks <- c(0, 0.5, 1.5, 2.5, 4.5, 8.5, 16.5, 32.5, 64.5, 128.5, 256.5, 1000)
breaks <- c(0, 1, 2*(1:256))
hh <- hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 57 &
                                     mov.data.diff$time.diff > 0 & 
                                     mov.data.diff$time.diff < 512],
									 main = 'Individual 57',xlab = 'Time between fixes (h)',
           plot = T, breaks = breaks)

#X58
X58=unstacked$X58
X58_df <- as(X58, "data.frame")
head(X58_df)
nrow(X58_df)
range(X58_df$date)
m <- get_map(bbox(extent(X58)*1.1), source="google", zoom=12, maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=X58_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
m <- get_map(bbox(extent(X58)*1.1), source="google", zoom=13, maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=X58_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
#' ## Display space-time cube. 
open3d()
# To get a bigger window than the default
par3d(windowRect = c(100, 100, 612, 612))
Sys.sleep(0.1) # Allow sluggish window managers to catch up
with(X58_df, plot3d(x,y,date, type="l", col=as.numeric(X58_df$Event_ID)))
(stcube<-with(X58_df, plot3d(x,y,date, type="l",col=as.numeric(cut(X58_df$dist,3)), alpha=0.4)))
#Timelag
head(timeLag(X58, units="mins"))
head(timeLag(X58, units="hours"))
median(timeLag(X58, units="hours"))
max(timeLag(X58, units="hours"))
max(timeLag(X58, units="days"))
plot(X58_df$individual.local.identifier..ID.~X58_df $date,xlab = "Sampling period",ylab = "Individuals ID Number")
#write.table(X58,file="c:/RWorkDir/jaguardatapaper/ID/X58.txt",row.names = F,quote=F,col.names=T,sep="\t")
#Histograms and plots
hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 58],
     main = 'Individual 58', xlab = "Time between fixes (h)")
plot(density.default(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 58],
     na.rm = T), main = 'Individual 58')
	 
breaks <- c(0, 0.5, 1.5, 2.5, 4.5, 8.5, 16.5, 32.5, 64.5, 128.5, 256.5, 1000)
breaks <- c(0, 1, 2*(1:256))
hh <- hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 58 &
                                     mov.data.diff$time.diff > 0 & 
                                     mov.data.diff$time.diff < 512],
									 main = 'Individual 58',xlab = 'Time between fixes (h)',
           plot = T, breaks = breaks)

#X59
X59=unstacked$X59
X59_df <- as(X59, "data.frame")
head(X59_df)
nrow(X59_df)
range(X59_df$date)
m <- get_map(bbox(extent(X59)*1.1), source="google", zoom=12, maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=X59_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
#' ## Display space-time cube. 
open3d()
# To get a bigger window than the default
par3d(windowRect = c(100, 100, 612, 612))
Sys.sleep(0.1) # Allow sluggish window managers to catch up
with(X59_df, plot3d(x,y,date, type="l", col=as.numeric(X59_df$Event_ID)))
(stcube<-with(X59_df, plot3d(x,y,date, type="l",col=as.numeric(cut(X59_df$dist,3)), alpha=0.4)))
#Timelag
head(timeLag(X59, units="mins"))
head(timeLag(X59, units="hours"))
median(timeLag(X59, units="hours"))
max(timeLag(X59, units="hours"))
max(timeLag(X59, units="days"))
plot(X59_df$individual.local.identifier..ID.~X59_df $date,xlab = "Sampling period",ylab = "Individuals ID Number")
#write.table(X59,file="c:/RWorkDir/jaguardatapaper/ID/X59.txt",row.names = F,quote=F,col.names=T,sep="\t")
#Histograms and plots
hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 59],
     main = 'Individual 59', xlab = "Time between fixes (h)")
plot(density.default(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 59],
     na.rm = T), main = 'Individual 59')
	 
breaks <- c(0, 0.5, 1.5, 2.5, 4.5, 8.5, 16.5, 32.5, 64.5, 128.5, 256.5, 1000)
breaks <- c(0, 1, 2*(1:256))
hh <- hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 59 &
                                     mov.data.diff$time.diff > 0 & 
                                     mov.data.diff$time.diff < 512],
									 main = 'Individual 59',xlab = 'Time between fixes (h)',
           plot = T, breaks = breaks)

#X60
X60=unstacked$X60
X60_df <- as(X60, "data.frame")
head(X60_df)
nrow(X60_df)
range(X60_df$date)
m <- get_map(bbox(extent(X60)*1.1), source="google", zoom=12, maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=X60_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
#' ## Display space-time cube. 
open3d()
# To get a bigger window than the default
par3d(windowRect = c(100, 100, 612, 612))
Sys.sleep(0.1) # Allow sluggish window managers to catch up
with(X60_df, plot3d(x,y,date, type="l", col=as.numeric(X60_df$Event_ID)))
(stcube<-with(X60_df, plot3d(x,y,date, type="l",col=as.numeric(cut(X60_df$dist,3)), alpha=0.4)))
#Timelag
head(timeLag(X60, units="mins"))
head(timeLag(X60, units="hours"))
median(timeLag(X60, units="hours"))
max(timeLag(X60, units="hours"))
max(timeLag(X60, units="days"))
plot(X60_df$individual.local.identifier..ID.~X60_df $date,xlab = "Sampling period",ylab = "Individuals ID Number")
#write.table(X60,file="c:/RWorkDir/jaguardatapaper/ID/X60.txt",row.names = F,quote=F,col.names=T,sep="\t")
#Histograms and plots
hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 60],
     main = 'Individual 60', xlab = "Time between fixes (h)")
plot(density.default(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 60],
     na.rm = T), main = 'Individual 60')
	 
breaks <- c(0, 0.5, 1.5, 2.5, 4.5, 8.5, 16.5, 32.5, 64.5, 128.5, 256.5, 1000)
breaks <- c(0, 1, 2*(1:256))
hh <- hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 60 &
                                     mov.data.diff$time.diff > 0 & 
                                     mov.data.diff$time.diff < 512],
									 main = 'Individual 60',xlab = 'Time between fixes (h)',
           plot = T, breaks = breaks)

#X61
X61=unstacked$X61
X61_df <- as(X61, "data.frame")
head(X61_df)
nrow(X61_df)
range(X61_df$date)
m <- get_map(bbox(extent(X61)*1.1), source="google", zoom=11, maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=X61_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
#' ## Display space-time cube. 
open3d()
# To get a bigger window than the default
par3d(windowRect = c(100, 100, 612, 612))
Sys.sleep(0.1) # Allow sluggish window managers to catch up
with(X61_df, plot3d(x,y,date, type="l", col=as.numeric(X61_df$Event_ID)))
(stcube<-with(X61_df, plot3d(x,y,date, type="l",col=as.numeric(cut(X61_df$dist,3)), alpha=0.4)))
#Timelag
head(timeLag(X61, units="mins"))
head(timeLag(X61, units="hours"))
median(timeLag(X61, units="hours"))
max(timeLag(X61, units="hours"))
max(timeLag(X61, units="days"))
plot(X61_df$individual.local.identifier..ID.~X61_df $date,xlab = "Sampling period",ylab = "Individuals ID Number")
#write.table(X61,file="c:/RWorkDir/jaguardatapaper/ID/X61.txt",row.names = F,quote=F,col.names=T,sep="\t")
#Histograms and plots
hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 61],
     main = 'Individual 61', xlab = "Time between fixes (h)")
plot(density.default(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 61],
     na.rm = T), main = 'Individual 61')
	 
breaks <- c(0, 0.5, 1.5, 2.5, 4.5, 8.5, 16.5, 32.5, 64.5, 128.5, 256.5, 1000)
breaks <- c(0, 1, 2*(1:256))
hh <- hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 61 &
                                     mov.data.diff$time.diff > 0 & 
                                     mov.data.diff$time.diff < 512],
									 main = 'Individual 61',xlab = 'Time between fixes (h)',
           plot = T, breaks = breaks)

#X62
X62=unstacked$X62
X62_df <- as(X62, "data.frame")
head(X62_df)
nrow(X62_df)
range(X62_df$date)
m <- get_map(bbox(extent(X62)*1.1), source="google", zoom=11, maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=X62_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
#' ## Display space-time cube. 
open3d()
# To get a bigger window than the default
par3d(windowRect = c(100, 100, 612, 612))
Sys.sleep(0.1) # Allow sluggish window managers to catch up
with(X62_df, plot3d(x,y,date, type="l", col=as.numeric(X62_df$Event_ID)))
(stcube<-with(X62_df, plot3d(x,y,date, type="l",col=as.numeric(cut(X62_df$dist,3)), alpha=0.4)))
#Timelag
head(timeLag(X62, units="mins"))
head(timeLag(X62, units="hours"))
median(timeLag(X62, units="hours"))
max(timeLag(X62, units="hours"))
max(timeLag(X62, units="days"))
plot(X62_df$individual.local.identifier..ID.~X62_df $date,xlab = "Sampling period",ylab = "Individuals ID Number")
#write.table(X62,file="c:/RWorkDir/jaguardatapaper/ID/X62.txt",row.names = F,quote=F,col.names=T,sep="\t")
#Histograms and plots
hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 62],
     main = 'Individual 62', xlab = "Time between fixes (h)")
plot(density.default(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 62],
     na.rm = T), main = 'Individual 62')
	 
breaks <- c(0, 0.5, 1.5, 2.5, 4.5, 8.5, 16.5, 32.5, 64.5, 128.5, 256.5, 1000)
breaks <- c(0, 1, 2*(1:256))
hh <- hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 62 &
                                     mov.data.diff$time.diff > 0 & 
                                     mov.data.diff$time.diff < 512],
									 main = 'Individual 62',xlab = 'Time between fixes (h)',
           plot = T, breaks = breaks)

#X63
X63=unstacked$X63
X63_df <- as(X63, "data.frame")
head(X63_df)
nrow(X63_df)
range(X63_df$date)
m <- get_map(bbox(extent(X63)*1.1), source="google", zoom=11, maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=X63_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
#' ## Display space-time cube. 
open3d()
# To get a bigger window than the default
par3d(windowRect = c(100, 100, 612, 612))
Sys.sleep(0.1) # Allow sluggish window managers to catch up
with(X63_df, plot3d(x,y,date, type="l", col=as.numeric(X63_df$Event_ID)))
(stcube<-with(X63_df, plot3d(x,y,date, type="l",col=as.numeric(cut(X63_df$dist,3)), alpha=0.4)))
#Timelag
head(timeLag(X63, units="mins"))
head(timeLag(X63, units="hours"))
median(timeLag(X63, units="hours"))
max(timeLag(X63, units="hours"))
max(timeLag(X63, units="days"))
plot(X63_df$individual.local.identifier..ID.~X63_df $date,xlab = "Sampling period",ylab = "Individuals ID Number")
#write.table(X63,file="c:/RWorkDir/jaguardatapaper/ID/X63.txt",row.names = F,quote=F,col.names=T,sep="\t")
#Histograms and plots
hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 63],
     main = 'Individual 63', xlab = "Time between fixes (h)")
plot(density.default(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 63],
     na.rm = T), main = 'Individual 63')
	 
breaks <- c(0, 0.5, 1.5, 2.5, 4.5, 8.5, 16.5, 32.5, 64.5, 128.5, 256.5, 1000)
breaks <- c(0, 1, 2*(1:256))
hh <- hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 63 &
                                     mov.data.diff$time.diff > 0 & 
                                     mov.data.diff$time.diff < 512],
									 main = 'Individual 63',xlab = 'Time between fixes (h)',
           plot = T, breaks = breaks)

#X64
X64=unstacked$X64
X64_df <- as(X64, "data.frame")
head(X64_df)
nrow(X64_df)
range(X64_df$date)
m <- get_map(bbox(extent(X64)*1.1), source="google", zoom=10, maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=X64_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
#' ## Display space-time cube. 
open3d()
# To get a bigger window than the default
par3d(windowRect = c(100, 100, 612, 612))
Sys.sleep(0.1) # Allow sluggish window managers to catch up
with(X64_df, plot3d(x,y,date, type="l", col=as.numeric(X64_df$Event_ID)))
(stcube<-with(X64_df, plot3d(x,y,date, type="l",col=as.numeric(cut(X64_df$dist,3)), alpha=0.4)))
#Timelag
head(timeLag(X64, units="mins"))
head(timeLag(X64, units="hours"))
median(timeLag(X64, units="hours"))
max(timeLag(X64, units="hours"))
max(timeLag(X64, units="days"))
plot(X64_df$individual.local.identifier..ID.~X64_df $date,xlab = "Sampling period",ylab = "Individuals ID Number")
#write.table(X64,file="c:/RWorkDir/jaguardatapaper/ID/X64.txt",row.names = F,quote=F,col.names=T,sep="\t")
#Histograms and plots
hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 64],
     main = 'Individual 64', xlab = "Time between fixes (h)")
plot(density.default(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 64],
     na.rm = T), main = 'Individual 64')
	 
breaks <- c(0, 0.5, 1.5, 2.5, 4.5, 8.5, 16.5, 32.5, 64.5, 128.5, 256.5, 1000)
breaks <- c(0, 1, 2*(1:256))
hh <- hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 64 &
                                     mov.data.diff$time.diff > 0 & 
                                     mov.data.diff$time.diff < 512],
									 main = 'Individual 64',xlab = 'Time between fixes (h)',
           plot = T, breaks = breaks)

#X65
X65=unstacked$X65
X65_df <- as(X65, "data.frame")
head(X65_df)
nrow(X65_df)
range(X65_df$date)
m <- get_map(bbox(extent(X65)*1.1), source="google", zoom=11, maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=X65_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
m <- get_map(bbox(extent(X65)*1.1), source="google", zoom=10, maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=X65_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
#' ## Display space-time cube. 
open3d()
# To get a bigger window than the default
par3d(windowRect = c(100, 100, 612, 612))
Sys.sleep(0.1) # Allow sluggish window managers to catch up
with(X65_df, plot3d(x,y,date, type="l", col=as.numeric(X65_df$Event_ID)))
(stcube<-with(X65_df, plot3d(x,y,date, type="l",col=as.numeric(cut(X65_df$dist,3)), alpha=0.4)))
#Timelag
head(timeLag(X65, units="mins"))
head(timeLag(X65, units="hours"))
median(timeLag(X65, units="hours"))
max(timeLag(X65, units="hours"))
max(timeLag(X65, units="days"))
plot(X65_df$individual.local.identifier..ID.~X65_df $date,xlab = "Sampling period",ylab = "Individuals ID Number")
#write.table(X65,file="c:/RWorkDir/jaguardatapaper/ID/X65.txt",row.names = F,quote=F,col.names=T,sep="\t")
#Histograms and plots
hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 65],
     main = 'Individual 65', xlab = "Time between fixes (h)")
plot(density.default(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 65],
     na.rm = T), main = 'Individual 65')
	 
breaks <- c(0, 0.5, 1.5, 2.5, 4.5, 8.5, 16.5, 32.5, 64.5, 128.5, 256.5, 1000)
breaks <- c(0, 1, 2*(1:256))
hh <- hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 65 &
                                     mov.data.diff$time.diff > 0 & 
                                     mov.data.diff$time.diff < 512],
									 main = 'Individual 65',xlab = 'Time between fixes (h)',
           plot = T, breaks = breaks)

#X66
X66=unstacked$X66
X66_df <- as(X66, "data.frame")
head(X66_df)
nrow(X66_df)
range(X66_df$date)
m <- get_map(bbox(extent(X66)*1.1), source="google", zoom=12, maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=X66_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
#' ## Display space-time cube. 
open3d()
# To get a bigger window than the default
par3d(windowRect = c(100, 100, 612, 612))
Sys.sleep(0.1) # Allow sluggish window managers to catch up
with(X66_df, plot3d(x,y,date, type="l", col=as.numeric(X66_df$Event_ID)))
(stcube<-with(X66_df, plot3d(x,y,date, type="l",col=as.numeric(cut(X66_df$dist,3)), alpha=0.4)))
#Timelag
head(timeLag(X66, units="mins"))
head(timeLag(X66, units="hours"))
median(timeLag(X66, units="hours"))
max(timeLag(X66, units="hours"))
max(timeLag(X66, units="days"))
plot(X66_df$individual.local.identifier..ID.~X66_df $date,xlab = "Sampling period",ylab = "Individuals ID Number")
#write.table(X66,file="c:/RWorkDir/jaguardatapaper/ID/X66.txt",row.names = F,quote=F,col.names=T,sep="\t")
#Histograms and plots
hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 66],
     main = 'Individual 66', xlab = "Time between fixes (h)")
plot(density.default(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 66],
     na.rm = T), main = 'Individual 66')
	 
breaks <- c(0, 0.5, 1.5, 2.5, 4.5, 8.5, 16.5, 32.5, 64.5, 128.5, 256.5, 1000)
breaks <- c(0, 1, 2*(1:256))
hh <- hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 66 &
                                     mov.data.diff$time.diff > 0 & 
                                     mov.data.diff$time.diff < 512],
									 main = 'Individual 66',xlab = 'Time between fixes (h)',
           plot = T, breaks = breaks)

#X67
X67=unstacked$X67
X67_df <- as(X67, "data.frame")
head(X67_df)
nrow(X67_df)
range(X67_df$date)
m <- get_map(bbox(extent(X67)*1.1), source="google", zoom=11, maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=X67_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
#' ## Display space-time cube. 
open3d()
# To get a bigger window than the default
par3d(windowRect = c(100, 100, 612, 612))
Sys.sleep(0.1) # Allow sluggish window managers to catch up
with(X67_df, plot3d(x,y,date, type="l", col=as.numeric(X67_df$Event_ID)))
(stcube<-with(X67_df, plot3d(x,y,date, type="l",col=as.numeric(cut(X67_df$dist,3)), alpha=0.4)))
#Timelag
head(timeLag(X67, units="mins"))
head(timeLag(X67, units="hours"))
median(timeLag(X67, units="hours"))
max(timeLag(X67, units="hours"))
max(timeLag(X67, units="days"))
plot(X67_df$individual.local.identifier..ID.~X67_df $date,xlab = "Sampling period",ylab = "Individuals ID Number")
#write.table(X67,file="c:/RWorkDir/jaguardatapaper/ID/X67.txt",row.names = F,quote=F,col.names=T,sep="\t")
#Histograms and plots
hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 67],
     main = 'Individual 67', xlab = "Time between fixes (h)")
plot(density.default(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 67],
     na.rm = T), main = 'Individual 67')
	 
breaks <- c(0, 0.5, 1.5, 2.5, 4.5, 8.5, 16.5, 32.5, 64.5, 128.5, 256.5, 1000)
breaks <- c(0, 1, 2*(1:256))
hh <- hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 67 &
                                     mov.data.diff$time.diff > 0 & 
                                     mov.data.diff$time.diff < 512],
									 main = 'Individual 67',xlab = 'Time between fixes (h)',
           plot = T, breaks = breaks)

#X68
X68=unstacked$X68
X68_df <- as(X68, "data.frame")
head(X68_df)
nrow(X68_df)
range(X68_df$date)
m <- get_map(bbox(extent(X68)*1.1), source="google", zoom=12, maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=X68_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
#' ## Display space-time cube. 
open3d()
# To get a bigger window than the default
par3d(windowRect = c(100, 100, 612, 612))
Sys.sleep(0.1) # Allow sluggish window managers to catch up
with(X68_df, plot3d(x,y,date, type="l", col=as.numeric(X68_df$Event_ID)))
(stcube<-with(X68_df, plot3d(x,y,date, type="l",col=as.numeric(cut(X68_df$dist,3)), alpha=0.4)))
#Timelag
head(timeLag(X68, units="mins"))
head(timeLag(X68, units="hours"))
median(timeLag(X68, units="hours"))
max(timeLag(X68, units="hours"))
max(timeLag(X68, units="days"))
plot(X68_df$individual.local.identifier..ID.~X68_df $date,xlab = "Sampling period",ylab = "Individuals ID Number")
#write.table(X68,file="c:/RWorkDir/jaguardatapaper/ID/X68.txt",row.names = F,quote=F,col.names=T,sep="\t")
#Histograms and plots
hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 68],
     main = 'Individual 68', xlab = "Time between fixes (h)")
plot(density.default(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 68],
     na.rm = T), main = 'Individual 68')
	 
breaks <- c(0, 0.5, 1.5, 2.5, 4.5, 8.5, 16.5, 32.5, 64.5, 128.5, 256.5, 1000)
breaks <- c(0, 1, 2*(1:256))
hh <- hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 68 &
                                     mov.data.diff$time.diff > 0 & 
                                     mov.data.diff$time.diff < 512],
									 main = 'Individual 68',xlab = 'Time between fixes (h)',
           plot = T, breaks = breaks)

#X69
X69=unstacked$X69
X69_df <- as(X69, "data.frame")
head(X69_df)
nrow(X69_df)
range(X69_df$date)
m <- get_map(bbox(extent(X69)*1.1), source="google", zoom=12, maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=X69_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
#' ## Display space-time cube. 
open3d()
# To get a bigger window than the default
par3d(windowRect = c(100, 100, 612, 612))
Sys.sleep(0.1) # Allow sluggish window managers to catch up
with(X69_df, plot3d(x,y,date, type="l", col=as.numeric(X69_df$Event_ID)))
(stcube<-with(X69_df, plot3d(x,y,date, type="l",col=as.numeric(cut(X69_df$dist,3)), alpha=0.4)))
#Timelag
head(timeLag(X69, units="mins"))
head(timeLag(X69, units="hours"))
median(timeLag(X69, units="hours"))
max(timeLag(X69, units="hours"))
max(timeLag(X69, units="days"))
plot(X69_df$individual.local.identifier..ID.~X69_df $date,xlab = "Sampling period",ylab = "Individuals ID Number")
#write.table(X69,file="c:/RWorkDir/jaguardatapaper/ID/X69.txt",row.names = F,quote=F,col.names=T,sep="\t")
#Histograms and plots
hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 69],
     main = 'Individual 69', xlab = "Time between fixes (h)")
plot(density.default(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 69],
     na.rm = T), main = 'Individual 69')
	 
breaks <- c(0, 0.5, 1.5, 2.5, 4.5, 8.5, 16.5, 32.5, 64.5, 128.5, 256.5, 1000)
breaks <- c(0, 1, 2*(1:256))
hh <- hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 69 &
                                     mov.data.diff$time.diff > 0 & 
                                     mov.data.diff$time.diff < 512],
									 main = 'Individual 69',xlab = 'Time between fixes (h)',
           plot = T, breaks = breaks)

#X70
X70=unstacked$X70
X70_df <- as(X70, "data.frame")
head(X70_df)
nrow(X70_df)
range(X70_df$date)
m <- get_map(bbox(extent(X70)*1.1), source="google", zoom=9, maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=X70_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
#' ## Display space-time cube. 
open3d()
# To get a bigger window than the default
par3d(windowRect = c(100, 100, 612, 612))
Sys.sleep(0.1) # Allow sluggish window managers to catch up
with(X70_df, plot3d(x,y,date, type="l", col=as.numeric(X70_df$Event_ID)))
(stcube<-with(X70_df, plot3d(x,y,date, type="l",col=as.numeric(cut(X70_df$dist,3)), alpha=0.4)))
#Timelag
head(timeLag(X70, units="mins"))
head(timeLag(X70, units="hours"))
median(timeLag(X70, units="hours"))
max(timeLag(X70, units="hours"))
max(timeLag(X70, units="days"))
plot(X70_df$individual.local.identifier..ID.~X70_df $date,xlab = "Sampling period",ylab = "Individuals ID Number")
#write.table(X70,file="c:/RWorkDir/jaguardatapaper/ID/X70.txt",row.names = F,quote=F,col.names=T,sep="\t")
#Histograms and plots
hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 70],
     main = 'Individual 70', xlab = "Time between fixes (h)")
plot(density.default(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 70],
     na.rm = T), main = 'Individual 70')
	 
breaks <- c(0, 0.5, 1.5, 2.5, 4.5, 8.5, 16.5, 32.5, 64.5, 128.5, 256.5, 1000)
breaks <- c(0, 1, 2*(1:256))
hh <- hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 70 &
                                     mov.data.diff$time.diff > 0 & 
                                     mov.data.diff$time.diff < 512],
									 main = 'Individual 70',xlab = 'Time between fixes (h)',
           plot = T, breaks = breaks)

#X71
X71=unstacked$X71
X71_df <- as(X71, "data.frame")
head(X71_df)
nrow(X71_df)
range(X71_df$date)
m <- get_map(bbox(extent(X71)*1.1), source="google", zoom=11, maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=X71_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
#' ## Display space-time cube. 
open3d()
# To get a bigger window than the default
par3d(windowRect = c(100, 100, 612, 612))
Sys.sleep(0.1) # Allow sluggish window managers to catch up
with(X71_df, plot3d(x,y,date, type="l", col=as.numeric(X71_df$Event_ID)))
(stcube<-with(X71_df, plot3d(x,y,date, type="l",col=as.numeric(cut(X71_df$dist,3)), alpha=0.4)))
#Timelag
head(timeLag(X71, units="mins"))
head(timeLag(X71, units="hours"))
median(timeLag(X71, units="hours"))
max(timeLag(X71, units="hours"))
max(timeLag(X71, units="days"))
plot(X71_df$individual.local.identifier..ID.~X71_df $date,xlab = "Sampling period",ylab = "Individuals ID Number")
#write.table(X71,file="c:/RWorkDir/jaguardatapaper/ID/X71.txt",row.names = F,quote=F,col.names=T,sep="\t")
#Histograms and plots
hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 71],
     main = 'Individual 71', xlab = "Time between fixes (h)")
plot(density.default(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 71],
     na.rm = T), main = 'Individual 71')
	 
breaks <- c(0, 0.5, 1.5, 2.5, 4.5, 8.5, 16.5, 32.5, 64.5, 128.5, 256.5, 1000)
breaks <- c(0, 1, 2*(1:256))
hh <- hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 71 &
                                     mov.data.diff$time.diff > 0 & 
                                     mov.data.diff$time.diff < 512],
									 main = 'Individual 71',xlab = 'Time between fixes (h)',
           plot = T, breaks = breaks)

#X72
X72=unstacked$X72
X72_df <- as(X72, "data.frame")
head(X72_df)
nrow(X72_df)
range(X72_df$date)
m <- get_map(bbox(extent(X72)*1.1), source="google", zoom=11, maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=X72_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
#' ## Display space-time cube. 
open3d()
# To get a bigger window than the default
par3d(windowRect = c(100, 100, 612, 612))
Sys.sleep(0.1) # Allow sluggish window managers to catch up
with(X72_df, plot3d(x,y,date, type="l", col=as.numeric(X72_df$Event_ID)))
(stcube<-with(X72_df, plot3d(x,y,date, type="l",col=as.numeric(cut(X72_df$dist,3)), alpha=0.4)))
#Timelag
head(timeLag(X72, units="mins"))
head(timeLag(X72, units="hours"))
median(timeLag(X72, units="hours"))
max(timeLag(X72, units="hours"))
max(timeLag(X72, units="days"))
plot(X72_df$individual.local.identifier..ID.~X72_df $date,xlab = "Sampling period",ylab = "Individuals ID Number")
#write.table(X72,file="c:/RWorkDir/jaguardatapaper/ID/X72.txt",row.names = F,quote=F,col.names=T,sep="\t")
#Histograms and plots
hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 72],
     main = 'Individual 72', xlab = "Time between fixes (h)")
plot(density.default(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 72],
     na.rm = T), main = 'Individual 72')
	 
breaks <- c(0, 0.5, 1.5, 2.5, 4.5, 8.5, 16.5, 32.5, 64.5, 128.5, 256.5, 1000)
breaks <- c(0, 1, 2*(1:256))
hh <- hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 72 &
                                     mov.data.diff$time.diff > 0 & 
                                     mov.data.diff$time.diff < 512],
									 main = 'Individual 72',xlab = 'Time between fixes (h)',
           plot = T, breaks = breaks)

#X73
X73=unstacked$X73
X73_df <- as(X73, "data.frame")
head(X73_df)
nrow(X73_df)
range(X73_df$date)
m <- get_map(bbox(extent(X73)*1.1), source="google", zoom=11, maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=X73_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
#' ## Display space-time cube. 
open3d()
# To get a bigger window than the default
par3d(windowRect = c(100, 100, 612, 612))
Sys.sleep(0.1) # Allow sluggish window managers to catch up
with(X73_df, plot3d(x,y,date, type="l", col=as.numeric(X73_df$Event_ID)))
(stcube<-with(X73_df, plot3d(x,y,date, type="l",col=as.numeric(cut(X73_df$dist,3)), alpha=0.4)))
#Timelag
head(timeLag(X73, units="mins"))
head(timeLag(X73, units="hours"))
median(timeLag(X73, units="hours"))
max(timeLag(X73, units="hours"))
max(timeLag(X73, units="days"))
plot(X73_df$individual.local.identifier..ID.~X73_df $date,xlab = "Sampling period",ylab = "Individuals ID Number")
#write.table(X73,file="c:/RWorkDir/jaguardatapaper/ID/X73.txt",row.names = F,quote=F,col.names=T,sep="\t")
#Histograms and plots
hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 73],
     main = 'Individual 73', xlab = "Time between fixes (h)")
plot(density.default(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 73],
     na.rm = T), main = 'Individual 73')
	 
breaks <- c(0, 0.5, 1.5, 2.5, 4.5, 8.5, 16.5, 32.5, 64.5, 128.5, 256.5, 1000)
breaks <- c(0, 1, 2*(1:256))
hh <- hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 73 &
                                     mov.data.diff$time.diff > 0 & 
                                     mov.data.diff$time.diff < 512],
									 main = 'Individual 73',xlab = 'Time between fixes (h)',
           plot = T, breaks = breaks)

#X74
X74=unstacked$X74
X74_df <- as(X74, "data.frame")
head(X74_df)
nrow(X74_df)
range(X74_df$date)
m <- get_map(bbox(extent(X74)*1.1), source="google", zoom=10, maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=X74_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
#' ## Display space-time cube. 
open3d()
# To get a bigger window than the default
par3d(windowRect = c(100, 100, 612, 612))
Sys.sleep(0.1) # Allow sluggish window managers to catch up
with(X74_df, plot3d(x,y,date, type="l", col=as.numeric(X74_df$Event_ID)))
(stcube<-with(X74_df, plot3d(x,y,date, type="l",col=as.numeric(cut(X74_df$dist,3)), alpha=0.4)))
#Timelag
head(timeLag(X74, units="mins"))
head(timeLag(X74, units="hours"))
median(timeLag(X74, units="hours"))
max(timeLag(X74, units="hours"))
max(timeLag(X74, units="days"))
plot(X74_df$individual.local.identifier..ID.~X74_df $date,xlab = "Sampling period",ylab = "Individuals ID Number")
#write.table(X74,file="c:/RWorkDir/jaguardatapaper/ID/X74.txt",row.names = F,quote=F,col.names=T,sep="\t")
#Histograms and plots
hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 74],
     main = 'Individual 74', xlab = "Time between fixes (h)")
plot(density.default(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 74],
     na.rm = T), main = 'Individual 74')
	 
breaks <- c(0, 0.5, 1.5, 2.5, 4.5, 8.5, 16.5, 32.5, 64.5, 128.5, 256.5, 1000)
breaks <- c(0, 1, 2*(1:256))
hh <- hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 74 &
                                     mov.data.diff$time.diff > 0 & 
                                     mov.data.diff$time.diff < 512],
									 main = 'Individual 74',xlab = 'Time between fixes (h)',
           plot = T, breaks = breaks)

#X75
X75=unstacked$X75
X75_df <- as(X75, "data.frame")
head(X75_df)
nrow(X75_df)
range(X75_df$date)
m <- get_map(bbox(extent(X75)*1.1), source="google", zoom=12, maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=X75_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
#' ## Display space-time cube. 
open3d()
# To get a bigger window than the default
par3d(windowRect = c(100, 100, 612, 612))
Sys.sleep(0.1) # Allow sluggish window managers to catch up
with(X75_df, plot3d(x,y,date, type="l", col=as.numeric(X75_df$Event_ID)))
(stcube<-with(X75_df, plot3d(x,y,date, type="l",col=as.numeric(cut(X75_df$dist,3)), alpha=0.4)))
#Timelag
head(timeLag(X75, units="mins"))
head(timeLag(X75, units="hours"))
median(timeLag(X75, units="hours"))
max(timeLag(X75, units="hours"))
max(timeLag(X75, units="days"))
plot(X75_df$individual.local.identifier..ID.~X75_df $date,xlab = "Sampling period",ylab = "Individuals ID Number")
#write.table(X75,file="c:/RWorkDir/jaguardatapaper/ID/X75.txt",row.names = F,quote=F,col.names=T,sep="\t")
#Histograms and plots
hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 75],
     main = 'Individual 75', xlab = "Time between fixes (h)")
plot(density.default(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 75],
     na.rm = T), main = 'Individual 75')
	 
breaks <- c(0, 0.5, 1.5, 2.5, 4.5, 8.5, 16.5, 32.5, 64.5, 128.5, 256.5, 1000)
breaks <- c(0, 1, 2*(1:256))
hh <- hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 75 &
                                     mov.data.diff$time.diff > 0 & 
                                     mov.data.diff$time.diff < 512],
									 main = 'Individual 75',xlab = 'Time between fixes (h)',
           plot = T, breaks = breaks)

#X76
X76=unstacked$X76
X76_df <- as(X76, "data.frame")
head(X76_df)
nrow(X76_df)
range(X76_df$date)
m <- get_map(bbox(extent(X76)*1.1), source="google", zoom=11, maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=X76_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
#' ## Display space-time cube. 
open3d()
# To get a bigger window than the default
par3d(windowRect = c(100, 100, 612, 612))
Sys.sleep(0.1) # Allow sluggish window managers to catch up
with(X76_df, plot3d(x,y,date, type="l", col=as.numeric(X76_df$Event_ID)))
(stcube<-with(X76_df, plot3d(x,y,date, type="l",col=as.numeric(cut(X76_df$dist,3)), alpha=0.4)))
#Timelag
head(timeLag(X76, units="mins"))
head(timeLag(X76, units="hours"))
median(timeLag(X76, units="hours"))
max(timeLag(X76, units="hours"))
max(timeLag(X76, units="days"))
plot(X76_df$individual.local.identifier..ID.~X76_df $date,xlab = "Sampling period",ylab = "Individuals ID Number")
#write.table(X76,file="c:/RWorkDir/jaguardatapaper/ID/X76.txt",row.names = F,quote=F,col.names=T,sep="\t")
#Histograms and plots
hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 76],
     main = 'Individual 76', xlab = "Time between fixes (h)")
plot(density.default(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 76],
     na.rm = T), main = 'Individual 76')
	 
breaks <- c(0, 0.5, 1.5, 2.5, 4.5, 8.5, 16.5, 32.5, 64.5, 128.5, 256.5, 1000)
breaks <- c(0, 1, 2*(1:256))
hh <- hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 76 &
                                     mov.data.diff$time.diff > 0 & 
                                     mov.data.diff$time.diff < 512],
									 main = 'Individual 76',xlab = 'Time between fixes (h)',
           plot = T, breaks = breaks)

#X77
X77=unstacked$X77
X77_df <- as(X77, "data.frame")
head(X77_df)
nrow(X77_df)
range(X77_df$date)
m <- get_map(bbox(extent(X77)*1.1), source="google", zoom=10, maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=X77_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
#' ## Display space-time cube. 
open3d()
# To get a bigger window than the default
par3d(windowRect = c(100, 100, 612, 612))
Sys.sleep(0.1) # Allow sluggish window managers to catch up
with(X77_df, plot3d(x,y,date, type="l", col=as.numeric(X77_df$Event_ID)))
(stcube<-with(X77_df, plot3d(x,y,date, type="l",col=as.numeric(cut(X77_df$dist,3)), alpha=0.4)))
#Timelag
head(timeLag(X77, units="mins"))
head(timeLag(X77, units="hours"))
median(timeLag(X77, units="hours"))
max(timeLag(X77, units="hours"))
max(timeLag(X77, units="days"))
#m <- get_map(bbox(extent(X77)*1.1), source="google", zoom=11, maptype="satellite") #define the box for map based on the trajectory coordinates
#ggmap(m)+geom_path(data=X77_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
plot(X77_df$individual.local.identifier..ID.~X77_df $date,xlab = "Sampling period",ylab = "Individuals ID Number")
#write.table(X77,file="c:/RWorkDir/jaguardatapaper/ID/X77.txt",row.names = F,quote=F,col.names=T,sep="\t")
#Histograms and plots
hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 77],
     main = 'Individual 77', xlab = "Time between fixes (h)")
plot(density.default(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 77],
     na.rm = T), main = 'Individual 77')
	 
breaks <- c(0, 0.5, 1.5, 2.5, 4.5, 8.5, 16.5, 32.5, 64.5, 128.5, 256.5, 1000)
breaks <- c(0, 1, 2*(1:256))
hh <- hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 77 &
                                     mov.data.diff$time.diff > 0 & 
                                     mov.data.diff$time.diff < 512],
									 main = 'Individual 77',xlab = 'Time between fixes (h)',
           plot = T, breaks = breaks)

#X78
X78=unstacked$X78
X78_df <- as(X78, "data.frame")
head(X78_df)
nrow(X78_df)
range(X78_df$date)
m <- get_map(bbox(extent(X78)*1.1), source="google", zoom=12, maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=X78_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
m <- get_map(bbox(extent(X78)*1.1), source="google", zoom=11, maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=X78_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
#' ## Display space-time cube. 
open3d()
# To get a bigger window than the default
par3d(windowRect = c(100, 100, 612, 612))
Sys.sleep(0.1) # Allow sluggish window managers to catch up
with(X78_df, plot3d(x,y,date, type="l", col=as.numeric(X78_df$Event_ID)))
(stcube<-with(X78_df, plot3d(x,y,date, type="l",col=as.numeric(cut(X78_df$dist,3)), alpha=0.4)))
#Timelag
head(timeLag(X78, units="mins"))
head(timeLag(X78, units="hours"))
median(timeLag(X78, units="hours"))
max(timeLag(X78, units="hours"))
max(timeLag(X78, units="days"))
plot(X78_df$individual.local.identifier..ID.~X78_df $date,xlab = "Sampling period",ylab = "Individuals ID Number")
#write.table(X78,file="c:/RWorkDir/jaguardatapaper/ID/X78.txt",row.names = F,quote=F,col.names=T,sep="\t")
#Histograms and plots
hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 78],
     main = 'Individual 78', xlab = "Time between fixes (h)")
plot(density.default(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 78],
     na.rm = T), main = 'Individual 78')
	 
breaks <- c(0, 0.5, 1.5, 2.5, 4.5, 8.5, 16.5, 32.5, 64.5, 128.5, 256.5, 1000)
breaks <- c(0, 1, 2*(1:256))
hh <- hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 78 &
                                     mov.data.diff$time.diff > 0 & 
                                     mov.data.diff$time.diff < 512],
									 main = 'Individual 78',xlab = 'Time between fixes (h)',
           plot = T, breaks = breaks)
		   
#X79
X79=unstacked$X79
X79_df <- as(X79, "data.frame")
head(X79_df)
nrow(X79_df)
range(X79_df$date)
m <- get_map(bbox(extent(X79)*1.1), source="google", zoom=12, maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=X79_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
#' ## Display space-time cube. 
open3d()
# To get a bigger window than the default
par3d(windowRect = c(100, 100, 612, 612))
Sys.sleep(0.1) # Allow sluggish window managers to catch up
with(X79_df, plot3d(x,y,date, type="l", col=as.numeric(X79_df$Event_ID)))
(stcube<-with(X79_df, plot3d(x,y,date, type="l",col=as.numeric(cut(X79_df$dist,3)), alpha=0.4)))
#Timelag
head(timeLag(X79, units="mins"))
head(timeLag(X79, units="hours"))
median(timeLag(X79, units="hours"))
max(timeLag(X79, units="hours"))
max(timeLag(X79, units="days"))
plot(X79_df$individual.local.identifier..ID.~X79_df $date,xlab = "Sampling period",ylab = "Individuals ID Number")
#write.table(X79,file="c:/RWorkDir/jaguardatapaper/ID/X79.txt",row.names = F,quote=F,col.names=T,sep="\t")
#Histograms and plots
hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 79],
     main = 'Individual 79', xlab = "Time between fixes (h)")
plot(density.default(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 79],
     na.rm = T), main = 'Individual 79')
	 
breaks <- c(0, 0.5, 1.5, 2.5, 4.5, 8.5, 16.5, 32.5, 64.5, 128.5, 256.5, 1000)
breaks <- c(0, 1, 2*(1:256))
hh <- hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 79 &
                                     mov.data.diff$time.diff > 0 & 
                                     mov.data.diff$time.diff < 512],
									 main = 'Individual 79',xlab = 'Time between fixes (h)',
           plot = T, breaks = breaks)

#X80
X80=unstacked$X80
X80_df <- as(X80, "data.frame")
head(X80_df)
nrow(X80_df)
range(X80_df$date)
m <- get_map(bbox(extent(X80)*1.1), source="google", zoom=12, maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=X80_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
#' ## Display space-time cube. 
open3d()
# To get a bigger window than the default
par3d(windowRect = c(100, 100, 612, 612))
Sys.sleep(0.1) # Allow sluggish window managers to catch up
with(X80_df, plot3d(x,y,date, type="l", col=as.numeric(X80_df$Event_ID)))
(stcube<-with(X80_df, plot3d(x,y,date, type="l",col=as.numeric(cut(X80_df$dist,3)), alpha=0.4)))
#Timelag
head(timeLag(X80, units="mins"))
head(timeLag(X80, units="hours"))
median(timeLag(X80, units="hours"))
max(timeLag(X80, units="hours"))
max(timeLag(X80, units="days"))
#m <- get_map(bbox(extent(X80)*1.1), source="google", zoom=13, maptype="satellite") #define the box for map based on the trajectory coordinates
#ggmap(m)+geom_path(data=X80_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
plot(X80_df$individual.local.identifier..ID.~X80_df $date,xlab = "Sampling period",ylab = "Individuals ID Number")
#write.table(X80,file="c:/RWorkDir/jaguardatapaper/ID/X80.txt",row.names = F,quote=F,col.names=T,sep="\t")
#Histograms and plots
hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 80],
     main = 'Individual 80', xlab = "Time between fixes (h)")
plot(density.default(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 80],
     na.rm = T), main = 'Individual 80')
	 
breaks <- c(0, 0.5, 1.5, 2.5, 4.5, 8.5, 16.5, 32.5, 64.5, 128.5, 256.5, 1000)
breaks <- c(0, 1, 2*(1:256))
hh <- hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 80 &
                                     mov.data.diff$time.diff > 0 & 
                                     mov.data.diff$time.diff < 512],
									 main = 'Individual 80',xlab = 'Time between fixes (h)',
           plot = T, breaks = breaks)

#X81
X81=unstacked$X81
X81_df <- as(X81, "data.frame")
head(X81_df)
nrow(X81_df)
range(X81_df$date)
m <- get_map(bbox(extent(X81)*1.1), source="google", zoom=11, maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=X81_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
#' ## Display space-time cube. 
open3d()
# To get a bigger window than the default
par3d(windowRect = c(100, 100, 612, 612))
Sys.sleep(0.1) # Allow sluggish window managers to catch up
with(X81_df, plot3d(x,y,date, type="l", col=as.numeric(X81_df$Event_ID)))
(stcube<-with(X81_df, plot3d(x,y,date, type="l",col=as.numeric(cut(X81_df$dist,3)), alpha=0.4)))
#Timelag
head(timeLag(X81, units="mins"))
head(timeLag(X81, units="hours"))
median(timeLag(X81, units="hours"))
max(timeLag(X81, units="hours"))
max(timeLag(X81, units="days"))
plot(X81_df$individual.local.identifier..ID.~X81_df $date,xlab = "Sampling period",ylab = "Individuals ID Number")
#write.table(X81,file="c:/RWorkDir/jaguardatapaper/ID/X81.txt",row.names = F,quote=F,col.names=T,sep="\t")
#Histograms and plots
hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 81],
     main = 'Individual 81', xlab = "Time between fixes (h)")
plot(density.default(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 81],
     na.rm = T), main = 'Individual 81')
	 
breaks <- c(0, 0.5, 1.5, 2.5, 4.5, 8.5, 16.5, 32.5, 64.5, 128.5, 256.5, 1000)
breaks <- c(0, 1, 2*(1:256))
hh <- hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 81 &
                                     mov.data.diff$time.diff > 0 & 
                                     mov.data.diff$time.diff < 512],
									 main = 'Individual 81',xlab = 'Time between fixes (h)',
           plot = T, breaks = breaks)

#X82
X82=unstacked$X82
X82_df <- as(X82, "data.frame")
head(X82_df)
nrow(X82_df)
range(X82_df$date)
m <- get_map(bbox(extent(X82)*1.1), source="google", zoom=11, maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=X82_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
#' ## Display space-time cube. 
open3d()
# To get a bigger window than the default
par3d(windowRect = c(100, 100, 612, 612))
Sys.sleep(0.1) # Allow sluggish window managers to catch up
#parent <- currentSubscene3d()
with(X82_df, plot3d(x,y,date, type="l", col=as.numeric(X82_df$Event_ID)))
(stcube<-with(X82_df, plot3d(x,y,date, type="l",col=as.numeric(cut(X82_df$dist,3)), alpha=0.4)))
#Timelag
head(timeLag(X82, units="mins"))
head(timeLag(X82, units="hours"))
median(timeLag(X82, units="hours"))
max(timeLag(X82, units="hours"))
max(timeLag(X82, units="days"))
plot(X82_df$individual.local.identifier..ID.~X82_df $date,xlab = "Sampling period",ylab = "Individuals ID Number")
#write.table(X82,file="c:/RWorkDir/jaguardatapaper/ID/X82.txt",row.names = F,quote=F,col.names=T,sep="\t")
#Histograms and plots
hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 82],
     main = 'Individual 82', xlab = "Time between fixes (h)")
plot(density.default(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 82],
     na.rm = T), main = 'Individual 82')
	 
breaks <- c(0, 0.5, 1.5, 2.5, 4.5, 8.5, 16.5, 32.5, 64.5, 128.5, 256.5, 1000)
breaks <- c(0, 1, 2*(1:256))
hh <- hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 82 &
                                     mov.data.diff$time.diff > 0 & 
                                     mov.data.diff$time.diff < 512],
									 main = 'Individual 82',xlab = 'Time between fixes (h)',
           plot = T, breaks = breaks)

#X83
X83=unstacked$X83
X83_df <- as(X83, "data.frame")
head(X83_df)
nrow(X83_df)
range(X83_df$date)
m <- get_map(bbox(extent(X83)*1.1), source="google", zoom=11, maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=X83_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
#' ## Display space-time cube. 
open3d()
# To get a bigger window than the default
par3d(windowRect = c(100, 100, 612, 612))
Sys.sleep(0.1) # Allow sluggish window managers to catch up
#parent <- currentSubscene3d()
with(X83_df, plot3d(x,y,date, type="l", col=as.numeric(X83_df$Event_ID)))
(stcube<-with(X83_df, plot3d(x,y,date, type="l",col=as.numeric(cut(X83_df$dist,3)), alpha=0.4)))
head(timeLag(X83, units="mins"))
head(timeLag(X83, units="hours"))
median(timeLag(X83, units="hours"))
max(timeLag(X83, units="hours"))
max(timeLag(X83, units="days"))
plot(X83_df$individual.local.identifier..ID.~X83_df $date,xlab = "Sampling period",ylab = "Individuals ID Number")
#write.table(X83,file="c:/RWorkDir/jaguardatapaper/ID/X83.txt",row.names = F,quote=F,col.names=T,sep="\t")
#Histograms and plots
hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 83],
     main = 'Individual 83', xlab = "Time between fixes (h)")
plot(density.default(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 83],
     na.rm = T), main = 'Individual 83')
	 
breaks <- c(0, 0.5, 1.5, 2.5, 4.5, 8.5, 16.5, 32.5, 64.5, 128.5, 256.5, 1000)
breaks <- c(0, 1, 2*(1:256))
hh <- hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 83 &
                                     mov.data.diff$time.diff > 0 & 
                                     mov.data.diff$time.diff < 512],
									 main = 'Individual 83',xlab = 'Time between fixes (h)',
           plot = T, breaks = breaks)

#X84
X84=unstacked$X84
X84_df <- as(X84, "data.frame")
head(X84_df)
nrow(X84_df)
range(X84_df$date)
m <- get_map(bbox(extent(X84)*1.1), source="google", zoom=12, maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=X84_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
#' ## Display space-time cube. 
open3d()
# To get a bigger window than the default
par3d(windowRect = c(100, 100, 612, 612))
Sys.sleep(0.1) # Allow sluggish window managers to catch up
with(X84_df, plot3d(x,y,date, type="l", col=as.numeric(X84_df$Event_ID)))
(stcube<-with(X84_df, plot3d(x,y,date, type="l",col=as.numeric(cut(X84_df$dist,3)), alpha=0.4)))
#Timelag
head(timeLag(X84, units="mins"))
head(timeLag(X84, units="hours"))
median(timeLag(X84, units="hours"))
max(timeLag(X84, units="hours"))
max(timeLag(X84, units="days"))
plot(X84_df$individual.local.identifier..ID.~X84_df $date,xlab = "Sampling period",ylab = "Individuals ID Number")
#write.table(X84,file="c:/RWorkDir/jaguardatapaper/ID/X84.txt",row.names = F,quote=F,col.names=T,sep="\t")
#Histograms and plots
hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 84],
     main = 'Individual 84', xlab = "Time between fixes (h)")
plot(density.default(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 84],
     na.rm = T), main = 'Individual 84')
	 
breaks <- c(0, 0.5, 1.5, 2.5, 4.5, 8.5, 16.5, 32.5, 64.5, 128.5, 256.5, 1000)
breaks <- c(0, 1, 2*(1:256))
hh <- hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 84 &
                                     mov.data.diff$time.diff > 0 & 
                                     mov.data.diff$time.diff < 512],
									 main = 'Individual 84',xlab = 'Time between fixes (h)',
           plot = T, breaks = breaks)

#X85
X85=unstacked$X85
X85_df <- as(X85, "data.frame")
head(X85_df)
nrow(X85_df)
range(X85_df$date)
m <- get_map(bbox(extent(X85)*1.1), source="google", zoom=10, maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=X85_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
#' ## Display space-time cube. 
open3d()
# To get a bigger window than the default
par3d(windowRect = c(100, 100, 612, 612))
Sys.sleep(0.1) # Allow sluggish window managers to catch up
with(X85_df, plot3d(x,y,date, type="l", col=as.numeric(X85_df$Event_ID)))
(stcube<-with(X85_df, plot3d(x,y,date, type="l",col=as.numeric(cut(X85_df$dist,3)), alpha=0.4)))
#Timelag
head(timeLag(X85, units="mins"))
head(timeLag(X85, units="hours"))
median(timeLag(X85, units="hours"))
max(timeLag(X85, units="hours"))
max(timeLag(X85, units="days"))
#m <- get_map(bbox(extent(X85)*1.1), source="google", zoom=11, maptype="satellite") #define the box for map based on the trajectory coordinates
#ggmap(m)+geom_path(data=X85_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
plot(X85_df$individual.local.identifier..ID.~X85_df $date,xlab = "Sampling period",ylab = "Individuals ID Number")
#write.table(X85,file="c:/RWorkDir/jaguardatapaper/ID/X85.txt",row.names = F,quote=F,col.names=T,sep="\t")
#Histograms and plots
hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 85],
     main = 'Individual 85', xlab = "Time between fixes (h)")
plot(density.default(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 85],
     na.rm = T), main = 'Individual 85')
	 
breaks <- c(0, 0.5, 1.5, 2.5, 4.5, 8.5, 16.5, 32.5, 64.5, 128.5, 256.5, 1000)
breaks <- c(0, 1, 2*(1:256))
hh <- hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 85 &
                                     mov.data.diff$time.diff > 0 & 
                                     mov.data.diff$time.diff < 512],
									 main = 'Individual 85',xlab = 'Time between fixes (h)',
           plot = T, breaks = breaks)

#X86
X86=unstacked$X86
X86_df <- as(X86, "data.frame")
head(X86_df)
nrow(X86_df)
range(X86_df$date)
m <- get_map(bbox(extent(X86)*1.1), source="google", zoom=12, maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=X86_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
#' ## Display space-time cube. 
open3d()
# To get a bigger window than the default
par3d(windowRect = c(100, 100, 612, 612))
Sys.sleep(0.1) # Allow sluggish window managers to catch up
with(X86_df, plot3d(x,y,date, type="l", col=as.numeric(X86_df$Event_ID)))
(stcube<-with(X86_df, plot3d(x,y,date, type="l",col=as.numeric(cut(X86_df$dist,3)), alpha=0.4)))
#Timelag
head(timeLag(X86, units="mins"))
head(timeLag(X86, units="hours"))
median(timeLag(X86, units="hours"))
max(timeLag(X86, units="hours"))
max(timeLag(X86, units="days"))
#m <- get_map(bbox(extent(X86)*1.1), source="google", zoom=13, maptype="satellite") #define the box for map based on the trajectory coordinates
#ggmap(m)+geom_path(data=X86_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
plot(X86_df$individual.local.identifier..ID.~X86_df $date,xlab = "Sampling period",ylab = "Individuals ID Number")
#write.table(X86,file="c:/RWorkDir/jaguardatapaper/ID/X86.txt",row.names = F,quote=F,col.names=T,sep="\t")
#Histograms and plots
hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 86],
     main = 'Individual 86', xlab = "Time between fixes (h)")
plot(density.default(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 86],
     na.rm = T), main = 'Individual 86')
	 
breaks <- c(0, 0.5, 1.5, 2.5, 4.5, 8.5, 16.5, 32.5, 64.5, 128.5, 256.5, 1000)
breaks <- c(0, 1, 2*(1:256))
hh <- hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 86 &
                                     mov.data.diff$time.diff > 0 & 
                                     mov.data.diff$time.diff < 512],
									 main = 'Individual 86',xlab = 'Time between fixes (h)',
           plot = T, breaks = breaks)

#X87
X87=unstacked$X87
X87_df <- as(X87, "data.frame")
head(X87_df)
nrow(X87_df)
range(X87_df$date)
m <- get_map(bbox(extent(X87)*1.1), source="google", zoom=13, maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=X87_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
#' ## Display space-time cube. 
open3d()
# To get a bigger window than the default
par3d(windowRect = c(100, 100, 612, 612))
Sys.sleep(0.1) # Allow sluggish window managers to catch up
with(X87_df, plot3d(x,y,date, type="l", col=as.numeric(X87_df$Event_ID)))
(stcube<-with(X87_df, plot3d(x,y,date, type="l",col=as.numeric(cut(X87_df$dist,3)), alpha=0.4)))
#Timelag
head(timeLag(X87, units="mins"))
head(timeLag(X87, units="hours"))
median(timeLag(X87, units="hours"))
max(timeLag(X87, units="hours"))
max(timeLag(X87, units="days"))
plot(X87_df$individual.local.identifier..ID.~X87_df $date,xlab = "Sampling period",ylab = "Individuals ID Number")
#write.table(X87,file="c:/RWorkDir/jaguardatapaper/ID/X87.txt",row.names = F,quote=F,col.names=T,sep="\t")
#Histograms and plots
hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 87],
     main = 'Individual 87', xlab = "Time between fixes (h)")
plot(density.default(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 87],
     na.rm = T), main = 'Individual 87')
	 
breaks <- c(0, 0.5, 1.5, 2.5, 4.5, 8.5, 16.5, 32.5, 64.5, 128.5, 256.5, 1000)
breaks <- c(0, 1, 2*(1:256))
hh <- hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 87 &
                                     mov.data.diff$time.diff > 0 & 
                                     mov.data.diff$time.diff < 512],
									 main = 'Individual 87',xlab = 'Time between fixes (h)',
           plot = T, breaks = breaks)

#X88
X88=unstacked$X88
X88_df <- as(X88, "data.frame")
head(X88_df)
nrow(X88_df)
range(X88_df$date)
m <- get_map(bbox(extent(X88)*1.1), source="google", zoom=13, maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=X88_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
#' ## Display space-time cube. 
open3d()
# To get a bigger window than the default
par3d(windowRect = c(100, 100, 612, 612))
Sys.sleep(0.1) # Allow sluggish window managers to catch up
with(X88_df, plot3d(x,y,date, type="l", col=as.numeric(X88_df$Event_ID)))
(stcube<-with(X88_df, plot3d(x,y,date, type="l",col=as.numeric(cut(X88_df$dist,3)), alpha=0.4)))
#Timelag
head(timeLag(X88, units="mins"))
head(timeLag(X88, units="hours"))
median(timeLag(X88, units="hours"))
max(timeLag(X88, units="hours"))
max(timeLag(X88, units="days"))
plot(X88_df$individual.local.identifier..ID.~X88_df $date,xlab = "Sampling period",ylab = "Individuals ID Number")
#write.table(X88,file="c:/RWorkDir/jaguardatapaper/ID/X88.txt",row.names = F,quote=F,col.names=T,sep="\t")
#Histograms and plots
hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 88],
     main = 'Individual 88', xlab = "Time between fixes (h)")
plot(density.default(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 88],
     na.rm = T), main = 'Individual 88')
	 
breaks <- c(0, 0.5, 1.5, 2.5, 4.5, 8.5, 16.5, 32.5, 64.5, 128.5, 256.5, 1000)
breaks <- c(0, 1, 2*(1:256))
hh <- hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 88 &
                                     mov.data.diff$time.diff > 0 & 
                                     mov.data.diff$time.diff < 512],
									 main = 'Individual 88',xlab = 'Time between fixes (h)',
           plot = T, breaks = breaks)


#X89 # Cerrado Fragmented
X89=unstacked$X89
X89_df <- as(X89, "data.frame")
head(X89_df)
nrow(X89_df)
range(X89_df$date)
m <- get_map(bbox(extent(X89)*1.1), source="google", zoom=11, maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=X89_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
#' ## Display space-time cube. 
open3d()
# To get a bigger window than the default
par3d(windowRect = c(100, 100, 612, 612))
Sys.sleep(0.1) # Allow sluggish window managers to catch up
with(X89_df, plot3d(x,y,date, type="l", col=as.numeric(X89_df$Event_ID)))
(stcube<-with(X89_df, plot3d(x,y,date, type="l",col=as.numeric(cut(X89_df$dist,3)), alpha=0.4)))
#Timelag
mhead(timeLag(X89, units="mins"))
head(timeLag(X89, units="hours"))
median(timeLag(X89, units="hours"))
max(timeLag(X89, units="hours"))
max(timeLag(X89, units="days"))
plot(X89_df$individual.local.identifier..ID.~X89_df $date,xlab = "Sampling period",ylab = "Individuals ID Number")
#write.table(X89,file="c:/RWorkDir/jaguardatapaper/ID/X89.txt",row.names = F,quote=F,col.names=T,sep="\t")
#Histograms and plots
hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 89],
     main = 'Individual 89', xlab = "Time between fixes (h)")
plot(density.default(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 89],
     na.rm = T), main = 'Individual 89')
	 
breaks <- c(0, 0.5, 1.5, 2.5, 4.5, 8.5, 16.5, 32.5, 64.5, 128.5, 256.5, 1000)
breaks <- c(0, 1, 2*(1:256))
hh <- hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 89 &
                                     mov.data.diff$time.diff > 0 & 
                                     mov.data.diff$time.diff < 512],
									 main = 'Individual 89',xlab = 'Time between fixes (h)',
           plot = T, breaks = breaks)

#X90
X90=unstacked$X90
X90_df <- as(X90, "data.frame")
head(X90_df)
nrow(X90_df)
range(X90_df$date)
m <- get_map(bbox(extent(X90)*1.1), source="google", zoom=12, maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=X90_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
#' ## Display space-time cube. 
open3d()
# To get a bigger window than the default
par3d(windowRect = c(100, 100, 612, 612))
Sys.sleep(0.1) # Allow sluggish window managers to catch up
with(X90_df, plot3d(x,y,date, type="l", col=as.numeric(X90_df$Event_ID)))
(stcube<-with(X90_df, plot3d(x,y,date, type="l",col=as.numeric(cut(X90_df$dist,3)), alpha=0.4)))
#Timelag
head(timeLag(X90, units="mins"))
head(timeLag(X90, units="hours"))
median(timeLag(X90, units="hours"))
max(timeLag(X90, units="hours"))
max(timeLag(X90, units="days"))
plot(X90_df$individual.local.identifier..ID.~X90_df $date,xlab = "Sampling period",ylab = "Individuals ID Number")
#write.table(X90,file="c:/RWorkDir/jaguardatapaper/ID/X90.txt",row.names = F,quote=F,col.names=T,sep="\t")
#Histograms and plots
hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 90],
     main = 'Individual 90', xlab = "Time between fixes (h)")
plot(density.default(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 90],
     na.rm = T), main = 'Individual 90')
	 
breaks <- c(0, 0.5, 1.5, 2.5, 4.5, 8.5, 16.5, 32.5, 64.5, 128.5, 256.5, 1000)
breaks <- c(0, 1, 2*(1:256))
hh <- hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 90 &
                                     mov.data.diff$time.diff > 0 & 
                                     mov.data.diff$time.diff < 512],
									 main = 'Individual 90',xlab = 'Time between fixes (h)',
           plot = T, breaks = breaks)
		   
#X91
X91=unstacked$X91
X91_df <- as(X91, "data.frame")
head(X91_df)
nrow(X91_df)
range(X91_df$date)
m <- get_map(bbox(extent(X91)*1.1), source="google", zoom=12, maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=X91_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
#' ## Display space-time cube. 
open3d()
# To get a bigger window than the default
par3d(windowRect = c(100, 100, 612, 612))
Sys.sleep(0.1) # Allow sluggish window managers to catch up
with(X91_df, plot3d(x,y,date, type="l", col=as.numeric(X91_df$Event_ID)))
(stcube<-with(X91_df, plot3d(x,y,date, type="l",col=as.numeric(cut(X91_df$dist,3)), alpha=0.4)))
#Timelag
head(timeLag(X91, units="mins"))
head(timeLag(X91, units="hours"))
median(timeLag(X91, units="hours"))
max(timeLag(X91, units="hours"))
max(timeLag(X91, units="days"))
plot(X91_df$individual.local.identifier..ID.~X91_df $date,xlab = "Sampling period",ylab = "Individuals ID Number")
#write.table(X91,file="c:/RWorkDir/jaguardatapaper/ID/X91.txt",row.names = F,quote=F,col.names=T,sep="\t")
#Histograms and plots
hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 91],
     main = 'Individual 91', xlab = "Time between fixes (h)")
plot(density.default(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 91],
     na.rm = T), main = 'Individual 91')
	 
breaks <- c(0, 0.5, 1.5, 2.5, 4.5, 8.5, 16.5, 32.5, 64.5, 128.5, 256.5, 1000)
breaks <- c(0, 1, 2*(1:256))
hh <- hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 91 &
                                     mov.data.diff$time.diff > 0 & 
                                     mov.data.diff$time.diff < 512],
									 main = 'Individual 91',xlab = 'Time between fixes (h)',
           plot = T, breaks = breaks)

#X92
X92=unstacked$X92
X92_df <- as(X92, "data.frame")
head(X92_df)
nrow(X92_df)
range(X92_df$date)
m <- get_map(bbox(extent(X92)*1.1), source="google", zoom=12, maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=X92_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
#' ## Display space-time cube. 
open3d()
# To get a bigger window than the default
par3d(windowRect = c(100, 100, 612, 612))
Sys.sleep(0.1) # Allow sluggish window managers to catch up
with(X92_df, plot3d(x,y,date, type="l", col=as.numeric(X92_df$Event_ID)))
(stcube<-with(X92_df, plot3d(x,y,date, type="l",col=as.numeric(cut(X92_df$dist,3)), alpha=0.4)))
#Timelag
head(timeLag(X92, units="mins"))
head(timeLag(X92, units="hours"))
median(timeLag(X92, units="hours"))
max(timeLag(X92, units="hours"))
max(timeLag(X92, units="days"))
plot(X92_df$individual.local.identifier..ID.~X92_df $date,xlab = "Sampling period",ylab = "Individuals ID Number")
#write.table(X92,file="c:/RWorkDir/jaguardatapaper/ID/X92.txt",row.names = F,quote=F,col.names=T,sep="\t")
#Histograms and plots
hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 92],
     main = 'Individual 92', xlab = "Time between fixes (h)")
plot(density.default(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 92],
     na.rm = T), main = 'Individual 92')
	 
breaks <- c(0, 0.5, 1.5, 2.5, 4.5, 8.5, 16.5, 32.5, 64.5, 128.5, 256.5, 1000)
breaks <- c(0, 1, 2*(1:256))
hh <- hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 92 &
                                     mov.data.diff$time.diff > 0 & 
                                     mov.data.diff$time.diff < 512],
									 main = 'Individual 92',xlab = 'Time between fixes (h)',
           plot = T, breaks = breaks)

#X93
X93=unstacked$X93
X93_df <- as(X93, "data.frame")
head(X93_df)
nrow(X93_df)
range(X93_df$date)
m <- get_map(bbox(extent(X93)*1.1), source="google", zoom=12, maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=X93_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
#' ## Display space-time cube. 
open3d()
# To get a bigger window than the default
par3d(windowRect = c(100, 100, 612, 612))
Sys.sleep(0.1) # Allow sluggish window managers to catch up
with(X93_df, plot3d(x,y,date, type="l", col=as.numeric(X93_df$Event_ID)))
(stcube<-with(X93_df, plot3d(x,y,date, type="l",col=as.numeric(cut(X93_df$dist,3)), alpha=0.4)))
#Timelag
head(timeLag(X93, units="mins"))
head(timeLag(X93, units="hours"))
median(timeLag(X93, units="hours"))
max(timeLag(X93, units="hours"))
max(timeLag(X93, units="days"))
plot(X93_df$individual.local.identifier..ID.~X93_df $date,xlab = "Sampling period",ylab = "Individuals ID Number")
#write.table(X93,file="c:/RWorkDir/jaguardatapaper/ID/X93.txt",row.names = F,quote=F,col.names=T,sep="\t")
#Histograms and plots
hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 93],
     main = 'Individual 93', xlab = "Time between fixes (h)")
plot(density.default(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 93],
     na.rm = T), main = 'Individual 93')
	 
breaks <- c(0, 0.5, 1.5, 2.5, 4.5, 8.5, 16.5, 32.5, 64.5, 128.5, 256.5, 1000)
breaks <- c(0, 1, 2*(1:256))
hh <- hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 93 &
                                     mov.data.diff$time.diff > 0 & 
                                     mov.data.diff$time.diff < 512],
									 main = 'Individual 93',xlab = 'Time between fixes (h)',
           plot = T, breaks = breaks)

#X94
X94=unstacked$X94
X94_df <- as(X94, "data.frame")
head(X94_df)
nrow(X94_df)
range(X94_df$date)
m <- get_map(bbox(extent(X94)*1.1), source="google", zoom=11, maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=X94_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
#m <- get_map(bbox(extent(X94)*1.1), source="google", zoom=12, maptype="satellite") #define the box for map based on the trajectory coordinates
#ggmap(m)+geom_path(data=X94_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
#' ## Display space-time cube. 
open3d()
# To get a bigger window than the default
par3d(windowRect = c(100, 100, 612, 612))
Sys.sleep(0.1) # Allow sluggish window managers to catch up
with(X94_df, plot3d(x,y,date, type="l", col=as.numeric(X94_df$Event_ID)))
(stcube<-with(X94_df, plot3d(x,y,date, type="l",col=as.numeric(cut(X94_df$dist,3)), alpha=0.4)))
#Timelag
head(timeLag(X94, units="mins"))
head(timeLag(X94, units="hours"))
median(timeLag(X94, units="hours"))
max(timeLag(X94, units="hours"))
max(timeLag(X94, units="days"))
plot(X94_df$individual.local.identifier..ID.~X94_df $date,xlab = "Sampling period",ylab = "Individuals ID Number")
#write.table(X94,file="c:/RWorkDir/jaguardatapaper/ID/X94.txt",row.names = F,quote=F,col.names=T,sep="\t")
#Histograms and plots
hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 94],
     main = 'Individual 94', xlab = "Time between fixes (h)")
plot(density.default(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 94],
     na.rm = T), main = 'Individual 94')
	 
breaks <- c(0, 0.5, 1.5, 2.5, 4.5, 8.5, 16.5, 32.5, 64.5, 128.5, 256.5, 1000)
breaks <- c(0, 1, 2*(1:256))
hh <- hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 94 &
                                     mov.data.diff$time.diff > 0 & 
                                     mov.data.diff$time.diff < 512],
									 main = 'Individual 94',xlab = 'Time between fixes (h)',
           plot = T, breaks = breaks)

#X95
X95=unstacked$X95
X95_df <- as(X95, "data.frame")
head(X95_df)
nrow(X95_df)
range(X95_df$date)
m <- get_map(bbox(extent(X95)*1.1), source="google", zoom=12, maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=X95_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
m <- get_map(bbox(extent(X95)*1.1), source="google", zoom=13, maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=X95_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
#' ## Display space-time cube. 
open3d()
# To get a bigger window than the default
par3d(windowRect = c(100, 100, 612, 612))
Sys.sleep(0.1) # Allow sluggish window managers to catch up
with(X95_df, plot3d(x,y,date, type="l", col=as.numeric(X95_df$Event_ID)))
(stcube<-with(X95_df, plot3d(x,y,date, type="l",col=as.numeric(cut(X95_df$dist,3)), alpha=0.4)))
#Timelag
head(timeLag(X95, units="mins"))
head(timeLag(X95, units="hours"))
median(timeLag(X95, units="hours"))
max(timeLag(X95, units="hours"))
max(timeLag(X95, units="days"))
plot(X95_df$individual.local.identifier..ID.~X95_df $date,xlab = "Sampling period",ylab = "Individuals ID Number")
#write.table(X95,file="c:/RWorkDir/jaguardatapaper/ID/X95.txt",row.names = F,quote=F,col.names=T,sep="\t")
#Histograms and plots
hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 95],
     main = 'Individual 95', xlab = "Time between fixes (h)")
plot(density.default(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 95],
     na.rm = T), main = 'Individual 95')
	 
breaks <- c(0, 0.5, 1.5, 2.5, 4.5, 8.5, 16.5, 32.5, 64.5, 128.5, 256.5, 1000)
breaks <- c(0, 1, 2*(1:256))
hh <- hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 95 &
                                     mov.data.diff$time.diff > 0 & 
                                     mov.data.diff$time.diff < 512],
									 main = 'Individual 95',xlab = 'Time between fixes (h)',
           plot = T, breaks = breaks)

#X96
X96=unstacked$X96
X96_df <- as(X96, "data.frame")
head(X96_df)
nrow(X96_df)
range(X96_df$date)
m <- get_map(bbox(extent(X96)*1.1), source="google", zoom=12, maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=X96_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
#m <- get_map(bbox(extent(X96)*1.1), source="google", zoom=13, maptype="satellite") #define the box for map based on the trajectory coordinates
#ggmap(m)+geom_path(data=X96_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
#' ## Display space-time cube. 
open3d()
# To get a bigger window than the default
par3d(windowRect = c(100, 100, 612, 612))
Sys.sleep(0.1) # Allow sluggish window managers to catch up
with(X96_df, plot3d(x,y,date, type="l", col=as.numeric(X96_df$Event_ID)))
(stcube<-with(X96_df, plot3d(x,y,date, type="l",col=as.numeric(cut(X96_df$dist,3)), alpha=0.4)))
#Timelag
head(timeLag(X96, units="mins"))
head(timeLag(X96, units="hours"))
median(timeLag(X96, units="hours"))
max(timeLag(X96, units="hours"))
max(timeLag(X96, units="days"))
plot(X96_df$individual.local.identifier..ID.~X96_df $date,xlab = "Sampling period",ylab = "Individuals ID Number")
#write.table(X96,file="c:/RWorkDir/jaguardatapaper/ID/X96.txt",row.names = F,quote=F,col.names=T,sep="\t")
#Histograms and plots
hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 96],
     main = 'Individual 96', xlab = "Time between fixes (h)")
plot(density.default(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 96],
     na.rm = T), main = 'Individual 96')
	 
breaks <- c(0, 0.5, 1.5, 2.5, 4.5, 8.5, 16.5, 32.5, 64.5, 128.5, 256.5, 1000)
breaks <- c(0, 1, 2*(1:256))
hh <- hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 96 &
                                     mov.data.diff$time.diff > 0 & 
                                     mov.data.diff$time.diff < 512],
									 main = 'Individual 96',xlab = 'Time between fixes (h)',
           plot = T, breaks = breaks)

#X97
X97=unstacked$X97
X97_df <- as(X97, "data.frame")
head(X97_df)
nrow(X97_df)
range(X97_df$date)
m <- get_map(bbox(extent(X97)*1.1), source="google", zoom=12, maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=X97_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
#' ## Display space-time cube. 
open3d()
# To get a bigger window than the default
par3d(windowRect = c(100, 100, 612, 612))
Sys.sleep(0.1) # Allow sluggish window managers to catch up
with(X97_df, plot3d(x,y,date, type="l", col=as.numeric(X97_df$Event_ID)))
(stcube<-with(X97_df, plot3d(x,y,date, type="l",col=as.numeric(cut(X97_df$dist,3)), alpha=0.4)))
#Timelag
head(timeLag(X97, units="mins"))
head(timeLag(X97, units="hours"))
median(timeLag(X97, units="hours"))
max(timeLag(X97, units="hours"))
max(timeLag(X97, units="days"))
plot(X97_df$individual.local.identifier..ID.~X97_df $date,xlab = "Sampling period",ylab = "Individuals ID Number")
#write.table(X97,file="c:/RWorkDir/jaguardatapaper/ID/X97.txt",row.names = F,quote=F,col.names=T,sep="\t")
#Histograms and plots
hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 97],
     main = 'Individual 97', xlab = "Time between fixes (h)")
plot(density.default(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 97],
     na.rm = T), main = 'Individual 97')
	 
breaks <- c(0, 0.5, 1.5, 2.5, 4.5, 8.5, 16.5, 32.5, 64.5, 128.5, 256.5, 1000)
breaks <- c(0, 1, 2*(1:256))
hh <- hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 97 &
                                     mov.data.diff$time.diff > 0 & 
                                     mov.data.diff$time.diff < 512],
									 main = 'Individual 97',xlab = 'Time between fixes (h)',
           plot = T, breaks = breaks)

#X98
X98=unstacked$X98
X98_df <- as(X98, "data.frame")
head(X98_df)
nrow(X98_df)
range(X98_df$date)
m <- get_map(bbox(extent(X98)*1.1), source="google", zoom=13, maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=X98_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
m <- get_map(bbox(extent(X98)*1.1), source="google", zoom=12, maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=X98_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
#' ## Display space-time cube. 
open3d()
# To get a bigger window than the default
par3d(windowRect = c(100, 100, 612, 612))
Sys.sleep(0.1) # Allow sluggish window managers to catch up
with(X98_df, plot3d(x,y,date, type="l", col=as.numeric(X98_df$Event_ID)))
(stcube<-with(X98_df, plot3d(x,y,date, type="l",col=as.numeric(cut(X98_df$dist,3)), alpha=0.4)))
#Timelag
head(timeLag(X98, units="mins"))
head(timeLag(X98, units="hours"))
median(timeLag(X98, units="hours"))
max(timeLag(X98, units="hours"))
max(timeLag(X98, units="days"))
plot(X98_df$individual.local.identifier..ID.~X98_df $date,xlab = "Sampling period",ylab = "Individuals ID Number")
#write.table(X98,file="c:/RWorkDir/jaguardatapaper/ID/X98.txt",row.names = F,quote=F,col.names=T,sep="\t")
#Histograms and plots
hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 98],
     main = 'Individual 98', xlab = "Time between fixes (h)")
plot(density.default(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 98],
     na.rm = T), main = 'Individual 98')
	 
breaks <- c(0, 0.5, 1.5, 2.5, 4.5, 8.5, 16.5, 32.5, 64.5, 128.5, 256.5, 1000)
breaks <- c(0, 1, 2*(1:256))
hh <- hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 98 &
                                     mov.data.diff$time.diff > 0 & 
                                     mov.data.diff$time.diff < 512],
									 main = 'Individual 98',xlab = 'Time between fixes (h)',
           plot = T, breaks = breaks)

#X99
X99=unstacked$X99
X99_df <- as(X99, "data.frame")
head(X99_df)
nrow(X99_df)
range(X99_df$date)
m <- get_map(bbox(extent(X99)*1.1), source="google", zoom=12, maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=X99_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
#' ## Display space-time cube. 
open3d()
# To get a bigger window than the default
par3d(windowRect = c(100, 100, 612, 612))
Sys.sleep(0.1) # Allow sluggish window managers to catch up
with(X99_df, plot3d(x,y,date, type="l", col=as.numeric(X99_df$Event_ID)))
(stcube<-with(X99_df, plot3d(x,y,date, type="l",col=as.numeric(cut(X99_df$dist,3)), alpha=0.4)))
#Timelag
head(timeLag(X99, units="mins"))
head(timeLag(X99, units="hours"))
median(timeLag(X99, units="hours"))
max(timeLag(X99, units="hours"))
max(timeLag(X99, units="days"))
plot(X99_df$individual.local.identifier..ID.~X99_df $date,xlab = "Sampling period",ylab = "Individuals ID Number")
#write.table(X99,file="c:/RWorkDir/jaguardatapaper/ID/X99.txt",row.names = F,quote=F,col.names=T,sep="\t")
#Histograms and plots
hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 99],
     main = 'Individual 99', xlab = "Time between fixes (h)")
plot(density.default(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 99],
     na.rm = T), main = 'Individual 99')
	 
breaks <- c(0, 0.5, 1.5, 2.5, 4.5, 8.5, 16.5, 32.5, 64.5, 128.5, 256.5, 1000)
breaks <- c(0, 1, 2*(1:256))
hh <- hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 99 &
                                     mov.data.diff$time.diff > 0 & 
                                     mov.data.diff$time.diff < 512],
									 main = 'Individual 99',xlab = 'Time between fixes (h)',
           plot = T, breaks = breaks)

#X100
X100=unstacked$X100
X100_df <- as(X100, "data.frame")
head(X100_df)
nrow(X100_df)
range(X100_df$date)
m <- get_map(bbox(extent(X100)*1.1), source="google", zoom=12, maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=X100_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
#' ## Display space-time cube. 
open3d()
# To get a bigger window than the default
par3d(windowRect = c(100, 100, 612, 612))
Sys.sleep(0.1) # Allow sluggish window managers to catch up
with(X100_df, plot3d(x,y,date, type="l", col=as.numeric(X100_df$Event_ID)))
(stcube<-with(X100_df, plot3d(x,y,date, type="l",col=as.numeric(cut(X100_df$dist,3)), alpha=0.4)))
#Timelag
head(timeLag(X100, units="mins"))
head(timeLag(X100, units="hours"))
median(timeLag(X100, units="hours"))
max(timeLag(X100, units="hours"))
max(timeLag(X100, units="days"))
plot(X100_df$individual.local.identifier..ID.~X100_df $date,xlab = "Sampling period",ylab = "Individuals ID Number")
#write.table(X100,file="c:/RWorkDir/jaguardatapaper/ID/X100.txt",row.names = F,quote=F,col.names=T,sep="\t")
#Histograms and plots
hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 100],
     main = 'Individual 100', xlab = "Time between fixes (h)")
plot(density.default(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 100],
     na.rm = T), main = 'Individual 100')
	 
breaks <- c(0, 0.5, 1.5, 2.5, 4.5, 8.5, 16.5, 32.5, 64.5, 128.5, 256.5, 1000)
breaks <- c(0, 1, 2*(1:256))
hh <- hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 100 &
                                     mov.data.diff$time.diff > 0 & 
                                     mov.data.diff$time.diff < 512],
									 main = 'Individual 100',xlab = 'Time between fixes (h)',
           plot = T, breaks = breaks)

#X101
X101=unstacked$X101
X101_df <- as(X101, "data.frame")
head(X101_df)
nrow(X101_df)
range(X101_df$date)
m <- get_map(bbox(extent(X101)*1.1), source="google", zoom=11, maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=X101_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
#m <- get_map(bbox(extent(X101)*1.1), source="google", zoom=12, maptype="satellite") #define the box for map based on the trajectory coordinates
#ggmap(m)+geom_path(data=X101_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
#' ## Display space-time cube. 
open3d()
# To get a bigger window than the default
par3d(windowRect = c(100, 100, 612, 612))
Sys.sleep(0.1) # Allow sluggish window managers to catch up
with(X101_df, plot3d(x,y,date, type="l", col=as.numeric(X101_df$Event_ID)))
(stcube<-with(X101_df, plot3d(x,y,date, type="l",col=as.numeric(cut(X101_df$dist,3)), alpha=0.4)))
#Timelag
head(timeLag(X101, units="mins"))
head(timeLag(X101, units="hours"))
median(timeLag(X101, units="hours"))
max(timeLag(X101, units="hours"))
max(timeLag(X101, units="days"))
plot(X101_df$individual.local.identifier..ID.~X101_df $date,xlab = "Sampling period",ylab = "Individuals ID Number")
#write.table(X101,file="c:/RWorkDir/jaguardatapaper/ID/X101.txt",row.names = F,quote=F,col.names=T,sep="\t")
#Histograms and plots
hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 101],
     main = 'Individual 101', xlab = "Time between fixes (h)")
plot(density.default(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 101],
     na.rm = T), main = 'Individual 101')
	 
breaks <- c(0, 0.5, 1.5, 2.5, 4.5, 8.5, 16.5, 32.5, 64.5, 128.5, 256.5, 1000)
breaks <- c(0, 1, 2*(1:256))
hh <- hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 101 &
                                     mov.data.diff$time.diff > 0 & 
                                     mov.data.diff$time.diff < 512],
									 main = 'Individual 101',xlab = 'Time between fixes (h)',
           plot = T, breaks = breaks)
#X102
X102=unstacked$X102
X102_df <- as(X102, "data.frame")
head(X102_df)
nrow(X102_df)
range(X102_df$date)
m <- get_map(bbox(extent(X102)*1.1), source="google", zoom=11, maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=X102_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
#' ## Display space-time cube. 
open3d()
# To get a bigger window than the default
par3d(windowRect = c(100, 100, 612, 612))
Sys.sleep(0.1) # Allow sluggish window managers to catch up
with(X102_df, plot3d(x,y,date, type="l", col=as.numeric(X102_df$Event_ID)))
(stcube<-with(X102_df, plot3d(x,y,date, type="l",col=as.numeric(cut(X102_df$dist,3)), alpha=0.4)))
#Timelag
head(timeLag(X102, units="mins"))
head(timeLag(X102, units="hours"))
median(timeLag(X102, units="hours"))
max(timeLag(X102, units="hours"))
max(timeLag(X102, units="days"))
plot(X102_df$individual.local.identifier..ID.~X102_df $date,xlab = "Sampling period",ylab = "Individuals ID Number")
#write.table(X102,file="c:/RWorkDir/jaguardatapaper/ID/X102.txt",row.names = F,quote=F,col.names=T,sep="\t")
#Histograms and plots
hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 102],
     main = 'Individual 102', xlab = "Time between fixes (h)")
plot(density.default(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 102],
     na.rm = T), main = 'Individual 102')
	 
breaks <- c(0, 0.5, 1.5, 2.5, 4.5, 8.5, 16.5, 32.5, 64.5, 128.5, 256.5, 1000)
breaks <- c(0, 1, 2*(1:256))
hh <- hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 102 &
                                     mov.data.diff$time.diff > 0 & 
                                     mov.data.diff$time.diff < 512],
									 main = 'Individual 102',xlab = 'Time between fixes (h)',
           plot = T, breaks = breaks)
#X103
X103=unstacked$X103
X103_df <- as(X103, "data.frame")
head(X103_df)
nrow(X103_df)
range(X103_df$date)
m <- get_map(bbox(extent(X103)*1.1), source="google", zoom=13, maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=X103_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
m <- get_map(bbox(extent(X103)*1.1), source="google", zoom=12, maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=X103_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
#' ## Display space-time cube. 
open3d()
# To get a bigger window than the default
par3d(windowRect = c(100, 100, 612, 612))
Sys.sleep(0.1) # Allow sluggish window managers to catch up
with(X103_df, plot3d(x,y,date, type="l", col=as.numeric(X103_df$Event_ID)))
(stcube<-with(X103_df, plot3d(x,y,date, type="l",col=as.numeric(cut(X103_df$dist,3)), alpha=0.4)))
#Timelag
head(timeLag(X103, units="mins"))
head(timeLag(X103, units="hours"))
median(timeLag(X103, units="hours"))
max(timeLag(X103, units="hours"))
max(timeLag(X103, units="days"))
plot(X103_df$individual.local.identifier..ID.~X103_df $date,xlab = "Sampling period",ylab = "Individuals ID Number")
#write.table(X103,file="c:/RWorkDir/jaguardatapaper/ID/X103.txt",row.names = F,quote=F,col.names=T,sep="\t")
#Histograms and plots
hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 103],
     main = 'Individual 103', xlab = "Time between fixes (h)")
plot(density.default(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 103],
     na.rm = T), main = 'Individual 103')
	 
breaks <- c(0, 0.5, 1.5, 2.5, 4.5, 8.5, 16.5, 32.5, 64.5, 128.5, 256.5, 1000)
breaks <- c(0, 1, 2*(1:256))
hh <- hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 103 &
                                     mov.data.diff$time.diff > 0 & 
                                     mov.data.diff$time.diff < 512],
									 main = 'Individual 103',xlab = 'Time between fixes (h)',
           plot = T, breaks = breaks)

#X104
X104=unstacked$X104
X104_df <- as(X104, "data.frame")
head(X104_df)
nrow(X104_df)
range(X104_df$date)
m <- get_map(bbox(extent(X104)*1.1), source="google", zoom=11, maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=X104_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
#' ## Display space-time cube. 
open3d()
# To get a bigger window than the default
par3d(windowRect = c(100, 100, 612, 612))
Sys.sleep(0.1) # Allow sluggish window managers to catch up
with(X104_df, plot3d(x,y,date, type="l", col=as.numeric(X104_df$Event_ID)))
(stcube<-with(X104_df, plot3d(x,y,date, type="l",col=as.numeric(cut(X104_df$dist,3)), alpha=0.4)))
#Timelag
head(timeLag(X104, units="mins"))
head(timeLag(X104, units="hours"))
median(timeLag(X104, units="hours"))
max(timeLag(X104, units="hours"))
max(timeLag(X104, units="days"))
plot(X104_df$individual.local.identifier..ID.~X104_df $date,xlab = "Sampling period",ylab = "Individuals ID Number")
#write.table(X104,file="c:/RWorkDir/jaguardatapaper/ID/X104.txt",row.names = F,quote=F,col.names=T,sep="\t")
#Histograms and plots
hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 104],
     main = 'Individual 104', xlab = "Time between fixes (h)")
plot(density.default(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 104],
     na.rm = T), main = 'Individual 104')
	 
breaks <- c(0, 0.5, 1.5, 2.5, 4.5, 8.5, 16.5, 32.5, 64.5, 128.5, 256.5, 1000)
breaks <- c(0, 1, 2*(1:256))
hh <- hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 104 &
                                     mov.data.diff$time.diff > 0 & 
                                     mov.data.diff$time.diff < 512],
									 main = 'Individual 104',xlab = 'Time between fixes (h)',
           plot = T, breaks = breaks)

#X105
X105=unstacked$X105
X105_df <- as(X105, "data.frame")
head(X105_df)
nrow(X105_df)
range(X105_df$date)
m <- get_map(bbox(extent(X105)*1.1), source="google", zoom=12, maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=X105_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
#' ## Display space-time cube. 
open3d()
# To get a bigger window than the default
par3d(windowRect = c(100, 100, 612, 612))
Sys.sleep(0.1) # Allow sluggish window managers to catch up
with(X105_df, plot3d(x,y,date, type="l", col=as.numeric(X105_df$Event_ID)))
(stcube<-with(X105_df, plot3d(x,y,date, type="l",col=as.numeric(cut(X105_df$dist,3)), alpha=0.4)))
#Timelag
head(timeLag(X105, units="mins"))
head(timeLag(X105, units="hours"))
median(timeLag(X105, units="hours"))
max(timeLag(X105, units="hours"))
max(timeLag(X105, units="days"))
plot(X105_df$individual.local.identifier..ID.~X105_df $date,xlab = "Sampling period",ylab = "Individuals ID Number")
#write.table(X105,file="c:/RWorkDir/jaguardatapaper/ID/X105.txt",row.names = F,quote=F,col.names=T,sep="\t")
#Histograms and plots
hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 105],
     main = 'Individual 105', xlab = "Time between fixes (h)")
plot(density.default(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 105],
     na.rm = T), main = 'Individual 105')
	 
breaks <- c(0, 0.5, 1.5, 2.5, 4.5, 8.5, 16.5, 32.5, 64.5, 128.5, 256.5, 1000)
breaks <- c(0, 1, 2*(1:256))
hh <- hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 105 &
                                     mov.data.diff$time.diff > 0 & 
                                     mov.data.diff$time.diff < 512],
									 main = 'Individual 105',xlab = 'Time between fixes (h)',
           plot = T, breaks = breaks)

#X106
X106=unstacked$X106
X106_df <- as(X106, "data.frame")
head(X106_df)
nrow(X106_df)
range(X106_df$date)
m <- get_map(bbox(extent(X106)*1.1), source="google", zoom=12, maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=X106_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
#' ## Display space-time cube. 
open3d()
# To get a bigger window than the default
par3d(windowRect = c(100, 100, 612, 612))
Sys.sleep(0.1) # Allow sluggish window managers to catch up
with(XX106_df, plot3d(x,y,date, type="l", col=as.numeric(XX106_df$Event_ID)))
(stcube<-with(XX106_df, plot3d(x,y,date, type="l",col=as.numeric(cut(XX106_df$dist,3)), alpha=0.4)))
#Timelag
head(timeLag(X106, units="mins"))
head(timeLag(X106, units="hours"))
median(timeLag(X106, units="hours"))
max(timeLag(X106, units="hours"))
max(timeLag(X106, units="days"))
plot(X106_df$individual.local.identifier..ID.~X106_df $date,xlab = "Sampling period",ylab = "Individuals ID Number")
#write.table(X106,file="c:/RWorkDir/jaguardatapaper/ID/X106.txt",row.names = F,quote=F,col.names=T,sep="\t")
#Histograms and plots
hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 106],
     main = 'Individual 106', xlab = "Time between fixes (h)")
plot(density.default(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 106],
     na.rm = T), main = 'Individual 106')
	 
breaks <- c(0, 0.5, 1.5, 2.5, 4.5, 8.5, 16.5, 32.5, 64.5, 128.5, 256.5, 1000)
breaks <- c(0, 1, 2*(1:256))
hh <- hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 106 &
                                     mov.data.diff$time.diff > 0 & 
                                     mov.data.diff$time.diff < 512],
									 main = 'Individual 106',xlab = 'Time between fixes (h)',
           plot = T, breaks = breaks)
#X107
X107=unstacked$X107
X107_df <- as(X107, "data.frame")
head(X107_df)
nrow(X107_df)
range(X107_df$date)
m <- get_map(bbox(extent(X107)*1.1), source="google", zoom=13, maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=X107_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
m <- get_map(bbox(extent(X107)*1.1), source="google", zoom=12, maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=X107_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
#' ## Display space-time cube. 
open3d()
# To get a bigger window than the default
par3d(windowRect = c(100, 100, 612, 612))
Sys.sleep(0.1) # Allow sluggish window managers to catch up
with(X107_df, plot3d(x,y,date, type="l", col=as.numeric(X107_df$Event_ID)))
(stcube<-with(X107_df, plot3d(x,y,date, type="l",col=as.numeric(cut(X107_df$dist,3)), alpha=0.4)))
#Timelag
head(timeLag(X107, units="mins"))
head(timeLag(X107, units="hours"))
median(timeLag(X107, units="hours"))
max(timeLag(X107, units="hours"))
max(timeLag(X107, units="days"))
plot(X107_df$individual.local.identifier..ID.~X107_df $date,xlab = "Sampling period",ylab = "Individuals ID Number")
#write.table(X107,file="c:/RWorkDir/jaguardatapaper/ID/X107.txt",row.names = F,quote=F,col.names=T,sep="\t")
#Histograms and plots
hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 107],
     main = 'Individual 107', xlab = "Time between fixes (h)")
plot(density.default(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 107],
     na.rm = T), main = 'Individual 107')
	 
breaks <- c(0, 0.5, 1.5, 2.5, 4.5, 8.5, 16.5, 32.5, 64.5, 128.5, 256.5, 1000)
breaks <- c(0, 1, 2*(1:256))
hh <- hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 107 &
                                     mov.data.diff$time.diff > 0 & 
                                     mov.data.diff$time.diff < 512],
									 main = 'Individual 107',xlab = 'Time between fixes (h)',
           plot = T, breaks = breaks)

#X108
X108=unstacked$X108
X108_df <- as(X108, "data.frame")
head(X108_df)
nrow(X108_df)
range(X108_df$date)
m <- get_map(bbox(extent(X108)*1.1), source="google", zoom=12, maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=X108_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
#' ## Display space-time cube. 
open3d()
# To get a bigger window than the default
par3d(windowRect = c(100, 100, 612, 612))
Sys.sleep(0.1) # Allow sluggish window managers to catch up
with(X108_df, plot3d(x,y,date, type="l", col=as.numeric(X108_df$Event_ID)))
(stcube<-with(X108_df, plot3d(x,y,date, type="l",col=as.numeric(cut(X108_df$dist,3)), alpha=0.4)))
#Timelag
head(timeLag(X108, units="mins"))
head(timeLag(X108, units="hours"))
median(timeLag(X108, units="hours"))
max(timeLag(X108, units="hours"))
max(timeLag(X108, units="days"))
plot(X108_df$individual.local.identifier..ID.~X108_df $date,xlab = "Sampling period",ylab = "Individuals ID Number")
#write.table(X108,file="c:/RWorkDir/jaguardatapaper/ID/X108.txt",row.names = F,quote=F,col.names=T,sep="\t")
#Histograms and plots
hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 108],
     main = 'Individual 108', xlab = "Time between fixes (h)")
plot(density.default(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 108],
     na.rm = T), main = 'Individual 108')
	 
breaks <- c(0, 0.5, 1.5, 2.5, 4.5, 8.5, 16.5, 32.5, 64.5, 128.5, 256.5, 1000)
breaks <- c(0, 1, 2*(1:256))
hh <- hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 108 &
                                     mov.data.diff$time.diff > 0 & 
                                     mov.data.diff$time.diff < 512],
									 main = 'Individual 108',xlab = 'Time between fixes (h)',
           plot = T, breaks = breaks)
#X109
X109=unstacked$X109
X109_df <- as(X109, "data.frame")
head(X109_df)
nrow(X109_df)
range(X109_df$date)
m <- get_map(bbox(extent(X109)*1.1), source="google", zoom=13, maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=X109_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
m <- get_map(bbox(extent(X109)*1.1), source="google", zoom=12, maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=X109_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
#' ## Display space-time cube. 
open3d()
# To get a bigger window than the default
par3d(windowRect = c(100, 100, 612, 612))
Sys.sleep(0.1) # Allow sluggish window managers to catch up
with(X109_df, plot3d(x,y,date, type="l", col=as.numeric(X109_df$Event_ID)))
(stcube<-with(X109_df, plot3d(x,y,date, type="l",col=as.numeric(cut(X109_df$dist,3)), alpha=0.4)))
#Timelag
head(timeLag(X109, units="mins"))
head(timeLag(X109, units="hours"))
median(timeLag(X109, units="hours"))
max(timeLag(X109, units="hours"))
max(timeLag(X109, units="days"))
plot(X109_df$individual.local.identifier..ID.~X109_df $date,xlab = "Sampling period",ylab = "Individuals ID Number")
#write.table(X109,file="c:/RWorkDir/jaguardatapaper/ID/X109.txt",row.names = F,quote=F,col.names=T,sep="\t")
#Histograms and plots
hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 109],
     main = 'Individual 109', xlab = "Time between fixes (h)")
plot(density.default(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 109],
     na.rm = T), main = 'Individual 109')
	 
breaks <- c(0, 0.5, 1.5, 2.5, 4.5, 8.5, 16.5, 32.5, 64.5, 128.5, 256.5, 1000)
breaks <- c(0, 1, 2*(1:256))
hh <- hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 109 &
                                     mov.data.diff$time.diff > 0 & 
                                     mov.data.diff$time.diff < 512],
									 main = 'Individual 109',xlab = 'Time between fixes (h)',
           plot = T, breaks = breaks)
#X110
X110=unstacked$X110
X110_df <- as(X110, "data.frame")
head(X110_df)
nrow(X110_df)
range(X110_df$date)
m <- get_map(bbox(extent(X110)*1.1), source="google", zoom=12, maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=X110_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
#' ## Display space-time cube. 
open3d()
# To get a bigger window than the default
par3d(windowRect = c(100, 100, 612, 612))
Sys.sleep(0.1) # Allow sluggish window managers to catch up
with(X110_df, plot3d(x,y,date, type="l", col=as.numeric(X110_df$Event_ID)))
(stcube<-with(X110_df, plot3d(x,y,date, type="l",col=as.numeric(cut(X110_df$dist,3)), alpha=0.4)))
#Timelag
head(timeLag(X110, units="mins"))
head(timeLag(X110, units="hours"))
median(timeLag(X110, units="hours"))
max(timeLag(X110, units="hours"))
max(timeLag(X110, units="days"))
plot(X110_df$individual.local.identifier..ID.~X110_df $date,xlab = "Sampling period",ylab = "Individuals ID Number")
#write.table(X110,file="c:/RWorkDir/jaguardatapaper/ID/X110.txt",row.names = F,quote=F,col.names=T,sep="\t")
#Histograms and plots
hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 110],
     main = 'Individual 110', xlab = "Time between fixes (h)")
plot(density.default(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 110],
     na.rm = T), main = 'Individual 110')
	 
breaks <- c(0, 0.5, 1.5, 2.5, 4.5, 8.5, 16.5, 32.5, 64.5, 128.5, 256.5, 1000)
breaks <- c(0, 1, 2*(1:256))
hh <- hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 110 &
                                     mov.data.diff$time.diff > 0 & 
                                     mov.data.diff$time.diff < 512],
									 main = 'Individual 110',xlab = 'Time between fixes (h)',
           plot = T, breaks = breaks)
#X111
X111=unstacked$X111
X111_df <- as(X111, "data.frame")
head(X111_df)
nrow(X111_df)
range(X111_df$date)
m <- get_map(bbox(extent(X111)*1.1), source="google", zoom=13, maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=X111_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
m <- get_map(bbox(extent(X111)*1.1), source="google", zoom=12, maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=X111_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
#' ## Display space-time cube. 
open3d()
# To get a bigger window than the default
par3d(windowRect = c(100, 100, 612, 612))
Sys.sleep(0.1) # Allow sluggish window managers to catch up
with(X111_df, plot3d(x,y,date, type="l", col=as.numeric(X111_df$Event_ID)))
(stcube<-with(X111_df, plot3d(x,y,date, type="l",col=as.numeric(cut(X111_df$dist,3)), alpha=0.4)))
#Timelag
head(timeLag(X111, units="mins"))
head(timeLag(X111, units="hours"))
median(timeLag(X111, units="hours"))
max(timeLag(X111, units="hours"))
max(timeLag(X111, units="days"))
plot(X111_df$individual.local.identifier..ID.~X111_df $date,xlab = "Sampling period",ylab = "Individuals ID Number")
#write.table(X111,file="c:/RWorkDir/jaguardatapaper/ID/X111.txt",row.names = F,quote=F,col.names=T,sep="\t")
#Histograms and plots
hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 111],
     main = 'Individual 111', xlab = "Time between fixes (h)")
plot(density.default(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 111],
     na.rm = T), main = 'Individual 111')
	 
breaks <- c(0, 0.5, 1.5, 2.5, 4.5, 8.5, 16.5, 32.5, 64.5, 128.5, 256.5, 1000)
breaks <- c(0, 1, 2*(1:256))
hh <- hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 111 &
                                     mov.data.diff$time.diff > 0 & 
                                     mov.data.diff$time.diff < 512],
									 main = 'Individual 111',xlab = 'Time between fixes (h)',
           plot = T, breaks = breaks)

#X112
X112=unstacked$X112
X112_df <- as(X112, "data.frame")
head(X112_df)
nrow(X112_df)
range(X112_df$date)
m <- get_map(bbox(extent(X112)*1.1), source="google", zoom=13, maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=X112_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
m <- get_map(bbox(extent(X112)*1.1), source="google", zoom=12, maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=X112_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
#' ## Display space-time cube. 
open3d()
# To get a bigger window than the default
par3d(windowRect = c(100, 100, 612, 612))
Sys.sleep(0.1) # Allow sluggish window managers to catch up
with(X112_df, plot3d(x,y,date, type="l", col=as.numeric(X112_df$Event_ID)))
(stcube<-with(X112_df, plot3d(x,y,date, type="l",col=as.numeric(cut(X112_df$dist,3)), alpha=0.4)))
#Timelag
head(timeLag(X112, units="mins"))
head(timeLag(X112, units="hours"))
median(timeLag(X112, units="hours"))
max(timeLag(X112, units="hours"))
max(timeLag(X112, units="days"))
plot(X112_df$individual.local.identifier..ID.~X112_df $date,xlab = "Sampling period",ylab = "Individuals ID Number")
#write.table(X112,file="c:/RWorkDir/jaguardatapaper/ID/X112.txt",row.names = F,quote=F,col.names=T,sep="\t")
#Histograms and plots
hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 112],
     main = 'Individual 112', xlab = "Time between fixes (h)")
plot(density.default(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 112],
     na.rm = T), main = 'Individual 112')
	 
breaks <- c(0, 0.5, 1.5, 2.5, 4.5, 8.5, 16.5, 32.5, 64.5, 128.5, 256.5, 1000)
breaks <- c(0, 1, 2*(1:256))
hh <- hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 112 &
                                     mov.data.diff$time.diff > 0 & 
                                     mov.data.diff$time.diff < 512],
									 main = 'Individual 112',xlab = 'Time between fixes (h)',
           plot = T, breaks = breaks)

#X113
X113=unstacked$X113
X113_df <- as(X113, "data.frame")
head(X113_df)
nrow(X113_df)
range(X113_df$date)
m <- get_map(bbox(extent(X113)*1.1), source="google", zoom=13, maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=X113_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
m <- get_map(bbox(extent(X113)*1.1), source="google", zoom=12, maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=X113_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
#' ## Display space-time cube. 
open3d()
# To get a bigger window than the default
par3d(windowRect = c(100, 100, 612, 612))
Sys.sleep(0.1) # Allow sluggish window managers to catch up
with(X113_df, plot3d(x,y,date, type="l", col=as.numeric(X113_df$Event_ID)))
(stcube<-with(X113_df, plot3d(x,y,date, type="l",col=as.numeric(cut(X113_df$dist,3)), alpha=0.4)))
#Timelag
head(timeLag(X113, units="mins"))
head(timeLag(X113, units="hours"))
median(timeLag(X113, units="hours"))
max(timeLag(X113, units="hours"))
max(timeLag(X113, units="days"))
plot(X113_df$individual.local.identifier..ID.~X113_df $date,xlab = "Sampling period",ylab = "Individuals ID Number")
#write.table(X113,file="c:/RWorkDir/jaguardatapaper/ID/X113.txt",row.names = F,quote=F,col.names=T,sep="\t")
#Histograms and plots
hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 113],
     main = 'Individual 113', xlab = "Time between fixes (h)")
plot(density.default(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 113],
     na.rm = T), main = 'Individual 113')
	 
breaks <- c(0, 0.5, 1.5, 2.5, 4.5, 8.5, 16.5, 32.5, 64.5, 128.5, 256.5, 1000)
breaks <- c(0, 1, 2*(1:256))
hh <- hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 113 &
                                     mov.data.diff$time.diff > 0 & 
                                     mov.data.diff$time.diff < 512],
									 main = 'Individual 113',xlab = 'Time between fixes (h)',
           plot = T, breaks = breaks)

#X114
X114=unstacked$X114
X114_df <- as(X114, "data.frame")
head(X114_df)
nrow(X114_df)
range(X114_df$date)
m <- get_map(bbox(extent(X114)*1.1), source="google", zoom=12, maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=X114_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
#' ## Display space-time cube. 
open3d()
# To get a bigger window than the default
par3d(windowRect = c(100, 100, 612, 612))
Sys.sleep(0.1) # Allow sluggish window managers to catch up
with(X114_df, plot3d(x,y,date, type="l", col=as.numeric(X114_df$Event_ID)))
(stcube<-with(X114_df, plot3d(x,y,date, type="l",col=as.numeric(cut(X114_df$dist,3)), alpha=0.4)))
#Timelag
head(timeLag(X114, units="mins"))
head(timeLag(X114, units="hours"))
median(timeLag(X114, units="hours"))
max(timeLag(X114, units="hours"))
max(timeLag(X114, units="days"))
plot(X114_df$individual.local.identifier..ID.~X114_df $date,xlab = "Sampling period",ylab = "Individuals ID Number")
#write.table(X114,file="c:/RWorkDir/jaguardatapaper/ID/X114.txt",row.names = F,quote=F,col.names=T,sep="\t")
#Histograms and plots
hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 114],
     main = 'Individual 114', xlab = "Time between fixes (h)")
plot(density.default(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 114],
     na.rm = T), main = 'Individual 114')
	 
breaks <- c(0, 0.5, 1.5, 2.5, 4.5, 8.5, 16.5, 32.5, 64.5, 128.5, 256.5, 1000)
breaks <- c(0, 1, 2*(1:256))
hh <- hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 114 &
                                     mov.data.diff$time.diff > 0 & 
                                     mov.data.diff$time.diff < 512],
									 main = 'Individual 114',xlab = 'Time between fixes (h)',
           plot = T, breaks = breaks)

#X115
X115=unstacked$X115
X115_df <- as(X115, "data.frame")
head(X115_df)
nrow(X115_df)
range(X115_df$date)
m <- get_map(bbox(extent(X115)*1.1), source="google", zoom=12, maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=X115_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
#' ## Display space-time cube. 
open3d()
# To get a bigger window than the default
par3d(windowRect = c(100, 100, 612, 612))
Sys.sleep(0.1) # Allow sluggish window managers to catch up
with(X115_df, plot3d(x,y,date, type="l", col=as.numeric(X115_df$Event_ID)))
(stcube<-with(X115_df, plot3d(x,y,date, type="l",col=as.numeric(cut(X115_df$dist,3)), alpha=0.4)))
#Timelag
head(timeLag(X115, units="mins"))
head(timeLag(X115, units="hours"))
median(timeLag(X115, units="hours"))
max(timeLag(X115, units="hours"))
max(timeLag(X115, units="days"))
plot(X115_df$individual.local.identifier..ID.~X115_df $date,xlab = "Sampling period",ylab = "Individuals ID Number")
#write.table(X115,file="c:/RWorkDir/jaguardatapaper/ID/X115.txt",row.names = F,quote=F,col.names=T,sep="\t")
#Histograms and plots
hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 115],
     main = 'Individual 115', xlab = "Time between fixes (h)")
plot(density.default(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 115],
     na.rm = T), main = 'Individual 115')
	 
breaks <- c(0, 0.5, 1.5, 2.5, 4.5, 8.5, 16.5, 32.5, 64.5, 128.5, 256.5, 1000)
breaks <- c(0, 1, 2*(1:256))
hh <- hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 115 &
                                     mov.data.diff$time.diff > 0 & 
                                     mov.data.diff$time.diff < 512],
									 main = 'Individual 115',xlab = 'Time between fixes (h)',
           plot = T, breaks = breaks)

#X116
X116=unstacked$X116
X116_df <- as(X116, "data.frame")
head(X116_df)
nrow(X116_df)
range(X116_df$date)
m <- get_map(bbox(extent(X116)*1.1), source="google", zoom=12, maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=X116_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
#' ## Display space-time cube. 
open3d()
# To get a bigger window than the default
par3d(windowRect = c(100, 100, 612, 612))
Sys.sleep(0.1) # Allow sluggish window managers to catch up
with(X116_df, plot3d(x,y,date, type="l", col=as.numeric(X116_df$Event_ID)))
(stcube<-with(X116_df, plot3d(x,y,date, type="l",col=as.numeric(cut(X116_df$dist,3)), alpha=0.4)))
#Timelag
head(timeLag(X116, units="mins"))
head(timeLag(X116, units="hours"))
median(timeLag(X116, units="hours"))
max(timeLag(X116, units="hours"))
max(timeLag(X116, units="days"))
plot(X116_df$individual.local.identifier..ID.~X116_df $date,xlab = "Sampling period",ylab = "Individuals ID Number")
#write.table(X116,file="c:/RWorkDir/jaguardatapaper/ID/X116.txt",row.names = F,quote=F,col.names=T,sep="\t")
#Histograms and plots
hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 116],
     main = 'Individual 116', xlab = "Time between fixes (h)")
plot(density.default(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 116],
     na.rm = T), main = 'Individual 116')
	 
breaks <- c(0, 0.5, 1.5, 2.5, 4.5, 8.5, 16.5, 32.5, 64.5, 128.5, 256.5, 1000)
breaks <- c(0, 1, 2*(1:256))
hh <- hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 116 &
                                     mov.data.diff$time.diff > 0 & 
                                     mov.data.diff$time.diff < 512],
									 main = 'Individual 116',xlab = 'Time between fixes (h)',
           plot = T, breaks = breaks)

#X117
X117=unstacked$X117
X117_df <- as(X117, "data.frame")
head(X117_df)
nrow(X117_df)
range(X117_df$date)
m <- get_map(bbox(extent(X117)*1.1), source="google", zoom=13, maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=X117_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
m <- get_map(bbox(extent(X117)*1.1), source="google", zoom=12, maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=X117_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
#' ## Display space-time cube. 
open3d()
# To get a bigger window than the default
par3d(windowRect = c(100, 100, 612, 612))
Sys.sleep(0.1) # Allow sluggish window managers to catch up
with(X117_df, plot3d(x,y,date, type="l", col=as.numeric(X117_df$Event_ID)))
(stcube<-with(X117_df, plot3d(x,y,date, type="l",col=as.numeric(cut(X117_df$dist,3)), alpha=0.4)))
#Timelag
head(timeLag(X117, units="mins"))
head(timeLag(X117, units="hours"))
median(timeLag(X117, units="hours"))
max(timeLag(X117, units="hours"))
max(timeLag(X117, units="days"))
plot(X117_df$individual.local.identifier..ID.~X117_df $date,xlab = "Sampling period",ylab = "Individuals ID Number")
#write.table(X117,file="c:/RWorkDir/jaguardatapaper/ID/X117.txt",row.names = F,quote=F,col.names=T,sep="\t")
#Histograms and plots
hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 117],
     main = 'Individual 117', xlab = "Time between fixes (h)")
plot(density.default(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 117],
     na.rm = T), main = 'Individual 117')
	 
breaks <- c(0, 0.5, 1.5, 2.5, 4.5, 8.5, 16.5, 32.5, 64.5, 128.5, 256.5, 1000)
breaks <- c(0, 1, 2*(1:256))
hh <- hist(mov.data.diff$time.diff[mov.data.diff$individual.local.identifier..ID. == 117 &
                                     mov.data.diff$time.diff > 0 & 
                                     mov.data.diff$time.diff < 512],
									 main = 'Individual 117',xlab = 'Time between fixes (h)',
           plot = T, breaks = breaks)

########################################################################################################
# Ploting all data
move::plot(move.data, pch = 19, xlab = 'Longitude', ylab = 'Latitude')
cols = rainbow(nlevels(move.data@trackId))[move.data@trackId]
m <- get_map(bbox(extent(move.data)*1.1), source = 'google', zoom = 3)
ggmap(m) + geom_point(data=jaguar_df, aes(x=x, y=y), col = cols)

########################  Ploting regions or projects   ######################################################
# 1)  Cerrado and Para (Parauapebas) IOP      Cerrado(17, 65, 67, 82, 85)   ### Parauapebas (24) translocatade!!!
unstacked <- split(move.data)
X24=unstacked$X24  #Para
X24_df <- as(X24, "data.frame")
#Cerrado IOP
X17=unstacked$X17  
X17_df <- as(X17, "data.frame")
X65=unstacked$X65
X65_df <- as(X65, "data.frame")
X67=unstacked$X67
X67_df <- as(X67, "data.frame")
X82=unstacked$X82
X82_df <- as(X82, "data.frame")
X85=unstacked$X85
X85_df <- as(X85, "data.frame")

#####  1) IOP Cerrado e Para
m <- get_map(bbox(extent(X24)*1.1), source="google", zoom=10, maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=X24_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
m <- get_map(bbox(extent(X17)*1.1), source="google", zoom=10, maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=X17_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
m <- get_map(bbox(extent(X65)*1.1), source="google", zoom=10, maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=X65_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
m <- get_map(bbox(extent(X67)*1.1), source="google", zoom=11, maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=X67_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
ggmap(m)+geom_path(data=X82_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
ggmap(m)+geom_path(data=X85_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory

animals <- as(move.data, 'data.frame')
jaguar_df <- as(move.data, 'data.frame')
iop=subset(jaguar_df,study.name=='IOP')
iopPA=subset(iop,id=='24')
iop17=subset(iop,id=='17')
iop65=subset(iop,id=='65')
iop67=subset(iop,id=='67')
iop82=subset(iop,id=='82')
iop85=subset(iop,id=='85')
# Sampling period
plot(iop$individual.local.identifier..ID.~iop$date,xlab = "Sampling period",ylab = "Individuals ID Number")

#IOP Cerrado
iopCerrado=rbind(iop17, iop65, iop67, iop82, iop85)
head(iopCerrado)
str(iopCerrado)
m <- get_map(bbox(extent(iopCerrado)*1.1), source="google", zoom=10,maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=iop65,aes(x=x,y=y),colour = "red")+geom_path(data=iop17,aes(x=x,y=y),colour = "yellow")+geom_path(data=iop67,aes(x=x, y=y),colour = "green")+geom_path(data=iop82,aes(x=x, y=y),colour = "orange")+geom_path(data=iop85,aes(x=x, y=y),colour = "purple")
#ggmap(m)+geom_point(data=iop65,aes(x=x,y=y),colour = "red")+geom_point(data=iop17,aes(x=x,y=y),colour = "yellow")+geom_point(data=iop67,aes(x=x, y=y),colour = "green")+geom_point(data=iop82,aes(x=x, y=y),colour = "orange")+geom_point(data=iop85,aes(x=x, y=y),colour = "purple")

# Sampling period
plot(iopCerrado$individual.local.identifier..ID.~iopCerrado$date,xlab = "Sampling period",ylab = "Individuals ID Number")

#or using the individuals id above (same thing)
#IOP Cerrado
iopCerrado=rbind(iop17, iop65, iop67, iop82, iop85)
m <- get_map(bbox(extent(iopCerrado)*1.1), source="google", zoom=10,maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=X65_df,aes(x=x,y=y),colour = "red")+geom_path(data=X17_df,aes(x=x,y=y),colour = "yellow")+geom_path(data=X67_df,aes(x=x, y=y),colour = "green")+geom_path(data=X82_df,aes(x=x, y=y),colour = "orange")+geom_path(data=X85_df,aes(x=x, y=y),colour = "purple")
#write.table(iopCerrado,file="c:/RWorkDir/jaguardatapaper/ID/iopCerrado.txt",row.names = F,quote=F,col.names=T,sep="\t")


#IOP Para
m <- get_map(bbox(extent(X24)*1.1), source="google", zoom=10, maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=X24_df,aes(x=x, y=y),colour = "blue") #plot the map and the trajectory

##########################################################################################################
#Cerrado Fragmented Landscape
#X89
X89=unstacked$X89
X89_df <- as(X89, "data.frame")
head(X89_df)
nrow(X89_df)
range(X89_df$date)
m <- get_map(bbox(extent(X89)*1.1), source="google", zoom=11, maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=X89_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
# Sampling period
plot(X89_df$individual.local.identifier..ID.~X89_df$date,xlab = "Sampling period",ylab = "Individuals ID Number")

#Cerrado Total  #IOP Cerrado + #Cerrado Fragmented Landscape
jaguar_df <- as(move.data, 'data.frame')
iop=subset(jaguar_df,study.name=='IOP')
jfl=subset(jaguar_df,tag.local.identifier=="2314")
iop17=subset(iop,id=='17')
iop65=subset(iop,id=='65')
iop67=subset(iop,id=='67')
iop82=subset(iop,id=='82')
iop85=subset(iop,id=='85')
jfl89=subset(jfl,id=='89')

CerradoTotal=rbind(iop17, iop65, iop67, iop82, iop85, jfl89)
m <- get_map(bbox(extent(CerradoTotal)*1.1), source="google", zoom=7,maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=iop65,aes(x=x,y=y),colour = "red")+geom_path(data=iop17,aes(x=x,y=y),colour = "yellow")+geom_path(data=iop67,aes(x=x, y=y),colour = "green")+geom_path(data=iop82,aes(x=x, y=y),colour = "orange")+geom_path(data=iop85,aes(x=x, y=y),colour = "purple")+geom_path(data=jfl89,aes(x=x,y=y),colour = "blue")
#ggmap(m)+geom_path(data=X65_df,aes(x=x,y=y),colour = "red")+geom_path(data=X17_df,aes(x=x,y=y),colour = "yellow")+geom_path(data=X67_df,aes(x=x, y=y),colour = "green")+geom_path(data=X82_df,aes(x=x, y=y),colour = "orange")+geom_path(data=X85_df,aes(x=x, y=y),colour = "purple")+geom_path(data=X89_df,aes(x=x,y=y),colour = "blue")
plot(CerradoTotal$individual.local.identifier..ID.~CerradoTotal$date,xlab = "Sampling period",ylab = "Individuals ID Number")

######################################################

#Humid Chaco Paraguay
jaguar_df <- as(move.data, 'data.frame')
hch=subset(jaguar_df,study.name=='Humid Chaco')
hch1=subset(hch,id=='1')
hch3=subset(hch,id=='3')
hch4=subset(hch,id=='4')
hch5=subset(hch,id=='5')
hch7=subset(hch,id=='7')
hch9=subset(hch,id=='9')
hch10=subset(hch,id=='10')
hch11=subset(hch,id=='11')

#Humid Chaco
#humidchaco=rbind(hch1, hch3, hch4, hch5, hch7, hch9, hch10, hch11)
m <- get_map(bbox(extent(hch)*1.1), source="google", zoom=10,maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=hch1,aes(x=x,y=y),colour = "red")+geom_path(data=hch3,aes(x=x,y=y),colour = "yellow")+geom_path(data=hch4,aes(x=x, y=y),colour = "green")+geom_path(data=hch5,aes(x=x, y=y),colour = "orange")+geom_path(data=hch7,aes(x=x, y=y),colour = "purple")+geom_path(data=hch9,aes(x=x, y=y),colour = "blue")+geom_path(data=hch10,aes(x=x, y=y),colour = "aquamarine")+geom_path(data=hch11,aes(x=x, y=y),colour = "cyan")
# Sampling period
plot(hch$individual.local.identifier..ID.~hch$date,xlab = "Sampling period",ylab = "Individuals ID Number")

#################################################################################
##Humid Chaco Individuals
#X1	
unstacked <- split(move.data)
X1=unstacked$X1
X1_df <- as(X1, "data.frame")
head(X1_df)
m <- get_map(bbox(extent(X1)*1.1), source="google", zoom=12, maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=X1_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory

#X3
X3=unstacked$X3
X3_df <- as(X3, "data.frame")
head(X3_df)
m <- get_map(bbox(extent(X3)*1.1), source="google", zoom=10, maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=X3_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory

#m <- get_map(bbox(extent(X3)*1.1), source="google", zoom=11, maptype="satellite") #define the box for map based on the trajectory coordinates
#ggmap(m)+geom_path(data=X3_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
#Miss 1 row

#X4
X4=unstacked$X4
X4_df <- as(X4, "data.frame")
head(X4_df)
m <- get_map(bbox(extent(X4)*1.1), source="google", zoom=12, maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=X4_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory

#X5
X5=unstacked$X5
X5_df <- as(X5, "data.frame")
head(X5_df)
m <- get_map(bbox(extent(X5)*1.1), source="google", zoom=13, maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=X5_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory

#X7
X7=unstacked$X7
X7_df <- as(X7, "data.frame")
head(X7_df)
m <- get_map(bbox(extent(X7)*1.1), source="google", zoom=12, maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=X7_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory

#X9
X9=unstacked$X9
X9_df <- as(X9, "data.frame")
head(X9_df)
m <- get_map(bbox(extent(X9)*1.1), source="google", zoom=13, maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=X9_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
m <- get_map(bbox(extent(X9)*1.1), source="google", zoom=12, maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=X9_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory

#X10	
unstacked <- split(move.data)
X10=unstacked$X10
X10_df <- as(X10, "data.frame")
head(X1_df)
m <- get_map(bbox(extent(X10)*1.1), source="google", zoom=12, maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=X10_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory

#X11	
unstacked <- split(move.data)
X11=unstacked$X11
X11_df <- as(X11, "data.frame")
head(X11_df)
m <- get_map(bbox(extent(X11)*1.1), source="google", zoom=12, maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=X11_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
#############################################################################################################

## Atlantic Forest - Paraguay  ###
jaguar_df <- as(move.data, 'data.frame')
Af=subset(jaguar_df,study.name=='Atlantic forest') #Nao inclui Iguacu e outros Projetos na AF
AfParaguay=subset(Af,country=='Paraguay')

Af2=subset(AfParaguay,id=='2')
Af6=subset(AfParaguay,id=='6')           ### AF junto aos do humid chaco
Af8=subset(AfParaguay,id=='8')
Af21=subset(AfParaguay,id=='21')
Af76=subset(AfParaguay,id=='76')         ### AF junto aos do dry chaco
Af78=subset(AfParaguay,id=='78')

#Atlantic Forest - Paraguay
m <- get_map(bbox(extent(AfParaguay)*1.1), source="google", zoom=7,maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=Af2,aes(x=x,y=y),colour = "red")+geom_path(data=Af6,aes(x=x,y=y),colour = "yellow")+geom_path(data=Af8,aes(x=x, y=y),colour = "green")+geom_path(data=Af21,aes(x=x, y=y),colour = "orange")+geom_path(data=Af76,aes(x=x, y=y),colour = "purple")+geom_path(data=Af78,aes(x=x, y=y),colour = "blue")
plot(AfParaguay$individual.local.identifier..ID.~AfParaguay$date,xlab = "Sampling period",ylab = "Individuals ID Number")


##Atlantic Forest - Paraguay Individuals
#X2 
X2=unstacked$X2
X2_df <- as(X2, "data.frame")
head(X2_df)
m <- get_map(bbox(extent(X2)*1.1), source="google", zoom=13, maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=X2_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
#m <- get_map(bbox(extent(X2)*1.1), source="google", zoom=12, maptype="satellite") #define the box for map based on the trajectory coordinates
#ggmap(m)+geom_path(data=X2_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory

#X6
X6=unstacked$X6
X6_df <- as(X6, "data.frame")
head(X6_df)
m <- get_map(bbox(extent(X6)*1.1), source="google", zoom=11, maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=X6_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
m <- get_map(bbox(extent(X6)*1.1), source="google", zoom=12, maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=X6_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory

#X8
X8=unstacked$X8
X8_df <- as(X8, "data.frame")
head(X8_df)
m <- get_map(bbox(extent(X8)*1.1), source="google", zoom=11, maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=X8_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
#m <- get_map(bbox(extent(X8)*1.1), source="google", zoom=12, maptype="satellite") #define the box for map based on the trajectory coordinates
#ggmap(m)+geom_path(data=X8_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
#6 row warning

#X21
unstacked <- split(move.data)
X21=unstacked$X21
X21_df <- as(X21, "data.frame")
head(X21_df)
nrow(X21_df)
m <- get_map(bbox(extent(X21)*1.1), source="google", zoom=13, maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=X21_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory

#X76
X76=unstacked$X76
X76_df <- as(X76, "data.frame")
head(X76_df)
m <- get_map(bbox(extent(X76)*1.1), source="google", zoom=10, maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=X76_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory

#X78
X78=unstacked$X78
X78_df <- as(X78, "data.frame")
head(X78_df)
m <- get_map(bbox(extent(X78)*1.1), source="google", zoom=10, maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=X78_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory

##############################################################################################################
#Atlantic Forest & Humid Chaco Paraguay
jaguar_df <- as(move.data, 'data.frame')
Paraguay=subset(jaguar_df,country=='Paraguay')
AfParaguay=subset(Af,country=='Paraguay')

Af2=subset(Paraguay,id=='2')
Af6=subset(Paraguay,id=='6')
Af8=subset(Paraguay,id=='8')
Af21=subset(Paraguay,id=='21')
Af76=subset(Paraguay,id=='76')
Af78=subset(Paraguay,id=='78')
hch1=subset(Paraguay,id=='1')
hch3=subset(Paraguay,id=='3')
hch4=subset(Paraguay,id=='4')
hch5=subset(Paraguay,id=='5')
hch7=subset(Paraguay,id=='7')
hch9=subset(Paraguay,id=='9')
hch10=subset(Paraguay,id=='10')
hch11=subset(Paraguay,id=='11')
hch=subset(jaguar_df,study.name=='Humid Chaco')
AfParaguay=subset(Af,country=='Paraguay')

hchAfParaguay=rbind(hch,AfParaguay)
m <- get_map(bbox(extent(hchAfParaguay)*1.1), source="google", zoom=7,maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=Af2,aes(x=x,y=y),colour = "red")+geom_path(data=Af6,aes(x=x,y=y),colour = "yellow")+geom_path(data=Af8,aes(x=x, y=y),colour = "green")+geom_path(data=Af21,aes(x=x, y=y),colour = "orange")+geom_path(data=Af76,aes(x=x, y=y),colour = "purple")+geom_path(data=Af78,aes(x=x, y=y),colour = "blue")+geom_path(data=hch1,aes(x=x,y=y),colour = "red")+geom_path(data=hch3,aes(x=x,y=y),colour = "yellow")+geom_path(data=hch4,aes(x=x, y=y),colour = "green")+geom_path(data=hch5,aes(x=x, y=y),colour = "orange")+geom_path(data=hch7,aes(x=x, y=y),colour = "purple")+geom_path(data=hch9,aes(x=x, y=y),colour = "blue")+geom_path(data=hch10,aes(x=x, y=y),colour = "aquamarine")+geom_path(data=hch11,aes(x=x, y=y),colour = "cyan")
# Sampling period
plot(hchAfParaguay$individual.local.identifier..ID.~hchAfParaguay$date,xlab = "Sampling period",ylab = "Individuals ID Number")


#################################################################
#Dry Chaco
jaguar_df <- as(move.data, 'data.frame')
Paraguay=subset(jaguar_df,country=='Paraguay')
drych16=subset(Paraguay,id=='16')
drych70=subset(Paraguay,id=='70')
drych71=subset(Paraguay,id=='71')
drych72=subset(Paraguay,id=='72')
drych73=subset(Paraguay,id=='73')
drych77=subset(Paraguay,id=='77')
drych=rbind(drych16, drych70, drych71, drych72, drych73,drych77)
m <- get_map(bbox(extent(drych)*1.1), source="google", zoom=7,maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=drych16,aes(x=x,y=y),colour = "red")+geom_path(data=drych70,aes(x=x,y=y),colour = "yellow")+geom_path(data=drych71,aes(x=x, y=y),colour = "green")+geom_path(data=drych72,aes(x=x, y=y),colour = "orange")+geom_path(data=drych73,aes(x=x, y=y),colour = "purple")+geom_path(data=drych77,aes(x=x, y=y),colour = "blue")
plot(drych$individual.local.identifier..ID.~drych$date,xlab = "Sampling period",ylab = "Individuals ID Number")

###############################################################

#Paraguayan Pantanal
jaguar_df <- as(move.data, 'data.frame')
Paraguay=subset(jaguar_df,country=='Paraguay')

pant51=subset(Paraguay,id=='51')
pant74=subset(Paraguay,id=='74')
pant75=subset(Paraguay,id=='75')
pantParaguay=rbind(pant51, pant74, pant75)
m <- get_map(bbox(extent(pantParaguay)*1.1), source="google", zoom=10,maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=pant51,aes(x=x,y=y),colour = "cyan1")+geom_path(data=pant75,aes(x=x, y=y),colour = "purple")+geom_path(data=pant74,aes(x=x,y=y),colour = "blue")
plot(pantParaguay$individual.local.identifier..ID.~pantParaguay$date,xlab = "Sampling period",ylab = "Individuals ID Number")

######################################################################
#Paraguay (Total)

drych=rbind(drych16, drych70, drych71, drych72, drych73,drych77)
hch=subset(jaguar_df,study.name=='Humid Chaco')
AfParaguay=subset(Af,country=='Paraguay')
pantParaguay=rbind(pant51, pant74, pant75)

m <- get_map(bbox(extent(Paraguay)*1.1), source="google", zoom=7,maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=AfParaguay,aes(x=x,y=y),colour = "green")+geom_path(data=hch,aes(x=x,y=y),colour = "cyan")+geom_path(data=drych,aes(x=x,y=y),colour = "orange")+geom_path(data=pantParaguay,aes(x=x,y=y),colour = "blue")
#ggmap(m)+geom_point(data=AfParaguay,aes(x=x,y=y),colour = "green")+geom_path(data=hch,aes(x=x,y=y),colour = "cyan")+geom_path(data=drych,aes(x=x,y=y),colour = "orange")+geom_path(data=pantParaguay,aes(x=x,y=y),colour = "blue")
plot(Paraguay$individual.local.identifier..ID.~Paraguay$date,xlab = "Sampling period",ylab = "Individuals ID Number")

##########################################################################################################
   ### Argentina (AF)  
jaguar_df <- as(move.data,'data.frame')
Argentina=subset(jaguar_df,country=='Argentina')
arg42=subset(Argentina,id=='42')
arg80=subset(Argentina,id=='80')
arg90=subset(Argentina,id=='90')
m <- get_map(bbox(extent(Argentina)*1.1), source="google", zoom=10,maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=arg42,aes(x=x,y=y),colour = "yellow")+geom_path(data=arg80,aes(x=x,y=y),colour = "red")+geom_path(data=arg90,aes(x=x,y=y),colour = "cyan")
plot(Argentina$individual.local.identifier..ID.~Argentina$date,xlab = "Sampling period",ylab = "Individuals ID Number")

###########################################################################################################
    ### Iguazu (Brazil) AF
jaguar_df <- as(move.data,'data.frame')
Brazil=subset(jaguar_df,country=='Brazil')
AFI66=subset(Brazil,id=='66')
AFI83=subset(Brazil,id=='83')
BrazilIguazu=rbind(AFI66,AFI83)
m <- get_map(bbox(extent(BrazilIguazu)*1.1), source="google", zoom=10,maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=AFI66,aes(x=x,y=y),colour = "orange")+geom_path(data=AFI83,aes(x=x,y=y),colour = "blue")	
plot(BrazilIguazu$individual.local.identifier..ID.~BrazilIguazu$date,xlab = "Sampling period",ylab = "Individuals ID Number")

################################################################################
################ AF Iguacu (Argentina - Brazil) ################################
jaguar_df <- as(move.data,'data.frame')
Argentina=subset(jaguar_df,country=='Argentina')
Brazil=subset(jaguar_df,country=='Brazil')
arg42=subset(Argentina,id=='42')
arg80=subset(Argentina,id=='80')
arg90=subset(Argentina,id=='90')
AFI66=subset(Brazil,id=='66')
AFI83=subset(Brazil,id=='83')
AFIguazuTotal=rbind(AFI66,AFI83,arg42,arg80,arg90)
m <- get_map(bbox(extent(AFIguazuTotal)*1.1), source="google", zoom=9,maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=AFI66,aes(x=x,y=y),colour = "orange")+geom_path(data=AFI83,aes(x=x,y=y),colour = "blue")	+geom_path(data=arg42,aes(x=x,y=y),colour = "yellow")+geom_path(data=arg80,aes(x=x,y=y),colour = "red")+geom_path(data=arg90,aes(x=x,y=y),colour = "cyan")
plot(AFIguazuTotal$individual.local.identifier..ID.~AFIguazuTotal$date,xlab = "Sampling period",ylab = "Individuals ID Number")

################################################################################
##########################################################################################################
#AF Brazil
jaguar_df <- as(move.data,'data.frame')
Brazil=subset(jaguar_df,country=='Brazil')
AF34=subset(Brazil,id=='34')
AF35=subset(Brazil,id=='35')
AF36=subset(Brazil,id=='36')
AF37=subset(Brazil,id=='37')
AF38=subset(Brazil,id=='38')
AF39=subset(Brazil,id=='39')
AF40=subset(Brazil,id=='40')
AF58=subset(Brazil,id=='58')
AF62=subset(Brazil,id=='62')
AF63=subset(Brazil,id=='63')
AFBrazilNW=rbind(AF34,AF35,AF36,AF37,AF38,AF39,AF40,AF58,AF62,AF63) #Laury
m <- get_map(bbox(extent(AFBrazilNW)*1.1), source="google", zoom=8,maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=AF36,aes(x=x,y=y),colour = "cyan")+geom_path(data=AF34,aes(x=x,y=y),colour = "yellow")+geom_path(data=AF35,aes(x=x,y=y),colour = "red")+geom_path(data=AF37,aes(x=x,y=y),colour = "deeppink")+geom_path(data=AF38,aes(x=x,y=y),colour = "orange")+geom_path(data=AF58,aes(x=x,y=y),colour = "white")+geom_path(data=AF62,aes(x=x,y=y),colour = "purple")+geom_path(data=AF40,aes(x=x,y=y),colour = "yellow")+geom_path(data=AF39,aes(x=x,y=y),colour = "blue")+geom_path(data=AF63,aes(x=x,y=y),colour = "aquamarine")
#m <- get_map(bbox(extent(AFBrazilNW)*1.1), source="google", zoom=9,maptype="satellite") #define the box for map based on the trajectory coordinates
#ggmap(m)+geom_path(data=AF36,aes(x=x,y=y),colour = "cyan")+geom_path(data=AF34,aes(x=x,y=y),colour = "yellow")+geom_path(data=AF35,aes(x=x,y=y),colour = "red")+geom_path(data=AF37,aes(x=x,y=y),colour = "deeppink")+geom_path(data=AF38,aes(x=x,y=y),colour = "orange")+geom_path(data=AF58,aes(x=x,y=y),colour = "white")+geom_path(data=AF62,aes(x=x,y=y),colour = "purple")+geom_path(data=AF40,aes(x=x,y=y),colour = "yellow")+geom_path(data=AF39,aes(x=x,y=y),colour = "blue")+geom_path(data=AF63,aes(x=x,y=y),colour = "aquamarine")
plot(AFBrazilNW$individual.local.identifier..ID.~AFBrazilNW$date,xlab = "Sampling period",ylab = "Individuals ID Number")

###################################################################################################
   ### Caatinga (Brazil)
jaguar_df <- as(move.data,'data.frame')
Brazil=subset(jaguar_df,country=='Brazil')
ca20=subset(Brazil,id=='20')
ca50=subset(Brazil,id=='50')
Caatinga=rbind(ca20,ca50)
m <- get_map(bbox(extent(Caatinga)*1.1), source="google", zoom=10,maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=ca20,aes(x=x,y=y),colour = "orange")+geom_path(data=ca50,aes(x=x,y=y),colour = "red")
plot(Caatinga$individual.local.identifier..ID.~Caatinga$date,xlab = "Sampling period",ylab = "Individuals ID Number")

###################################################################################################
   ### Amazonia Mamiraua (Brazil)
jaguar_df <- as(move.data,'data.frame')
Brazil=subset(jaguar_df,country=='Brazil')
mam93=subset(Brazil,id=='93')
mam94=subset(Brazil,id=='94')
mam95=subset(Brazil,id=='95')
mam96=subset(Brazil,id=='96')
mam97=subset(Brazil,id=='97')
mam98=subset(Brazil,id=='98')
mam99=subset(Brazil,id=='99')
mam100=subset(Brazil,id=='100')
Mamiraua=rbind(mam93,mam94,mam95,mam96,mam97,mam98,mam99,mam100)
m <- get_map(bbox(extent(Mamiraua)*1.1), source="google", zoom=11,maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=mam93,aes(x=x,y=y),colour = "yellow")+geom_path(data=mam94,aes(x=x,y=y),colour = "orange")+geom_path(data=mam99,aes(x=x,y=y),colour = "blue")+geom_path(data=mam96,aes(x=x,y=y),colour = "deeppink")+geom_path(data=mam97,aes(x=x,y=y),colour = "red")+geom_path(data=mam98,aes(x=x,y=y),colour = "cyan")+geom_path(data=mam95,aes(x=x,y=y),colour = "green")+geom_path(data=mam100,aes(x=x,y=y),colour = "white")
plot(Mamiraua$individual.local.identifier..ID.~Mamiraua$date,xlab = "Sampling period",ylab = "Individuals ID Number")

###################################################################################################
   ### Pantanal (Brazil)
   ###Oncafari
jaguar_df <- as(move.data,'data.frame')
Brazil=subset(jaguar_df,country=='Brazil')
oncafari14=subset(Brazil,id=='14')
oncafari15=subset(Brazil,id=='15')
oncafari19=subset(Brazil,id=='19')
oncafari25=subset(Brazil,id=='25')
oncafari68=subset(Brazil,id=='68')
oncafari69=subset(Brazil,id=='69')
oncafari79=subset(Brazil,id=='79')
oncafari84=subset(Brazil,id=='84')
oncafari86=subset(Brazil,id=='86')
oncafari87=subset(Brazil,id=='87')
Oncafari=rbind(oncafari14,oncafari15,oncafari19,oncafari25,oncafari68,oncafari69,oncafari79,oncafari84,oncafari86,oncafari87)
m <- get_map(bbox(extent(Oncafari)*1.1), source="google", zoom=11,maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=oncafari84,aes(x=x,y=y),colour = "white")+geom_path(data=oncafari15,aes(x=x,y=y),colour = "orange")+geom_path(data=oncafari79,aes(x=x,y=y),colour = "green")+geom_path(data=oncafari19,aes(x=x,y=y),colour = "blue")+geom_path(data=oncafari25,aes(x=x,y=y),colour = "deeppink")+geom_path(data=oncafari68,aes(x=x,y=y),colour = "red")+geom_path(data=oncafari69,aes(x=x,y=y),colour = "cyan")+geom_path(data=oncafari86,aes(x=x,y=y),colour = "purple")+geom_path(data=oncafari87,aes(x=x,y=y),colour = "aquamarine")+geom_path(data=oncafari14,aes(x=x,y=y),colour = "yellow")
plot(Oncafari$individual.local.identifier..ID.~Oncafari$date,xlab = "Sampling period",ylab = "Individuals ID Number")


###Taiama
jaguar_df <- as(move.data,'data.frame')
Brazil=subset(jaguar_df,country=='Brazil')
taiama12=subset(Brazil,id=='12')
taiama13=subset(Brazil,id=='13')
taiama18=subset(Brazil,id=='18')
taiama22=subset(Brazil,id=='22')
taiama23=subset(Brazil,id=='23')
taiama41=subset(Brazil,id=='41')
taiama52=subset(Brazil,id=='52')
taiama81=subset(Brazil,id=='81')
taiama88=subset(Brazil,id=='88')
taiama91=subset(Brazil,id=='91')
taiama92=subset(Brazil,id=='92')
taiama116=subset(Brazil,id=='116')
taiama117=subset(Brazil,id=='117')
Taiama=rbind(taiama12,taiama13,taiama18,taiama22,taiama23,taiama41,taiama52,taiama81,taiama88,taiama91,taiama92,taiama116,taiama117)
m <- get_map(bbox(extent(Taiama)*1.1), source="google", zoom=11,maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=taiama12,aes(x=x,y=y),colour = "white")+geom_path(data=taiama13,aes(x=x,y=y),colour = "orange")+geom_path(data=taiama18,aes(x=x,y=y),colour = "green")+geom_path(data=taiama22,aes(x=x,y=y),colour = "blue")+geom_path(data=taiama23,aes(x=x,y=y),colour = "deeppink")+geom_path(data=taiama41,aes(x=x,y=y),colour = "red")+geom_path(data=taiama52,aes(x=x,y=y),colour = "cyan")+geom_path(data=taiama81,aes(x=x,y=y),colour = "purple")+geom_path(data=taiama88,aes(x=x,y=y),colour = "aquamarine")+geom_path(data=taiama91,aes(x=x,y=y),colour = "yellow")+geom_path(data=taiama92,aes(x=x,y=y),colour = "beige")+geom_path(data=taiama116,aes(x=x,y=y),colour = "coral1")+geom_path(data=taiama117,aes(x=x,y=y),colour = "chartreuse")
m <- get_map(bbox(extent(Taiama)*1.1), source="google", zoom=10,maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=taiama12,aes(x=x,y=y),colour = "white")+geom_path(data=taiama13,aes(x=x,y=y),colour = "orange")+geom_path(data=taiama18,aes(x=x,y=y),colour = "green")+geom_path(data=taiama22,aes(x=x,y=y),colour = "blue")+geom_path(data=taiama23,aes(x=x,y=y),colour = "deeppink")+geom_path(data=taiama41,aes(x=x,y=y),colour = "red")+geom_path(data=taiama52,aes(x=x,y=y),colour = "cyan")+geom_path(data=taiama81,aes(x=x,y=y),colour = "purple")+geom_path(data=taiama88,aes(x=x,y=y),colour = "aquamarine")+geom_path(data=taiama91,aes(x=x,y=y),colour = "yellow")+geom_path(data=taiama92,aes(x=x,y=y),colour = "beige")+geom_path(data=taiama116,aes(x=x,y=y),colour = "coral1")+geom_path(data=taiama117,aes(x=x,y=y),colour = "chartreuse")
plot(Taiama$individual.local.identifier..ID.~Taiama$date,xlab = "Sampling period",ylab = "Individuals ID Number")

###Panthera
jaguar_df <- as(move.data,'data.frame')
Brazil=subset(jaguar_df,country=='Brazil')
panthera27=subset(Brazil,id=='27')
panthera28=subset(Brazil,id=='28')
panthera29=subset(Brazil,id=='29')
panthera30=subset(Brazil,id=='30')
panthera31=subset(Brazil,id=='31')
panthera32=subset(Brazil,id=='32')
panthera33=subset(Brazil,id=='33')
panthera53=subset(Brazil,id=='53')
panthera54=subset(Brazil,id=='54')
panthera55=subset(Brazil,id=='55')
panthera56=subset(Brazil,id=='56')
panthera57=subset(Brazil,id=='57')
panthera59=subset(Brazil,id=='59')
panthera60=subset(Brazil,id=='60')
panthera61=subset(Brazil,id=='61')
#### Panthera All
Panthera=rbind(panthera27,panthera28,panthera29,panthera30,panthera31,panthera32,panthera33,panthera53,panthera54,panthera55,panthera56,panthera57,panthera59,panthera60,panthera61)
m <- get_map(bbox(extent(Panthera)*1.1), source="google", zoom=9,maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=panthera27,aes(x=x,y=y),colour = "white")+geom_path(data=panthera30,aes(x=x,y=y),colour = "blue")+geom_path(data=panthera31,aes(x=x,y=y),colour = "deeppink")+geom_path(data=panthera32,aes(x=x,y=y),colour = "red")+geom_path(data=panthera33,aes(x=x,y=y),colour = "cyan")+geom_path(data=panthera53,aes(x=x,y=y),colour = "purple")+geom_path(data=panthera55,aes(x=x,y=y),colour = "yellow")+geom_path(data=panthera59,aes(x=x,y=y),colour = "chartreuse")+geom_path(data=panthera60,aes(x=x,y=y),colour = "deepskyblue")+geom_path(data=panthera61,aes(x=x,y=y),colour = "chocolate")+geom_path(data=panthera28,aes(x=x,y=y),colour = "orange")+geom_path(data=panthera54,aes(x=x,y=y),colour = "aquamarine")+geom_path(data=panthera29,aes(x=x,y=y),colour = "green")+geom_path(data=panthera56,aes(x=x,y=y),colour = "beige")+geom_path(data=panthera57,aes(x=x,y=y),colour = "red")
#m <- get_map(bbox(extent(Panthera)*1.1), source="google", zoom=10,maptype="satellite") #define the box for map based on the trajectory coordinates
#ggmap(m)+geom_path(data=panthera27,aes(x=x,y=y),colour = "white")+geom_path(data=panthera30,aes(x=x,y=y),colour = "blue")+geom_path(data=panthera31,aes(x=x,y=y),colour = "deeppink")+geom_path(data=panthera32,aes(x=x,y=y),colour = "red")+geom_path(data=panthera33,aes(x=x,y=y),colour = "cyan")+geom_path(data=panthera53,aes(x=x,y=y),colour = "purple")+geom_path(data=panthera55,aes(x=x,y=y),colour = "yellow")+geom_path(data=panthera57,aes(x=x,y=y),colour = "coral1")+geom_path(data=panthera59,aes(x=x,y=y),colour = "chartreuse")+geom_path(data=panthera60,aes(x=x,y=y),colour = "deepskyblue")+geom_path(data=panthera61,aes(x=x,y=y),colour = "chocolate")+geom_path(data=panthera28,aes(x=x,y=y),colour = "orange")+geom_path(data=panthera54,aes(x=x,y=y),colour = "aquamarine")+geom_path(data=panthera29,aes(x=x,y=y),colour = "green")+geom_path(data=panthera56,aes(x=x,y=y),colour = "beige")
plot(Panthera$individual.local.identifier..ID.~Panthera$date,xlab = "Sampling period",ylab = "Individuals ID Number")


Panthera1=rbind(panthera28,panthera54,panthera29,panthera56,panthera57)
m <- get_map(bbox(extent(Panthera1)*1.1), source="google", zoom=11,maptype="satellite")
ggmap(m)+geom_path(data=panthera28,aes(x=x,y=y),colour = "orange")+geom_path(data=panthera54,aes(x=x,y=y),colour = "aquamarine")+geom_path(data=panthera29,aes(x=x,y=y),colour = "green")+geom_path(data=panthera56,aes(x=x,y=y),colour = "beige")+geom_path(data=panthera57,aes(x=x,y=y),colour = "red")
plot(Panthera1$individual.local.identifier..ID.~Panthera1$date,xlab = "Sampling period",ylab = "Individuals ID Number")


Panthera2=rbind(panthera27,panthera30,panthera31,panthera32,panthera33,panthera53,panthera55,panthera59,panthera60,panthera61)
m <- get_map(bbox(extent(Panthera2)*1.1), source="google", zoom=11,maptype="satellite")
ggmap(m)+geom_path(data=panthera27,aes(x=x,y=y),colour = "white")+geom_path(data=panthera30,aes(x=x,y=y),colour = "blue")+geom_path(data=panthera31,aes(x=x,y=y),colour = "deeppink")+geom_path(data=panthera32,aes(x=x,y=y),colour = "red")+geom_path(data=panthera33,aes(x=x,y=y),colour = "cyan")+geom_path(data=panthera53,aes(x=x,y=y),colour = "purple")+geom_path(data=panthera55,aes(x=x,y=y),colour = "yellow")+geom_path(data=panthera59,aes(x=x,y=y),colour = "chartreuse")+geom_path(data=panthera60,aes(x=x,y=y),colour = "deepskyblue")+geom_path(data=panthera61,aes(x=x,y=y),colour = "chocolate")
plot(Panthera2$individual.local.identifier..ID.~Panthera2$date,xlab = "Sampling period",ylab = "Individuals ID Number")

###Rio Negro - Pantanal Brazil
jaguar_df <- as(move.data,'data.frame')
Brazil=subset(jaguar_df,country=='Brazil')
rionegro101=subset(Brazil,id=='101')
rionegro102=subset(Brazil,id=='102')
rionegro103=subset(Brazil,id=='103')
rionegro104=subset(Brazil,id=='104')
RioNegro=rbind(rionegro101,rionegro102,rionegro103,rionegro104)
m <- get_map(bbox(extent(RioNegro)*1.1), source="google", zoom=10,maptype="satellite")
ggmap(m)+geom_path(data=rionegro101,aes(x=x,y=y),colour = "orange")+geom_path(data=rionegro102,aes(x=x,y=y),colour = "aquamarine")+geom_path(data=rionegro103,aes(x=x,y=y),colour = "green")+geom_path(data=rionegro104,aes(x=x,y=y),colour = "red")
#m <- get_map(bbox(extent(RioNegro)*1.1), source="google", zoom=11,maptype="satellite")
#ggmap(m)+geom_path(data=rionegro101,aes(x=x,y=y),colour = "orange")+geom_path(data=rionegro102,aes(x=x,y=y),colour = "aquamarine")+geom_path(data=rionegro103,aes(x=x,y=y),colour = "green")+geom_path(data=rionegro104,aes(x=x,y=y),colour = "red")
plot(RioNegro$individual.local.identifier..ID.~RioNegro$date,xlab = "Sampling period",ylab = "Individuals ID Number")


###Sao Bento - Pantanal Brazil
jaguar_df <- as(move.data,'data.frame')
Brazil=subset(jaguar_df,country=='Brazil')

saobento105=subset(Brazil,id=='105')
saobento106=subset(Brazil,id=='106')
saobento107=subset(Brazil,id=='107')
saobento108=subset(Brazil,id=='108')
saobento109=subset(Brazil,id=='109')
saobento110=subset(Brazil,id=='110')
saobento111=subset(Brazil,id=='111')
saobento112=subset(Brazil,id=='112')
saobento113=subset(Brazil,id=='113')
saobento114=subset(Brazil,id=='114')
saobento115=subset(Brazil,id=='115')
SaoBento=rbind(saobento105,saobento106,saobento107,saobento108,saobento109,saobento110,saobento111,saobento112,saobento113,saobento114,saobento115)
m <- get_map(bbox(extent(SaoBento)*1.1), source="google", zoom=11,maptype="satellite")
ggmap(m)+geom_path(data=saobento114,aes(x=x,y=y),colour = "purple")+geom_path(data=saobento105,aes(x=x,y=y),colour = "orange")+geom_path(data=saobento106,aes(x=x,y=y),colour = "aquamarine")+geom_path(data=saobento107,aes(x=x,y=y),colour = "green")+geom_path(data=saobento108,aes(x=x,y=y),colour = "coral1")+geom_path(data=saobento109,aes(x=x,y=y),colour = "white")+geom_path(data=saobento110,aes(x=x,y=y),colour = "blue")+geom_path(data=saobento111,aes(x=x,y=y),colour = "deeppink")+geom_path(data=saobento112,aes(x=x,y=y),colour = "red")+geom_path(data=saobento113,aes(x=x,y=y),colour = "cyan")+geom_path(data=saobento115,aes(x=x,y=y),colour = "yellow")
plot(SaoBento$individual.local.identifier..ID.~SaoBento$date,xlab = "Sampling period",ylab = "Individuals ID Number")


########### Brazilian Pantanal  (All)

BrPantanal=rbind(SaoBento,RioNegro,Panthera,Taiama,Oncafari)
m <- get_map(bbox(extent(BrPantanal)*1.1), source="google", zoom=7,maptype="satellite")
ggmap(m)+geom_path(data=saobento114,aes(x=x,y=y),colour = "purple")+geom_path(data=saobento105,aes(x=x,y=y),colour = "orange")+geom_path(data=saobento106,aes(x=x,y=y),colour = "aquamarine")+geom_path(data=saobento107,aes(x=x,y=y),colour = "green")+geom_path(data=saobento108,aes(x=x,y=y),colour = "coral1")+geom_path(data=saobento109,aes(x=x,y=y),colour = "white")+geom_path(data=saobento110,aes(x=x,y=y),colour = "blue")+geom_path(data=saobento111,aes(x=x,y=y),colour = "deeppink")+geom_path(data=saobento112,aes(x=x,y=y),colour = "red")+geom_path(data=saobento113,aes(x=x,y=y),colour = "cyan")+geom_path(data=saobento115,aes(x=x,y=y),colour = "yellow")+geom_path(data=rionegro101,aes(x=x,y=y),colour = "orange")+geom_path(data=rionegro102,aes(x=x,y=y),colour = "aquamarine")+geom_path(data=rionegro103,aes(x=x,y=y),colour = "green")+geom_path(data=rionegro104,aes(x=x,y=y),colour = "red")+geom_path(data=panthera27,aes(x=x,y=y),colour = "white")+geom_path(data=panthera30,aes(x=x,y=y),colour = "blue")+geom_path(data=panthera31,aes(x=x,y=y),colour = "deeppink")+geom_path(data=panthera32,aes(x=x,y=y),colour = "red")+geom_path(data=panthera33,aes(x=x,y=y),colour = "cyan")+geom_path(data=panthera53,aes(x=x,y=y),colour = "purple")+geom_path(data=panthera55,aes(x=x,y=y),colour = "yellow")+geom_path(data=panthera59,aes(x=x,y=y),colour = "chartreuse")+geom_path(data=panthera60,aes(x=x,y=y),colour = "deepskyblue")+geom_path(data=panthera61,aes(x=x,y=y),colour = "chocolate")+geom_path(data=panthera28,aes(x=x,y=y),colour = "orange")+geom_path(data=panthera54,aes(x=x,y=y),colour = "aquamarine")+geom_path(data=panthera29,aes(x=x,y=y),colour = "green")+geom_path(data=panthera56,aes(x=x,y=y),colour = "beige")+geom_path(data=panthera57,aes(x=x,y=y),colour = "red")+geom_path(data=taiama12,aes(x=x,y=y),colour = "white")+geom_path(data=taiama13,aes(x=x,y=y),colour = "orange")+geom_path(data=taiama18,aes(x=x,y=y),colour = "green")+geom_path(data=taiama22,aes(x=x,y=y),colour = "blue")+geom_path(data=taiama23,aes(x=x,y=y),colour = "deeppink")+geom_path(data=taiama41,aes(x=x,y=y),colour = "red")+geom_path(data=taiama52,aes(x=x,y=y),colour = "cyan")+geom_path(data=taiama81,aes(x=x,y=y),colour = "purple")+geom_path(data=taiama88,aes(x=x,y=y),colour = "aquamarine")+geom_path(data=taiama91,aes(x=x,y=y),colour = "yellow")+geom_path(data=taiama92,aes(x=x,y=y),colour = "beige")+geom_path(data=taiama116,aes(x=x,y=y),colour = "coral1")+geom_path(data=taiama117,aes(x=x,y=y),colour = "chartreuse")+geom_path(data=oncafari84,aes(x=x,y=y),colour = "white")+geom_path(data=oncafari15,aes(x=x,y=y),colour = "orange")+geom_path(data=oncafari79,aes(x=x,y=y),colour = "green")+geom_path(data=oncafari19,aes(x=x,y=y),colour = "blue")+geom_path(data=oncafari25,aes(x=x,y=y),colour = "deeppink")+geom_path(data=oncafari68,aes(x=x,y=y),colour = "red")+geom_path(data=oncafari69,aes(x=x,y=y),colour = "cyan")+geom_path(data=oncafari86,aes(x=x,y=y),colour = "purple")+geom_path(data=oncafari87,aes(x=x,y=y),colour = "aquamarine")+geom_path(data=oncafari14,aes(x=x,y=y),colour = "yellow")
#m <- get_map(bbox(extent(BrPantanal)*1.1), source="google", zoom=8,maptype="satellite")
#ggmap(m)+geom_path(data=saobento114,aes(x=x,y=y),colour = "purple")+geom_path(data=saobento105,aes(x=x,y=y),colour = "orange")+geom_path(data=saobento106,aes(x=x,y=y),colour = "aquamarine")+geom_path(data=saobento107,aes(x=x,y=y),colour = "green")+geom_path(data=saobento108,aes(x=x,y=y),colour = "coral1")+geom_path(data=saobento109,aes(x=x,y=y),colour = "white")+geom_path(data=saobento110,aes(x=x,y=y),colour = "blue")+geom_path(data=saobento111,aes(x=x,y=y),colour = "deeppink")+geom_path(data=saobento112,aes(x=x,y=y),colour = "red")+geom_path(data=saobento113,aes(x=x,y=y),colour = "cyan")+geom_path(data=saobento115,aes(x=x,y=y),colour = "yellow")+geom_path(data=rionegro101,aes(x=x,y=y),colour = "orange")+geom_path(data=rionegro102,aes(x=x,y=y),colour = "aquamarine")+geom_path(data=rionegro103,aes(x=x,y=y),colour = "green")+geom_path(data=rionegro104,aes(x=x,y=y),colour = "red")+geom_path(data=panthera27,aes(x=x,y=y),colour = "white")+geom_path(data=panthera30,aes(x=x,y=y),colour = "blue")+geom_path(data=panthera31,aes(x=x,y=y),colour = "deeppink")+geom_path(data=panthera32,aes(x=x,y=y),colour = "red")+geom_path(data=panthera33,aes(x=x,y=y),colour = "cyan")+geom_path(data=panthera53,aes(x=x,y=y),colour = "purple")+geom_path(data=panthera55,aes(x=x,y=y),colour = "yellow")+geom_path(data=panthera59,aes(x=x,y=y),colour = "chartreuse")+geom_path(data=panthera60,aes(x=x,y=y),colour = "deepskyblue")+geom_path(data=panthera61,aes(x=x,y=y),colour = "chocolate")+geom_path(data=panthera28,aes(x=x,y=y),colour = "orange")+geom_path(data=panthera54,aes(x=x,y=y),colour = "aquamarine")+geom_path(data=panthera29,aes(x=x,y=y),colour = "green")+geom_path(data=panthera56,aes(x=x,y=y),colour = "beige")+geom_path(data=panthera57,aes(x=x,y=y),colour = "red")+geom_path(data=taiama12,aes(x=x,y=y),colour = "white")+geom_path(data=taiama13,aes(x=x,y=y),colour = "orange")+geom_path(data=taiama18,aes(x=x,y=y),colour = "green")+geom_path(data=taiama22,aes(x=x,y=y),colour = "blue")+geom_path(data=taiama23,aes(x=x,y=y),colour = "deeppink")+geom_path(data=taiama41,aes(x=x,y=y),colour = "red")+geom_path(data=taiama52,aes(x=x,y=y),colour = "cyan")+geom_path(data=taiama81,aes(x=x,y=y),colour = "purple")+geom_path(data=taiama88,aes(x=x,y=y),colour = "aquamarine")+geom_path(data=taiama91,aes(x=x,y=y),colour = "yellow")+geom_path(data=taiama92,aes(x=x,y=y),colour = "beige")+geom_path(data=taiama116,aes(x=x,y=y),colour = "coral1")+geom_path(data=taiama117,aes(x=x,y=y),colour = "chartreuse")+geom_path(data=oncafari84,aes(x=x,y=y),colour = "white")+geom_path(data=oncafari15,aes(x=x,y=y),colour = "orange")+geom_path(data=oncafari79,aes(x=x,y=y),colour = "green")+geom_path(data=oncafari19,aes(x=x,y=y),colour = "blue")+geom_path(data=oncafari25,aes(x=x,y=y),colour = "deeppink")+geom_path(data=oncafari68,aes(x=x,y=y),colour = "red")+geom_path(data=oncafari69,aes(x=x,y=y),colour = "cyan")+geom_path(data=oncafari86,aes(x=x,y=y),colour = "purple")+geom_path(data=oncafari87,aes(x=x,y=y),colour = "aquamarine")+geom_path(data=oncafari14,aes(x=x,y=y),colour = "yellow")
plot(BrPantanal$individual.local.identifier..ID.~BrPantanal$date,xlab = "Sampling period",ylab = "Individuals ID Number")

###################################################################################################
   ### Greater Lacandona Ecosystem(Mex)   44 to 48
jaguar_df <- as(move.data,'data.frame')
Mexico=subset(jaguar_df,country=='Mexico')
glMex44=subset(Mexico,id=='44')
glMex45=subset(Mexico,id=='45')
glMex46=subset(Mexico,id=='46')
glMex47=subset(Mexico,id=='47')
glMex48=subset(Mexico,id=='48')
GLMex=rbind(glMex44,glMex45,glMex46,glMex47,glMex48)
m <- get_map(bbox(extent(GLMex)*1.1), source="google", zoom=10,maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=glMex44,aes(x=x,y=y),colour = "orange")+geom_path(data=glMex46,aes(x=x,y=y),colour = "blue")+geom_path(data=glMex47,aes(x=x,y=y),colour = "white")+geom_path(data=glMex48,aes(x=x,y=y),colour = "purple")+geom_path(data=glMex45,aes(x=x,y=y),colour = "red")
plot(GLMex$individual.local.identifier..ID.~GLMex$date,xlab = "Sampling period",ylab = "Individuals ID Number")

################
  #### Mexico East
#X49
X49=unstacked$X49
X49_df <- as(X49, "data.frame")
head(X49_df)
nrow(X49_df)
range(X49_df$date)
m <- get_map(bbox(extent(X49)*1.1), source="google", zoom=11, maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=X49_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
plot(X49_df$individual.local.identifier..ID.~X49_df $date,xlab = "Sampling period",ylab = "Individuals ID Number")

####################### Mexico Sonora #################
jaguar_df <- as(move.data,'data.frame')
Mexico=subset(jaguar_df,country=='Mexico')
Sonora43=subset(Mexico,id=='43')
Sonora64=subset(Mexico,id=='64')
Sonora=rbind(Sonora43,Sonora64)
m <- get_map(bbox(extent(Sonora)*1.1), source="google", zoom=9,maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=Sonora64,aes(x=x,y=y),colour = "orange")+geom_path(data=Sonora43,aes(x=x,y=y),colour = "yellow")
plot(Sonora$individual.local.identifier..ID.~Sonora$date,xlab = "Sampling period",ylab = "Individuals ID Number")

####################### Mexico East & Greater Lacandona#################
jaguar_df <- as(move.data,'data.frame')
Mexico=subset(jaguar_df,country=='Mexico')
EastMex49=subset(Mexico,id=='49')
## Mexico East & Greater Lacandona ##
GL_e_EMex=rbind(glMex44,glMex45,glMex46,glMex47,glMex48,EastMex49)
m <- get_map(bbox(extent(GL_e_EMex)*1.1), source="google", zoom=7,maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=glMex44,aes(x=x,y=y),colour = "orange")+geom_path(data=glMex46,aes(x=x,y=y),colour = "blue")+geom_path(data=glMex47,aes(x=x,y=y),colour = "white")+geom_path(data=glMex48,aes(x=x,y=y),colour = "purple")+geom_path(data=glMex45,aes(x=x,y=y),colour = "red")+geom_path(data=X49_df,aes(x=x, y=y),colour = "red")
plot(GL_e_EMex$individual.local.identifier..ID.~GL_e_EMex$date,xlab = "Sampling period",ylab = "Individuals ID Number")

########## Mexico East & Greater Lacandona & Sonora ######################################
Sonora_GL_e_EMex=rbind(glMex44,glMex45,glMex46,glMex47,glMex48,EastMex49,Sonora)
m <- get_map(bbox(extent(Sonora_GL_e_EMex)*1.1), source="google", zoom=5,maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=glMex44,aes(x=x,y=y),colour = "orange")+geom_path(data=glMex46,aes(x=x,y=y),colour = "blue")+geom_path(data=glMex47,aes(x=x,y=y),colour = "white")+geom_path(data=glMex48,aes(x=x,y=y),colour = "purple")+geom_path(data=glMex45,aes(x=x,y=y),colour = "red")+geom_path(data=X49_df,aes(x=x, y=y),colour = "red")+geom_path(data=Sonora64,aes(x=x,y=y),colour = "orange")+geom_path(data=Sonora43,aes(x=x,y=y),colour = "yellow")
plot(Sonora_GL_e_EMex$individual.local.identifier..ID.~Sonora_GL_e_EMex$date,xlab = "Sampling period",ylab = "Individuals ID Number")

  
#####################################  Costa Rica ##############
#X26                     
unstacked <- split(move.data)
X26=unstacked$X26
X26_df <- as(X26, "data.frame")
head(X26_df)
nrow(X26_df)
range(X26_df$date)
m <- get_map(bbox(extent(X26)*1.1), source="google", zoom=12, maptype="satellite") #define the box for map based on the trajectory coordinates
ggmap(m)+geom_path(data=X26_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
#m <- get_map(bbox(extent(X26)*1.1), source="google", zoom=13, maptype="satellite") #define the box for map based on the trajectory coordinates
#ggmap(m)+geom_path(data=X26_df,aes(x=x, y=y),colour = "red") #plot the map and the trajectory
plot(X26_df$individual.local.identifier..ID.~X26_df$date,xlab = "Sampling period",ylab = "Individuals ID Number")


##########################################################################################################
##########################################################################################################

### Loop to produce frames for Gifs #################################################
unstacked <- split(move.data)
head(unstacked)
										   										   									   
#Steps
for(i in 1:length(unstacked))
{
	print(i)
	oc=unstacked[[i]]
	oc_df <- as(oc, "data.frame")

	print(oc_df$id[1])
	
	Oc<- move(x=oc_df$x, y=oc_df$y,time=as.POSIXct(oc_df$date,
	format = "%m/%d/%Y %H:%M", tz = 'GMT'),data=oc_df,
	proj=CRS('+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'),animal=oc_df$id[1],sensor="GPS")

	m <- get_map(bbox(extent(Oc)*1.1), source="google", zoom=14, maptype="satellite") #define the box for map based on the trajectory coordinates
	dir.create(paste0("oc","-",i))
	setwd(paste0("oc","-",i))
	for(j in 1:nrow(oc_df)){
	png(paste0("oc","-",i,"-",j,".png"))
	p <- ggmap(m)+geom_path(data=oc_df[1:j,],aes(x=x, y=y),colour = "red") + ggtitle(oc_df$id[1]) #plot the map and the trajectory
	print(p)
	dev.off()
	
	}
	setwd("..")
	
}

#End
### This only save the frames 
## To build the gifs it requires external use of other program, internet or further link it to ImageMagick
