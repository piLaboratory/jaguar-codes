#' # Fitting RSFs to individual animals
#'
#' Goals: 
#' 
#' - Illustrate data development for fitting RSF models
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
library(amt) 
library(spatstat)
library(maptools)

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

#' Lets use amt to prepare data for RSF fitting via logistic regression
trk<-with(SHd,track(x=LONGITUDE , y=LATITUDE , t=TIMESTAMP, id=IND))  


#' ## RSF prep
#' 
#' Generate random points within MCPs for a single individual using amt functions.
#' Notes:
#' 
#' - It is common to generate points randomly, but other options are possible. 
#' - In particular, it can beneficial to generate a systematically placed sample
#' - Samples can also be generated using the *spsample* function in the sp library or
#' using a GIS (note: amt using the spsample function in random_points)
#' - Other home range polygons could be used (e.g., kernel density, local convex hull
#' etc.)
#' 
#' #### Random points:
Apollo<-filter(trk, id=="Apollo")
rnd_pts <- random_points(Apollo, factor = 20)

#' Plot the mcp and random points
#+ fig.width=5, fig.height=5
plot(rnd_pts)

#' Illustrate systematic points (to do this, we need to create the mcp first). This 
#' should be read as "filter the data (choose only data from Apollo), then determine
#' the mcp, then generate systematically placed points, then plot!
#+ fig.width=5, fig.height=5
filter(trk, id=="Apollo")%>%hr_mcp()%>%random_points( factor = 20, type="regular")%>%plot() 

#' Now, lets generate points for all individuals using a loop.  We will append
#' the name of the individual to the points, then continually bind together
#' data from each individual (using rbind) as we loop throug data from each
#' individual. 
avail.pts<-NULL # set up object to hold data
uid<-unique(trk$id) # individual identifiers
luid<-length(uid) # number of unique individuals
for(i in 1:luid){
  # random_points will generate random points within mcp
  # Add on the individual id and combine all data
  temp<-cbind(id=uid[i],trk%>%filter(id==uid[i])%>%random_points)
  avail.pts<-rbind(avail.pts, temp)
}
avail.pts<-as_tibble(avail.pts)

#' Turn the data back into a track object so we can easily extract covariates.
avail.pts<- with(avail.pts, track(x=x_, y=y_, id=id, case_=case_))

#' Now, use extract_covariates to determine the spatial covariates for the rasters
#' (ndvi and veg category) for both used and available points.
#' This takes some time!!!
rsfdat<-avail.pts %>% extract_covariates(env)

#' Now create "distance to" covariates (also takes time!).  To do this
#' we will make use of the nncross function in the spatstat library.
rsf.ppp<-as.ppp(as_sp(rsfdat)) 
rsfdat$distroad<- nncross(rsf.ppp, as.psp(roads))$dist
rsfdat$distcamps<-  nncross(rsf.ppp, as.ppp(camps))$dist
rsfdat$distrivers<- nncross(rsf.ppp, as.psp(rivers))$dist
rsfdat$distcattlep<-nncross(rsf.ppp, as.ppp(cattleP))$dist

#' Scale and center all covariates by subtracting their mean and dividing by
#' their standard deviation.  Regression coefficients will give the relative 
#' risk of use per 1 sd change in each covariate.
rsfdat<-mutate(rsfdat, 
               distroad.sc=scale(distroad),
               distcamps.sc=scale(distcamps),
               distrivers.sc=scale(distrivers), 
               distcattlep.sc=scale(distcattlep),
               ndvi.sc=scale(NDVI))

#' Have to be careful with categorical variables, make sure 
#' all animals are found in all categories.  Note: veg categories 1 and 5 are
#' problematic for some individuals (e.g., Gin, Fly, Bumble, Tori).
with(rsfdat, table(id, vegetation, case_))

#' ## Explore the data:

#' Look Distribution of variables for used and available
#' ### Vegetation
#+fig.width=8, fig.height=8
ggplot(rsfdat, aes(as.factor(vegetation), y=..prop..,group=case_, colour=case_))+
  geom_bar(position="dodge", aes(fill=case_))+facet_wrap(~id, scales="free")+
  xlab("Vegatation")

#' ### NDVI
#+fig.width=8, fig.height=8
ggplot(rsfdat, aes(ndvi.sc, colour=case_))+geom_density()+facet_wrap(~id, scales="free")+
  xlab("NDVI")

#' ### Distance to nearest Road
#+fig.width=8, fig.height=8
ggplot(rsfdat, aes(distroad.sc, colour=case_))+geom_density()+facet_wrap(~id, scales="free")+
  xlab("Distance to Road")

#' ### Distance to nearest Camp
#+fig.width=8, fig.height=8
ggplot(rsfdat, aes(distcamps.sc, colour=case_))+geom_density()+facet_wrap(~id, scales="free")+
  xlab("Distance to Camp")

#' ### Distance to nearest River
#+fig.width=8, fig.height=8
ggplot(rsfdat, aes(distrivers.sc, colour=case_))+geom_density()+facet_wrap(~id, scales="free")+
  xlab("Distance to River")

#' ### Distance to nearest Cattle Post
#+fig.width=8, fig.height=8
ggplot(rsfdat, aes(distcattlep.sc, colour=case_))+geom_density()+facet_wrap(~id, scales="free")+
  xlab("Distance to Cattle Posts")
 

#' ## RSF fitting

#' Weight available data 
rsfdat$w<-ifelse(rsfdat$case_==1, 1, 5000)

#' We can fit an RSF model to a single animal using logistic regression.  Note, 
#' however, that the SEs and p-values assume all data are independent, when 
#' they are clearly not.  Rather than make inference for any particular individual,
#' we will fit models separately to all individuals & treat the regression 
#' coefficients as "data" for later population-level inference.
summary(glm(case_ ~ as.factor(vegetation)+ndvi.sc+distroad.sc+distcamps.sc+
              distrivers.sc+distcattlep.sc, data = subset(rsfdat, id=="Apollo"),
            weight=w,family = binomial))


#' Now, fit an RSF model to data from each animal. To do this, we will write
#' a function that fits the model, then use a nested data frame, purrr,
#' and tidy/broom to fit the model to each individual and summarize output.
fit_rsf <- function(data){
  mod <- glm(case_ ~ as.factor(vegetation)+ndvi.sc+distroad.sc+distcamps.sc+
               distrivers.sc+distcattlep.sc, data = data, weight=w,family = binomial)
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
rsf_coefs2 <- rsffits %>%
  tidyr::unnest(tidy) %>%
  dplyr::select(-(std.error:p.value))  
rsf_coefs2

#' ## Summary of Individual-level Models

#' Summarize mean and SE of the coefficients (treating them as data)
pop.coef<-rsf_coefs2 %>% group_by(term) %>% 
  summarize(mean.beta=mean(estimate, na.rm=TRUE),
            se.mean.beta=sd(estimate, na.rm=TRUE)/sqrt(n()))

#' Plot coefficients
#+ fig.height=16, fig.width=16
ggplot(filter(rsf_coefs2, term!="(Intercept)"), aes(x=term, y=estimate, col=id)) +
  geom_dotplot(binaxis="y", stackdir="center")+geom_hline(yintercept=0) +
  theme(text = element_text(size = 25))+ theme(axis.text.x = element_text(angle = 90))+
  scale_x_discrete(breaks= unique(pop.coef$term)[-1],
                   labels=c("Veg2", "Veg3", "veg4", "Veg5", "Veg5", "Dist Camp",
                            "Dist CattleP","Dist Rivers","Dist Roads","NDVI"))

#' ## Population-level Summaries
#' 
#' Plot of mean coefficient in the population
#+ fig.height=12, fig.width=16
ggplot(filter(pop.coef, term!="(Intercept)"), aes(x=term, y=mean.beta)) +
    geom_errorbar(aes(ymin=mean.beta-1.96*se.mean.beta, ymax=mean.beta+1.96*se.mean.beta), width=0.1)+
    geom_point()+xlab("Predictor")+ylab(expression(paste("Mean  ", hat(beta))))+
    geom_hline(yintercept=0) + theme(text = element_text(size = 25))+
  theme(axis.text.x = element_text(angle = 90))+
  scale_x_discrete(breaks= unique(pop.coef$term)[-1],
                   labels=c("Veg2", "Veg3", "veg4", "Veg5", "Veg5", "Dist Camp",
                            "Dist CattleP","Dist Rivers","Dist Roads","NDVI"))
    
#' ## Document Footer	
#' 	
#' 	
#' Session Information:	
#' 	
devtools::session_info()  
