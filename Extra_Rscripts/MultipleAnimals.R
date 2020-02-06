#' ---
#' title: "Multiple RSF's with Fisher data"
#' author: "John Fieberg"
#' date: ""
#' ---

#' **Purpose**:  demonstrate methods for analyzing data from multiple animals

#' #### Preamble
#' 
#' Load libraries
#+warning=FALSE, message=FALSE
library(knitr)
library(lubridate)
library(raster)
library(move)
library(amt) 
library(broom)
library(nlme)
library(lme4)
library(tidyverse)
library(geepack)
library(tictoc)
options(width=165)
opts_chunk$set(fig.width=12,fig.height=4.5, error=TRUE,cache = FALSE)


#' Record time for running all code
tic("total")

#' Read in annotated available data for RSF modeling
rsfdat<-read.csv("data/FisherRSF2018-EnvDATA-results.csv")

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


#' Weight available data 
rsfdat$w<-ifelse(rsfdat$case_==1, 1, 5000)


#' Now, fit an RSF model to data from each animal.  Since not all animals experience
#' all habitat types, lets just explore forest versus non-forest
rsfdat$forest<-ifelse(rsfdat$landC=="forest", 1, 0)

#' # Population-level patterns
#' 
#' ## Cluster-level bootstrap
#' 
#' Use 100 bootstraps for illustrative purposes
#' (ideally, you would want many more bootstrap resamples)
tic("bootstrap")
nboot<-500 
beta.hat<-matrix(NA,nboot, 3)
uids<-unique(rsfdat$id)
n.uids<-length(uids)
for(i in 1:nboot){
   # reasample individuals
   ids.boot<-data.frame(id=sample(uids, n.uids, replace=T))
   
#   # Take all obs from these individuals
   bootdat<-merge(ids.boot,rsfdat)
   
#   # Now fit lm and pull off coeficients
   boot.fit<-glm(case_~elev+popD+forest, weight=w,family=binomial(), data=bootdat)
   beta.hat[i,]<-coef(boot.fit)[-1]
}
 colnames(beta.hat)<-c("elev", "popD", "forest")
 apply(beta.hat, 2, mean)
 apply(beta.hat, 2,sd)
 
 #' Compare to glm
 summary(glm(case_~elev+popD+forest,weight=w,family=binomial(), data=rsfdat))
 toc()
 
#' ## GEE
#'
#' Here is the call, but in my first attempt, it caused R to crash
#geefit<-geeglm(case_~elev+popD+forest,family=binomial(), id=id, corstr="independence", data=bootdat)


#' Mixed models, random intercept only. fixed intercepts, random slopes. This takes a long time to fit
tic("glmerfit.ri")
rsfdat$id<-as.factor(rsfdat$id)
glmerfit.ri<-glmer(case_~elev+popD+forest+ (1|id), weight=w,
                   family=binomial(), data=rsfdat)
summary(glmerfit.ri)
toc()

#' A better model:  fixed intercepts and random slopes. This takes a long time to fit
#' and doesn't converge
tic("glmerfit.rc")
glmerfit.rc<-glmer(case_~as.factor(id)+elev+popD+forest+(0+forest+elev+popD|id),
                family=binomial(),  weight=w, data=rsfdat)
summary(glmerfit.rc)
toc()


#' Read in coefficients from fit to individuals   
load("data/rsfcoefs.Rdata")
rsf_coefs<-rsf_coefs %>% tidyr::spread(term, estimate)
   
#' Compare fixed  effects summaries:
rsf_coefs<-rsf_coefs %>% select("elev", "forest", "popD")
apply(rsf_coefs,2,mean) # from fit to individuals
apply(rsf_coefs,2,sd) # from fit to individuals

#' Compare estimates for each animal
rsf_coefs
fixef(glmerfit.rc)[9:11]+ranef(glmerfit.rc)$id
  
#' Total Elapsed time      
toc()

#' ## Document Footer	
   
#' Session Information:	
#' 	
sessionInfo()	  
   