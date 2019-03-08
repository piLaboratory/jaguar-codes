#' ---
#' title: "Fisher RSF"
#' author: "John Fieberg"
#' date: ""
#' ---

#' ### Preamble
#' 
#' Load  libraries
#+warning=FALSE, message=FALSE
library(ezknitr)
library(knitr)
library(lubridate)
library(raster)
library(move)
library(amt) 
library(broom)
library(nlme)
library(lme4)
library(tidyverse)
options(width=165)
opts_chunk$set(fig.width=12,fig.height=4.5, error=TRUE,cache = F)

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
 