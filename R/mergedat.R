#' ---       
#' title: "Merge Data"
#' author: "John Fieberg"
#' date: ""
#' ---
#' 
#' 
#' ### Preamble
#rm(list= ls())                                                  ### For a fresh start
setwd("C:/RWorkDir/AnimoveRscripts")                           ### Set directory
#' 
#' **Purpose**: Merge together original Use/available data with annotated Env-Data.
#'  
#' Load libraries
#+warning=FALSE, message=FALSE
library(ezknitr)
library(knitr)
library(lubridate)
library(raster)
library(move)
library(amt) 
library(tidyverse)
options(width=150)

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


## SSF data

#' Read in original data (used and available points) and merge on environmental
#' data.
#ssfdattemp<-read.csv("data/AllStepsFisher2018.csv")
#annotated<-read.csv("data/FisherSSFannotate.csv-1402112999909362686.csv")
ssfdattemp<-read.csv("c:/RWorkDir/AnimoveRscripts/AllStepsFisher2018.csv")
annotated<-read.csv("c:/RWorkDir/AnimoveRscripts/FisherSSFannotate.csv-1402112999909362686.csv")

#' Now, merge these by "timestamp", "location-long", "location-lat")
ssfdat<-merge(ssfdattemp,annotated)

#Gambiarra
#ssfdat<-merge(annotated,ssfdattemp,by = c("location.lat","location.long"))
#SSF <- rbind(data frameA, data frameB)
#ssfdat<- na.omit(SSF[16]) 
#SSF1=subset(SSF,timestamp=='NA')
#iop=subset(animals,study.name=='IOP')
#SSF<-merge(annotated2,ssfdattemp,by = c("timestamp","location.long","location.lat"))
#SSF1 <- na.omit(SSF[16]) 
#ssfdat
#SSF<-merge(ssfdattemp, annotated, all=TRUE)
#SSF<-rbind.fill(annotated,ssfdattemp)
#aggregate(x=carros.numeros,by=list(carros.marcas),FUN=mean)
#Subs1 <- na.omit(DATA[2:3]) 
#library(dplyr)
#SSF<-match(ssfdattemp,annotated, by = c("timestamp", "location.long ", "location.lat"))
#cols<- intersect(col(ssfdattemp), col(annotated))
#annotatedSSF<-merge(ssfdattemp[,cols], annotated[,cols])
#library(prodlim)
#SSF<-row.match(ssfdat, annotated, nomatch = NA)
#SSF<-merge(annotated,ssfdat)
#SSF<-match(annotated,ssfdat)
#cols <- intersect(colnames(df1), colnames(df2))
#rbind(df1[,cols], df2[,cols])
#annotatedSSF<-plyr::rbind.fill(annotated, annotated2)


#' Write out data for use in FisherSSF2018.R
#write.csv(ssfdat, file="data/AllStepsFisher2018-EnvDATA-results.csv", row.names=FALSE)
write.csv(ssfdat, file="c:/RWorkDir/AnimoveRscripts/AllStepsFisher2018-EnvDATA-results.csv", row.names=FALSE)
