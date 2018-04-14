#!/usr/bin/Rscript
#Author: "Emma Fox (eafox@ucla.edu)"

#This tutorial takes you through making a gradient forest for the plant ITS shrub/scrub data as the response variable and the 19 standard bioclim layers and elevation as the predictors.

#This tutorial takes you through making a deicision tree for the plant ITS shrub/scrub data as the response variable and the 19 standard bioclim layers and elevation as the predictors.

#############
# LIBRARIES #
#############
library(raster)
library(elevatr)
library(rgdal)

####################################
# LOAD IN COORDINATES AND METADATA #
####################################
metaDF_base<-read.csv("Base_MetaData.csv", header = TRUE)
metaDF<-metaDF_base[,c(3:4)]
metaDF<-cbind(metaDF,metaDF_base[,c(7:ncol(metaDF_base))])
row.names(metaDF)<-metaDF_base$MatchName

#########################################
# LOAD IN SPECIES DATA TABLE AND SUBSET #
#########################################
asvDF<-read.csv("CombinedASV_Tables.csv",header = TRUE)
row.names(asvDF)<-asvDF$X
asvDF<-asvDF[,2:length(asvDF)]
algaeNames<-grep("Chlorophyta*",colnames(asvDF),value = TRUE)
respVar<-subset(asvDF,select=algaeNames)

##########################
# COMBINE BOTH DATA SETS #
##########################
datDF<-merge(metaDF,respVar,by="row.names")

###########################
# GRADIENT FOREST EXAMPLE #
###########################
library(gradientForest)
eDNA_gradFor<-gradientForest(datDF,
                   predictor.vars = colnames(metaDF), 
                   response.vars = colnames(respVar),
                   ntree = 500, transform = NULL, compact = T,
                   nbin = 201, corr.threshold = 0.5)
eDNA_gradFor
plot(eDNA_gradFor,plot.type="Overall.Importance")
