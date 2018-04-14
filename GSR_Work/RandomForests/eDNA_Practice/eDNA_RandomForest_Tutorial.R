#!/usr/bin/Rscript
#Author: "Emma Fox (eafox@ucla.edu)"

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
metaDF<-read.csv("Base_MetaData.csv", header = TRUE)
row.names(metaDF)<-metaDF$MatchName
metaDF<-metaDF[,3:ncol(metaDF)]

#########################################
# LOAD IN SPECIES DATA TABLE AND SUBSET #
#########################################
asvDF<-read.csv("CombinedASV_Tables.csv",header = TRUE)
row.names(asvDF)<-asvDF$X
asvDF<-asvDF[,2:length(asvDF)]
#Tree works with one response variable so we'll use a particularly prevalent green algae
respVar<-subset(asvDF,select="Chlorophyta.Chlorophyceae.Chlamydomonadales.Chlamydomonadaceae.Chlamydomonas.")
colnames(respVar)<-"Abundance"

##########################
# COMBINE BOTH DATA SETS #
##########################
datDF<-merge(metaDF,respVar,by="row.names")
datDF<-na.omit(datDF)

#########################
# RANDOM FOREST EXAMPLE #
#########################
library(randomForest)
#Regression tree
eDNA_ranFor<-randomForest(x=datDF[,c(2:25)],y=datDF$Abundance, importance=TRUE)
eDNA_ranFor
class(eDNA_ranFor)
str(eDNA_ranFor)
names(eDNA_ranFor)
## Inspect the confusion matrix of the OOB error assessment
eDNA_ranFor$confusion
varImpPlot(eDNA_ranFor)

#Classification tree
#Transform to presence absence
datDF$Abundance[datDF$Abundance > 0] <- "Present"
datDF$Abundance[datDF$Abundance == 0] <- "Absent"
datDF$Abundance<-as.factor(datDF$Abundance)

eDNA_ranFor<-randomForest(x=datDF[,c(2:25)],y=datDF$Abundance, importance=TRUE)
eDNA_ranFor
class(eDNA_ranFor)
str(eDNA_ranFor)
names(eDNA_ranFor)
## Inspect the confusion matrix of the OOB error assessment
eDNA_ranFor$confusion
varImpPlot(eDNA_ranFor)
