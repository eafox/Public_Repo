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
asvDF<-read.csv("04b_fullDF_naRM.csv",header = TRUE)
row.names(asvDF)<-asvDF$X
asvDF<-asvDF[,2:length(asvDF)]
#Tree works with one response variable so we'll use a particularly prevalent green algae
respVar<-subset(asvDF,select="Chlorophyta.Chlorophyceae.Chlamydomonadales.Chlamydomonadaceae.Chlamydomonas.")
colnames(respVar)<-"Abundance"

##########################
# COMBINE BOTH DATA SETS #
##########################
datDF<-merge(metaDF,respVar,by="row.names")

################
# TREE EXAMPLE #
################
library(tree)
#Regression tree
eDNA_tree<-tree(Abundance~Latitude+Longitude+elev+
                  bio1+bio2+bio3+bio4+bio5+bio6+bio7+
                  bio8+bio9+bio10+bio11+bio12+bio13+
                  bio14+bio15+bio16+bio17+bio18+bio19,
                data = datDF)

summary(eDNA_tree)
plot(eDNA_tree)
text(eDNA_tree, all=T)

#Classification tree
#Transform to presence absence
datDF$Abundance[datDF$Abundance > 0] <- "Present"
datDF$Abundance[datDF$Abundance == 0] <- "Absent"
datDF$Abundance<-as.factor(datDF$Abundance)
eDNA_tree<-tree(Abundance~Latitude+Longitude+elev+
                  bio1+bio2+bio3+bio4+bio5+bio6+bio7+
                  bio8+bio9+bio10+bio11+bio12+bio13+
                  bio14+bio15+bio16+bio17+bio18+bio19,
                data = datDF)
summary(eDNA_tree)
plot(eDNA_tree)
text(eDNA_tree, all=T)
#See if you should prune down the tree (see how much each parameter matters)
cv.tree(eDNA_tree)
#The best tree size appears to be 6 nodes. Why?
ecoli.tree.2<-prune.misclass(ecoli.tree.base, best = 6)
summary(ecoli.tree.2)
plot(ecoli.tree.2)
text(ecoli.tree.2, all=T)
