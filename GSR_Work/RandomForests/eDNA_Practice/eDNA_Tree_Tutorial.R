#!/usr/bin/Rscript
#Author: "Emma Fox (eafox@ucla.edu)"
#Code for extracting bioclim layers from Sabrina Shiraz

#MAKE SURE SYSTEMS MATCH UP
#DEGREES, MINUTES/SECONDS, 1 KM PIXELS

#This tutorial takes you through making a gradient forest for the plant ITS shrub/scrub data as the response variable and the 19 standard bioclim layers and elevation as the predictors.

#############
# LIBRARIES #
#############
library(raster)
library(elevatr)
library(rgdal)

#######################
# LOAD IN COORDINATES #
#######################
coords<-read.csv("~/Code/UCLA/Starter_Packs/RandomForests/eDNA_Practice/FCSs_Coords.csv")
#Change coords names to match with eventual anacapa data
newNamesA<-sub("L",".",coords$Name) #replace L with period
coords$MatchName<-sub("-S","",newNamesA) #remove everything after S
#Make list of coords to extract GIS data with and put in spatial points format
ext_coord<-cbind(coords$Longitude,coords$Latitude)
pts<-SpatialPoints(ext_coord,proj4string=CRS("+init=epsg:4326"))

#######################################################
# LOAD GIS DATA AND EXTRACT INFO FOR YOUR COORDINATES #
#######################################################
#BIOCLIM DATA
r <- getData("worldclim",var="bio",res=10)                                                 
#Use all 19 variables
r <-r[[c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19)]]
#Extract values for your points
values<-extract(r,pts)
#ELEV DATA
prj_dd <- "+proj=longlat +datum=WGS84 +no_defs"
elev_df<-get_elev_point(pts, prj = prj_dd, src = "epqs")
#Put values and points into 1 data frame
metaDF<-cbind.data.frame(coords$MatchName,coordinates(ext_coord),values,elev_df$elevation)
colnames(metaDF)[1:3]<-c("MatchName","Longitude","Latitude")
colnames(metaDF)[23]<-"elev"
metaDF<-metaDF[order(coords$MatchName),]

#############################
# LOAD IN PITS ANACAPA DATA #
#############################
pruneCSV<-function(inFile){
  readIn<-read.csv(inFile,stringsAsFactors = FALSE)
  #remove row where otu name is NA
  remNA<-subset(readIn, !is.na(readIn$sum.taxonomy))
  #change row 1 to the column name
  Orig <- remNA[,-1]
  rownames(Orig) <- remNA[,1]
  #remove asvs not present in 15% of samples
  Filt<-Orig[rowSums(Orig == 0) <= (1-0.15)*length(Orig), ]
  #transpose so each row is a sample
  tr.Filt<-t(Filt)
  #Match sample names
  oldNames<-row.names(tr.Filt)
  newNames<-sub(".*_", "",oldNames) #remove primer tag
  row.names(tr.Filt)<-regmatches(newNames,regexpr("\\w*.{3}",newNames)) #keep part of sample name that matches meta name
  #remove any blanks still present
  finalDF1<-data.frame(tr.Filt[!grepl("*neg*",rownames(tr.Filt)),])
  finalDF2<-data.frame(finalDF1[!grepl("*lank*",rownames(finalDF1)),])
  #Data frame to add in info by match and sample
  return(finalDF2)
}
# INSERT STEP SOMEWHERE TO MAKE THIS TAXON INSTEAD OF ASV ########################################

sDF<-pruneCSV("~/Code/UCLA/Starter_Packs/RandomForests/eDNA_Practice/PITS_ShrubScrub.csv")
#FIX SAMPLE 177
cDF<-pruneCSV("~/Code/UCLA/Starter_Packs/RandomForests/eDNA_Practice/PITS_Coastal.csv")
fDF<-pruneCSV("~/Code/UCLA/Starter_Packs/RandomForests/eDNA_Practice/PITS_Forest.csv")
#REMOVE DOUBLES FROM ONE OR THE OTHER?
#WOULD TAXON FIX?

################
# TREE EXAMPLE #
################
library(tree)
names(sort(colSums(sDF), decreasing = TRUE))
names(sort(colSums(cDF), decreasing = TRUE))
names(sort(colSums(fDF), decreasing = TRUE))
anaDF<-data.frame(c(sDF$Chlorophyta.Chlorophyceae.Chlamydomonadales.Chlamydomonadaceae.Chlamydomonas.,
                 cDF$Chlorophyta.Chlorophyceae.Chlamydomonadales.Chlamydomonadaceae.Chlamydomonas.,
                 fDF$Chlorophyta.Chlorophyceae.Chlamydomonadales.Chlamydomonadaceae.Chlamydomonas.))
names(anaDF)<-"grass_abund"
anaDF$MatchName<-c(row.names(sDF),row.names(cDF),row.names(fDF))
anaDF$SetName<-c(rep("ShrubScrub",length(sDF[,1])),rep("Coastal",length(cDF[,1])),rep("Forest",length(fDF[,1])))
fullDF<-merge(metaDF,anaDF,by="MatchName")

treeDF<-cbind(anaDF_sub,metaDF)
eDNA_tree<-tree(grass_abund~Latitude+Longitude+elev+
                  bio1+bio2+bio3+bio4+bio5+bio6+bio7+
                  bio8+bio9+bio10+bio11+bio12+bio13+
                  bio14+bio15+bio16+bio17+bio18+bio19,
                data = treeDF)




summary(eDNA_tree)
plot(eDNA_tree)
text(eDNA_tree, all=T)