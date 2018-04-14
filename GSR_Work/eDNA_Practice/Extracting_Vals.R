#!/usr/bin/Rscript
#Author: "Emma Fox (eafox@ucla.edu)"

# How to get your GIS data by coordinate #
# Based off answers provided here:
# https://gis.stackexchange.com/questions/200417/very-basic-question-on-extracting-data-from-tif-raster-layer-in-r-projection-n?utm_medium=organic&utm_source=google_rich_qa&utm_campaign=google_rich_qa
# User manual
# http://rspatial.org/spatial/rst/8-rastermanip.html

library(raster)
library(elevatr)
library(rgdal)
library(reshape2)

#Import a raster layer
str_name<-"~/Downloads/wc2.0_30s_prec/wc2.0_30s_prec_01.tif" 
imported_raster<-raster(str_name)
#Import points to use
coords<-read.csv("~/Code/UCLA/Starter_Packs/RandomForests/eDNA_Practice/Fox_SampleCoords.csv")
#Change coords names to match with eventual anacapa data
newNamesA<-sub("L",".",coords$Name)
coords$MatchName<-sub("-S","",newNamesA)
#Extract data at those points
pts<-cbind(coords$Longitude,coords$Latitude)
coords$Val<-extract(imported_raster,pts)
#Crop to California
e<-extent(-126,-114,32,44)
calmap<-crop(imported_raster,e)
plot(calmap)
pts2<-SpatialPoints(pts,proj4string=CRS("+init=epsg:4326")) #puts points in a format than can be 
plot(pts2,add=TRUE)

elev_df = get_elev_point(pts2, prj = NULL, src = "epqs")

#Bring in ANACAPA output to merge with data frame
readIn<-read.csv("~/Code/UCLA/Starter_Packs/SpiecEasi/eDNA_Practice/16s.csv",stringsAsFactors = FALSE)
#remove row where otu name is NA
remNA<-subset(readIn, !is.na(readIn$sum.taxonomy))
#change row 1 to the column name
Orig <- remNA[,-1]
rownames(Orig) <- remNA[,1]
#remove asvs not present in x percent of samples
Filt<-Orig[rowSums(Orig == 0) <= (1-0.15)*length(Orig), ]
tr.Filt<-t(Filt)
#Match sample names
oldNames<-row.names(tr.Filt)
newNames<-sub(".*_", "",oldNames)
newNames2<-regmatches(newNames,regexpr("\\w*.{3}",newNames))
row.names(tr.Filt)<-newNames2
anaDat<-tr.Filt[!grepl("*Blank*",rownames(tr.Filt)),]
#Data frame to add in info by match and sample
anaDF<-melt(anaDat)
colnames(anaDF)<-c("MatchName","ASV","Abund")

#Merge the two data sets
coords_sub<-subset(coords,Set=="ShrubScrub",select=c(Set,MatchName,Val))
fullDF<-merge(anaDF,coords_sub,by="MatchName")