#!/usr/bin/Rscript
#Editor: "Emma Fox (eafox@ucla.edu)"
#Code Provided by: Sabrina Shiraz

#MAKE SURE SYSTEMS MATCH UP
#DEGREES, MINUTES/SECONDS, 1 KM PIXELS

library(raster)
library(elevatr)
library(rgdal)
#make matrix of coordinates
#can also import a matrix of coordinates longitude, latitude
coord = read.csv("02_datCoords.csv")
ext_coord=cbind(coord$Longitude,coord$Latitude)
pts<-SpatialPoints(ext_coord,proj4string=CRS("+init=epsg:4326"))

#BIOCLIM DATA
#matrix to extract the 19 bioclimateic variables... can shorten this if just interested in a few. 
r <- getData("worldclim",var="bio",res=10)                                                 
r <-r[[c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19)]]
bc_values <- extract(r,pts)

#ELEV DATA
prj_dd <- "+proj=longlat +datum=WGS84 +no_defs"
elev_df<-get_elev_point(pts, prj = prj_dd, src = "epqs")

#Put values and points into 1 data frame
metaDF<-cbind.data.frame(coord$MatchName,coordinates(ext_coord),bc_values,elev_df$elevation)
colnames(metaDF)[1:3]<-c("MatchName","Longitude","Latitude")
colnames(metaDF)[23]<-"elev"
metaDF<-metaDF[order(coord$MatchName),]

#Write CSV
write.csv(metaDF,"Base_MetaData.csv")
