#!/usr/bin/Rscript
#Editor: "Emma Fox (eafox@ucla.edu)"
#Code Provided by: Sabrina Shiraz

#MAKE SURE SYSTEMS MATCH UP
#DEGREES, MINUTES/SECONDS, 1 KM PIXELS

library(raster)
#matrix to extract the 19 bioclimateic variables... can shorten this if just interested in a few. 
r <- getData("worldclim",var="bio",res=10)                                                 
r <-r[[c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19)]]
#make matrix of coordinates
#can also import a matrix of coordinates longitude, latitude
coord = read.csv("~/Code/UCLA/Starter_Packs/RandomForests/eDNA_Practice/Fox_SampleCoords.csv")
ext_coord=cbind(coord$Longitude,coord$Latitude)
values = extract(r,ext_coord)
df= cbind.data.frame(coordinates(ext_coord),values)

