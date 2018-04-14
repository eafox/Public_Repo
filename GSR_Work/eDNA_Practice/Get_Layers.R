#!/usr/bin/Rscript
#Editor: "Emma Fox (eafox@ucla.edu)"
#Code Provided by: Sabrina Shiraz

#MAKE SURE SYSTEMS MATCH UP
#DEGREES, MINUTES/SECONDS, 1 KM PIXELS

library(raster)
library(elevatr)
library(rgdal)

################
# BASIC LAYERS #
################
#make matrix of coordinates
#can also import a matrix of coordinates longitude, latitude
coord<-read.csv("02_datCoords.csv")
#add in metadata from eDNA collection (soil or sediment)
SoSDF<-read.csv("SoS_base.csv")
SoSDF2<-unique(merge(coord,SoSDF,by="MatchName",all.x=TRUE))


# FORMAT COORDINATE SET IN LONGITUDE THEN LATITUDE
ext_coord<-cbind(coord$Longitude,coord$Latitude)
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
metaDF<-cbind.data.frame(coord$MatchName,coordinates(ext_coord),coord$Set,SoSDF2$SoS,bc_values,elev_df$elevation)
colnames(metaDF)[1:5]<-c("MatchName","Longitude","Latitude","Set","SoS")
colnames(metaDF)[25]<-"elev"
metaDF<-metaDF[order(coord$MatchName),]

#Write CSV
write.csv(metaDF,"Base_MetaData.csv")

##############################
# GETTING VALUES FROM LAYERS #
##############################
#Have a landcover tiff
lanCoDF<- raster("LC_5min_global_2012.tif")
#View the California section
Cal_LandCo<-crop(lanCoDF, extent(-126,-114,32,44))
plot(Cal_LandCo)
#Transform coordinates to spatial object and then match format to raster layer
spaPTs<-SpatialPoints(ext_coord,proj4string=CRS("+init=epsg:4326"))
ptsLC<-spTransform(spaPTs, projection(lanCoDF))
#Plot points to california map
plot(ptsLC, add=TRUE)
#Extract points and add to meta-data-frame
LanCo<-extract(lanCoDF, ptsLC)
metaDF_full<-cbind(metaDF,LanCo)
#TRANSFORM ME TO A FACTOR###############################################################################


#NDVI FOR APRIL 2017#
#Read in
NDVIdf<-raster("MOD_NDVI_M_2017-12-01_rgb_3600x1800.TIFF")
#Plot 
CalCropped<-crop(NDVIdf, extent(-126,-114,32,44))
plot(CalCropped)
#Extract points
ptsND<-spTransform(spaPTs, projection(NDVIdf))
NDVI<-extract(NDVIdf,ptsND)
#Add to data frame
metaDF_full<-cbind(metaDF_full,NDVI)

#SOIL TYPE#
soilTypeLayer<-readOGR(dsn = "USDA_Soil/", layer = toString("wss_gsmsoil_CA_[2016-10-13]_spatial_gsmsoilmu_a_ca.shp"))









#############################################
# FUNCTION FOR EXTRACTION WHEN SINGLE LAYER #
#############################################
extractVals<-function(coords,layer){
  #Coords must be a spatial object
  newLayer<-raster(layer)
  newPts<-spTransform(coords,projection(newLayer))
  vals<-extract(newLayer,newPts)
  return(vals)
}

layerList<-c()
for (i in 1:length(layerList)){
  metaDF_full<-cbind(metaDF_full,extractVals(spaPTs,layerList[[i]]))
}
colnames(metaDF_full)[25:length(metaDF_full)]<-c("","")









