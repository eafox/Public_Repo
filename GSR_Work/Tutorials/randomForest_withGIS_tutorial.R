#!/usr/bin/Rscript
#Author: "Emma Fox (eafox@ucla.edu)"

##############
# TUTORIAL 1 #
##############
#tutorial: https://geoscripting-wur.github.io/AdvancedRasterAnalysis/
#data files found here: https://github.com/GeoScripting-WUR/AdvancedRasterAnalysis.git
## ALL CODE COPIED FROM THE TUTORIAL CITED ABOVE ##

#Load necessary libraries
library(raster)
library(randomForest)

#Load data
load("GewataB2.rda")
load("GewataB3.rda")
load("GewataB4.rda")

#Examining the loaded data
GewataB2 #attributes
summary(GewataB2)

#Combine into one brick
g_datBrick<-brick(GewataB2, GewataB3, GewataB4)
hist(g_datBrick, xlim = c(0, 5000), ylim = c(0, 750000), breaks = seq(0, 5000, by = 100)) #make sure bins are same size
pairs(g_datBrick) #compare layers

#Calc normalized difference vegetation index (NDVI)
ndvi <- overlay(GewataB4, GewataB3, fun=function(x,y){(x-y)/(x+y)})
par(mfrow = c(1, 1))
plot(ndvi)

#Load in Vegetation Continuous Field (VCF) 
load("vcfGewata.rda")
vcfGewata
plot(vcfGewata)
#Replace values above 100 with NA since these are clouds or water
vcfGewata[vcfGewata>100]<-NA
plot(vcfGewata)

#Add new values to existing brick
g_datBrick<-calc(g_datBrick, fun = function(x) x/10000) #scale current brick
covs<-addLayer(g_datBrick,ndvi,vcfGewata)
names(covs)<-c("band2", "band3", "band4", "NDVI", "VCF")
plot(covs)

#Load training data sections
load("trainingPoly.rda")
plot(ndvi)
plot(trainingPoly,add=TRUE)
#Inspect training object
trainingPoly@data
trainingPoly@data$Class
str(trainingPoly@data$Class)
#convert to numeric
trainingPoly@data$Code <- as.numeric(trainingPoly@data$Class)
trainingPoly@data
#Add those values to a raster layer
classes <- rasterize(trainingPoly, ndvi, field='Code') #ignore warnings
#View type raster layer
cols <- c("orange", "dark green", "light blue")
plot(classes, col=cols, legend=FALSE)
legend("topright", legend=c("cropland", "forest", "wetland"), fill=cols, bg="white")
#View other raster layers that overlap with that type
covmasked <- mask(covs, classes)
plot(covmasked)

#make training brick with predictor and response variable
names(classes) <- "class"
trainingbrick <- addLayer(covmasked, classes)
plot(trainingbrick)
#put values into a matrix
valuetable <- getValues(trainingbrick)
#get rid of na rows
valuetable <- na.omit(valuetable)
#convert to data frame
valuetable <- as.data.frame(valuetable)
head(valuetable, n = 10)
tail(valuetable, n = 10)
#convert class back to a factor
valuetable$class <- factor(valuetable$class, levels = c(1:3))

#Okay up to here
#Run random forest model fitting
modelRF <- randomForest(x=valuetable[ ,c(1:5)], y=valuetable$class,
                        importance = TRUE)
#Inspect Results
modelRF
class(modelRF)
str(modelRF)
names(modelRF)
## Inspect the confusion matrix of the OOB error assessment
modelRF$confusion
# to make the confusion matrix more readable
colnames(modelRF$confusion) <- c("cropland", "forest", "wetland", "class.error")
rownames(modelRF$confusion) <- c("cropland", "forest", "wetland")
modelRF$confusion
#examine the importance of each variable in making the determination
varImpPlot(modelRF)

#Use model for prediction
#check input and output names match
names(covs)
names(valuetable)
#Run prediction for full raster layer
predLC <- predict(covs, model=modelRF, na.rm=TRUE)
#Plot prediction
cols <- c("orange", "dark green", "light blue")
plot(predLC, col=cols, legend=FALSE)
legend("bottomright", 
       legend=c("cropland", "forest", "wetland"), 
       fill=cols, bg="white")

