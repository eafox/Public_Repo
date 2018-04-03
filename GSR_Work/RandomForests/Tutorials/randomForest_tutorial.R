#!/usr/bin/Rscript
#Author: "Emma Fox (eafox@ucla.edu)"

library(randomForest)

##############
# TUTORIAL 1 #
##############
# short
# https://www.r-bloggers.com/random-forests-in-r/

#get data set
library(MASS)
attach(Boston)
set.seed(101)
dim(Boston)
#Separate out training from test data
train<-sample(1:nrow(Boston),300)
#Fit random forest with all predictors (.)
Boston.rf<-randomForest(medv ~ . , data = Boston , subset = train)
Boston.rf
#error vs number of trees considered
plot(Boston.rf)
#Looking at out of bag sample erros and error on test set for all 13 predictors
oob.err=double(13)
test.err=double(13)
#mtry is no of Variables randomly chosen at each split
for(mtry in 1:13) {
  rf=randomForest(medv ~ . , data = Boston , subset = train,mtry=mtry,ntree=400) 
  oob.err[mtry] = rf$mse[400] #Error of all Trees fitted
  pred<-predict(rf,Boston[-train,]) #Predictions on Test Set for each Tree
  test.err[mtry]= with(Boston[-train,], mean( (medv - pred)^2)) #Mean Squared Test Error
  cat(mtry," ") #printing the output to the console
}
test.err ; oob.err
#plot the errors for comparison
matplot(1:mtry , cbind(oob.err,test.err), pch=19 , col=c("red","blue"),type="b",ylab="Mean Squared Error",xlab="Number of Predictors Considered at each Split")
legend("topright",legend=c("Out of Bag Error","Test Error"),pch=19, col=c("red","blue"))

##############
# TUTORIAL 2 #
##############
# http://www.listendata.com/2014/11/random-forest-with-r.html

#DATA PREP
mydata=read.csv("https://sites.google.com/site/pocketecoworld/german_credit.csv")
# Check types of variables
str(mydata)
# Check number of rows and columns
dim(mydata)
# Make dependent variable as a factor (categorical)
mydata$Creditability = as.factor(mydata$Creditability)
#RUN RANDOM FOREST
set.seed(71) 
rf <-randomForest(Creditability~.,data=mydata, ntree=500) 
print(rf)
#OPTIMAL mtry VALUE
# -1 because dependent variable is in data set
mtry <- tuneRF(mydata[-1],mydata$Creditability, ntreeTry=500,
               stepFactor=1.5,improve=0.01, trace=TRUE, plot=TRUE)
best.m <- mtry[mtry[, 2] == min(mtry[, 2]), 1]
print(mtry)
print(best.m)
#BUILD WITH OPTIMAL mtry
set.seed(71)
rf <-randomForest(Creditability~.,data=mydata, mtry=best.m, importance=TRUE,ntree=500)
print(rf)
#EXAMINE VARIABLE IMPORTANCE
importance(rf)
varImpPlot(rf)

##############
# TUTORIAL 3 #
##############
#Load necessary libraries
library(raster)
library(tree)
library(rgeos)
library(rgdal)
library(maptools)