#!/usr/bin/Rscript
#Author: "Emma Fox (eafox@ucla.edu)"

library(tree)

##############
# TUTORIAL 1 #
##############
# CODE AND DATA COPIED FROM:
# https://www.r-bloggers.com/classification-trees/

#Load in data set
ecoli.df<-read.table("ecoli.txt")
colnames(ecoli.df)<-c("Sequence","mcv","gvh","lip","chg","aac","a1m1","a1m2","class")
head(ecoli.df)

#View number of types in each column
xtabs(~class, data=ecoli.df)

#Make a tree with all variables
ecoli.tree.base<-tree(class~mcv+gvh+lip+chg+aac+a1m1+a1m2, data = ecoli.df)
summary(ecoli.tree.base)
plot(ecoli.tree.base)
text(ecoli.tree.base, all=T)

#Prune tree using cross validation to see how much each parameter matters
cv.tree(ecoli.tree.base)
#The best tree size appears to be 6 nodes. Why?
ecoli.tree.2<-prune.misclass(ecoli.tree.base, best = 6)
summary(ecoli.tree.2)
plot(ecoli.tree.2)
text(ecoli.tree.2, all=T)

##############
# TUTORIAL 2 #
##############
# http://plantecology.syr.edu/fridley/bio793/cart.html

#Import data and get a red oak subset.
#Reformat data to get a dataframe of red oak presence and absence
dat<-read.csv("../Tutorials/treedata.csv")
dat2<-subset(dat,dat$species=="Quercus rubra") 
ro.plots<- unique(as.character(dat2$plotID))  #plots containing red oak
u.plots<-unique(as.character(dat$plotID))  #list of all plots
nob.plots<-u.plots[is.element(u.plots,ro.plots)==F] #plots with no red oak
dat3<-subset(dat,is.element(as.character(dat$plotID),nob.plots)) 
dat4<-subset(dat3,duplicated(dat3$plotID)==F)  #one row is now a plot
dat4$cover<-rep(0,nrow(dat4))     #cover of red oak is zero in these plots
oak.df<-rbind(dat2,dat4)     #new dataframe of presences and absences

#Regression tree of elevation to explain presence absence
library(tree)
rt1=tree(cover~elev,data=oak.df)
rt1 ; summary(rt1)
#summary shows the amount of deviance explained by each terminal node
#verify the root node deviance is the total deviance
sum(sapply(oak.df$cover,function(x)(x-mean(oak.df$cover))^2)) 
#residual deviance (remaining unexplained deviance)
sum(sapply(resid(rt1),function(x)(x-mean(resid(rt1)))^2)) 
#plot the tree
plot(rt1); text(rt1)
#plot the predicted relationship
with(oak.df,plot(elev,cover))
x<-seq(0,2000)
lines(x,predict(rt1,newdata=list(elev=x)),col="lightblue",lwd=3)
#Examining how tree splits when the percent deviance that must be explained threshold is lowered
#Plots deviance explained vs number of terminal nodes
rt2<-tree(cover~elev,data=oak.df,control=tree.control(1039,mindev=.003))
plot(prune.tree(rt2))
abline(v=3,col="red")
#test data is used to determine where to put the threshold of number of terminal nodes
#Plot training tree
plot(cv.tree(rt2))
#Compare to a GLM
glm1<-glm(cover~elev+I(elev^2),data=oak.df,family=poisson)
with(oak.df,plot(elev,cover))
x<-seq(0,2000)
lines(x,predict(rt1,newdata=list(elev=x)),col="lightblue",lwd=3)
lines(x,predict(glm1,newdata=list(elev=x),type="response"),col="orange",lwd=3)
#Compare psuedo r^2
1-(deviance(rt1)/7405)
1-(glm1$dev/glm1$null)

#A categorical predictor variable
rt3<-tree(cover~disturb,data=oak.df)
rt3
plot(rt3);text(rt3)

#More complex model
ct1<-tree(factor(cover>0)~utme+utmn+elev+tci+streamdist+beers,data=oak.df)
ct1
plot(ct1); text(ct1)
#Decide if you should prune
plot(cv.tree(ct1))
btree<-prune.tree(ct1,best=7)
btree
summary(btree)