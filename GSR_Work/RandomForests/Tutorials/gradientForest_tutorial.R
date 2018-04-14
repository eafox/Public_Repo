#!/usr/bin/Rscript
#Author: "Emma Fox (eafox@ucla.edu)"

library(gradientForest)
# doc page:
# http://gradientforest.r-forge.r-project.org/biodiversity-survey.pdf
# https://rdrr.io/rforge/gradientForest/man/gradientForest.html
# https://r-forge.r-project.org/projects/gradientforest


##############
# TUTORIAL 1 #
##############
# short
# https://www.r-bloggers.com/random-forests-in-r/


# long
#http://gradientforest.r-forge.r-project.org/biodiversity-survey.pdf
#http://esapubs.org/archive/ecol/E093/015/suppl-2.htm#filelist
load("~/Downloads/GZ.sps.mat.Rdata")
dim(Sp_mat)

load("~/Downloads/GZ.phys.site.Rdata")
dim(Phys_site)

nSites<-dim(Sp_mat)[1]
nSpecs<-dim(Sp_mat)[2]
lev<-floor(log2(nSites*0.368/2))

gf<-gradientForest(cbind(Phys_site,Sp_mat),
                   predictor.vars = colnames(Phys_site), 
                   response.vars = colnames(Sp_mat),
                   ntree = 500, transform = NULL, compact = T,
                   nbin = 201, maxLevel = lev, corr.threshold = 0.5)

gf

#predictor overall importance plot
plot(gf,plot.type="Overall.Importance")

#split density plot
most_important<-names(importance(gf))[1:25]
par(mgp=c(2,0.75,0))

plot(gf, plot.type="Split.Density", imp.vars=most_important, leg.posn="topright", 
     cex.legend=0.4, cex.axis=0.6, cex.lab=0.7, line.ylab=0.9, 
     par.args=list(mgp=c(1.5,0.5,0), mar=c(3.1,1.5,0.1,1)))

#species cumulative plot
plot(gf, plot.type="Cumulative.Importance", imp.vars=most_important,
     show.overall=F, legend=T, leg.posn="topleft", leg.nspecies=5, cex.lab=0.7,
     cex.legend=0.4, cex.axis=0.6, line.ylab=0.9, 
     par.args=list(mgp=c(1.5,0.5,0),mar=c(2.5,1,0.1,0.5),
                   omi=c(0,0.3,0,0)))


#predictor cumulative plot
plot(gf, plot.type = "Cumulative.Importance", imp.vars = most_important,
     show.species = F, common.scale = T, cex.axis = 0.6, cex.lab = 0.7, line.ylab = 0.9, 
     par.args = list(mgp = c(1.5,0.5, 0), mar = c(2.5, 1, 0.1, 0.5), 
                     omi = c(0,0.3, 0, 0)))

#R^2 fit for each species
plot(gf, plot.type = "Performance", show.names = F, horizontal = F,
     cex.axis = 1, cex.labels = 0.7, line = 2.5)

###############
# PREDICTIONS #
###############
#Transforming data
#originally is in a grid system
load("~/Downloads/GZ.phys.grid.Rdata")
dim(Phys_grid)
names(Phys_grid)
imp.vars<-names(importance(gf))
Trns_grid<-cbind(Phys_grid[,c("EAST","NORTH")],
                 predict(gf,Phys_grid[,imp.vars]))
Trns_site<-predict(gf)
# GO BACK OVER THIS TO UNDESTAND ###############################################################################

#Principle Component Analysis
PCs<-prcomp(Trns_grid[,imp.vars])
a1<-PCs$x[,1]
a2<-PCs$x[,2]
a3<-PCs$x[,3]
r<-a1+a2
g<- -a2
b<- a3+a2-a1

r<-(r - min(r))/(max(r) - min(r)) * 255
g<-(g - min(g))/(max(g) - min(g)) * 255
b<-(b - min(b))/(max(b) - min(b)) * 255

par(mfrow=c(1,1))
#Only a few of the most important predictors are used for pca
nvs <- dim(PCs$rotation)[1]
vec <- c("BSTRESS", "MUD", "SST_AV", "T_AV", "CHLA_AV", "SAND", "CRBNT", "GRAVEL")
lv <- length(vec)
vind <- rownames(PCs$rotation) %in% vec
scal <- 40
xrng <- range(PCs$x[, 1], PCs$rotation[, 1]/scal) * 1.1
yrng <- range(PCs$x[, 2], PCs$rotation[, 2]/scal) * 1.1
plot((PCs$x[, 1:2]), xlim = xrng, ylim = yrng, pch = ".", cex = 4, col = rgb(r, g, b, max = 255), asp = 1)
points(PCs$rotation[!vind, 1:2]/scal, pch = "+")
arrows(rep(0, lv), rep(0, lv), PCs$rotation[vec,1]/scal, PCs$rotation[vec, 2]/scal, length = 0.0625)
jit <- 0.0015
text(PCs$rotation[vec, 1]/scal + jit * sign(PCs$rotation[vec,1]), PCs$rotation[vec, 2]/scal + jit * sign(PCs$rotation[vec,2]), labels = vec)

#add location of plots
PCsites <- predict(PCs, Trns_site[, imp.vars])
points(PCsites[, 1:2])
SpsWtd <- sweep(gf$Y, 2, apply(gf$Y, 2, min), "-")
SpsWtdPCs <- (t(SpsWtd) %*% (PCsites[, 1:2]))/colSums(SpsWtd)
points(SpsWtdPCs, col = "red", pch = "+")

#plot abundance of specific species
sp <- colnames(SpsWtd)[1]
points(PCsites[, 1:2], col = "blue", cex = SpsWtd[,sp]/2)

#Map prediction in geographic space
plot(Trns_grid[, c("EAST", "NORTH")], pch = ".", cex = 3, asp = 1, col = rgb(r, g, b, max = 255))


#Clustered for inferred assemblages rather than continuous biodiv composition
library(cluster)
ncl <- 8
clPCs <- clara(PCs$x, ncl, sampsize = 1000)
medcolR <- r[clPCs$i.med]
medcolG <- g[clPCs$i.med]
medcolB <- b[clPCs$i.med]
plot((PCs$x[, 1:2]), xlim = xrng, ylim = yrng, pch = ".", cex = 4, 
     col = rgb(medcolR[clPCs$clustering],medcolG[clPCs$clustering],medcolB[clPCs$clustering], max = 255), asp = 1)
points(PCs$rotation[!vind, 1:2]/scal, pch = "+")
arrows(rep(0, lv), rep(0, lv), PCs$rotation[vec, 1]/scal, PCs$rotation[vec, 2]/scal, length = 0.0625)
text(PCs$rotation[vec, 1]/scal + jit * sign(PCs$rotation[vec, 1]), PCs$rotation[vec, 2]/scal + jit * sign(PCs$rotation[vec, 2]), labels = vec)
text(clPCs$medoids[, 1:2], labels = seq(1, ncl))
legend("bottomleft", as.character(seq(1, ncl)), pch = 15, cex = 1, col = rgb(medcolR, medcolG, medcolB, max = 255))


#Plot clusters on map
plot(Trns_grid[, c("EAST", "NORTH")], pch = ".", cex = 3, asp = 1, 
     col = rgb(medcolR[clPCs$clustering], medcolG[clPCs$clustering], medcolB[clPCs$clustering],max = 255))
points(Trns_grid[clPCs$i.med, c("EAST", "NORTH")],pch = as.character(seq(1, ncl)))
legend("bottomleft", as.character(seq(1, ncl)),pch = 15, cex = 1, 
       col = rgb(medcolR, medcolG, medcolB, max = 255))
