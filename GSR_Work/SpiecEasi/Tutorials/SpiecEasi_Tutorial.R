#!/usr/bin/Rscript

#Author: "Emma Fox (eafox@ucla.edu)"

# INSTALLATION #
source('http://bioconductor.org/biocLite.R')
biocLite('phyloseq')
library(devtools)
install_github("ltipton/SpiecEasi")

#####################
# AMERICAN GUT DATA #
#####################
library(SpiecEasi)
library(phyloseq)
load('phydata.RData')
#Number of cores
nc<-4
#Control vs IBD
se.mb.amgut <- spiec.easi(amgut1.filt, method='mb', lambda.min.ratio=1e-2, 
                          nlambda=20, icov.select.params=list(rep.num=50))
se.gl.amgut <- spiec.easi(amgut1.filt, method='glasso', lambda.min.ratio=1e-2,
                          nlambda=20, icov.select.params=list(rep.num=50))
sparcc.amgut <- sparcc(amgut1.filt)

## Define arbitrary threshold for SparCC correlation matrix for the graph
sparcc.graph <- abs(sparcc.amgut$Cor) >= 0.3
diag(sparcc.graph) <- 0
sparcc.graph <- Matrix(sparcc.graph, sparse=TRUE)

## Create igraph objects
library(igraph)
ig.mb <- adj2igraph(se.mb.amgut$refit)
ig.gl <- adj2igraph(se.gl.amgut$refit)
ig.sparcc <- adj2igraph(sparcc.graph)

## set size of vertex proportional to clr-mean
vsize <- rowMeans(clr(amgut1.filt, 1))+6
am.coord <- layout.fruchterman.reingold(ig.mb)

par(mfrow=c(1,3))
plot(ig.mb, layout=am.coord, vertex.size=vsize, vertex.label=NA, main="MB")
plot(ig.gl, layout=am.coord, vertex.size=vsize, vertex.label=NA, main="glasso")
plot(ig.sparcc, layout=am.coord, vertex.size=vsize, vertex.label=NA, main="sparcc")

#Exporting graphs with Matrix
library(Matrix)
elist.unweight<-summary(as(se.amgut.ctl$refit,"symmetricMatrix"))
elist.weighted<-summary(symBeta(getOptBeta(se.amgut.ctl),mode='maxabs'))

#Network degree statistics
dd.ctl<-degree_distribution(ig.ctl,cumulative=FALSE)
dd.ibd<-degree_distribution(ig.ibd,cumulative=FALSE)
## ave degree
sum(seq_along(dd.ibd)*dd.ibd)-1
## [1] 2.340741
sum(seq_along(dd.ctl)*dd.ctl)-1
## [1] 3.9
## plot degree distributions
plot(seq_along(dd.ibd)-1, dd.ibd,type='b',xlim=c(0,20),
     ylab="Frequency",xlab="Degree",col='red')
points(seq_along(dd.ctl)-1, dd.ctl ,type='b')
legend("topright",c("Control","IBD"),col=c("black","red"),pch=1,lty=1)

#Testing stability
natcon<-function(ig) {
  N<-vcount(ig)
  adj<-get.adjacency(ig)
  evals<-eigen(adj)$value
  nc<-log(mean(exp(evals)))
  nc/(N-log(N))
}
nc.attack<-function(ig) {
  hubord<-order(rank(betweenness(ig)),rank(degree(ig)),decreasing=TRUE)
  sapply(1:round(vcount(ig)*.8),function(i) {
    ind<-hubord[1:i]
    tmp<-delete_vertices(ig,V(ig)$name[ind])
    natcon(tmp)
  })
}
nc.ibd<-nc.attack(ig.ibd)
nc.ctl<-nc.attack(ig.ctl)
plot(seq(0,.8,len=length(nc.ctl)), nc.ctl,type='l',ylim=c(0,max(nc.ctl)),
     xlab="Proportion of removed nodes",ylab="natural connectivity")
points(seq(0,.8,len=length(nc.ibd)), nc.ibd,type='l',col='red')
hist(diff(nc.ibd),breaks=30,col='red',ylim=c(0,80),main="Fragility rates")
hist(diff(nc.ctl),breaks=30,add=TRUE)

#Cross-sectional datasets to compare condition-dependent networks
selist<-lapply(list(ccfa.ctl, ccfa.uc, ccfa.cd, hmp),function(data) {
  spiec.easi(data,method='mb',lambda.min.ratio=1e-2,nlambda=20,
             sel.criterion='stars',icov.select.params=list(rep.num=20,ncores=nc))
})
selist<-c(list(se.amgut.ctl, se.amgut.ibd), selist)

#Graphlet correlation distance
library(orca)
orbits<-c(0,2,5,7,8,10,11,6,9,4,1)+1
## Compute graphlet frequency distribution
graphlets<-lapply(selist,function(x)
  count4(summary(Matrix::triu(x$refit))[,1:2])[,orbits])
## Compute graphlet correlation matrix
gcormat<-sapply(graphlets,function(x)
  SpiecEasi::triu(cor(rbind(x,1),method="spearman")))
## Compute graphlet correlation distances
gcordist<-dist(t(gcormat))
## Compute the MDS
points<-cmdscale(gcordist)
map<-cbind(c("CTL","IBD","CTL","UC","CD","CTL"),
           c('red','red','blue','blue','blue','black'),
           c('AGP','AGP','RISK','RISK','RISK','HMP'))
legend<-unique(map[,-1])
plot(points,type='n',xlab="MDS1",ylab="MDS2")
text(points,labels=map[,1],col=map[,2])
legend("topleft", legend[,2],text.col=legend[,1])

##################
# SYNTHETIC DATA #
##################
#data parameters
depths<-colSums(amgut.ctl@otu_table@.Data)
minDepth<-summary(depths)[3]
amgut.filt<-prune_samples(depths>=minDepth, amgut.ctl)
amgut.f<-t(apply(amgut.filt@otu_table@.Data,2, norm_to_total))
## Common-scale normalized data
amgut.cs<-round(amgut.f*minDepth)
d<-ncol(amgut.cs)
n<-nrow(amgut.cs)
e<-d*2

#synthesize data
set.seed(10010)
graph<-SpiecEasi::make_graph('cluster', d, e)
Prec<-graph2prec(graph)
Cor<-cov2cor(prec2cov(Prec))
X<-synth_comm_from_counts(amgut.cs,mar=2,distr='zipois',Sigma=Cor,n=n)

#Analysis with different methods
se.mb<-spiec.easi(X,method='mb',lambda.min.ratio=1e-2,nlambda=20,
                  icov.select.params=list(rep.num=20,ncores=nc))
se.gl<-spiec.easi(X,method='glasso',lambda.min.ratio=1e-2,nlambda=20,
                  icov.select.params=list(rep.num=20,ncores=nc))
sp.est<-sparcc(X)
## define arbitrary correlation threshold, see also sparccboot
sp.net<-abs(sp.est$Cor)>=0.15
diag(sp.net)<-0

#visualize network estimated and inferred network
ig.tru<-adj2igraph(Cor*graph)
ig.mb<-adj2igraph(symBeta(getOptBeta(se.mb),mode='maxabs'))
ig.gl<-adj2igraph(se.gl$opt.cov*se.gl$refit)
ig.sp<-adj2igraph(sp.est$Cor*sp.net)
## save layout for side-by-side plotting
agcoord<-layout.fruchterman.reingold(ig.tru)
plotnet<-function(ig,coord=agcoord,...)
  plot(ig,layout=coord,vertex.size=2,vertex.label=NA,
       edge.color=ifelse(E(ig)$weight>0,'green','red'), ...)
## uncommenct for side-by-side plotting
# par(mfrow=c(2,2))
plotnet(ig.tru,main="True")
plotnet(ig.mb,main="MB")
plotnet(ig.gl,main="glasso")
plotnet(ig.sp,main="SparCC")

#Quantify accuracy
## Max F1
max(stars.roc(getOptMerge(se.mb), graph,plot=FALSE,verbose=FALSE)$F1)
## [1] 0.9416334
max(stars.roc(getOptMerge(se.gl), graph,plot=FALSE,verbose=FALSE)$F1)
## [1] 0.9209065
max(stars.roc(abs(sp.est$Cor*sp.net), graph,plot=FALSE,verbose=FALSE)$F1)
## [1] 0.7559874
## Hamming distances
sum(se.mb$refit!=as.matrix(graph))/2
## [1] 181
sum(se.gl$refit!=as.matrix(graph))/2
## [1] 251
sum(sp.net!=as.matrix(graph))/2
## [1] 309
