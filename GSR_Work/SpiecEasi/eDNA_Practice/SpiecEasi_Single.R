#!/usr/bin/Rscript
#Author: "Emma Fox (eafox@ucla.edu)" and "Rachel Meyer (rsmeyer@ucla.edu)"
#Assistance provided by "Zachary Kurtz (zdkurtz@gmail.com)"

#########
# NOTES #
#########
# Convert files from .txt to .csv before using.
# (Unix) Commandline way:
# cat foo.txt | tr -s "\t" "," >> foo.csv 

#####################
# REQUIRED PACKAGES #
#####################
library(SpiecEasi)
library(igraph)
library(Matrix)

##########################################
# READ IN FILES AND CORRECT SAMPLE NAMES #
##########################################
minPres<-0.15
n_lambda<-50
rep_num<-50
num_cores<-2

#read in file
readIn<-read.csv("16s.csv", stringsAsFactors = FALSE)
#remove row where otu name is NA
remNA<-subset(readIn, !is.na(readIn$sum.taxonomy))
#change row 1 to the column name
Orig <- remNA[,-1]
rownames(Orig) <- remNA[,1]
#remove asvs not present in x percent of samples
Filt<-Orig[rowSums(Orig == 0) <= (1-minPres)*length(Orig), ]
#transpose
tr.Filt<-t(Filt)
#fix ANACAPA sample names to match across primers and remove blanks
oldNames<-row.names(tr.Filt)
newNames<-sub(".*_", "",oldNames)
row.names(tr.Filt)<-newNames
tr.Filt<-tr.Filt[!grepl("*Blank*",rownames(tr.Filt)),]

####################
# RUNNING THE DATA #
####################
se.mb.result<-spiec.easi(tr.Filt, method='mb', lambda.min.ratio=1e-2, 
                         nlambda=n_lambda, icov.select.params=list(rep.num=rep_num,
                                                                   ncores=num_cores))
se.glasso.result<-spiec.easi(tr.Filt, method='glasso', lambda.min.ratio=1e-2, 
                             nlambda=n_lambda, icov.select.params=list(rep.num=rep_num,
                                                                       ncores=num_cores))

########
# PLOT #
########
#import networks into iGraph
ig.mb.PS <- adj2igraph(se.mb.result$refit)
ig.gl.PS <- adj2igraph(se.glasso.result$refit)  
#set size of vertex proportional to clr-mean
vsize <- rowMeans(clr(tr.Filt, 1))+6
am.coord <- layout.fruchterman.reingold(ig.mb.PS)
plot(ig.mb.PS, layout=am.coord, vertex.size=vsize, vertex.label=NA, main="MB")
plot(ig.gl.PS, layout=am.coord, vertex.size=vsize, vertex.label=NA, main="glasso")

######################
# SAVE OUTPUT TABLES #
######################
#glasso
gl_result_table<-formatSparseM(se.glasso.result$refit, zero.print = ".", align = c("fancy", "right"),asLogical=NULL, uniDiag=NULL, digits=NULL)
asv_names<-colnames(tr.Filt)
colnames(gl_result_table)<-asv_names
rownames(gl_result_table)<-asv_names
write.csv(gl_result_table, file="gl_single.csv")

#mb
mb_result_table<-formatSparseM(se.mb.result$refit, zero.print = ".", align = c("fancy", "right"),asLogical=NULL, uniDiag=NULL, digits=NULL)
colnames(mb_result_table)<-asv_names
rownames(mb_result_table)<-asv_names
write.csv(mb_result_table, file="mb_single.csv")