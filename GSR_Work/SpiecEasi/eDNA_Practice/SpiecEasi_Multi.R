#!/usr/bin/Rscript
#Author: "Emma Fox (eafox@ucla.edu)" and "Rachel Meyer (rsmeyer@ucla.edu)"
#Assistance provided by "Zachary Kurtz (zdkurtz@gmail.com)"

#########
# NOTES #
#########
# Convert files from .txt to .csv before using.
# (Unix) Commandline way:
# cat foo.txt | tr -s "\t" "," >> foo.csv 

# To make list of files for multi-spiec easi:
# ls *.txt > inList

#####################
# REQUIRED PACKAGES #
#####################
library(SpiecEasi)
library(network)

##########################################
# READ IN FILES AND CORRECT SAMPLE NAMES #
##########################################
fileList<-c("16s.csv",  "CO1.csv",  "PITS.csv")
minPres<-0.15
n_lambda<-50
rep_num<-50
stars_thresh<-0.05
num_cores<-2
  
datFrames<-list()
for (i in 1:length(fileList)){
  #read in file
  readIn<-read.csv(fileList[i], stringsAsFactors = FALSE)
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
  oldNames<-row.names(tr.Filt) # comment out these lines to keep original sample names
  newNames<-sub(".*_", "",oldNames) # comment out these lines to keep original sample names
  row.names(tr.Filt)<-newNames # comment out these lines to keep original sample names
  tr.Filt<-tr.Filt[!grepl("*Blank*",rownames(tr.Filt)),] # comment out these lines to keep original sample names
  #return filtered and formatted file
  datFrames[[i]]<-tr.Filt
}

####################
# RUNNING THE DATA #
####################
se.mb.result<-multi.spiec.easi(datFrames, method='mb', sel.criterion='stars', 
                               lambda.min.ratio=1e-3, nlambda=n_lambda, 
                               icov.select.params=list(rep.num=rep_num, 
                                                       stars.thresh=stars_thresh,
                                                       ncores=num_cores))
se.glasso.result<-multi.spiec.easi(datFrames, method='glasso', sel.criterion='stars', 
                                   lambda.min.ratio=1e-3, nlambda=n_lambda, 
                                   icov.select.params=list(rep.num=rep_num, 
                                                           stars.thresh=stars_thresh,
                                                           ncores=num_cores))

########
# PLOT #
########
taxtype <- rep(1:3, times=sapply(datFrames, ncol))
plot(network::network(se.mb.result$refit), usearrows=FALSE, vertex.col=taxtype, main="MB")
plot(network::network(se.glasso.result$refit), usearrows=FALSE, vertex.col=taxtype, main="Glasso")
