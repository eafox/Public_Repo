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
library(optparse)
library(tools)
library(network)
library(igraph)

#######################
# COMMANDLINE OPTIONS #
#######################
option_list<-list(
  make_option(c("--input_file"), default=NULL,
              help = "Name of the file(s) containing OTUs. ANACAPA output"),
  make_option(c("--name_correct"), default=FALSE,
              help="If files are ANACAPA output, specify TRUE to standardize sample names across files"),
  make_option(c("--multi_se"), default = FALSE, help = "Specify '--multi TRUE' to use multi-spiec easi"),
  make_option(c("--min_presence"), default=0.15,
              help="Proportion of samples OTU must be present in to be considered. 0<x<1. Default=0.25"),
  make_option(c("--n_lambda"), default=50,
              help=""),
  make_option(c("--stars_thresh"), default=0.05,
              help=""),
  make_option(c("--rep_num"), default=50,
              help=""),
  make_option(c("--plot_output"), type="character",
              help="name of file to output the network plots to")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

####################
# SCRIPT FUNCTIONS #
####################
pruneFiles<-function(inFileName,minPres){
  #read in file
  readIn<-read.csv(inFileName, stringsAsFactors = FALSE)
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
  if (opt$name_correct==TRUE){
    oldNames<-row.names(tr.Filt)
    newNames<-sub(".*_", "",oldNames)
    row.names(tr.Filt)<-newNames
    tr.Filt<-tr.Filt[!grepl("*Blank*",rownames(tr.Filt)),]
  } 
  #return filtered and formatted file
  return(tr.Filt)
}


SE_mb<-function(pruneFileName,multi){
  #Runs mb method of network fitting
  if (multi==TRUE){ 
    se.mb.result<-multi.spiec.easi(pruneFileName, method='mb', sel.criterion='stars', 
                                   lambda.min.ratio=1e-3, nlambda=opt$n_lambda, 
                                   icov.select.params=list(rep.num=opt$rep_num, 
                                                           stars.thresh=opt$stars_thresh))
  } else{
    se.mb.result<-spiec.easi(pruneFileName, method='mb', lambda.min.ratio=1e-2, 
                             nlambda=opt$n_lambda, icov.select.params=list(rep.num=opt$rep_num))
  }
  return(se.mb.result)
}

SE_glasso<-function(pruneFileName,multi){
  #Runs glasso method of network fitting
  if (multi==TRUE){ 
    se.glasso.result<-multi.spiec.easi(pruneFileName, method='glasso', sel.criterion='stars', 
                                       lambda.min.ratio=1e-3, nlambda=opt$n_lambda, 
                                       icov.select.params=list(rep.num=opt$rep_num, 
                                                               stars.thresh=opt$stars_thresh))
  } else{
    se.glasso.result<-spiec.easi(pruneFileName, method='glasso', lambda.min.ratio=1e-2, 
                             nlambda=opt$n_lambda, icov.select.params=list(rep.num=opt$rep_num))
  }
  return(se.glasso.result)
}

####################
# RUNNING THE DATA #
####################
#single file
if (opt$multi==FALSE){
  datFile<-pruneFiles(opt$input_file,opt$min_presence)
  mb.output<-SE_mb(datFile,opt$multi)
  glasso.output<-SE_glasso(datFile,opt$multi)
  #import networks into iGraph
  ig.mb.PS <- adj2igraph(mb.output$refit)
  ig.gl.PS <- adj2igraph(glasso.output$refit)  
  #set size of vertex proportional to clr-mean
  vsize <- rowMeans(clr(datFile, 1))+6
  am.coord <- layout.fruchterman.reingold(ig.mb.PS)
  #Save plots to pdf
  if (is.null(opt$plot_output)){
    BaseName<-basename(file_path_sans_ext(file_path_sans_ext(opt$input_file)))
    PlotName<-paste(BaseName,".pdf", sep = "")
  } else {
    PlotName<-paste(opt$plot_output,".pdf", sep = "")  
  }
  pdf(PlotName)
  par(mfrow=c(1,2))
  plot(ig.mb.PS, layout=am.coord, vertex.size=vsize, vertex.label=NA, main="MB")
  plot(ig.gl.PS, layout=am.coord, vertex.size=vsize, vertex.label=NA, main="glasso")
  dev.off()
}

#multiple files
if (opt$multi==TRUE){
  inFiles<-read.table(opt$input_file)
  inFiles$V1<-as.character(inFiles$V1)
  datFrames<-list()
  for (i in 1:nrow(inFiles)){
    dat<-pruneFiles(inFiles[i,1],opt$min_presence)
    datFrames[[i]]<-dat
  }
  mb.output<-SE_mb(datFrames,opt$multi)
  glasso.output<-SE_glasso(datFrames,opt$multi)
  #plot
  taxtype <- rep(1:length(datFrames), times=sapply(datFrames, ncol))
  if (is.null(opt$plot_output)){
    BaseName<-basename(file_path_sans_ext(file_path_sans_ext(opt$input_file)))
    PlotName<-paste(BaseName,".pdf", sep = "")
  } else {
    PlotName<-paste(opt$plot_output,".pdf", sep = "")  
  }
  pdf(PlotName)
  par(mfrow=c(1,2))
  plot(network::network(mb.output$refit), usearrows=FALSE, vertex.col=taxtype, main="MB")
  plot(network::network(glasso.output$refit), usearrows=FALSE, vertex.col=taxtype, main="Glasso")
  dev.off()
}
