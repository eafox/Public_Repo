#!/usr/bin/Rscript
#Author: "Emma Fox (eafox@ucla.edu)"

library(reshape2)

#Right now, we have 3 data sets (shrub scrub, forest, coastal) with 5 primers each

##########################
# READ IN ALL DATA FILES #
##########################
#List all data files
inList<-c("~/Downloads/ShrubScrub/ShrubScrub_16S_ASV_sum_by_taxonomy_60.csv",
          "~/Downloads/ShrubScrub/ShrubScrub_18S_ASV_sum_by_taxonomy_60.csv",
          "~/Downloads/ShrubScrub/ShrubScrub_CO1_ASV_sum_by_taxonomy_60.csv",
          "~/Downloads/ShrubScrub/ShrubScrub_FITS_ASV_sum_by_taxonomy_60.csv",
          "~/Downloads/ShrubScrub/ShrubScrub_PITS_ASV_sum_by_taxonomy_60.csv",
          "~/Downloads/Forest/Forest_16S_ASV_sum_by_taxonomy_60.csv",
          "~/Downloads/Forest/Forest_18S_ASV_sum_by_taxonomy_60.csv",
          "~/Downloads/Forest/Forest_CO1_ASV_sum_by_taxonomy_60.csv",
          "~/Downloads/Forest/Forest_FITS_ASV_sum_by_taxonomy_60.csv",
          "~/Downloads/Forest/Forest_PITS_ASV_sum_by_taxonomy_60.csv",
          "~/Downloads/Coastal/Coastal_12S_ASV_sum_by_taxonomy_60.csv",
          "~/Downloads/Coastal/Coastal_16S_ASV_sum_by_taxonomy_60.csv",
          "~/Downloads/Coastal/Coastal_18S_ASV_sum_by_taxonomy_60.csv",
          "~/Downloads/Coastal/Coastal_CO1_ASV_sum_by_taxonomy_60.csv",
          "~/Downloads/Coastal/Coastal_FITS_ASV_sum_by_taxonomy_60.csv",
          "~/Downloads/Coastal/Coastal_PITS_ASV_sum_by_taxonomy_60.csv")
datList<-list()

#Read in each of the files to a list of data frames and adjust sample names to match across sets
for (i in 1:length(inList)){
  readIn<-read.csv(inList[i],header = TRUE)
  remNA<-subset(readIn, !is.na(readIn$sum.taxonomy))
  #change row 1 to the column name
  Orig <- remNA[,-1]
  rownames(Orig) <- remNA[,1]
  #remove asvs not present in 15% of samples
  Filt<-Orig[rowSums(Orig == 0) <= (1-0.15)*length(Orig), ]
  #transpose so each row is a sample
  tr.Filt<-t(Filt)
  #Match sample names
  oldNames<-row.names(tr.Filt)
  newNames<-sub(".*_", "",oldNames) #remove primer tag
  newNames2<-regmatches(newNames,regexpr("\\w*.{3}",newNames))
  newNames3<-sub("[[:punct:]]","",newNames2)
  newNames4<-regmatches(newNames3,regexpr("\\w\\d*\\w\\d",newNames3))
  row.names(tr.Filt)<-regmatches(newNames4,regexpr("\\w*",newNames4))#keep part of sample name that matches meta name
  #remove any blanks or negatives still present
  remDF<-tr.Filt[!grepl("*eg*",rownames(tr.Filt)),]
  anaDF<-remDF[!grepl("*nk*",rownames(remDF)),]
  anaDF2<-anaDF[!grepl("*gS*",rownames(anaDF)),]
  finalDF<-melt(anaDF2)
  colnames(finalDF)<-c("MatchName","ASV","Abund")
  if (i<6){
    finalDF$Set<-rep("ShrubScrub",length(finalDF$MatchName))
  } else if (5<i & i<11){
    finalDF$Set<-rep("Forest",length(finalDF$MatchName))
  } else if (10<i & i<17){
    finalDF$Set<-rep("Coastal",length(finalDF$MatchName))
  }
  
  datList[[i]]<-finalDF
}

#Put all of the data frames in one data frame
combDatList<-datList[[1]]
for (i in 2:length(datList)){
  combDatList<-rbind(combDatList,datList[[i]])
}
write.csv(combDatList,"01_readIn_data.csv",row.names = TRUE)

#get unique coordinates from the data read in
nameDF<-subset(combDatList,select=c("MatchName","Set"))
nameDF_filt<-subset(nameDF, !duplicated(MatchName))
#278 unique sites from our available data

########################################
# GET COORDINATES FOR THE SAMPLE SITES #
########################################
coords<-read.csv("FCSs_Coords.csv")
#Make data frame with sample name and coordinates
coordsDF<-subset(coords,select = c("MatchName","Set","Latitude","Longitude"))
combDF<-merge(coordsDF,nameDF_filt,by=c("MatchName","Set"))

write.csv(combDF,file = "02_datCoords.csv",row.names = TRUE)
#seem to have lost a site. I am aware and working on it#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#
#########################################
# FILTER DATA FOR COORDINATES AVAILABLE #
#########################################
newDatList<-list()
for (i in 1:length(datList)){
  filtDF<-datList[[i]][datList[[i]]$MatchName %in% combDF$MatchName,]
  newDatList[[i]]<-subset(filtDF,select = c("MatchName","ASV","Abund"))
}

#Merge to DF that has all data points for each species by sample
abunDF<-newDatList[[1]]
for (i in 2:length(newDatList)){
  abunDF<-rbind(abunDF,newDatList[[i]])
}
abunDF$Abund<-as.numeric(abunDF$Abund)
write.csv(abunDF,"03_fullMelt_abund.csv",row.names = TRUE)

###########################################
# CAST THE DATA FRAME TO THE RIGHT FORMAT #
###########################################
#Put all of the data into the form of match name as rows and ASV as columns with abundance as the full values
castDF<-dcast(abunDF,MatchName~ASV,value.var = "Abund",sum)

#Last column is all NA
completeDF<-castDF[,2:2543]
row.names(completeDF)<-castDF[,1]
#Change NAs to 0
completeDF_narm<-completeDF
completeDF_narm[is.na(completeDF_narm)]<-0

write.csv(completeDF,"04_fullDF.csv",row.names = TRUE)
write.csv(completeDF_narm,"04b_fullDF_naRM.csv",row.names = TRUE)

###########################################
# FINAL SET OF COORDS TO GET METADATA FOR #
###########################################
#This step is just checking we have the same set from earlier
final_coords<-data.frame(row.names(completeDF))
colnames(final_coords)<-c("MatchName")
coords_match<-merge(final_coords, combDF,by="MatchName")
#If they do match, this should equal 277
length(coords_match$MatchName)
#THEY DO







