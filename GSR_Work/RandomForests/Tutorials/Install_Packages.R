#!/usr/bin/Rscript
#Author: "Emma Fox (eafox@ucla.edu)"

#tree
install.packages("tree")

#randomForest
install.packages("randomForest")

#gradientForest
# may need R 3.3.3 or earlier
install.packages("gradientForest", repos="http://R-Forge.R-project.org")
library(gradientForest)
# if having problems, go to terminal and type:
# xcode-select --install
# brew install gcc
# and comment out FLIBS in /Library/Frameworks/R.framework/Resources/etc/Makeconf

#LFMM from LEA
install.packages("devtools")
devtools::install_github("bcm-uga/LEA")
#get from commandline instead

#GDM
#real old version