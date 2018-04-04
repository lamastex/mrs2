##UniformMixtures.R
##Purpose:
##To plot uniform mixtures based on given leaf levels 
##Date updated: 15 December 2015
###########################################################
rm(list = ls(all = TRUE))
dev.off() #close all plots
cat("\014") #clear console

#set your parameters here
UM=1

#set your directory here
mydir = "~/mrs/branches/gat41/examples/StatsSubPav/IAECalculations/UM" 
setwd(mydir)
source("plotLeafLevel1D.R") 
#source("plotLeafLevel2D.R") 
################################################
## actual leaf level
if (shape==1) {
  leafLevel = 0;
  filled = 1;
} else { #if (shape==2)
  leafLevel = c(3,4,4,2,2,3,3);
  filled = rep(1,7);
}


e