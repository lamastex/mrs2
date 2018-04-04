##MDEHistogram.R
##Purpose:
#To Plot the MDE histogram using CompCount PQ for different critLeaves, n and d
##Date updated: 11 December 2015

#Output:
#1. MDE histogram, IAE and delta plots
#2. Table of min delta (theta), IAE (theta)
###########################################################
rm(list = ls(all = TRUE))
dev.off() #close all plots
cat("\014") #clear console

#set your directory here
mydir = "~/mrs/branches/gat41/examples/StatsSubPav/IAECalculations/" 
setwd(mydir)
source("subpavhist.R") #to use user-defined function "subpavhist" 
source("subpavhist2D.R") 
#read documentation later on how to organized user-defined functions

#change the method and directory to the folder of the density that you want 
#(e.g. FM1, FM2, UM1, etc.)
mydensity = "UM4D1"
mydensitydir = paste(mydir, mydensity, sep="")
setwd(mydensitydir)

#method= "FinMix"
method = "Unif"

#log scale or not?
boolLog = 1

##Start plotting...
#########################################
L = c(10, 100, 500)
N=c("100", "1000", "10000", "100000", "1000000")
#N=c("1000000")

minDeltaIAE = matrix(0, length(N), length(L))
minIAE = matrix(0, length(N), length(L))
thetaDelta = matrix(0, length(N), length(L))
thetaIAE = matrix(0, length(N), length(L))

for(i in 1:length(N)) {
    n=N[i]
   
    setEPS()
    postscript(paste(mydensity, "n", n, ".eps", sep = ""))
    par(mfrow = c(length(L),3))
    
    for(j in 1:length(L)){
      leaves = L[j]
      theta = 0:(leaves-1)
      
    #read in the files
    thisdir = paste("n", n, "L", leaves, sep = "")
    #print(thisdir)
    histpq = read.table(paste(thisdir, "/HistPQ1.txt", sep=""), header = FALSE)
    deltamax = read.table(paste(thisdir, "/", method, "Method1DeltaMax1.txt",sep=""), header = FALSE)
    IAE = read.table(paste(thisdir, "/", method, "Method1IAEandTrueDelta1.txt",sep=""), header = FALSE)
    
    #plots
    subpavhist(histpq,1); title(paste("n=", n, "L = ", leaves, sep=""))
    #subpavhist2D(histpq,1); title(paste("n=", n, "L = ", leaves, sep=""))
    plot(theta, log(deltamax$V1), col = "red", "b", lwd = 2, ylab = "DeltaMax (log)", 
          main = "Delta Max")
    plot(theta, log(IAE$V1), col = "blue", "b", lwd = 2, ylab = "IAE (log)", 
          main = "IAE")
    
    #store in vector
    m = which.min(deltamax$V1)
    thetaDelta[i,j] = which.min(deltamax$V1) -1
    minDeltaIAE[i,j] = IAE$V1[m]
    minIAE[i,j] = min(IAE$V1)
    thetaIAE[i,j] = which.min(IAE$V1) -1
    
    }
    dev.off()
    
}

#for L10
L10 = matrix(0, length(N),4)
L10[,1] = minDeltaIAE[,1]
L10[,2] =(thetaDelta[,1])
L10[,3] =minIAE[,1]
L10[,4] =thetaIAE[,1]
L10

#for L100
L100 = matrix(0, length(N),4)
L100[,1] = minDeltaIAE[,2]
L100[,2] =(thetaDelta[,2])
L100[,3] =minIAE[,2]
L100[,4] =thetaIAE[,2]
L100

#for L500
L100 = matrix(0, length(N),4)
L100[,1] = minDeltaIAE[,3]
L100[,2] =(thetaDelta[,3])
L100[,3] =minIAE[,3]
L100[,4] =thetaIAE[,3]
L100
