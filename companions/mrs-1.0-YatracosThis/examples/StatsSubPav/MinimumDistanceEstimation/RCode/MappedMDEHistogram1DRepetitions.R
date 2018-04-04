
##MappedMDEHistogram1DRepetition.R
##Purpose:
#To Plot the MDE histogram using CompCount PQ for different critLeaves, n and d (for repeated simulations)
#Date updated: 5 August 2016

###########################################################
#rm(list = ls(all = TRUE))
dev.off() #close all plots
cat("\014") #clear console

#to get the needed functions
mydir = "~/mrs/branches/gat41/examples/StatsSubPav/IAECalculations" 
setwd(mydir)
source("subpavhist2D.R") #to use user-defined function "subpavhist" 
source("plotMap2D.R") #to use user-defined function "subpavhist" 
#read documentation later on how to organized user-defined functions

#set your directory here
mydir = "~/mrs/branches/gat41/examples/StatsSubPav/IAECalculations/RepMG2D/" 
setwd(mydir)

#change the method and directory to the folder of the density that you want 
#(e.g. MappedGaussian1D, MappedGaussian2D.)
rep = 1
lambda = 10000
L = 1000
n = 10000
mydensity = paste("RepMappedGaussian1DL10000/n10000", "L", L, sep = "")
mydensitydir = paste(mydir, mydensity, sep="")
mydensitydir
setwd(mydensitydir)

#read in IAE
#IAEdensitydir = "~/mrs/branches/gat41/examples/StatsSubPav/IAECalculations/"
#density = "MG1DIAELambda100/"
#setwd(paste(IAEdensitydir,density,sep = ""))

#myIAE = matrix(0, L, 1)

#for(i in 1:L){
#  val = read.table(paste("n", n, "L", L, "/","MappedIAE", i, ".txt",sep=""), header = FALSE)
#  myIAE[i] = val$V1
#}

#########################################
#optional
#import the PCF histogram
#par(mfrow = c(1,1))

pcf = read.table("PCF.txt")
plotMap1D(pcf, 1)
#x=seq(-3,3,by=10^-3)
#lines(x, dnorm(x), col="blue")

#########################################
#L = critLeaves, N = sample size
results = matrix(0, rep, 2)
numleaves = 1:L

for(i in 1:rep){
    i 
    #histpq = read.table(paste(thisdir, "/HistPQ1.txt", sep=""), header = FALSE)
    deltamax = read.table(paste(mydensitydir, "/", "MappedDeltaMax", i, ".txt",sep=""), header = FALSE)
    #IAEdensitydir = "~/mrs/branches/gat41/examples/StatsSubPav/IAECalculations/MappedIAE/MappedGaussian1D/"
    #IAE = read.table(paste(IAEdensitydir,"n",n, "L", L, "/","MappedIAE", i, ".txt",sep=""), header = FALSE)
    #IAE = myIAE
    
    #plots
    #subpavhist(histpq,1); title(paste("n=", n, ", L = ", leaves, sep=""))
    #lines(x, dnorm(x), col="blue")
    #par(mfrow=c(2,1))
    plot(numleaves, log(deltamax$V1), col = "black", "b", lwd = 2, ylab = "DeltaMax (log)", 
         pch = 20, main = "Delta Max")
    #plot(numleaves, log(IAE), col = "blue", "b", lwd = 2, ylab = "IAE (log)", 
    #     pch = 20, main = "IAE")
    
    #store in vector
    #IAE from minimum delta and corresponding number of leaves
    m = which.min(deltamax$V1)
    #results[i,1] = IAE[m]
    results[i,2] = m 
    #minimum IAE and corresponding number of leaves
    #results[i,3] = min(IAE) 
    #results[i,4] = which.min(IAE) 
  }


results  
#ave_results = apply(results, 2, mean)
#median_results = apply(results, 2, median)
#diff_leaves = results[, 2] - results[, 4]
#par(mfrow = c(1,1))
#plot(results[,2], results[,4], xlab = "Number of leaves (min delta)", ylab = "Number of leaves (min IAE)")
#abline(0,1)
