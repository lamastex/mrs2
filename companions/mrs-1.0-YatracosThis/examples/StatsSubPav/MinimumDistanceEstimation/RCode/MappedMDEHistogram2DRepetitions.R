##MappedMDEHistogram2DRepetitions.R
##Purpose:
#To Plot the MDE histogram using CompCount PQ for different critLeaves, n and d (for repeated simulations)
#Date updated: 6 August 2016

###########################################################
rm(list = ls(all = TRUE))
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
rep = 5
lambda = 100
L = 1000
n = 10000
mydensity = paste("RepMappedGaussian2D/n", n, "L", L, sep = "")
mydensitydir = paste(mydir, mydensity, sep="")
mydensitydir
setwd(mydensitydir)

#########################################
#optional
#import the PCF histogram
par(mfrow = c(1,1))
pcf = read.table("PCF.txt")
#plotMap2D(pcf, 1) #double check the 1 and 2
#subpavhist2D(pcf,2) #double check the 1 and 2

#########################################
#L = critLeaves, N = sample size
results = matrix(0, rep, 4)
numleaves = 1:200

for(i in 1:rep){
  i 
  #histpq = read.table(paste(thisdir, "/HistPQ1.txt", sep=""), header = FALSE)
  deltamax = read.table(paste(mydensitydir, "/", "MappedDeltaMax", i, ".txt",sep=""), header = FALSE)
  IAEdensitydir = "~/mrs/branches/gat41/examples/StatsSubPav/IAECalculations/MappedIAE/MappedGaussian2D/"
  IAE = read.table(paste(IAEdensitydir,"n",n, "L", L, "/","MappedIAE", i, ".txt",sep=""), header = FALSE)
  
  #plots
  #subpavhist(histpq,1); title(paste("n=", n, ", L = ", leaves, sep=""))
  #lines(x, dnorm(x), col="blue")
  par(mfrow=c(2,1))
  plot(numleaves, log(deltamax$V1), col = "black", "b", lwd = 2, ylab = "DeltaMax (log)", 
       pch = 20, main = paste("Delta Max", i))
  plot(numleaves, log(IAE$V1), col = "blue", "b", lwd = 2, ylab = "IAE (log)", 
       pch = 20, main = paste("IAE", i))
  
  #store in vector
  #IAE from minimum delta and corresponding number of leaves
  m = which.min(deltamax$V1)
  results[i,1] = IAE$V1[m]
  results[i,2] = m 
  #minimum IAE and corresponding number of leaves
  results[i,3] = min(IAE$V1) 
  results[i,4] = which.min(IAE$V1) 
}

results  
ave_results = apply(results, 2, mean)
median_results = apply(results, 2, median)
diff_leaves = results[, 2] - results[, 4]
par(mfrow = c(1,1))
plot(results[,2], results[,4], xlab = "Number of leaves (min delta)", ylab = "Number of leaves (min IAE)")
abline(0,1)
