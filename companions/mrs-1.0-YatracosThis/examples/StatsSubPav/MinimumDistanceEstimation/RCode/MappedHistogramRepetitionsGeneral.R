##MappedMDERepetitions.R
##Purpose:
#To Plot the MDE histogram using CompCount PQ for different critLeaves, n and d (for repeated simulations)
#Date updated: 17 August 2016

###########################################################
rm(list = ls(all = TRUE))
dev.off() #close all plots
cat("\014") #clear console

#set your directory here
mydir = "~/mrs/branches/gat41/examples/StatsSubPav/IAECalculations/" 
#change the method and directory to the folder of the density that you want 
#(e.g. MappedGaussian1D, MappedGaussian2D.)
rep = 1
k = 1
kk = 0
#lambda = 10000
#L = 1000
#n = 10000
#mydensity = paste("RepMG1D/Lambda100/n10000","L200", sep = "")
#mydensitydir = paste(mydir, mydensity, sep="")
#mydensitydir
#setwd(mydensitydir)
setwd(mydir)
#########################################
#L = critLeaves, N = sample size
results = matrix(0, rep, 4)

for(i in 1:rep){
  i 
  par(mfrow=c(2,1))
  
  deltamax = read.table(paste("MappedDeltaMax", i, ".txt",sep=""), header = FALSE)
  plot(1:NROW(deltamax), log(deltamax$V1), col = "black", "b", lwd = 2, ylab = "DeltaMax (log)", 
       pch = 20, main = paste("Delta Max", i))
  
  iaes =  read.table(paste("MappedIAE", i, ".txt",sep=""), header = FALSE)
  plot(1:NROW(iaes), log(iaes$V1), col = "blue", "b", lwd = 2, ylab = "IAE (log)", 
       pch = 20, main = paste("IAE", i))

  #store in vector
  #IAE from minimum delta and corresponding number of leaves
  m = which.min(deltamax$V1)
  results[i,1] = iaes$V1[(m-1)*k + kk]
  results[i,2] = (m-1) * k + kk
  #minimum IAE and corresponding number of leaves
  results[i,3] = min(iaes$V1) 
  results[i,4] = which.min(iaes$V1) 
}

results  
ave_results = apply(results, 2, mean)
std_results = apply(results, 2, sd)

ave_results
std_results

#diff_leaves = results[, 2] - results[, 4]
#par(mfrow = c(1,1))
#plot(results[,2], results[,4], xlab = "Number of leaves (min delta)", ylab = "Number of leaves (min IAE)")
#abline(0,1)
