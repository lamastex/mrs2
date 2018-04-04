##ObserveEachSequence.R
##Purpose:
#To plot the IAEs for histograms to determine the minimum IAE.
#Date updated: 1 September 2016

###########################################################
rm(list = ls(all = TRUE))
dev.off() #close all plots
cat("\014") #clear console

mydir = "~/mrs/branches/gat41/examples/StatsSubPav/IAECalculations/output/" 
mydensity = "Rosen5D/Lambda100000/"
mydensitydir = paste(mydir, mydensity, sep="")
mydensitydir
setwd(mydensitydir)
tot = 1

#read in tables
IAE = read.table(paste("MappedIAE", tot, ".txt", sep=""), header = FALSE)
attach(IAE)
par(mfrow=c(1,1))
plot(1:length(V1), log(V1), col = "blue", "b", lwd = 2, xlab = "histogram", ylab = "IAE (log)", 
     pch = 20, main = "IAE")
m = which.min(V1)
points(m, log(V1[m]), col = "red", pch = 20)
text(m, log(V1[m]) + 0.2, m, pos = 4)

par(mfrow=c(2,1))
iaes = matrix(0, tot+1, 1)
leaves = matrix(0, tot+1, 1)
leaves[1] = m
iaes[1] = IAE$V1[m]


for ( i in 0:tot) {
  
  deltamax0 = read.table(paste("MappedDeltaMax", i, ".txt", sep=""), header = FALSE)
  sequence0 = read.table(paste("Sequence", i, ".txt", sep = ""), header=FALSE)
  plot(sequence0$V1, log(deltamax0$V1), col = "blue", "b", lwd = 2, xlab = "histogram", ylab = "deltamax (log)", 
     pch = 20, main = paste("Iteration", i+1))
  m0 = which.min(deltamax0$V1)
  points(sequence0$V1[m0], log(deltamax0$V1[m0]), col = "red", pch = 20)
  text(sequence0$V1[m0], log(deltamax0$V1[m0]) +0.5, sequence0$V1[m0])

  leaves[i+2] = sequence0$V1[m0]
  iaes[i+2] = IAE$V1[sequence0$V1[m0]]
}


iaes
leaves

detach(IAE)
