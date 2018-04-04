##MappedMDEHistogram1D.R
##Purpose:
#To Plot the MDE histogram using CompCount PQ for different critLeaves, n and d
##Date updated: 11 December 2015

###########################################################
rm(list = ls(all = TRUE))
dev.off() #close all plots
cat("\014") #clear console

#set your directory here
mydir = "~/mrs/branches/gat41/examples/StatsSubPav/IAECalculations/" 
setwd(mydir)
source("subpavhist.R") #to use user-defined function "subpavhist" 
source("plotMap1D.R") #to use user-defined function "subpavhist" 
#read documentation later on how to organized user-defined functions

#change the method and directory to the folder of the density that you want 
#(e.g. MappedGaussian1D, MappedGaussian2D.)
mydensity = "MappedGaussian1D"
mydensitydir = paste(mydir, mydensity, sep="")
setwd(mydensitydir)

#########################################
#import the PCF histogram
par(mfrow = c(1,1))
pcf = read.table("PCF.txt")
plotMap1D(pcf, 1)
x=seq(-3,3,by=10^-3)
lines(x, dnorm(x), col="blue")

#########################################

L = c(10,100)
N=c("100", "1000", "10000", "100000")

minDeltaIAE = matrix(0, length(N), length(L))
minIAE = matrix(0, length(N), length(L))
thetaDelta = matrix(0, length(N), length(L))
thetaIAE = matrix(0, length(N), length(L))

for(i in 1:4){
  n=N[i]
  
  setEPS()
  postscript(paste(mydensity, "n", n, ".eps", sep = ""))
  par(mfrow = c(2,3))
  
  for(j in 1:2){
    leaves = L[j]
    theta = 0:(leaves-1)
    
    #read in the files
    thisdir = paste("n", n, "L", leaves, sep = "")
    #print(thisdir)
    histpq = read.table(paste(thisdir, "/HistPQ1.txt", sep=""), header = FALSE)
    deltamax = read.table(paste(thisdir, "/", "MappedDeltaMax1.txt",sep=""), header = FALSE)
    IAE = read.table(paste(thisdir, "/", "MappedIAEandTrueDelta1.txt",sep=""), header = FALSE)
    
    #plots
    subpavhist(histpq,1); title(paste("n=", n, "L = ", leaves, sep=""))
    lines(x, dnorm(x), col="blue")
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
L10 = matrix(0, 4,4)
L10[,1] = minDeltaIAE[,1]
L10[,2] =(thetaDelta[,1])
L10[,3] =minIAE[,1]
L10[,4] =thetaIAE[,1]
L10

#for L100
L100 = matrix(0, 4,4)
L100[,1] = minDeltaIAE[,2]
L100[,2] =(thetaDelta[,2])
L100[,3] =minIAE[,2]
L100[,4] =thetaIAE[,2]
L100


# #n=100, L = 10
# par(mfrow = c(2,3))
# L = 10
# theta = 0:(L-1)
# histpq = read.table("n100L10/HistPQ1.txt", header = FALSE)
# deltamax = read.table(paste("n100L10/MappedDeltaMax1.txt",sep=""), header = FALSE)
# IAE = read.table(paste("n100L10/MappedIAEandTrueDelta1.txt",sep=""), header = FALSE)
# subpavhist(histpq,1);title("n=100, L = 10")
# lines(x, dnorm(x), col="blue")
# 
# #delta and IAE
# plot(theta, deltamax$V1, col = "red", "b", lwd = 2, ylab = "DeltaMax", 
#      main = "Delta Max (n=100, L = 10)")
# plot(theta, IAE$V1, col = "blue", "b", lwd = 2, ylab = "IAE", 
#      main = "IAE (n=100, L = 10)")
# 
# ##########################
# #n=100, L = 100
# L = 100
# theta = 0:(L-1)
# histpq = read.table("n100L100/HistPQ1.txt", header = FALSE)
# deltamax = read.table(paste("n100L100/MappedDeltaMax1.txt",sep=""), header = FALSE)
# IAE = read.table(paste("n100L100/MappedIAEandTrueDelta1.txt",sep=""), header = FALSE)
# subpavhist(histpq,1);title("n=100, L = 100")
# lines(x, dnorm(x), col="blue")
# 
# #delta and IAE
# plot(theta, deltamax$V1, col = "red", "b", lwd = 2, ylab = "DeltaMax", 
#      main = "Delta Max (n=100, L = 100)")
# plot(theta, IAE$V1, col = "blue", "b", lwd = 2, ylab = "IAE", 
#      main = "IAE (n=100, L = 100)")
# 
# ######################################################
# par(mfrow = c(2,3))
# #n=1000, L = 10
# L = 10
# theta = 0:(L-1)
# histpq = read.table("n1000L10/HistPQ1.txt", header = FALSE)
# deltamax = read.table(paste("n1000L10/MappedDeltaMax1.txt",sep=""), header = FALSE)
# IAE = read.table(paste("n1000L10/MappedIAEandTrueDelta1.txt",sep=""), header = FALSE)
# subpavhist(histpq,1);title("n=1000, L = 10")
# lines(x, dnorm(x), col="blue")
# 
# #delta and IAE
# plot(theta, deltamax$V1, col = "red", "b", lwd = 2, ylab = "DeltaMax", 
#      main = "Delta Max (n=1000, L = 10)")
# plot(theta, IAE$V1, col = "blue", "b", lwd = 2, ylab = "IAE", 
#      main = "IAE (n=1000, L = 10)")
# 
# ##########################
# #n=1000, L = 100
# L = 100
# theta = 0:(L-1)
# histpq = read.table("n1000L100/HistPQ1.txt", header = FALSE)
# deltamax = read.table(paste("n1000L100/MappedDeltaMax1.txt",sep=""), header = FALSE)
# IAE = read.table(paste("n1000L100/MappedIAEandTrueDelta1.txt",sep=""), header = FALSE)
# subpavhist(histpq,1);title("n=1000, L = 100")
# lines(x, dnorm(x), col="blue")
# 
# #delta and IAE
# plot(theta, deltamax$V1, col = "red", "b", lwd = 2, ylab = "DeltaMax", 
#      main = "Delta Max (n=1000, L = 100)")
# plot(theta, IAE$V1, col = "blue", "b", lwd = 2, ylab = "IAE", 
#      main = "IAE (n=1000, L = 100)")
# 
# ######################################################
# 
# ######################################################
# par(mfrow = c(2,3))
# #n=10000, L = 10
# L = 10
# theta = 0:(L-1)
# histpq = read.table("n10000L10/HistPQ1.txt", header = FALSE)
# deltamax = read.table(paste("n10000L10/MappedDeltaMax1.txt",sep=""), header = FALSE)
# IAE = read.table(paste("n10000L10/MappedIAEandTrueDelta1.txt",sep=""), header = FALSE)
# subpavhist(histpq,1);title("n=10000, L = 10")
# lines(x, dnorm(x), col="blue")
# 
# #delta and IAE
# plot(theta, deltamax$V1, col = "red", "b", lwd = 2, ylab = "DeltaMax", 
#      main = "Delta Max (n=10000, L = 10)")
# plot(theta, IAE$V1, col = "blue", "b", lwd = 2, ylab = "IAE", 
#      main = "IAE (n=10000, L = 10)")
# 
# ##########################
# #n=10000, L = 100
# L = 100
# theta = 0:(L-1)
# histpq = read.table("n10000L100/HistPQ1.txt", header = FALSE)
# deltamax = read.table(paste("n10000L100/MappedDeltaMax1.txt",sep=""), header = FALSE)
# IAE = read.table(paste("n10000L100/MappedIAEandTrueDelta1.txt",sep=""), header = FALSE)
# subpavhist(histpq,1);title("n=10000, L = 100")
# lines(x, dnorm(x), col="blue")
# 
# #delta and IAE
# plot(theta, deltamax$V1, col = "red", "b", lwd = 2, ylab = "DeltaMax", 
#      main = "Delta Max (n=10000, L = 100)")
# plot(theta, IAE$V1, col = "blue", "b", lwd = 2, ylab = "IAE", 
#      main = "IAE (n=10000, L = 100)")
# 
# ######################################################
# 
# ######################################################
# par(mfrow = c(2,3))
# #n=100000, L = 10
# L = 10
# theta = 0:(L-1)
# histpq = read.table("n100000L10/HistPQ1.txt", header = FALSE)
# deltamax = read.table(paste("n100000L10/MappedDeltaMax1.txt",sep=""), header = FALSE)
# IAE = read.table(paste("n100000L10/MappedIAEandTrueDelta1.txt",sep=""), header = FALSE)
# subpavhist(histpq,1);title("n=100000, L = 10")
# lines(x, dnorm(x), col="blue")
# 
# #delta and IAE
# plot(theta, deltamax$V1, col = "red", "b", lwd = 2, ylab = "DeltaMax", 
#      main = "Delta Max (n=100000, L = 10)")
# plot(theta, IAE$V1, col = "blue", "b", lwd = 2, ylab = "IAE", 
#      main = "IAE (n=100000, L = 10)")
# 
# ##########################
# #n=100000, L = 100
# L = 100
# theta = 0:(L-1)
# histpq = read.table("n100000L100/HistPQ1.txt", header = FALSE)
# deltamax = read.table(paste("n100000L100/MappedDeltaMax1.txt",sep=""), header = FALSE)
# IAE = read.table(paste("n100000L100/MappedIAEandTrueDelta1.txt",sep=""), header = FALSE)
# subpavhist(histpq,1);title("n=100000, L = 100")
# lines(x, dnorm(x), col="blue")
# 
# #delta and IAE
# plot(theta, deltamax$V1, col = "red", "b", lwd = 2, ylab = "DeltaMax", 
#      main = "Delta Max (n=100000, L = 100)")
# plot(theta, IAE$V1, col = "blue", "b", lwd = 2, ylab = "IAE", 
#      main = "IAE (n=100000, L = 100)")
# 
# ######################################################
# 
# 
