##MappedIAE.R
##Purpose:
#To plot the IAEs for histograms to determine the minimum IAE.
#Date updated: 1 September 2016

###########################################################
rm(list = ls(all = TRUE))
dev.off() #close all plots
cat("\014") #clear console

mydir = "~/mrs/branches/gat41/examples/StatsSubPav/IAECalculations/output/" 
setwd(mydir)

#change the method and directory to the folder of the density that you want 
mydensity = paste("2DExample/")
mydensitydir = paste(mydir, mydensity, sep="")
mydensitydir
setwd(mydensitydir)

#read in tables
IAE = read.table("MappedIAE1.txt", header = FALSE)
attach(IAE)
deltamax = read.table("MappedDeltaMax1.txt", sep="", header = FALSE)
sequence = read.table("Sequence1.txt", sep = "", header=FALSE)

#plots
plot(1:length(V1), log(V1), col = "blue", "b", lwd = 2, xlab = "histogram", ylab = "IAE (log)", 
     pch = 20, main = "IAE")
m = which.min(V1)
points(m, log(V1[m]), col = "red", pch = 20)
text(m, log(V1[m]) + 0.5, m)

#Delta Max
plot(sequence$V1, log(deltamax$V1), col = "blue", "b", lwd = 2, xlab = "histogram", ylab = "IAE (log)", 
     pch = 20, main = "Delta Max")
m = which.min(deltamax$V1)
points(sequence$V1[m], log(deltamax$V1[m]), col = "red", pch = 20)
text(sequence$V1[m], log(deltamax$V1[m]) + 0.4, sequence$V1[m])

deltamin = sequence$V1[m]
deltamin
IAE$V1[deltamin]

IAEmin = which.min(IAE$V1)
IAEmin
IAE$V1[IAEmin]

detach(IAE)

