rm(list = ls(all = TRUE))
setwd("~/git2019/mrs2/companions/mrs-1.0-YatracosThis/examples/StatsSubPav/Trajectory")
####################################################
par(mfrow = c(2,2))

start = 0
end = 3

for (k in start:end){
  dat = read.table(paste("spaceColl", k, ".txt", sep=""))

  # V4 = xlow, V5 = xupp, V6 = ylow, V7 = yupp 
  xx = matrix(c(dat$V4, dat$V5, dat$V5, dat$V4), ncol = 4)
  yy = matrix(c(dat$V6, dat$V6, dat$V7, dat$V7), ncol = 4)
  
  plot(xx, yy, type = "n")
  points(0, 0, col = "red", pch = 16)
  
  for (i in 1:NROW(dat)) {
    polygon(xx[i, ], yy[i, ])
  }
}