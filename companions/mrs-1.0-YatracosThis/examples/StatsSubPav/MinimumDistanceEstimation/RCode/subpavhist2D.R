subpavhist2D = function(mydata, plotWhich){

  #Purpose: To plot the histogram using subpaving boxes.
  #
  #Input argument list:
  #       mydata     : data set
  #       plotWhich  : 1 if PQ hist (from ADH Validation), 2 if PCF
  #
  #########################################################################

n = nrow(mydata)
xlow = rep(0, n)
xupp = rep(0, n)
ylow = rep(0, n)
yupp = rep(0, n)

if(plotWhich == 1){
  xlow=mydata[,5];
  xupp=mydata[,6];
  ylow=mydata[,7];
  yupp=mydata[,8];
} else{
  xlow=mydata[,4];
  xupp=mydata[,5];
  ylow=mydata[,6];
  yupp=mydata[,7];
}

plot(c(min(xlow), max(xupp)), c(min(ylow), max(yupp)), type = "n", 
     xlab = expression("x"[1]), ylab = expression("x"[2]))
rect(xlow, ylow, xupp, yupp)

#rect(xleft, ybottom, xright, ytop, density = NULL, angle = 45,
#col = NA, border = NULL, lty = par("lty"), lwd = par("lwd"),
#...)

#later on can try to incorporate heat map (if i have energy)

} #end of function subpavhist2D
