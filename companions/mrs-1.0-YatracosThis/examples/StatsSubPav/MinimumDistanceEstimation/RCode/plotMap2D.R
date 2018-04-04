plotMap2D = function(mydata, plotWhat){

##Purpose: to plot mapped functions in 2D  
#Input arguments:
# mymydata    : the mydata (in mydata frame)
# hz        :   position of subpaving plot
# plotWhat  : 1 for height; 2 to compute fhat

##############################################################
n = NROW(xlow)
fhat = traincount/volBin/n

if(plotWhat == 1) {
  volBin=mydata[,2];
  fhat = mydata[,3]; 
  xlow=mydata[,4];
  xupp=mydata[,5];
  ylow=mydata[,6];
  yupp=mydata[,7];
} else {
  mydata = histpq
  volBin=mydata[,2];
  traincount = mydata[,3]
  xlow=mydata[,5];
  xupp=mydata[,6];
  ylow=mydata[,7];
  yupp=mydata[,8];
  fhat=trainCount/volBin/n;
}

box3D(x0 = xlow, y0 = ylow, z0 = rep(0, NROW(xlow)),
x1 = xupp, y1 = yupp, z1 = fhat,col = "pink", alpha = 0.3,  phi = 10, theta = 135,
border = "black", lwd = 1)

} #end of function plotMap2D