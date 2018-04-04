plotMap1D = function(mydata, str){
##Purpose: to plot mapped functions in 1D  
#Input arguments:
# mydata: the data (in data frame)
# str: if str = 1, a plot will be output
################################################
    #bin width
  vol=mydata[,2];
 
  #Density f
  range = mydata[,3];
  
  #Ub and Lb of boxes (output of subpaving should already be in order)
  Lb=mydata[,4]
  Ub=mydata[,5]
  
  #Number of boxex
  n_box=NROW(Lb);
  
  #plot histogram
  if (str==1) {
    YY=matrix(c(0, 1, 1, 0)) %*% t(matrix(range)); YY = c(YY);
    XX=matrix(c(1, 1, 0, 0)) %*% t(matrix(Lb)) + matrix(c(0, 0, 1, 1))%*% t(matrix(Ub)); 
    XX=c(XX);
    plot(XX, YY, xlab = "x", ylab = "f", "l")
    #polygon(XX, YY, col = "white", border = "white")
  } #end of plotting histogram
  
} #end of function plotMap1D