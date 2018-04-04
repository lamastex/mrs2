subpavhist = function(mydata, str){
  
  ##this is from matlab - can modify later
## function [fhat, Lb, Ub, vol, count]=subpavhist(data, str)
#Calling Function:
  #           [fhat, Lb, Ub, vol, count]=subpavhist(data, str)
  #           
  #Purpose: To plot the histogram using subpaving boxes.
  #
  #Input argument list:
    #       data: from the AdaptiveHistogram output
  #       str: plot histogram if 1
  #
  #Output argument list:
    #       fhat: density estimate
  #       Lb  : lower bound
  #       Ub  : upper bound
  #       vol  : the Lebesgue measure
  #       count: number of points in each bin
  ###########################################################################
  
  #Ub and Lb of boxes (output of subpaving should already be in order)
  Lb=mydata[,5]
  Ub=mydata[,6]
  
  #Number of points in each box
  traincount=mydata[,3]
  validcount = mydata[,4]
  
  #Number of boxex
  n_box=NROW(Lb);
  
  #Total number of points
  N=sum(traincount);
  
  #bin width
  vol=mydata[,2];
  
  #Density estimate
  fhat=traincount/vol/N;
  
  #plot histogram
  if (str==1) {
    YY=matrix(c(0, 1, 1, 0)) %*% t(matrix(fhat)); YY = c(YY);
    XX=matrix(c(1, 1, 0, 0)) %*% t(matrix(Lb)) + matrix(c(0, 0, 1, 1))%*% t(matrix(Ub)); 
    XX=c(XX);
    plot(XX, YY, xlab = "x", ylab = expression("f"[n]), "l")
  } #end of plotting histogram
  
} #end of function subpavhist