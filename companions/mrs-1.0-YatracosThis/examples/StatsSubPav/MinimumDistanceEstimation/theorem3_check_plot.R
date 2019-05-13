rm(list = ls(all = TRUE))
cat("\014")

#acer
setwd("C:/Users/Office/Dropbox/ResearchRaaz/theorems")

#campus
setwd("C:/Users/Office/Dropbox/ResearchRaaz/theorems")

############################
#read txt file and store 
d = c(2, 3) 
n = c(150, 1500)
num_n = NROW(n)
num_d = NROW(d)
theta = 10
phi = 0.3
rhs_term_1 = (1 + 2*phi/(1-phi) + 8*sqrt(phi))

avex = matrix(0, num_n, num_d)
avey = matrix(0, num_n, num_d)

for (i in 1:num_d){
  for (j in 1:num_n){
    filename = paste("uniform_", d[i], "d_", n[j], "n_theorem3_iaes_5sims.txt", sep = "")
    results = read.table(filename, header = FALSE, "\t")
    avex[j, i] = mean(results$V1)
    avey[j, i] = mean(results$V2) * rhs_term_1 + 8*sqrt((log(2*theta*(theta + 1)))/(phi * n))
  }
}

##########################
#build an empty plot

plot(NULL, xlab = "", ylab = "", cex = 1.5, 
     axes = FALSE, xlim = c(0, max(avey)), ylim = c(0, max(avey)))
box(lwd = 2)
abline(a = 0, b = 1, lty = 2)

xlabels = seq(0, max(avey), 1)
axis(1, at = xlabels, labels = NA, lwd = 2)
axis(1, at = xlabels, labels = xlabels, lwd = 0, line = -0.4, hadj = 0.4, cex.axis = 1.2, font.axis = 2)


ylabels = seq(0, max(avey), 1)
axis(2, at = ylabels, labels = NA)
axis(2, at = ylabels, labels = ylabels, lwd = 0, las = 2, line = 0,
     cex.axis = 1.2, font.axis = 2)

#Plot title and axis titles
mtext(text="LHS", side=1, line=2.2, cex = 1.7, font = 2)
mtext(text="RHS", side=2, line=2, cex = 1.7, font = 2)
mtext(text= "Theorem 3", side=3, line = 1, cex = 1.7, font = 2)

###############
mycol = c("blue", "black", "red")
#mycol =  c("#00AFBB", "#E7B800", "#FC4E07", "#d55e00", "#cc79a7", "#0072b2", "#f0e442", "#009e73")
mypch = c(3, 4, 15, 16, 17, 18)

for (i in 1:num_d){
  for (j in 1:num_n){
    points(avex[j,i], avey[j,i], pch = mypch[j], cex = 1.5, col = mycol[i])
  }
}


#####################
legend("bottomright", legend = c("2D, 150", "2D, 1500", "3D, 150", "3D, 1500"), 
       col = c(mycol[1], mycol[1], mycol[2], mycol[2]), 
       pch = c(mypch[1], mypch[2], mypch[1], mypch[2]), cex = 1)


#####################
