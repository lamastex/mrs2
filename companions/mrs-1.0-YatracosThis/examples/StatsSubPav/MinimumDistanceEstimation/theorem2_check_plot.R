rm(list = ls(all = TRUE))
cat("\014")
if(!is.null(dev.list())) dev.off()

#acer
setwd("C:/Users/Office/Dropbox/ResearchRaaz/theorems")

#campus
setwd("C:/Users/Office/Dropbox/ResearchRaaz/theorems")

##############################
#read txt file and store 
d = c(1, 2, 5, 10, 100, 1000) 
n = c("150", "1500", "15000", "150000", "1500000", "15000000")
num_n = NROW(n)
num_d = NROW(d)

avex = matrix(0, num_n, num_d)
avey = matrix(0, num_n, num_d)

for (i in 1:num_d){
  for (j in 1:num_n){
    filename = paste("uniform_", d[i], "d_", n[j], "n_theorem2_check_5sims.txt", sep = "")
    results = read.table(filename, header = FALSE, "\t")
    avex[j, i] = mean(results$V1)
    avey[j, i] = mean(results$V2)
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
mtext(text= "Theorem 2", side=3, line = 1, cex = 1.7, font = 2)

###############
mycol =  c("#00AFBB", "#E7B800", "#FC4E07", "#d55e00", "#cc79a7", "#0072b2", "#f0e442", "#009e73")
mypch = c(15, 16, 17, 18, 7, 9)
for (i in 1:num_d){
  for (j in 1:num_n){
    points(avex[j,i], avey[j,i], pch = mypch[j], cex = 1.5, col = mycol[i])
  }
}


#####################
vals = c(expression(paste('1D, ', "10"^"2")), expression(paste('1D, ', "10"^"3")), 
         expression(paste('1D, ', "10"^"4")), expression(paste('1D, ', "10"^"5")),
         expression(paste('1D, ', "10"^"6")), expression(paste('1D, ', "10"^"7")),
         expression(paste('2D, ', "10"^"2")), expression(paste('2D, ', "10"^"3")), 
         expression(paste('2D, ', "10"^"4")), expression(paste('2D, ', "10"^"5")),
         expression(paste('2D, ', "10"^"6")), expression(paste('2D, ', "10"^"7")),
         expression(paste('5D, ', "10"^"2")), expression(paste('5D, ', "10"^"3")), 
         expression(paste('5D, ', "10"^"4")), expression(paste('5D, ', "10"^"5")),
         expression(paste('5D, ', "10"^"6")), expression(paste('5D, ', "10"^"7")),
         expression(paste('10D, ', "10"^"2")), expression(paste('10D, ', "10"^"3")), 
         expression(paste('10D, ', "10"^"4")), expression(paste('10D, ', "10"^"5")),
         expression(paste('10D, ', "10"^"6")), expression(paste('10D, ', "10"^"7")),
         expression(paste('100D, ', "10"^"2")), expression(paste('100D, ', "10"^"3")), 
         expression(paste('100D, ', "10"^"4")), expression(paste('100D, ', "10"^"5")),
         expression(paste('100D, ', "10"^"6")), expression(paste('100D, ', "10"^"7")),
         expression(paste('1000D, ', "10"^"2")), expression(paste('1000D, ', "10"^"3")), 
         expression(paste('1000D, ', "10"^"4")), expression(paste('1000D, ', "10"^"5")),
         expression(paste('1000D, ', "10"^"6")), expression(paste('1000D, ', "10"^"7")))

c6 = c(rep(mycol[1], 6), rep(mycol[2], 6), rep(mycol[3], 6),
       rep(mycol[4], 6), rep(mycol[5], 6), rep(mycol[6], 6))


legend("bottomright", legend =vals,
       col = c6,
       ncol = 6, 
       cex = 0.8, 
       pch = mypch,
       text.font = 1, 
       text.col = c6)


#####################
