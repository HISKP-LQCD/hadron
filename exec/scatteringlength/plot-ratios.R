require(tikzDevice)

tikz(paste("ratios.tex", sep=""), standAlone = TRUE, width=6, height=5)
par(cex=.7, cex.lab=1.5, cex.axis=1.5)

plot(NA, xlim=c(6,32), ylim=c(1.6,1.95), ylab=c("$R(t)$"), xlab=c("$t/a$"))

datafilelist = c("A40.32/A40.32.14-26.ratio.dat", "A40.24/A40.24.12-18.ratio.dat", "A60.24/A60.24.12-18.ratio.dat", "A80.24/A80.24.11-23.ratio.dat", "B35.32/B35.32.17-24.ratio.dat")
fitfilelist = c("A40.32/A40.32.16-31.ratiofit.dat", "A40.24/A40.24.14-23.ratiofit.dat", "A60.24/A60.24.14-23.ratiofit.dat", "A80.24/A80.24.14-23.ratiofit.dat", "B35.32/B35.32.18-31.ratiofit.dat")
colours = c("red", "blue", "darkgreen", "navy", "orange")
pchlist =c(21, 22, 23, 24, 25)
ii <- c(1:5)

for(i in c(1:length(datafilelist))) {
  data <- read.table(datafilelist[i])
  plotwitherror(x=data$V1, y=data$V2, dy=data$V3, pch=pchlist[i], col=colours[i], bg=colours[i], rep=TRUE)
  fit <- read.table(fitfilelist[i])
  lines(fit$V1, fit$V2, col=colours[i])
}
legend("topright", legend=c("A40.32", "A40.24", "A60.24", "A80.24", "B35.32"), pch=pchlist[ii], col=colours[ii], pt.bg=colours[ii], bty="n", cex=1.5)

dev.off()
tools::texi2dvi(paste("ratios.tex", sep=""), pdf=T)
