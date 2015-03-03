require(tikzDevice)

types <- c("-efm", "-ratio", "")
odata <- read.table("scatteringlength.dat")

for(j in c(1:length(types))) {
  load(paste("res", types[j], "-allens.Rdata", sep=""))
  
  x <- odata$V4/odata$V5
  
  tikz(paste("mpia0", types[j], ".tex", sep=""), standAlone = TRUE, width=5, height=5)
  
  ii <- c(1,2,5:7)
  plotwitherror(x=x[ii], y=res[ii,9], dy=res[ii,10], xlab=c("$M_{\\pi}/f_{\\pi}$"), ylab=c("$M_{\\pi}a_0$"), col=c("red"), pch=21, xlim=c(0,3), ylim=c(-0.35,0))
  ii <- c(8,10,11)
  plotwitherror(x=x[ii], y=res[ii,9], dy=res[ii,10], xlab=c("$M_{\\pi}/f_{\\pi}$"), ylab=c("$M_{\\pi}a_0$"), col=c("blue"), pch=22, xlim=c(0,4), rep=TRUE)
  ii <- c(9)
  plotwitherror(x=x[ii], y=res[ii,9], dy=res[ii,10], xlab=c("$M_{\\pi}/f_{\\pi}$"), ylab=c("$M_{\\pi}a_0$"), col=c("darkgreen"), pch=23, xlim=c(0,4), rep=TRUE)
  
  legend(x="bottomleft", legend=c("A ensembles", "B ensembles", "D ensembles"), pt.bg=c("white", "white", "white"), col=c("red", "blue", "darkgreen"), pch=c(21:23))
  dev.off()
  tools::texi2dvi(paste("mpia0", types[j], ".tex", sep=""), pdf=T)
  
  tikz(paste("mpia0-fs", types[j], ".tex", sep=""), standAlone = TRUE, width=5, height=5)
  
  x <- odata$V4/odata$V10/odata$V5*odata$V12
  L <- odata$V2
  a0 <- rep(0, times=length(res[,1]))
  jj <- c(1:length(res[,1]))
  for(i in jj) {
    a0[i] <- fs.a0(a0=res[i,5], L=L[i], mps=odata$V4[i])
  }
  res[,9] <- res[,9]/odata$V10[jj]/res[,5]*a0
  ii <- c(1,2,5:7)
  plotwitherror(x=x[ii], y=res[ii,9], dy=res[ii,10], xlab=c("$M_{\\pi}/f_{\\pi}$"), ylab=c("$M_{\\pi}a_0$"), col=c("red"), pch=21, xlim=c(0,3), ylim=c(-0.35,0))
  ii <- c(8,10,11)
  plotwitherror(x=x[ii], y=res[ii,9], dy=res[ii,10], xlab=c("$M_{\\pi}/f_{\\pi}$"), ylab=c("$M_{\\pi}a_0$"), col=c("blue"), pch=22, xlim=c(0,4), rep=TRUE)
  ii <- c(9)
  plotwitherror(x=x[ii], y=res[ii,9], dy=res[ii,10], xlab=c("$M_{\\pi}/f_{\\pi}$"), ylab=c("$M_{\\pi}a_0$"), col=c("darkgreen"), pch=23, xlim=c(0,4), rep=TRUE)
  
  legend(x="bottomleft", legend=c("A ensembles", "B ensembles", "D ensembles"), pt.bg=c("white", "white", "white"), col=c("red", "blue", "darkgreen"), pch=c(21:23))
  ##dev.copy2pdf(file="qcotdelta.pdf")
  dev.off()
  tools::texi2dvi(paste("mpia0-fs", types[j], ".tex", sep=""), pdf=T)
}

## produce empty plot with correct limits
## plot(x,y, type="n")
## add shaded band
## x <- c(1:10)
## polygon(x=c(x, rex(x)), y=c(x^2+1, rev(x^2)-1), col="gray", lty=0, lwd=0.001, border="red")
## plot points on top of it
## points(x,y)
