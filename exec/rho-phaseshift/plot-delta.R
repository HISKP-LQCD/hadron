source("parameters.R")

cframes <- toupper(frames)
pch <- c(21, 22, 23, 24, 25)
col <- c("red", "darkblue", "darkgreen", "purple", "orange", "black")

load("pion.Rdata")
Mpi <-pion.matrixfit$opt.res$par[1]

tikzfiles <- tikz.init(basename=paste("delta-", ens[1], sep=""), width=5., height=4.)

plot(NA, xlim=c(2*Mpi,4*Mpi), ylim=c(0,pi), xlab=c("$aE_\\mathrm{CM}$"), ylab=c("$\\delta\\,[\\mathrm{rad}]$"))

for(e in ens) {
  for(f in frames) {
    if(f == "mf4") {
      file <- paste(toupper(f), "/res.", e, "mf1.Rdata", sep="")
    }
    else {
      file <- paste(toupper(f), "/res.", e, f, ".Rdata", sep="")
    }
    cat("loading file", file, "\n")
    load(file)
    i <- which(frames == f)
    arrows(res.all.pc1$Ecm[1], res.all.pc1$delta[1]+res.all.pc1$delta[3],
           res.all.pc1$Ecm[1], res.all.pc1$delta[1]+res.all.pc1$delta[4], length=0.01,angle=90,code=3, col="black", lwd=0.5)
    arrows(res.all.pc1$Ecm[1]+res.all.pc1$Ecm[3], res.all.pc1$delta[1],
           res.all.pc1$Ecm[1]+res.all.pc1$Ecm[4], res.all.pc1$delta[1], length=0.01,angle=90,code=3, col="black", lwd=0.5)
    plotwitherror(x=res.all.pc1$Ecm[1], dx=res.all.pc1$Ecm[2], y=res.all.pc1$delta[1], dy=res.all.pc1$delta[2], rep=TRUE, col=col[i], bg=col[i], pch=pch[1])
    cat("Ecm1\t", res.all.pc1$Ecm[1], res.all.pc1$Ecm[2], res.all.pc1$Ecm[3], res.all.pc1$Ecm[4], "\n")
    cat("delta1\t", res.all.pc1$delta[1], res.all.pc1$delta[2], res.all.pc1$delta[3], res.all.pc1$delta[4], "\n")

    arrows(res.all.pc2$Ecm[1], res.all.pc2$delta[1]+res.all.pc2$delta[3],
           res.all.pc2$Ecm[1], res.all.pc2$delta[1]+res.all.pc2$delta[4], length=0.01,angle=90,code=3, col="black", lwd=0.5)
    arrows(res.all.pc2$Ecm[1]+res.all.pc2$Ecm[3], res.all.pc2$delta[1],
           res.all.pc2$Ecm[1]+res.all.pc2$Ecm[4], res.all.pc2$delta[1], length=0.01,angle=90,code=3, col="black", lwd=0.5)
    plotwitherror(x=res.all.pc2$Ecm[1], dx=res.all.pc2$Ecm[2], y=res.all.pc2$delta[1], dy=res.all.pc2$delta[2], rep=TRUE, col=col[i], bg=col[i], pch=pch[2])
    cat("Ecm2\t", res.all.pc2$Ecm[1], res.all.pc2$Ecm[2], res.all.pc2$Ecm[3], res.all.pc2$Ecm[4], "\n")
    cat("delta2\t", res.all.pc2$delta[1], res.all.pc2$delta[2], res.all.pc2$delta[3], res.all.pc2$delta[4], "\n")
  }
}

legend("topleft", legend=toupper(frames), bty="n", text.col=col)
legend("bottomright", legend=c("ground state","1st excited state"), bty="n", pch=pch[1:2], col="black", pt.bg="black")

tikz.finalize(tikzfiles=tikzfiles)
