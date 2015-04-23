source("parameters.R")

##source("~/daten/workdir/hadron/exec/scatteringlength/ratio.R")
source("ratio.R")

ens <- "A40.32"

t1 <- 14
t2 <- 24
tr1 <- 10
tr2 <- 20
ind.vector <- c(2,2)
##srcpath <- path
p1 <- c(0,0,0)
p2 <- c(0,0,0)

irreplist <- c("A1", "E", "T2")
tplist <- list(c("TP0", "TP1", "TP2", "TP3"), c("TP0"), c("TP0"))
sizelist <- list(c(5,4,5,3), c(3), c(2))
idlist <- list(c(2,2,2,2), c(2), c(1))

for(i in c(1:length(irreplist))) {
  irrep <- irreplist[i]
  for(t in c(1:length(tplist[[i]]))) {
    tp <- tplist[[i]][t]
    if(tp != "TP0") mf <- TRUE
    else mf <- FALSE
    
    for(t1 in seq(10,26,2)) for(t2 in c(20,26,32)) {
      ##for(t1 in c(14)) {
      for(tr1 in seq(7,23,2)) {
        for(tr2 in seq(15,31,2)) {
          ##for(tr1 in c(14)) {
          ##for(tr2 in c(23)) {
          if(tr2 - tr1 > 5 && t2-t1 > 5) { 
            for(pc.id in c(1:idlist[[i]][t])) {
              if(tp == "TP0") {
                  p1 <- c(0,0,0)
                  p2 <- c(0,0,0)
              }
              if(tp == "TP1") {
                if(pc.id==1) {
                  p1 <- c(1,0,0)
                  p2 <- c(0,0,0)
                }
                if(pc.id==2) {
                  p1 <- c(1,0,0)
                  p2 <- c(-1,0,1)
                }
              }
              if(tp == "TP2") {
                if(pc.id==1) {
                  p1 <- c(1,1,0)
                  p2 <- c(0,0,0)
                }
                if(pc.id==2) {
                  p1 <- c(0,1,0)
                  p2 <- c(0,0,1)
                }
              }
              if(tp == "TP3") {
                if(pc.id==1) {
                  p1 <- c(1,1,1)
                  p2 <- c(0,0,0)
                }
                if(pc.id==2) {
                  p1 <- c(1,0,0)
                  p2 <- c(0,1,1)
                }
              }
              sp1 <- paste("p", sum(p1), sep="")
              sp2 <- paste("p", sum(p2), sep="")
              matrix.size <- sizelist[[i]][t]
              
              cat(ens, irrep, tp, pc.id, t1, t2, tr1, tr2, "\n")
              if(!file.exists(paste(path, "Cmatrix.gevp.", tp, ".", irrep, ".", ens, ".t", t1, "-", t2, ".Rdata", sep=""))) {
                pion.cor <- readtextcf(paste("pi_corr_p0.dat", sep=""), T=T, check.t=1, path=srcpath, ind.vector=ind.vector)
                pion.cor1 <- readtextcf(paste("pi_corr_", sp1, ".dat", sep=""), T=T, check.t=1, path=srcpath, ind.vector=ind.vector)
                pion.cor2 <- readtextcf(paste("pi_corr_", sp2, ".dat", sep=""), T=T, check.t=1, path=srcpath, ind.vector=ind.vector)
                Cmatrix <- getMatrix.pipi(N=matrix.size, path=srcpath, tp=tp, irrep=irrep, ens=ens, T=T)
                Cmatrix.gevp <- solveGEVP.pipi(Cmatrix, pion.cor=pion.cor, boot.R=boot.R, boot.l=boot.l, seed=seed, matrix.size=matrix.size,
                                               L=L, t1=t1, t2=t2, p1=p1, p2=p2)
                save(Cmatrix, Cmatrix.gevp, pion.cor, pion.cor1, pion.cor2, p1, p2, file=paste(path, "Cmatrix.gevp.", tp, ".", irrep, ".", ens, ".t", t1, "-", t2, ".Rdata", sep=""))
              }
              load(file=paste(path, "Cmatrix.gevp.", tp, ".", irrep, ".", ens, ".t", t1, "-", t2, ".Rdata", sep=""))
              
              dEres <- deltaEfromRatio(pipi.cor = Cmatrix.gevp$Cmatrix.gevp, pion = Cmatrix.gevp$pion, boot.R = boot.R, boot.l=boot.l,
                                       Thalf = T/2, tr1=tr1, tr2=tr2, pc.id=pc.id, pipimethod="gevp", mf=mf)
              
              save(dEres, file=paste("dEres.", tp, ".", irrep, ".", ens, ".", pc.id, ".tr1.", tr1, ".tr2.", tr2, ".t1.", t1, ".t2.", t2, ".Rdata", sep=""))
              
              if(interactive()) X11()
              tt <- c(2:(T/2))
              shift <- 0.5
              if(mf) shift <- 0.
              plotwitherror(x=tt-shift, y=dEres$Rpipi.tsboot[1,tt], dy=dEres$dRpipi[tt], main="Ratio Plot", xlab="t/a", ylab="R", log="y")
              lines(x=dEres$x, y=dEres$rat, col="red")
              legend("topright", legend=c(paste("Ensemble:", ens), paste("Irrep:", irrep), paste("frame:", tp), paste("state:", pc.id), paste("p-value =", format(dEres$Qval, digits=3)),
                                     paste("deltaE=", format(dEres$opt.tsboot[1,2], digits=4), "+-", format(sd(dEres$opt.tsboot[,2]), digits=2)),
                                     paste("E=", format(dEres$opt.tsboot[1,5], digits=4), "+-", format(sd(dEres$opt.tsboot[,5]), digits=2)),
                                     paste("Mpi=", format(dEres$opt.tsboot[1,4], digits=4), "+-", format(sd(dEres$opt.tsboot[,4]), digits=2)),
                                     paste("p-value (pion) =", format(Cmatrix.gevp$pion$Qval, digits=3))
                                     ),
                     bty="n")
              legend("bottomleft", legend=c(paste("fitrange ratio t=[", tr1-shift, ",", tr2-shift, "]"), paste("fitrange pion t=[", t1, ",", t2, "]")), bty="n")
              
              if(interactive()) X11()
              rc <- rchisq(n=boot.R, df=dEres$dof)
              qqplot(dEres$opt.tsboot[,3], rc, xlab=c("Sample Quantiels"), ylab=c("Theoretical Quantiles"), main=c("QQ-Plot Chisq"), xlim=c(0, 4*range(rc)[2]))
              
              cat("deltaE = ", dEres$opt.tsboot[1,2], "+-", sd(dEres$opt.tsboot[,2]), "\n")
              cat("E      = ", dEres$opt.tsboot[1,5], "+-", sd(dEres$opt.tsboot[,5]), "\n")
              cat("Mpi    = ", dEres$opt.tsboot[1,4], "+-", sd(dEres$opt.tsboot[,4]), "\n")
              cat("chisq =", dEres$opt.tsboot[1,3], "\n")
              cat("dof = ", dEres$dof, "\n")
              cat("Qval = ", dEres$Qval, "\n")
              cat("pionQval = ", dEres$pionQval, "\n\n")
            }
          }
        }
      }
    }
  }
}
