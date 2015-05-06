source("ratio.R")

t1 <- 14
t2 <- 24
tr1 <- 10
tr2 <- 20
ind.vector <- c(2,2)
##srcpath <- path
p1 <- c(0,0,0)
p2 <- c(0,0,0)

t0 <- 1
tlower <- seq(10,26,2)
tupper <- c(20,26,32)

trlower <- seq(6,20,2)
truppertp0 <- seq(15,31,2)
truppertpx <- seq(15,25,1)

irreplist <- c("A1", "E", "T2")
tplist <- list(c("TP0", "TP1", "TP2", "TP3"), c("TP0"), c("TP0"))
sizelist <- list(c(5,4,5,3), c(3), c(2))
idlist <- list(c(2,2,2,2), c(2), c(1))

source("parameters.R")

ld <- "lat"
if(!lat.disp) ld <- "cont"

for(i in c(1:length(irreplist))) {
  irrep <- irreplist[i]
  for(t in c(1:length(tplist[[i]]))) {
    tp <- tplist[[i]][t]
    if(tp != "TP0") mf <- TRUE
    else mf <- FALSE
    for(t1 in tlower) for(t2 in tupper) {
      if(!interactive()) pdf(file=paste("dEres-plots", tp, ".", irrep, ".", ens, ".t", t1, "-", t2, ".", ld, ".pdf", sep=""))    
      for(tr1 in trlower) {
        trupper <- truppertpx
        if(tp == "TP0") trupper <- truppertp0
        for(tr2 in trupper) {
          if(tr2 - tr1 > 5 && t2-t1 > 5) { 

            for(pc.id in c(1:idlist[[i]][t])) {
              if(tp == "TP0") {
                  p1 <- c(0,0,0)
                  p2 <- c(0,0,0)
                  ratiomethod <- "zeromom"
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
              
              cat(ens, irrep, tp, pc.id, t1, t2, tr1, tr2, p1, p2, "\n")

              filename <- paste(path, "Cmatrix.gevp.", tp, ".", irrep, ".", ens, ".t", t1, "-", t2, ".", ld, ".Rdata", sep="")
              if(t0 != 1) filename <- paste(path, "Cmatrix.gevp.", tp, ".", irrep, ".", ens, ".t0", t0, ".t", t1, "-", t2, ".", ld, ".Rdata", sep="")
              
              if(!file.exists(filename)) {
                pion.cor <- readtextcf(paste("pi_corr_p0.dat", sep=""), T=T, check.t=1, path=srcpath, ind.vector=ind.vector)
                pion.cor1 <- bootstrap.cf(readtextcf(paste("pi_corr_", sp1, ".dat", sep=""), T=T, check.t=1, path=srcpath, ind.vector=ind.vector), boot.R=boot.R, boot.l=boot.l, seed=seed)
                pion.cor2 <- bootstrap.cf(readtextcf(paste("pi_corr_", sp2, ".dat", sep=""), T=T, check.t=1, path=srcpath, ind.vector=ind.vector), boot.R=boot.R, boot.l=boot.l, seed=seed)
                Cmatrix <- getMatrix.pipi(N=matrix.size, path=srcpath, tp=tp, irrep=irrep, ens=ens, T=T)
                Cmatrix.gevp <- solveGEVP.pipi(Cmatrix, pion.cor=pion.cor, boot.R=boot.R, boot.l=boot.l, seed=seed, matrix.size=matrix.size,
                                               L=L, t1=t1, t2=t2, p1=p1, p2=p2, lat.disp=lat.disp, t0=t0)
                save(Cmatrix, Cmatrix.gevp, pion.cor, pion.cor1, pion.cor2, p1, p2, file=filename)
              }
              load(file=filename)

              dEres <- deltaEfromRatio(pipi.cor = Cmatrix.gevp$Cmatrix.gevp, pion = Cmatrix.gevp$pion, pion1 = pion.cor1, pion2 = pion.cor2,
                                       boot.R = boot.R, boot.l=boot.l, lat.disp=lat.disp, p1=p1, p2=p2, L=L,
                                       Thalf = T/2, tr1=tr1, tr2=tr2, pc.id=pc.id, pipimethod="gevp", ratiomethod="zeromom", mf=mf)
              effmass <- t(apply(dEres$Rpipi.tsboot, 1, effectivemass.cf, Thalf=T/2-1, type="log"))
              
              if(tp != "TP0") {
                dEresR2 <- deltaEfromRatio(pipi.cor = Cmatrix.gevp$Cmatrix.gevp, pion = Cmatrix.gevp$pion, pion1 = pion.cor1, pion2 = pion.cor2,
                                           boot.R = boot.R, boot.l=boot.l, lat.disp=lat.disp, p1=p1, p2=p2, L=L,
                                           Thalf = T/2, tr1=tr1, tr2=tr2, pc.id=pc.id, pipimethod="gevp", ratiomethod="mom", mf=mf)
                effmassR2 <- t(apply(dEresR2$Rpipi.tsboot, 1, effectivemass.cf, Thalf=T/2-1, type="log"))
              }
              else {
                dEresR2 <- NULL
                effmassR2 <- NULL
              }
              
              save(dEres, dEresR2, effmass, effmassR2, file=paste("dEres.", tp, ".", irrep, ".", ens, ".", pc.id, ".t0", t0, ".tr1.", tr1, ".tr2.", tr2, ".t1.", t1, ".t2.", t2, ".Rdata", sep=""))
              
              if(interactive()) X11()
              tt <- c((t0+1):(T/2))
              shift <- 0.5
              if(mf) shift <- 0.
              plotwitherror(x=tt-shift, y=dEres$Rpipi.tsboot[1,tt], dy=dEres$dRpipi[tt], main="Ratio Plot", xlab="t/a", ylab="R", log="y")
              if(tp != "TP0") plotwitherror(x=tt-shift, y=dEresR2$Rpipi.tsboot[1,tt], dy=dEresR2$dRpipi[tt], rep=TRUE, col="darkgreen")
              lines(x=dEres$x, y=dEres$rat, col="red")
              if(tp != "TP0") lines(x=dEres$x, y=dEresR2$rat, col="red")
              legend("topright", legend=c(paste("Ensemble:", ens), paste("Irrep:", irrep), paste("frame:", tp), paste("state:", pc.id), paste("p-value =", format(dEres$Qval, digits=3)),
                                     paste("deltaE=", format(dEres$opt.tsboot[1,2], digits=4), "+-", format(sd(dEres$opt.tsboot[,2]), digits=2)),
                                     paste("E=", format(dEres$opt.tsboot[1,5], digits=4), "+-", format(sd(dEres$opt.tsboot[,5]), digits=2)),
                                     paste("Mpi=", format(dEres$opt.tsboot[1,4], digits=4), "+-", format(sd(dEres$opt.tsboot[,4]), digits=2)),
                                     paste("p-value (pion) =", format(Cmatrix.gevp$pion$Qval, digits=3))
                                     ),
                     bty="n")
              legend("bottomleft", legend=c(paste("fitrange ratio t=[", tr1-shift, ",", tr2-shift, "]"), paste("fitrange pion t=[", t1, ",", t2, "]")), bty="n")
              if(tp != "TP0") legend("topleft", legend=c("R-Carsten", "R-Jinlong"), col=c("black", "darkgreen"), pch=c(21,21), bty="n")
              
              if(interactive()) X11()
              rc <- rchisq(n=boot.R, df=dEres$dof)
              qqplot(dEres$opt.tsboot[,3], rc, xlab=c("Sample Quantiels"), ylab=c("Theoretical Quantiles"), main=c("QQ-Plot Chisq"), xlim=c(0, 4*range(rc)[2]))

              cat("Carstens definition\n")
              cat("deltaE = ", dEres$opt.tsboot[1,2], "+-", sd(dEres$opt.tsboot[,2]), "\n")
              cat("E      = ", dEres$opt.tsboot[1,5], "+-", sd(dEres$opt.tsboot[,5]), "\n")
              cat("Mpi    = ", dEres$opt.tsboot[1,4], "+-", sd(dEres$opt.tsboot[,4]), "\n")
              cat("chisq =", dEres$opt.tsboot[1,3], "\n")
              cat("dof = ", dEres$dof, "\n")
              cat("Qval = ", dEres$Qval, "\n")
              cat("pionQval = ", dEres$pionQval, "\n\n")
              if(tp != "TP0") {
                cat("Jinlongs definition\n")
                cat("deltaE = ", dEresR2$opt.tsboot[1,2], "+-", sd(dEresR2$opt.tsboot[,2]), "\n")
                cat("E      = ", dEresR2$opt.tsboot[1,5], "+-", sd(dEresR2$opt.tsboot[,5]), "\n")
                cat("Mpi    = ", dEresR2$opt.tsboot[1,4], "+-", sd(dEresR2$opt.tsboot[,4]), "\n")
                cat("chisq =", dEresR2$opt.tsboot[1,3], "\n")
                cat("dof = ", dEresR2$dof, "\n")
                cat("Qval = ", dEresR2$Qval, "\n")
              }
            }
          }
        }
      }
      if(!interactive()) dev.off()
    }
  }
}


