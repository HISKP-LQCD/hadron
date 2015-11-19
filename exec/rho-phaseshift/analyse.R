source("/home/urbach/daten/workdir/hadron/exec/rho-phaseshift/phaseshift.rho.R")
source("/home/urbach/daten/workdir/hadron/exec/rho-phaseshift/summarise.R")
reread <- TRUE
n <- 1
t10 <- 7
t11 <- 15
source("parameters.R")

Thalf <- T/2
load("../dataA40.24/pion.Rdata")

if(reread || !file.exists(paste("./","Cmatrix.", ens, frame, ".Rdata", sep=""))) {
  Cmatrix <- cf()
  tmp <- readtextcf(paste("rho_corr_", tp, "_00.dat", sep=""), T=T, check.t=0, path=path)
  Cmatrix <- c(Cmatrix, tmp)
  tmp <- readtextcf(paste("rho_corr_", tp, "_01.dat", sep=""), T=T, check.t=0, path=path)
  Cmatrix <- c(Cmatrix, tmp)
  tmp <- readtextcf(paste("rho_corr_", tp, "_11.dat", sep=""), T=T, check.t=0, path=path)
  
  Cmatrix <- c(Cmatrix, tmp)
  Cmatrix <- bootstrap.cf(Cmatrix, boot.R=boot.R, boot.l=boot.l, seed=1234)
  Cmatrix.rt <- removeTemporal.cf(Cmatrix, p1=p, p2=c(0,0,0), single.cf1=pion.matrixfit, L=L, lat.disp=TRUE)

  save(Cmatrix, Cmatrix.rt, file=paste("./","Cmatrix.", ens, frame, ".Rdata", sep=""))
}
load(paste("./","Cmatrix.", ens, frame, ".Rdata", sep=""))

Cmatrix.bootstrap.gevp.n <- bootstrap.gevp(Cmatrix, matrix.size=2, element.order=c(1,2,2,3))
pc.n <- gevp2cf(Cmatrix.bootstrap.gevp.n, id=1)
pc2.n <- gevp2cf(Cmatrix.bootstrap.gevp.n, id=2)

Cmatrix.bootstrap.gevp <- bootstrap.gevp(Cmatrix.rt, matrix.size=2, element.order=c(1,2,2,3))
pc <- gevp2cf(Cmatrix.bootstrap.gevp, id=1)
pc2 <- gevp2cf(Cmatrix.bootstrap.gevp, id=2)

pc.effectivemass.n <- bootstrap.effectivemass(cf=pc.n, type="solve", boot.R=boot.R, boot.l=boot.l)
pc2.effectivemass.n <- bootstrap.effectivemass(cf=pc2.n, type="solve", boot.R=boot.R, boot.l=boot.l)

pc.effectivemass <- bootstrap.effectivemass(cf=pc, type=type, boot.R=boot.R, boot.l=boot.l)
pc2.effectivemass <- bootstrap.effectivemass(cf=pc2, type=type, boot.R=boot.R, boot.l=boot.l)

for(t1 in seq(t10,t11,1)) {
  for(t2 in seq(t11,T/2,1)) {
    if(t2-t1>5) {
      cat(t1, t2, "\n")
      pc.matrixfit.n <- matrixfit(pc.n, t1=t1, t2=t2, useCov=TRUE, parlist=array(c(1,1), dim=c(2,1)), sym.vec=c("cosh"), fit.method="lm")
      pc2.matrixfit.n <- matrixfit(pc2.n, t1=t1, t2=t2, useCov=TRUE, parlist=array(c(1,1), dim=c(2,1)), sym.vec=c("cosh"), fit.method="lm")
      
      
      ## 10 to 14
      pc.matrixfit <- matrixfit(pc, t1=t1, t2=t2, useCov=TRUE, parlist=array(c(1,1), dim=c(2,1)), sym.vec=c("cosh"), fit.method="lm")
      ## 10 to 15
      pc2.matrixfit <- matrixfit(pc2, t1=t1, t2=t2, useCov=TRUE, parlist=array(c(1,1), dim=c(2,1)), sym.vec=c("cosh"), fit.method="lm")

      Rpipi.tsboot <- array(NA, dim=c(boot.R+1, Thalf))
      if(TRUE || frame == "cmf") {
        ## this is R(t+1/2)
        Rpipi.tsboot[1,] <- compRpipi2(c4=pc$cf0, c2=pion.cor$cf0, Thalf=Thalf)

        for(i in c(1:boot.R)) {
          Rpipi.tsboot[i+1,] <- compRpipi2(c4=pc$cf.tsboot$t[i,], c2=pion.cor$cf.tsboot$t[i,], Thalf=Thalf)
        }
      }
      else {
        Mps <- rep(0, times=boot.R+1)
        Mps[1] <- pion.matrixfit$opt.res$par[1]
        Mps[c(1:boot.R)+1] <- pion.matrixfit$opt.tsboot[1,]
          
        dE <- rep(0, times=boot.R+1)
        p2 <- c(0,0,0)
        if(disp=="lat") {
          p1shift <- 2*sum(sin(pi*p/L)^2)
          p2shift <- 2*sum(sin(pi*p2/L)^2)
          dE <- abs(acosh( cosh(Mps) + p1shift) - acosh( cosh(Mps) + p2shift))
        }
        else {
          p1shift <- sum((2*pi*p/L)^2)
          p2shift <- sum((2*pi*p2/L)^2)
          dE <- abs(sqrt( Mps^2 + p1shift ) - sqrt( Mps^2 + p2shift ))
        }

        ## this is R(t)
        Rpipi.tsboot[1,] <- compRpipi3(c4=pc$cf0, c21=pion1$cf0, c22=pion2$cf0, Thalf=Thalf, dE=dE[1])
        
        for(i in c(1:boot.R)) {
          Rpipi.tsboot[i+1,] <- compRpipi3(c4=pc$cf.tsboot$t[i,], c21=pion1$cf.tsboot$t[i,], c22=pion2$cf.tsboot$t[i,], Thalf=Thalf, dE=dE[i+1])
        }
      }
      
      gs <- phaseshift.rho(pcfit =pc.matrixfit, L=L, Mpi=pion.matrixfit$opt.res$par[1], Mpiboot=pion.matrixfit$opt.tsboot[1,], frame=frame, disp=disp, n=n)
      fes <- phaseshift.rho(pcfit =pc2.matrixfit, L=L, Mpi=pion.matrixfit$opt.res$par[1], Mpiboot=pion.matrixfit$opt.tsboot[1,], frame=frame, disp=disp, n=n)

      save(Cmatrix.bootstrap.gevp, Cmatrix.bootstrap.gevp.n, pc.matrixfit, pc.matrixfit.n,
           pc.effectivemass, pc.effectivemass.n, pc2.matrixfit, pc2.matrixfit.n,
           pc2.effectivemass, pc2.effectivemass.n, type, Rpipi.tsboot,
           L, T, pion.matrixfit, gs, fes, frame, ens, disp, file=paste("rhoana.", t1, "-", t2, ".", ens, frame, ".Rdata", sep=""))
      
      
      cat(t1, t2, gs$Ecm, sd(gs$Ecmboot, na.rm=TRUE), gs$delta, sd(gs$deltaboot, na.rm=TRUE), sin(gs$delta)^2, sd(sin(gs$deltaboot)^2, na.rm=TRUE),
          gs$tandelta, sd(gs$tandeltaboot, na.rm=TRUE), "Qval=", pc.matrixfit$Qval, "\n")
      cat(t1, t2, fes$Ecm, sd(fes$Ecmboot, na.rm=TRUE), fes$delta, sd(fes$deltaboot, na.rm=TRUE), sin(fes$delta)^2, sd(sin(fes$deltaboot)^2, na.rm=TRUE),
          fes$tandelta, sd(fes$tandeltaboot, na.rm=TRUE), "Qval=", pc2.matrixfit$Qval, "\n")
    }
  }
}


res <- summarise.rho(ens, frame)

res.all <- compute.error.rho(res)

res.boot <- array(0, dim=c(boot.R+1, 4))
for(i in c(1:4)) {
  res.boot[,i] <- compute.boots(res, index=i)
}

save(res, res.all, res.boot, file=paste("res.", ens, frame, ".Rdata", sep=""))
