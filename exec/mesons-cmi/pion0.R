pion0 <- function(cmicor, mu=0.1, kappa=0.156, t1, t2, S=1.5, pl=FALSE, skip=0,
                variational=list(ta=3, tb=4, N=6), ind.vec=c(1,3,4,5),
                no.masses=1, matrix.size=2, boot.R=99, boot.l=10, tsboot.sim="geom",
                method="uwerr", fit.routine="optim", mass.guess, par.guess, nrep) {
  
  if(missing(cmicor)) {
    stop("Error! Data is missing!")
  }
  if(missing(t1) || missing(t2)) {
    stop("Error! t1 and t2 must be specified!")
  }
  if(missing(mass.guess)) {
    mass.guess <- c(0.2, 1., 3.)
  }
  else {
    if(length(mass.guess) < no.masses) {
      stop("mass.guess has not the correct length!")
    }
  }
  if(missing(par.guess)) {
    par.guess <- c(1.,0.8,0.1,0.1,0.1,0.1, 0.1,0.1,0.1,0.1,0.1,0.1, 0.1,0.1,0.1,0.1,0.1,0.1) 
  }
  else{
    if(length(par.guess) < no.masses*matrix.size) {
      stop("par.guess has not the correct length!")
    }
  }
  par <- numeric()
  length(par) <- no.masses*(matrix.size+1)
  for(i in 1:no.masses) {
    par[i*(matrix.size+1)] <- mass.guess[i]
    par[(c(1:matrix.size)+(i-1)*(matrix.size+1))] = par.guess[(c(1:matrix.size)+(i-1)*(matrix.size))]
  }
  
  Time <-  2*max(cmicor[,ind.vec[2]])
  Thalf <- max(cmicor[,ind.vec[2]])
  T1 <- Thalf+1
  t1p1 <- (t1+1)
  t2p1 <- (t2+1)
  ##nrObs <- max(cmicor[,ind.vec[1]])
  ## number of gamma matrix combinations
  ## in this case only 1 for the moment (only scalar)
  nrObs <- 1
  Skip <- (skip*(T1)*nrObs*4+1)
  Length <- length(cmicor[,ind.vec[3]])
  if(missing(nrep)) {
    nrep <- c(length(cmicor[((Skip):Length),ind.vec[3]])/(nrObs*(T1)*4))
  }
  else {
    skip <- 0
    if(sum(nrep) != length(cmicor[((Skip):Length),ind.vec[3]])/(nrObs*(T1)*4)) {
      stop("sum of replica differs from total no of measurements!")
    }
  }
  cat(nrObs, Skip, Length, "\n")

  W <- array(cmicor[((Skip):Length),ind.vec[3]], 
             dim=c(nrObs*(T1)*4,(length(cmicor[((Skip):Length),ind.vec[3]])/(nrObs*(T1)*4))))

  rm(cmicor)

  pion.eff.ll <- effectivemass(from=(t1+1), to=(t2+1), Time, W[1:T1,] , pl=FALSE, S=1.5, nrep=nrep)
  pion.eff.lf <- effectivemass(from=(t1+1), to=(t2+1), Time, W[(T1+1):(2*T1),] , pl=FALSE, S=1.5, nrep=nrep)
  pion.eff.ff <- effectivemass(from=(t1+1), to=(t2+1), Time, W[(3*T1+1):(4*T1),] , pl=FALSE, S=1.5, nrep=nrep)

  pion.eff <- data.frame(t=pion.eff.ll$t, mll=pion.eff.ll$mass, dmll=pion.eff.ll$dmass,
                         mlf=pion.eff.lf$mass, dmlf=pion.eff.lf$dmass,
                         mff=pion.eff.ff$mass, dmff=pion.eff.ff$dmass)

  Cor <- rep(0., times=4*T1)
  E <- rep(0., times=4*T1)
  
  for(i in 1:(4*T1)) {
    Cor[i] <- mean(W[(i),])
    tmpe <- try(uwerrprimary(W[(i),], pl=F, nrep=nrep)$dvalue, TRUE)
    if(!inherits(tmpe, "try-error")) E[i] = tmpe
    else {
      warning("error of correlator replaced by naive estimate!\n", call.=F)
      E[i] = sd(W[(i),])/sqrt(length(W[(i),]))
    }
  }

  ## Index vector of data to be used in the analysis
  ii <- c((t1p1):(t2p1), (t1p1+T1):(t2p1+T1), (t1p1+3*T1):(t2p1+3*T1))

  pionfit <- optim(par, ChiSqr.1mass, method="BFGS", control=list(trace=50),Thalf=Thalf,
                   x=c((t1):(t2)), y=Cor[ii], err=E[ii], tr = (t2-t1+1), N=2)
  pionfit <- optim(pionfit$par, ChiSqr.1mass, method="BFGS", control=list(trace=50, parscale=1./pionfit$par),Thalf=Thalf,
                   x=c((t1):(t2)), y=Cor[ii], err=E[ii], tr = (t2-t1+1), N=2)

  fit.mass <- abs(pionfit$par[3])
  fit.dof <- (t2-t1+1)*3-length(pionfit$par)
  fit.chisqr <- pionfit$value

  if(TRUE) {
    plot.effmass(m=fit.mass, ll=pion.eff.ll, lf=pion.eff.lf, ff=pion.eff.ff)
  }

  fit.uwerrm <- NULL
  fit.boot <- NULL
  fit.tsboot <- NULL

  if(method == "uwerr" || method == "all") {
    fit.uwerrm <- uwerr(f=fitmasses.pion, data=W[ii,], S=S, pl=pl, nrep=nrep,
                        Time=Time, t1=t1, t2=t2, Err=E[ii], par=par, N=2,
                        no.masses=1, fit.routine=fit.routine)
  }

  if(method == "boot" || method == "all") {
    fit.boot <- boot(data=t(W[ii,]), statistic=fit.pion.boot, R=boot.R, stype="i",
                     Time=Time, t1=t1, t2=t2, Err=E[ii], par=par, N=2, no.masses=1,
                     kappa=kappa, mu=mu, fit.routine=fit.routine)

    fit.tsboot <- tsboot(tseries=t(W[ii,]), statistic=fit.pion.boot, R=boot.R, l=boot.l, sim=tsboot.sim,
                         Time=Time, t1=t1, t2=t2, Err=E[ii], par=par, N=2, no.masses=1,
                         kappa=kappa, mu=mu, fit.routine=fit.routine)
  }

  Chi <- rep(0., times=4*T1)
  Fit <- rep(0., times=4*T1)

  jj <-  c(t1p1:t2p1)
  Fit[jj] <- pionfit$par[1]^2*CExp(m=fit.mass[1], Time=2*Thalf, x=jj-1)
  Fit[jj+T1] <- pionfit$par[1]*pionfit$par[2]*CExp(m=fit.mass[1], Time=2*Thalf, x=jj-1)
  Fit[jj+2*T1] <- pionfit$par[1]*pionfit$par[2]*CExp(m=fit.mass[1], Time=2*Thalf, x=jj-1)
  Fit[jj+3*T1] <- pionfit$par[2]*pionfit$par[2]*CExp(m=fit.mass[1], Time=2*Thalf, x=jj-1)

  Chi[ii] <- (Fit[ii]-Cor[ii])/E[ii]
  
  res <- list(fitresult=pionfit, t1=t1, t2=t2, N=length(W[1,]), Time=Time,
              fitdata=data.frame(t=(jj-1), Fit=Fit[ii], Cor=Cor[ii], Err=E[ii], Chi=Chi[ii]),
              uwerrresultmps=fit.uwerrm, uwerrresultmps2=NULL, uwerrresultmps3=NULL,
              uwerrresultfps=NULL, uwerrresultmpcac=NULL, uwerrresultzv=NULL,
              boot=fit.boot, tsboot=fit.tsboot, method=method,
              effmass=pion.eff, kappa=kappa, mu=mu, fit.routine=fit.routine,
              no.masses=1,
              matrix.size = 2, nrep=nrep)
  attr(res, "class") <- c("pionfit", "list")  
  return(invisible(res))
}
