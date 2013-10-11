onlinemeas <- function(data, t1, t2, S=1.5, pl=FALSE, skip=0,
                       iobs=1, ind.vec=c(1,3,4,5), mu=0.1, kappa=0.125,
                       boot.R=99, boot.l=10, tsboot.sim="geom",
                       method="uwerr", fit.routine="optim", nrep,
                       oldnorm = FALSE) {
  
  if(missing(data)) {
    stop("Error! Data is missing!")
  }
  if(missing(t1) || missing(t2)) {
    stop("Error! t1 and t2 must be specified!")
  }
  par <- numeric(2)
  sign <- +1.
  # in the data which was read with readcmicor, load only the data
  # for PP and PA (V1 = 1 or 2)
  data <- data[data$V1<3,]
  # and gamma matrix combination iobs (by default iobs=1 <-> gamma5)
  data <- data[data$V2==iobs,]

  Time <-  2*max(data[,ind.vec[2]])
  Thalf <- max(data[,ind.vec[2]])
  T1 <- Thalf+1
  t1p1 <- (t1+1)
  t2p1 <- (t2+1)
  # number of observables in data
  nrObs <- 2
  # number of instances of certain observables (i.e.: smeared-smeared, local-smeared, smeared-local, local-local would be 4) 
  nrType <- 1
  Skip <- (skip*(T1)*nrType*nrObs+1)
  Length <- length(data[,ind.vec[3]])
  cat("time =", Time, "Thalf =", Thalf, "\n")
  if(missing(nrep)) {
    nrep <- c(length(data[((Skip):Length),ind.vec[3]])/(nrObs*(T1)*nrType))
  }
  else {
    skip <- 0
    if(sum(nrep) != length(data[((Skip):Length),ind.vec[3]])/(nrObs*(T1)*nrType)) {
      stop("sum of replica differs from total no of measurements!")
    }
  }

  # extract the negative and positive correlators from the data
  Z <- array(data[((Skip):Length),ind.vec[3]], 
             dim=c(nrObs*(T1)*nrType,(length(data[((Skip):Length),ind.vec[3]])/(nrObs*(T1)*nrType))))
  # negative times
  W <- array(data[((Skip):Length),ind.vec[4]], 
             dim=c(nrObs*(T1)*nrType,(length(data[((Skip):Length),ind.vec[4]])/(nrObs*(T1)*nrType))))

  rm(data)
  # extract the correlator from the data, averaging negative and positive times
  # and store it back to W
  W <- getCor(T1=T1, W=W, Z=Z, type=c("cosh","sinh"))
  if(oldnorm) {
    W <- W*(2*kappa)^2
  }
  rm(Z)

  # some effective pion masses
  eff <- effectivemass(from=(t1+1), to=(t2+1), Time, W[1:T1,] , pl=FALSE, S=1.5, nrep=nrep)
  options(show.error.messages = TRUE)

  mass.eff <- data.frame(t=eff$t, m=eff$mass, dm=eff$dmass)

  Cor <- rep(0., times=nrObs*T1)
  E <- rep(0., times=nrObs*T1)
  
  # W contains length(W[(i),]) samples of the PA and PP correlators
  # rows correspond to timeslices and columns to samples
  # Cor will therefore store averages of the two correlators, P from 1 to T1
  # and A from T1+1 to 2T1 
  for(i in 1:(T1*nrObs)) {
    Cor[i] <- mean(W[(i),])
    tmpe <- try(uwerrprimary(W[(i),], pl=F, nrep=nrep)$dvalue, TRUE)
    if(!inherits(tmpe, "try-error")) E[i] = tmpe
    else {
      warning("error of correlator replaced by naive estimate!\n", call.=F)
      E[i] = sd(W[(i),])/sqrt(length(W[(i),]))
    }
  }

  # now fit the PP correlator first _only_ 
  par[2] <- eff$mass[1]
  par[1] <- sqrt(abs(Cor[(t1+1)]*exp(par[2]*(t1+1))))
  
  # Index vector of data to be used in the analysis
  ii <- c((t1p1):(t2p1))

  #BFGS
  if(fit.routine != "gsl") {
    massfit <- optim(par, ChiSqr.singleCor, method="BFGS", control=list(trace=0),Thalf=Thalf,
                     x=c((t1):(t2)), y=Cor[ii], err=E[ii], tr = (t2-t1+1), sign=sign)
  }
  else {
    massfit <- gsl_fit_correlator(par, Thalf=Thalf,
                                  x=c((t1):(t2)), y=Cor[ii], err=E[ii], tr = (t2-t1+1))
  }
  sfit.mass <- abs(massfit$par[2])
  sfit.fpi <- 2*kappa*2*mu/sqrt(2)*abs(massfit$par[1])/sqrt(sfit.mass^3)
  if(fit.routine != "gsl" && massfit$convergence!=0) {
    warning("optim did not converge for massfit! ", massfit$convergence)
  }
  sfit.fpi <- 2*kappa*2*mu/sqrt(2)*abs(massfit$par[1])/sqrt(sfit.mass^3)
  cat("mpi =", sfit.mass, " fpi =",sfit.fpi, "\n")
  
  sfit.dof <- (t2-t1+1)-length(massfit$par)
  sfit.chisqr <- massfit$value

  sfit.uwerrm <- NULL
  sfit.boot <- NULL
  sfit.tsboot <- NULL
  if(method == "uwerr" || method == "all") {
    sfit.uwerrm <- uwerr(f=fitmass, data=W[ii,], S=S, pl=pl, nrep=nrep,
                        Time=Time, t1=t1, t2=t2, Err=E[ii], par=par, sign=sign,
                        fit.routine=fit.routine)
  }
  if(method == "boot" || method == "all") {
    sfit.boot <- boot(data=t(W[ii,]), statistic=fit.boot, R=boot.R, stype="i",
                     Time=Time, t1=t1, t2=t2, Err=E[ii], par=par, sign=sign,
                     fit.routine=fit.routine)

    sfit.tsboot <- tsboot(tseries=t(W[ii,]), statistic=fit.boot, R=boot.R, l=boot.l, sim=tsboot.sim,
                         Time=Time, t1=t1, t2=t2, Err=E[ii], par=par, sign=sign,
                         fit.routine=fit.routine)
  }

  
  # now do pcac mass, i.e. PP and PA
  # Index vector of data to be used in the analysis
  ii <- c((t1p1):(t2p1), (t1p1+T1):(t2p1+T1))

  par2 <- c(1, 0.1, 0.1)
  par2[3] <- eff$mass[1]
  par2[1] <- sqrt(abs(Cor[(t1+1)]*exp(par[2]*(t1+1))))
  par2[2] <- 0.
                                        #BFGS
  if(fit.routine != "gsl75") {
    pcacfit <- optim(par2, ChiSqr.pcac, method="BFGS", control=list(trace=0),Thalf=Thalf,
                     x=c((t1):(t2)), y=Cor[ii], err=E[ii], tr = (t2-t1+1))
  }
  else {
    pcacfit <- gsl_fit_correlator(par, Thalf=Thalf,
                                  x=c((t1):(t2)), y=Cor[ii], err=E[ii], tr = (t2-t1+1))
  }
  fit.pcac <- 0.5*abs(pcacfit$par[3])*pcacfit$par[2]/pcacfit$par[1]

  if(fit.routine != "gsl" && pcacfit$convergence!=0) {
    warning("optim did not converge for pcacfit! ", massfit$convergence)
  }
  fit.dof <- (t2-t1+1)-length(massfit$par)
  fit.chisqr <- massfit$value

  fit.uwerrm <- NULL
  fit.uwerrpcac <- NULL
  fit.uwerrfpi <- NULL
  fit.boot <- NULL
  fit.tsboot <- NULL
  if(method == "uwerr" || method == "all") {
    fit.uwerrm <- uwerr(f=fitmass.online, data=W[ii,], S=S, pl=pl, nrep=nrep,
                      Time=Time, t1=t1, t2=t2, Err=E[ii], par=par2,
                      fit.routine=fit.routine)
    fit.uwerrpcac <- uwerr(f=fitmpcac.online, data=W[ii,], S=S, pl=pl, nrep=nrep,
                           Time=Time, t1=t1, t2=t2, Err=E[ii], par=par2,
                           fit.routine=fit.routine)
    fit.uwerrfpi <- uwerr(f=fitf.online, data=W[ii,], S=S, pl=pl, nrep=nrep,
                      Time=Time, t1=t1, t2=t2, Err=E[ii], par=par2,
                      fit.routine=fit.routine)
  }
  if(method == "boot" || method == "all") {
    fit.boot <- boot(data=t(W[ii,]), statistic=fit.online.boot, R=boot.R, stype="i",
                     Time=Time, t1=t1, t2=t2, Err=E[ii], par=par2,
                     fit.routine=fit.routine)

    fit.tsboot <- tsboot(tseries=t(W[ii,]), statistic=fit.online.boot, R=boot.R, l=boot.l, sim=tsboot.sim,
                         Time=Time, t1=t1, t2=t2, Err=E[ii], par=par2, 
                         fit.routine=fit.routine)
  }

  Chi <- rep(0., times=2*T1)
  Fit <- rep(0., times=2*T1)

  jj <-  c(t1p1:t2p1)
  Fit[jj] <- pcacfit$par[1]^2*CExp(m=pcacfit$par[3], Time=2*Thalf, x=jj-1, sign=+1.)
  Fit[jj+T1] <- pcacfit$par[1]*pcacfit$par[2]*CExp(m=pcacfit$par[3], Time=2*Thalf, x=jj-1, sign=-1.)
  Chi[ii] <- (Fit[ii]-Cor[ii])/E[ii]

  dpaopp <- data.frame(t = array(0.,dim=c((t2-t1+1))), mass = array(0.,dim=c((t2-t1+1))),
                        dmass = array(0.,dim=c((t2-t1+1))), ddmass = array(0.,dim=c((t2-t1+1))),
                        tauint = array(0.,dim=c((t2-t1+1))), dtauint = array(0.,dim=c((t2-t1+1))))
  i <- 1
  for(t in t1p1:t2p1) {
    mass <- uwerrderived(pcacsym.online, data=W, S=S, pl=FALSE, t=t, T1=T1)
    dpaopp$t[i] <- t-1
    dpaopp$mass[i] <- mass$value[1]
    dpaopp$dmass[i] <- mass$dvalue[1]
    dpaopp$ddmass[i] <- mass$ddvalue[1]
    dpaopp$tauint[i] <- mass$tauint[1]
    dpaopp$dtauint[i] <- mass$dtauint[1]
    i=i+1
  }

  MChist.dpaopp <- rep(0., times=length(W[1,]))
  for(i in 1:length(W[1,])) {
    MChist.dpaopp[i] <- 0.
    for(t in t1p1:t2p1) {
      MChist.dpaopp[i] <- MChist.dpaopp[i] + pcacsym.online(data=W[,i], t=t, T1=T1)
    }
    MChist.dpaopp[i] <- MChist.dpaopp[i]/(t2p1-t1p1)
  }
  
  res <- list(fitresult=pcacfit, fitresultpp=massfit, t1=t1, t2=t2, N=length(W[1,]), Time=Time,
              fitdata=data.frame(t=(jj-1), Fit=Fit[ii], Cor=Cor[ii], Err=E[ii], Chi=Chi[ii]),
              uwerrresultmps=fit.uwerrm, uwerrresultmpcac=fit.uwerrpcac, uwerrresultfps=fit.uwerrfpi, 
              boot=fit.boot, tsboot=fit.tsboot, method=method, skip=skip,
              effmass=mass.eff, fit.routine=fit.routine, dpaopp=dpaopp, MChist.dpaopp=MChist.dpaopp,
              iobs=iobs, mu=mu, kappa=kappa,
              nrep=nrep, matrix.size=2)
  attr(res, "class") <- c("ofit", "list")  
  return(invisible(res))
}

fitmpcac.online <- function(Cor, Err, t1, t2, Time, par=c(1.,0.1,0.12),
                            N=2, no.masses=1, no=1, kappa, mu,
                            fit.routine="gsl") {
  Thalf <- Time/2
  T1 <- Thalf+1
  t1p1 <- (t1+1)
  t2p1 <- (t2+1)
  tr <- (t2-t1+1)

  if(no.masses == 1) {
    if(fit.routine != "gsl75") {
      fit <- optim(par, ChiSqr.pcac, method="BFGS", Thalf=Thalf,
                   x=c((t1):(t2)), y=Cor, err=Err, tr=tr)
    }
    else {
      fit <- gsl_fit_correlator_matrix(par, Thalf=Thalf,
                                       x=c((t1):(t2)), y=Cor, err=Err, tr = tr, N=N)
    }

    sort.ind <- c(1)
  }
  return(0.5*abs(fit$par[3])*fit$par[2]/fit$par[1])
}

fitmass.online <- function(Cor, Err, t1, t2, Time, par=c(1.,0.1,0.12),
                            N=2, no.masses=1, no=1, kappa, mu,
                            fit.routine="gsl") {
  Thalf <- Time/2
  T1 <- Thalf+1
  t1p1 <- (t1+1)
  t2p1 <- (t2+1)
  tr <- (t2-t1+1)

  if(no.masses == 1) {
    if(fit.routine != "gsl75") {
      fit <- optim(par, ChiSqr.pcac, method="BFGS", Thalf=Thalf,
                   x=c((t1):(t2)), y=Cor, err=Err, tr=tr)
    }
    else {
      fit <- gsl_fit_correlator_matrix(par, Thalf=Thalf,
                                       x=c((t1):(t2)), y=Cor, err=Err, tr = tr, N=N)
    }

    sort.ind <- c(1)
  }
  return(abs(fit$par[3]))
}

fitf.online <- function(Cor, Err, t1, t2, Time, par=c(1.,0.1,0.12),
                        N=2, no.masses=1, no=1, kappa, mu,
                        fit.routine="gsl") {
  Thalf <- Time/2
  T1 <- Thalf+1
  t1p1 <- (t1+1)
  t2p1 <- (t2+1)
  tr <- (t2-t1+1)

  if(no.masses == 1) {
    if(fit.routine != "gsl75") {
      fit <- optim(par, ChiSqr.pcac, method="BFGS", Thalf=Thalf,
                   x=c((t1):(t2)), y=Cor, err=Err, tr=tr)
    }
    else {
      fit <- gsl_fit_correlator_matrix(par, Thalf=Thalf,
                                       x=c((t1):(t2)), y=Cor, err=Err, tr = tr, N=N)
    }

    sort.ind <- c(1)
  }
  return(abs(fit$par[1])/sqrt(abs(fit$par[3]^3)))
}


pcacsym.online <- function(data, t, T1) {
  # PA[t+1]-PA[t-1]/4PP[t]
  return((data[t+1+T1]-data[t-1+T1])/(4.*data[t]))
}

fit.online.boot <- function(Z, d, Err, t1, t2, Time, par=c(1.,0.1,0.12),
                          kappa, mu, fit.routine="gsl") {
  Thalf <- Time/2
  T1 <- Thalf+1
  t1p1 <- (t1+1)
  t2p1 <- (t2+1)
  tr <- (t2-t1+1)
  Cor <- rep(0., times=length(Z[1,]))
  if(!missing(d)) {
    for(i in 1:length(Z[1,])) {
      Cor[i] = mean(Z[d,(i)])
    }
  }
  else {
    for(i in 1:length(Z[1,])) {
      Cor[i] = mean(Z[,(i)])
    }
  }

  if(fit.routine != "gsl75") {
    fit <- optim(par, ChiSqr.pcac, method="BFGS", Thalf=Thalf,
                 x=c((t1):(t2)), y=Cor, err=Err, tr=tr)
  }
  else {
    fit <- gsl_fit_correlator_matrix(par, Thalf=Thalf,
                                     x=c((t1):(t2)), y=Cor, err=Err, tr = tr)
  }
  #fit.fpi <- 2*kappa*2*mu/sqrt(2)*abs(fit$par[(sort.ind[1]-1)*(N+1)+1])/sqrt(abs(fit$par[sort.ind[1]*(N+1)])^3)
  fit.pcac <- abs(fit$par[3])*fit$par[2]/fit$par[1]/2.
  #fit.zv <- 2.*mu/abs(fit$par[sort.ind[1]*(N+1)])*fit$par[(sort.ind[1]-1)*(N+1)+1]/fit$par[(sort.ind[1]-1)*(N+1)+5]
  return(c(abs(fit$par[3]), fit.pcac, fit$par[1:2]))
}
