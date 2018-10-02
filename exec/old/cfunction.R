cfunction <- function(data, t1, t2, S=1.5, pl=FALSE, skip=0, cformat="cmi",
                      itype=1, iobs=1, type="cosh", ind.vec=c(1,3,4,5),
                      boot.R=99, boot.l=10, tsboot.sim="geom",
                      method="uwerr", fit.routine="optim", nrep) {
  
  if(missing(data)) {
    stop("Error! Data is missing!")
  }
  if(missing(t1) || missing(t2)) {
    stop("Error! t1 and t2 must be specified!")
  }
  par <- numeric(2)

  if(!any(data$V1 == itype & data$V2 == iobs)) {
    stop("The particular correlation function is missing!")
  }
  data <- data[(data$V1==itype & data$V2==iobs),]
  sign <- +1.
  if(type == "sinh") sign <- -1.
  
  Time <-  2*max(data[,ind.vec[2]])
#  Time <- Time + 1
  Thalf <- max(data[,ind.vec[2]])
  T1 <- Thalf+1
  t1p1 <- (t1+1)
  t2p1 <- (t2+1)
#  nrObs <- max(data[,ind.vec[1]])
  nrObs <- 1
  nrType <- 1
  Skip <- (skip*(T1)*nrType*nrObs+1)
  Length <- length(data[,ind.vec[3]])
  if(missing(nrep)) {
    nrep <- c(length(data[((Skip):Length),ind.vec[3]])/(nrObs*(T1)*nrType))
  }
  else {
    skip <- 0
    if(sum(nrep) != length(data[((Skip):Length),ind.vec[3]])/(nrObs*(T1)*nrType)) {
      stop("sum of replica differs from total no of measurements!")
    }
  }

  Z <- array(data[((Skip):Length),ind.vec[3]], 
             dim=c(nrObs*(T1)*nrType,(length(data[((Skip):Length),ind.vec[3]])/(nrObs*(T1)*nrType))))
  # negative times
  W <- array(data[((Skip):Length),ind.vec[4]], 
             dim=c(nrObs*(T1)*nrType,(length(data[((Skip):Length),ind.vec[4]])/(nrObs*(T1)*nrType))))

  rm(data)
  W <- getCor(T1=T1, W=W, Z=Z, type=type)
  rm(Z)
                                        #  options(show.error.messages = FALSE)
  eff <- effectivemass(from=(t1+1), to=(t2+1), Time, W[1:T1,] , pl=FALSE, S=1.5, nrep=nrep)
  options(show.error.messages = TRUE)

  mass.eff <- data.frame(t=eff$t, m=eff$mass, dm=eff$dmass)

  Cor <- rep(0., times=T1)
  E <- rep(0., times=T1)
  
  for(i in 1:(T1)) {
    Cor[i] <- mean(W[(i),])
    tmpe <- try(uwerrprimary(W[(i),], pl=F, nrep=nrep)$dvalue, TRUE)
    if(!inherits(tmpe, "try-error")) E[i] = tmpe
    else {
      warning("error of correlator replaced by naive estimate!\n", call.=F)
      E[i] = sd(W[(i),])/sqrt(length(W[(i),]))
    }
  }
  par[2] <- eff$mass[1]
  par[1] <- Cor[(t1+1)]*exp(par[2]*(t1+1))
  
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
  fit.mass <- abs(massfit$par[2])
  if(fit.routine != "gsl" && massfit$convergence!=0) {
    warning("optim did not converge for massfit! ", massfit$convergence)
  }

  fit.dof <- (t2-t1+1)*3-length(massfit$par)
  fit.chisqr <- massfit$value

  fit.uwerrm <- NULL
  fit.boot <- NULL
  fit.tsboot <- NULL
  if(method == "uwerr" || method == "all") {
    fit.uwerrm <- uwerr(f=fitmass, data=t(W[ii,]), S=S, pl=pl, nrep=nrep,
                        Time=Time, t1=t1, t2=t2, Err=E[ii], par=par, sign=sign,
                        fit.routine=fit.routine)
  }
  if(method == "boot" || method == "all") {
    fit.boot <- boot(data=t(W[ii,]), statistic=getfit.boot, R=boot.R, stype="i",
                     Time=Time, t1=t1, t2=t2, Err=E[ii], par=par, sign=sign,
                     fit.routine=fit.routine)

    fit.tsboot <- tsboot(tseries=t(W[ii,]), statistic=getfit.boot, R=boot.R, l=boot.l, sim=tsboot.sim,
                         Time=Time, t1=t1, t2=t2, Err=E[ii], par=par, sign=sign,
                         fit.routine=fit.routine)
  }

  
  Chi <- rep(0., times=T1)
  Fit <- rep(0., times=T1)

  jj <-  c(t1p1:t2p1)
  Fit[jj] <- massfit$par[1]*CExp(m=fit.mass[1], Time=2*Thalf, x=jj-1)
  
  Chi[ii] <- (Fit[ii]-Cor[ii])/E[ii]
  
  res <- list(fitresult=massfit, t1=t1, t2=t2, N=length(W[1,]), Time=Time,
              fitdata=data.frame(t=(jj-1), Fit=Fit[ii], Cor=Cor[ii], Err=E[ii], Chi=Chi[ii]),
              uwerrresultmps=fit.uwerrm,
              boot=fit.boot, tsboot=fit.tsboot, method=method,
              effmass=mass.eff, fit.routine=fit.routine,
              itype=itype, iobs=iobs, sign=sign,
              nrep=nrep, matrix.size=1)
  attr(res, "class") <- c("cfit", "list")  
  return(invisible(res))
}

getCor <- function(T1, W, Z, type=c("cosh")) {

  # iobs enumerating the gamma matrix combination
  # ityp enumeratiog the smearing level
  N <- length(type)
  sign = rep(+1., times=N)
  for(i in 1:N) {
    if(type[i]=="sinh") {
      sign[i] = -1.
    }
  }
  
  for(j in 1:N) {
    for(i in 1:(T1)) {
      two <- 2.
      if(i==1 || i==(T1)) {
                                        # Take care of zeros in the correlators when summing t and T-t+1
        two <- 1.
      }
      
      W[(i+(j-1)*T1),] <- (W[(i+(j-1)*T1),]
                + sign[j]*Z[(i+(j-1)*T1),])/two
    }
  }
  return(invisible(W))
}

fitmass <- function(Cor, Err, t1, t2, Time, par=c(1.,0.12), sign,
                    fit.routine="optim") {
  Thalf <- Time/2
  T1 <- Thalf+1
  t1p1 <- (t1+1)
  t2p1 <- (t2+1)
  tr <- (t2-t1+1)

  if(fit.routine != "gsl") {
    fit <- optim(par, ChiSqr.singleCor, method="BFGS", Thalf=Thalf,
                 x=c((t1):(t2)), y=Cor, err=Err, tr=tr, sign=sign)
  }
  else {
    fit <- gsl_fit_correlator_matrix(par, Thalf=Thalf,
                                     x=c((t1):(t2)), y=Cor, err=Err, tr = tr, sign=sign)
  }
  
  return(abs(fit$par[2]))
}

getfit.boot <- function(Z, d, Err, t1, t2, Time, par=c(1.,0.12),
                     fit.routine="optim", sign) {
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

  if(fit.routine != "gsl") {
    fit <- optim(par, ChiSqr.singleCor, method="BFGS", Thalf=Thalf,
                 x=c((t1):(t2)), y=Cor, err=Err, tr=tr, sign=sign)
  }
  else {
    fit <- gsl_fit_correlator_matrix(par, Thalf=Thalf,
                                     x=c((t1):(t2)), y=Cor, err=Err, tr = tr, sign=sign)
  }
  sort.ind <- c(1)
  return(c(abs(fit$par[2]), fit$par[1],
           fit$value))
}


