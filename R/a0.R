a0 <- function(cmicor, mu=0.1, kappa=0.156, t1, t2, S=1.5, pl=FALSE, skip=0,
               variational=list(ta=3, tb=4, N=2), ind.vec=c(1,3,4,5),
               no.masses=1, boot.R=99, boot.l=10, tsboot.sim="geom",
               method="uwerr", mass.guess, par.guess, nrep, fit.routine="gsl") {

  matrix.size=2
  if(missing(cmicor)) {
    stop("Error! Data is missing!")
  }
  if(missing(t1) || missing(t2)) {
    stop("Error! t1 and t2 must be specified!")
  }
  if(missing(mass.guess)) {
    mass.guess <- c(0.8, 1., 3.)
  }
  else {
    if(length(mass.guess) < no.masses) {
      stop("mass.guess has not the correct length!")
    }
  }
  if(missing(par.guess)) {
    par.guess <- c(1.,0.8, 0.1,0.1, 0.1,0.1,0.1) 
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
  nrObs <- max(cmicor[,ind.vec[1]])
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

  Z <- array(cmicor[((Skip):Length),ind.vec[3]], 
             dim=c(nrObs*(T1)*4,(length(cmicor[((Skip):Length),ind.vec[3]])/(nrObs*(T1)*4))))
  # negative times
  W <- array(cmicor[((Skip):Length),ind.vec[4]], 
             dim=c(nrObs*(T1)*4,(length(cmicor[((Skip):Length),ind.vec[4]])/(nrObs*(T1)*4))))

  rm(cmicor)
  W <- arrangeCor.a0(T1=T1, W=W, Z=Z)
  rm(Z)

  options(show.error.messages = FALSE)
  a0.eff.ll <- effectivemass(from=(t1+1), to=(t2+1), Time, W[1:T1,] , pl=FALSE, S=1.5, nrep=nrep)
  a0.eff.lf <- effectivemass(from=(t1+1), to=(t2+1), Time, W[(T1+1):(2*T1),] , pl=FALSE, S=1.5, nrep=nrep)
  a0.eff.ff <- effectivemass(from=(t1+1), to=(t2+1), Time, W[(3*T1+1):(4*T1),] , pl=FALSE, S=1.5, nrep=nrep)
  options(show.error.messages = TRUE)
  
  a0.eff <- data.frame(t=a0.eff.ll$t, mll=a0.eff.ll$mass, dmll=a0.eff.ll$dmass,
                         mlf=a0.eff.lf$mass, dmlf=a0.eff.lf$dmass,
                         mff=a0.eff.ff$mass, dmff=a0.eff.ff$dmass)
  
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

  N <- max(matrix.size,variational$N)
  ta <- (variational$ta+1)
  tb <- (variational$tb+1)
  C1 <- getNxNmatrix(Cor=Cor, t=ta, T1=T1, N=N)
  C2 <- getNxNmatrix(Cor=Cor, t=tb, T1=T1, N=N)
  C3 <- getNxNmatrix(Cor=Cor, t=tb, T1=T1, N=N)

  # first index is rows, second columns
  for(i in 1:N) {
    C3[,i] <- solve(C1, C2[,i])
  }
  variational.solve <- eigen(C3, symmetric=FALSE, only.values = FALSE, EISPACK=FALSE)
  # get the left eigenvectors, the eigenvectors have unit length
  for(i in 1:N) {
    if(abs(variational.solve$values[i]) > 0.95/(tb-ta)) {
      variational.solve$values[i] <- 0.0001
    }
  }
  # this does not quite work for the pion ?
  variational.sortindex <- order(-log(abs(variational.solve$values)*(tb-ta)))
  if(TRUE) {

    left.vectors <- array(0., dim=c(N,N))
    left.vectors <- crossprod(C1, variational.solve$vectors)
    X <- crossprod(left.vectors,variational.solve$vectors)
    for(i in 1:no.masses) {
      j <-  variational.sortindex[i]
      left.vectors[,j] <- left.vectors[,i]/sqrt(X[j,j])
      variational.solve$vectors[,j] <- variational.solve$vectors[,j]/sqrt(X[j,j])
    }
    
    
    par <- c(2*left.vectors[(1:matrix.size),1],
             -log(abs(variational.solve$values[variational.sortindex[1]]))/(tb-ta))
    if(no.masses > 1) {
      for(i in 2:(no.masses)) {
        par <- c(par,
                 left.vectors[(1:matrix.size),i],
                 -log(abs(variational.solve$values[variational.sortindex[i]]))/(tb-ta)) 
      }
    }
#    print(par)
  }
  
  variational.masses <-  -log(abs(variational.solve$values[variational.sortindex]))/(tb-ta)
  rm(C1, C2, C3, ta, tb, N)

  # Index vector of data to be used in the analysis
  ii <- c((t1p1):(t2p1), (t1p1+T1):(t2p1+T1), (t1p1+3*T1):(t2p1+3*T1))

  #BFGS
  if(fit.routine != "gsl") {
    if(no.masses == 1) {
      a0fit <- optim(par, ChiSqr.1mass, method="BFGS", control=list(trace=0),Thalf=Thalf,
                     x=c((t1):(t2)), y=Cor[ii], err=E[ii], tr = (t2-t1+1), N=matrix.size)
      fit.mass <- abs(a0fit$par[matrix.size+1])
    }
    else if(no.masses == 2) {
      a0fit <- optim(par, ChiSqr.2mass, method="BFGS", control=list(trace=0),Thalf=Thalf,
                     x=c((t1):(t2)), y=Cor[ii], err=E[ii], tr = (t2-t1+1), N=matrix.size)
      fit.mass <- sort(abs(a0fit$par[c((matrix.size+1),(2*matrix.size+2))]))
    }
    else if(no.masses > 2) {
      a0fit <- optim(par, ChiSqr.3mass, method="BFGS", control=list(trace=0),Thalf=Thalf,
                     x=c((t1):(t2)), y=Cor[ii], err=E[ii], tr = (t2-t1+1), N=matrix.size)
      fit.mass <- sort(abs(a0fit$par[c((matrix.size+1),(2*matrix.size+2),(3*matrix.size+3))]))
    }
    if(a0fit$convergence != 0) {
      warning("optim did not converge for a0fit!", call.=F)
    }
  }
  else {
    a0fit <- gsl_fit_correlator_matrix(par, Thalf=Thalf, x=c((t1):(t2)), y=Cor[ii],
                                       err=E[ii], tr = (t2-t1+1), N=matrix.size, no_masses=no.masses,
                                       prec=c(1.e-3,1.e-3))

    if(no.masses == 1) fit.mass <- abs(a0fit$par[matrix.size+1])
    if(no.masses == 2) fit.mass <- sort(abs(a0fit$par[c((matrix.size+1),(2*matrix.size+2))]))
    if(no.masses > 2) fit.mass <- sort(abs(a0fit$par[c((matrix.size+1),(2*matrix.size+2),(3*matrix.size+3))]))
    if(a0fit$convergence < 0) {
      warning("gsl multifit did not converge for a0fit!", call.=F)
    }

  }

  fit.dof <- (t2-t1+1)*3-length(a0fit$par)
  fit.chisqr <- a0fit$value

  if(pl) {
    plot.effmass(m=fit.mass, ll=a0.eff.ll, lf=a0.eff.lf, ff=a0.eff.ff)
  }

  fit.uwerrm <- NULL
  fit.uwerrf <- NULL
  fit.uwerrpcac <- NULL
  fit.uwerrm2 <- NULL
  fit.uwerrm3 <- NULL
  fit.boot <- NULL
  fit.tsboot <- NULL
  if(method == "uwerr" || method == "all") {
    fit.uwerrm <- uwerr(f=fitmasses.a0, data=W[ii,], S=S, pl=pl, nrep=nrep,
                        Time=Time, t1=t1, t2=t2, Err=E[ii], par=par, N=matrix.size,
                        no.masses=no.masses, fit.routine=fit.routine)

    if(no.masses == 2) {
      fit.uwerrm2 <- uwerr(f=fitmasses.a0, data=W[ii,], S=S, pl=pl, nrep=nrep,
                           Time=Time, t1=t1, t2=t2, Err=E[ii], par=par, N=matrix.size,
                           no.masses=no.masses, no=2, fit.routine=fit.routine)
    }
    if(no.masses > 2) {
      fit.uwerrm3 <- uwerr(f=fitmasses.a0, data=W[ii,], S=S, pl=pl, nrep=nrep,
                           Time=Time, t1=t1, t2=t2, Err=E[ii], par=par, N=matrix.size,
                           no.masses=no.masses, no=3, fit.routine=fit.routine)
    }
  }
  if(method == "boot" || method == "all") {
    fit.boot <- boot(data=t(W[ii,]), statistic=fit.a0.boot, R=boot.R, stype="i",
                     Time=Time, t1=t1, t2=t2, Err=E[ii], par=par, N=matrix.size, no.masses=no.masses,
                     kappa=kappa, mu=mu, fit.routine=fit.routine)

    fit.tsboot <- tsboot(tseries=t(W[ii,]), statistic=fit.a0.boot, R=boot.R, l=boot.l, sim=tsboot.sim,
                         Time=Time, t1=t1, t2=t2, Err=E[ii], par=par, N=matrix.size, no.masses=no.masses,
                         kappa=kappa, mu=mu, fit.routine=fit.routine)
  }

  Chi <- rep(0., times=4*T1)
  Fit <- rep(0., times=4*T1)
  N <- matrix.size
  jj <-  c(t1p1:t2p1)
  for(i in 1:no.masses) {
    Fit[jj] <- Fit[jj] +
      a0fit$par[1+(i-1)*(N+1)]^2*CExp(m=a0fit$par[i*(N+1)], Time=2*Thalf, x=jj-1)
    Fit[jj+T1] <- Fit[jj+T1] +
      a0fit$par[1+(i-1)*(N+1)]*a0fit$par[2+(i-1)*(N+1)]*CExp(m=a0fit$par[i*(N+1)], Time=2*Thalf, x=jj-1)
    Fit[jj+2*T1] <- Fit[jj+2*T1] +
      a0fit$par[1+(i-1)*(N+1)]*a0fit$par[2+(i-1)*(N+1)]*CExp(m=a0fit$par[i*(N+1)], Time=2*Thalf, x=jj-1)
    Fit[jj+3*T1] <- Fit[jj+3*T1] +
      a0fit$par[2+(i-1)*(N+1)]^2*CExp(m=a0fit$par[i*(N+1)], Time=2*Thalf, x=jj-1)
  }

  Chi[ii] <- (Fit[ii]-Cor[ii])/E[ii]
  
  res <- list(fitresult=a0fit, t1=t1, t2=t2, N=length(W[1,]), Time=Time,
              fitdata=data.frame(t=(jj-1), Fit=Fit[ii], Cor=Cor[ii], Err=E[ii], Chi=Chi[ii]),
              uwerrresultma0=fit.uwerrm, uwerrresultma02=fit.uwerrm2, uwerrresultma03=fit.uwerrm3,
              boot=fit.boot, tsboot=fit.tsboot,
              effmass=a0.eff, kappa=kappa, mu=mu,
              variational.masses=variational.masses, no.masses=no.masses,
              matrix.size = matrix.size, nrep=nrep)
  attr(res, "class") <- c("a0fit", "cfit", "list")  
  return(invisible(res))
}

arrangeCor.a0 <- function(T1, W, Z) {
  j <- 0
  for(i in (18*4*T1+1):(18*4*T1+T1)) {
    j <- j+1
    two <- 2.
    if(j==1 || j==(T1)) {
      # Take care of zeros in the correlators when summing t and T-t+1
      two <- 1.
    }
    # 11 for LL, (LF + FL)/2, FF -> symmetric cosh
    W[j,] <- -(W[i,] + Z[i,])/two
    W[(j+T1),] <- -(W[(i+T1),] + W[(i+2*T1),] + Z[(i+T1),] + Z[(i+2*T1),])/2./two
    W[(j+2*T1),] <- W[(j+T1),]
    W[(j+3*T1),] <- -(W[(i+3*T1),] + Z[(i+3*T1),])/two
  }
  return(invisible(W))
}

fitmasses.a0 <- function(Cor, Err, t1, t2, Time, par=c(1.,0.1,0.12),
                         N=2, no.masses=1, no=1, kludge=FALSE,
                         fit.routine = "gsl") {
  Thalf <- Time/2
  T1 <- Thalf+1
  t1p1 <- (t1+1)
  t2p1 <- (t2+1)
  tr <- (t2-t1+1)
  if(fit.routine != "gsl") {
    if(no.masses == 1) {
      fit <- optim(par, ChiSqr.1mass, method="BFGS", Thalf=Thalf,
                   x=c((t1):(t2)), y=Cor, err=Err, tr=tr, N=N)
      return(abs(fit$par[N+1]))
    }
    else if (no.masses == 2) {
      fit <- optim(par, ChiSqr.2mass, method="BFGS", Thalf=Thalf,
                   x=c((t1):(t2)), y=Cor, err=Err, tr=tr, N=N, kludge=kludge)
      return(sort(abs(fit$par[c((N+1),(2*N+2))]))[no])
    }
    else if (no.masses == 3) {
      fit <- optim(par, ChiSqr.3mass, method="BFGS", Thalf=Thalf,
                   x=c((t1):(t2)), y=Cor, err=Err, tr=tr, N=N, kludge=kludge)
      return(sort(abs(fit$par[c((N+1),(2*N+2),(3*N+3))]))[no])
    }
  }
  else {
    fit <- gsl_fit_correlator_matrix(par, Thalf=Thalf, x=c((t1):(t2)), y=Cor,
                                     err=Err, tr = tr, N=N, no_masses = no.masses,
                                     prec=c(1.e-3,1.e-3))
    if(no.masses == 1) return(abs(fit$par[N+1]))
    if(no.masses == 2) return(sort(abs(fit$par[c((N+1),(2*N+2))]))[no])
    if(no.masses == 3) return(sort(abs(fit$par[c((N+1),(2*N+2),(3*N+3))]))[no])
  }
}

fit.a0.boot <- function(Z, d, Err, t1, t2, Time, par=c(1.,0.1,0.12),
                        N=2, no.masses=1, kludge=FALSE, kappa, mu,
                        fit.routine = "gsl") {
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

  if(no.masses == 1) {
    if(fit.routine != "gsl") {
      fit <- optim(par, ChiSqr.1mass, method="BFGS", Thalf=Thalf,
                   x=c((t1):(t2)), y=Cor, err=Err, tr=tr, N=N)
    }
    else {
      fit <- gsl_fit_correlator_matrix(par, Thalf=Thalf,
                                       x=c((t1):(t2)), y=Cor, err=Err, tr = tr, N=N,
                                       prec=c(1.e-3,1.e-3))
    }
    sort.ind <- c(1)
    return(c(abs(fit$par[N+1]), fit$par[c(1:N)],
             fit$value))
  }
  else if (no.masses == 2) {
    if(fit.routine != "gsl") {
      fit <- optim(par, ChiSqr.2mass, method="BFGS", Thalf=Thalf,
                   x=c((t1):(t2)), y=Cor, err=Err, tr=tr, N=N, kludge=kludge)
    }
    else {
      fit <- gsl_fit_correlator_matrix(par, Thalf=Thalf, x=c((t1):(t2)), y=Cor,
                                       err=Err, tr = tr, N=N, no_masses = 2,
                                       prec=c(1.e-3,1.e-3))
    }
    sort.ind <- order(fit$par[c((N+1),(2*N+2))])
    return(c(abs(fit$par[sort.ind[1]*(N+1)]),
             fit$par[c(((sort.ind[1]-1)*(N+1)+1):((sort.ind[1])*(N+1)-1))],
             abs(fit$par[sort.ind[2]*(N+1)]),
             fit$par[c(((sort.ind[2]-1)*(N+1)+1):((sort.ind[2])*(N+1)-1))],
             fit$value))
  }
  else if (no.masses == 3) {
    if(fit.routine != "gsl") {
      fit <- optim(par, ChiSqr.3mass, method="BFGS", Thalf=Thalf,
                   x=c((t1):(t2)), y=Cor, err=Err, tr=tr, N=N, kludge=kludge)
    }
    else {
      fit <- gsl_fit_correlator_matrix(par, Thalf=Thalf, x=c((t1):(t2)), y=Cor,
                                       err=Err, tr = tr, N=N, no_masses = 3,
                                       prec=c(1.e-3,1.e-3))
    }
    sort.ind <- order(fit$par[c((N+1),(2*N+2),(3*N+3))])
    return(c(abs(fit$par[sort.ind[1]*(N+1)]),
             fit$par[c(((sort.ind[1]-1)*(N+1)+1):((sort.ind[1])*(N+1)-1))],
             abs(fit$par[sort.ind[2]*(N+1)]),
             fit$par[c(((sort.ind[2]-1)*(N+1)+1):((sort.ind[2])*(N+1)-1))],
             abs(fit$par[sort.ind[3]*(N+1)]),
             fit$par[c(((sort.ind[3]-1)*(N+1)+1):((sort.ind[3])*(N+1)-1))],
             fit$value))
  }
}


