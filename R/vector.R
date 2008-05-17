rho <- function(cmicor, mu=0.1, kappa=0.156, t1, t2, S=1.5, pl=FALSE, skip=0,
                variational=list(ta=4, tb=5, N=6), ind.vec=c(1,3,4,5),
                no.masses=1, matrix.size=2, boot.R=99, boot.l=10, tsboot.sim="geom",
                method="uwerr", nrep) {
  
  if(missing(cmicor)) {
    stop("Error! Data is missing!")
  }
  if(missing(t1) || missing(t2)) {
    stop("Error! t1 and t2 must be specified!")
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
  W <- arrangeCor.vector(T1=T1, W=W, Z=Z)
  rm(Z)

#  options(show.error.messages = FALSE)
  rho.eff.ll <- effectivemass(from=(t1+1), to=(t2+1), Time, W[1:T1,] , pl=FALSE, S=1.5, nrep=nrep)
  rho.eff.lf <- effectivemass(from=(t1+1), to=(t2+1), Time, W[(T1+1):(2*T1),] , pl=FALSE, S=1.5, nrep=nrep)
  rho.eff.ff <- effectivemass(from=(t1+1), to=(t2+1), Time, W[(3*T1+1):(4*T1),] , pl=FALSE, S=1.5, nrep=nrep)
  options(show.error.messages = TRUE)
  
  rho.eff <- data.frame(t=rho.eff.ll$t, mll=rho.eff.ll$mass, dmll=rho.eff.ll$dmass,
                         mlf=rho.eff.lf$mass, dmlf=rho.eff.lf$dmass,
                         mff=rho.eff.ff$mass, dmff=rho.eff.ff$dmass)
  
  Cor <- rep(0., times=9*4*T1)
  E <- rep(0., times=9*4*T1)
  
  for(i in 1:(9*4*T1)) {
    Cor[i] = mean(W[(i),])
    E[i] = uwerrprimary(W[(i),], pl=F, nrep=nrep)$dvalue
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
  variational.sortindex <- order(-log(abs(variational.solve$values)*(tb-ta)))
  left.vectors <- array(0., dim=c(N,N))
  left.vectors <- crossprod(C1, variational.solve$vectors)
  X <- crossprod(left.vectors,variational.solve$vectors)
  for(i in 1:N) {
    left.vectors[,i] <- left.vectors[,i]/sqrt(X[i,i])
    variational.solve$vectors[,i] <- variational.solve$vectors[,i]/sqrt(X[i,i])
  }

  par <- c(2*left.vectors[(1:matrix.size),1],
           -log(abs(variational.solve$values[variational.sortindex[1]]))/(tb-ta))
  if(no.masses > 1) {
    for(i in 2:(no.masses)) {
      par <- c(par,
               2.*left.vectors[(1:matrix.size),i],
               -log(abs(variational.solve$values[variational.sortindex[2]]))/(tb-ta)) 
    }
  }
                                        #  print(par)
  
  variational.masses <-  -log(abs(variational.solve$values[variational.sortindex]))/(tb-ta)
  rm(C1, C2, C3, ta, tb, X, left.vectors)

  # Index vector of data to be used in the analysis
  ii <- c((t1p1):(t2p1), (t1p1+T1):(t2p1+T1), (t1p1+3*T1):(t2p1+3*T1))
  if(matrix.size > 2) {
    ii <- c(ii, (t1p1+20*T1):(t2p1+20*T1), (t1p1+21*T1):(t2p1+21*T1),
            (t1p1+22*T1):(t2p1+22*T1), (t1p1+23*T1):(t2p1+23*T1),
            (t1p1+16*T1):(t2p1+16*T1), (t1p1+17*T1):(t2p1+17*T1),
            (t1p1+19*T1):(t2p1+19*T1))
  }
  if(matrix.size > 4) {
    ii <- c(ii, (t1p1+12*T1):(t2p1+12*T1), (t1p1+13*T1):(t2p1+13*T1),
            (t1p1+15*T1):(t2p1+15*T1),
            (t1p1+4*T1):(t2p1+4*T1), (t1p1+5*T1):(t2p1+5*T1),
            (t1p1+6*T1):(t2p1+6*T1), (t1p1+7*T1):(t2p1+7*T1),
            (t1p1+28*T1):(t2p1+28*T1), (t1p1+29*T1):(t2p1+29*T1),
            (t1p1+30*T1):(t2p1+30*T1), (t1p1+31*T1):(t2p1+31*T1))
  }
  #BFGS
  if(no.masses == 1) {
#    rhofit <- optim(par, ChiSqr.1mass, method="BFGS", control=list(trace=0),Thalf=Thalf,
#                    x=c((t1):(t2)), y=Cor[ii], err=E[ii], tr = (t2-t1+1), N=matrix.size)
    rhofit <- gsl_fit_correlator_matrix(par, Thalf=Thalf,
                    x=c((t1):(t2)), y=Cor[ii], err=E[ii], tr = (t2-t1+1), N=matrix.size)
    fit.mass <- abs(rhofit$par[matrix.size+1])
  }
  else if(no.masses == 2) {
#    rhofit <- optim(par, ChiSqr.2mass, method="BFGS", control=list(trace=0),Thalf=Thalf,
#                    x=c((t1):(t2)), y=Cor[ii], err=E[ii], tr = (t2-t1+1), N=matrix.size)
    rhofit <- gsl_fit_correlator_matrix(par, Thalf=Thalf,
                    x=c((t1):(t2)), y=Cor[ii], err=E[ii], tr = (t2-t1+1), N=matrix.size, no_masses = no.masses)

    fit.mass <- sort(abs(rhofit$par[c((matrix.size+1),(2*matrix.size+2))]))
  }
  else if(no.masses > 2) {
    rhofit <- optim(par, ChiSqr.3mass, method="BFGS", control=list(trace=0),Thalf=Thalf,
                    x=c((t1):(t2)), y=Cor[ii], err=E[ii], tr = (t2-t1+1), N=matrix.size)
    fit.mass <- sort(abs(rhofit$par[c((matrix.size+1),(2*matrix.size+2),(3*matrix.size+3))]))
  }
  if(rhofit$convergence < 0) {
    warning("optim did not converge for rhofit!", call.=F)
  }
                                        #  print(rhofit)

  fit.dof <- (t2-t1+1)*3-length(rhofit$par)
  fit.chisqr <- rhofit$value

  if(pl) {
    plot.effmass(m=fit.mass, ll=rho.eff.ll, lf=rho.eff.lf, ff=rho.eff.ff)
  }

  fit.uwerrm <- NULL
  fit.uwerrm2 <- NULL
  fit.uwerrm3 <- NULL
  fit.boot <- NULL
  fit.boot.ci <- NULL
  fit.tsboot <- NULL
  fit.tsboot.ci <- NULL
  if(method == "uwerr" || method == "all") {
    fit.uwerrm <- uwerr(f=fitmasses.vector, data=W[ii,], S=S, pl=pl, nrep=nrep,
                        Time=Time, t1=t1, t2=t2, Err=E[ii], par=par, N=matrix.size, no.masses=no.masses)

    if(no.masses == 2) {
      fit.uwerrm2 <- uwerr(f=fitmasses.vector, data=W[ii,], S=S, pl=pl, nrep=nrep,
                           Time=Time, t1=t1, t2=t2, Err=E[ii], par=par, N=matrix.size, no.masses=no.masses, no=2)
    }
    if(no.masses > 2) {
      fit.uwerrm3 <- uwerr(f=fitmasses.vector, data=W[ii,], S=S, pl=pl, nrep=nrep,
                           Time=Time, t1=t1, t2=t2, Err=E[ii], par=par, N=matrix.size, no.masses=no.masses, no=3)
    }
  }
  if(method == "boot" || method == "all") {
    fit.boot <- boot(data=t(W[ii,]), statistic=fitmasses.vector.boot, R=boot.R, stype="i",
                     Time=Time, t1=t1, t2=t2, Err=E[ii], par=par, N=matrix.size, no.masses=no.masses)

    fit.tsboot <- tsboot(tseries=t(W[ii,]), statistic=fitmasses.vector.boot, R=boot.R, l=boot.l, sim=tsboot.sim,
                         Time=Time, t1=t1, t2=t2, Err=E[ii], par=par, N=matrix.size, no.masses=no.masses)
  }

  
  Chi <- rep(0., times=9*4*T1)
  Fit <- rep(0., times=9*4*T1)

  jj <-  c(t1p1:t2p1)
  Fit[jj] <- rhofit$par[1]^2*CExp(m=fit.mass, Time=2*Thalf, x=jj-1)
  Fit[jj+T1] <- rhofit$par[1]*rhofit$par[2]*CExp(m=fit.mass, Time=2*Thalf, x=jj-1)
  Fit[jj+2*T1] <- rhofit$par[1]*rhofit$par[2]*CExp(m=fit.mass, Time=2*Thalf, x=jj-1)
  Fit[jj+3*T1] <- rhofit$par[2]*rhofit$par[2]*CExp(m=fit.mass, Time=2*Thalf, x=jj-1)
  if(matrix.size > 2) {
    Fit[jj+20*T1] <- rhofit$par[1]*rhofit$par[3]*CExp(m=fit.mass, Time=2*Thalf, x=jj-1, sign=-1.)
    Fit[jj+21*T1] <- rhofit$par[1]*rhofit$par[4]*CExp(m=fit.mass, Time=2*Thalf, x=jj-1, sign=-1.)
    Fit[jj+22*T1] <- rhofit$par[2]*rhofit$par[3]*CExp(m=fit.mass, Time=2*Thalf, x=jj-1, sign=-1.)
    Fit[jj+23*T1] <- rhofit$par[2]*rhofit$par[4]*CExp(m=fit.mass, Time=2*Thalf, x=jj-1, sign=-1.)

    Fit[jj+16*T1] <- rhofit$par[3]*rhofit$par[3]*CExp(m=fit.mass, Time=2*Thalf, x=jj-1)
    Fit[jj+17*T1] <- rhofit$par[3]*rhofit$par[4]*CExp(m=fit.mass, Time=2*Thalf, x=jj-1)
    Fit[jj+18*T1] <- rhofit$par[3]*rhofit$par[4]*CExp(m=fit.mass, Time=2*Thalf, x=jj-1)
    Fit[jj+19*T1] <- rhofit$par[4]*rhofit$par[4]*CExp(m=fit.mass, Time=2*Thalf, x=jj-1)
  }
  if(matrix.size > 4) {
    Fit[jj+12*T1] <- rhofit$par[5]*rhofit$par[5]*CExp(m=fit.mass, Time=2*Thalf, x=jj-1)
    Fit[jj+13*T1] <- rhofit$par[5]*rhofit$par[6]*CExp(m=fit.mass, Time=2*Thalf, x=jj-1)
    Fit[jj+14*T1] <- rhofit$par[5]*rhofit$par[6]*CExp(m=fit.mass, Time=2*Thalf, x=jj-1)
    Fit[jj+15*T1] <- rhofit$par[6]*rhofit$par[6]*CExp(m=fit.mass, Time=2*Thalf, x=jj-1)

    Fit[jj+4*T1] <- rhofit$par[1]*rhofit$par[5]*CExp(m=fit.mass, Time=2*Thalf, x=jj-1, sign=-1.)
    Fit[jj+5*T1] <- rhofit$par[1]*rhofit$par[6]*CExp(m=fit.mass, Time=2*Thalf, x=jj-1, sign=-1.)
    Fit[jj+6*T1] <- rhofit$par[2]*rhofit$par[5]*CExp(m=fit.mass, Time=2*Thalf, x=jj-1, sign=-1.)
    Fit[jj+7*T1] <- rhofit$par[2]*rhofit$par[6]*CExp(m=fit.mass, Time=2*Thalf, x=jj-1, sign=-1.)
    
    Fit[jj+28*T1] <- rhofit$par[3]*rhofit$par[5]*CExp(m=fit.mass, Time=2*Thalf, x=jj-1)
    Fit[jj+29*T1] <- rhofit$par[3]*rhofit$par[6]*CExp(m=fit.mass, Time=2*Thalf, x=jj-1)
    Fit[jj+30*T1] <- rhofit$par[4]*rhofit$par[5]*CExp(m=fit.mass, Time=2*Thalf, x=jj-1)
    Fit[jj+31*T1] <- rhofit$par[4]*rhofit$par[6]*CExp(m=fit.mass, Time=2*Thalf, x=jj-1)
  }
  
  Chi[ii] <- (Fit[ii]-Cor[ii])/E[ii]
  
  res <- list(fitresult=rhofit, t1=t1, t2=t2, N=length(W[1,]), Time=Time,
              fitdata=data.frame(t=(jj-1), Fit=Fit[ii], Cor=Cor[ii], Err=E[ii], Chi=Chi[ii]),
              uwerrresultmv=fit.uwerrm, uwerrresultmv2=fit.uwerrm2, uwerrresultmv3=fit.uwerrm3,
              mv.boot=fit.boot, mv.tsboot=fit.tsboot,
              effmass=rho.eff, kappa=kappa, mu=mu, nrep=nrep,
              variational.masses=variational.masses, no.masses=no.masses,
              matrix.size = matrix.size)
  attr(res, "class") <- c("rhofit", "cfit", "list")  
  return(invisible(res))
}


arrangeCor.vector <- function(T1, W, Z) {
  # vector starts with the 10-th correlator in the files
  # order in the file is 44 4V V4 VV AA 4A A4 VA AV
  # with LL, LF, FL, FF for each
  # for rho and a1: 4=gig4 V=gi A=gig5
  #
  # we order it as 44, 4V, V4, VV, AA, 4A, A4, VA, AV
  #
  # W contains t=0 to T1/2, Z t=T1-1 to T/2
  j <- 0
  for(i in (9*4*T1+1):(9*4*T1+T1)) {
    j <- j+1
    two <- 2.
    if(j==1 || j==(T1)) {
      # Take care of zeros in the correlators when summing t and T-t+1
      two <- 1.
    }
    # gamma_i gamma_4 (44) for LL, (LF + FL)/2, FF -> symmetric cosh
    W[j,] <- (W[i,] + Z[i,])/two
    W[(j+T1),] <- (W[(i+T1),] + W[(i+2*T1),] + Z[(i+T1),] + Z[(i+2*T1),])/2./two
    W[(j+2*T1),] <- W[(j+T1),]
    W[(j+3*T1),] <- (W[(i+3*T1),] + Z[(i+3*T1),])/two
    # (4V + V4)/2 for LL, LF, FL, FF -> anti-symmetric sinh
    W[(j+4*T1),] <- (W[(i+12*T1),] - Z[(i+12*T1),] + W[(i+16*T1),] - Z[(i+16*T1),])/2./two
    W[(j+5*T1),] <- (W[(i+13*T1),] - Z[(i+13*T1),] + W[(i+18*T1),] - Z[(i+18*T1),])/2./two
    W[(j+6*T1),] <- (W[(i+14*T1),] - Z[(i+14*T1),] + W[(i+17*T1),] - Z[(i+17*T1),])/2./two
    W[(j+7*T1),] <- (W[(i+15*T1),] - Z[(i+15*T1),] + W[(i+19*T1),] - Z[(i+19*T1),])/2./two
    # (4V + V4)/2 for LL, LF, FL, FF -> anti-symmetric sinh
    W[(j+8*T1),]  <- W[(j+4*T1),]
    W[(j+9*T1),]  <- W[(j+5*T1),]
    W[(j+10*T1),] <- W[(j+6*T1),]
    W[(j+11*T1),] <- W[(j+7*T1),]
    # gamma_i (VV) for LL, LF, FL, FF -> symmetric cosh
    W[(j+12*T1),] <- (W[(i+4*T1),] + Z[(i+4*T1),])/two
    W[(j+13*T1),] <- (W[(i+5*T1),] + W[(i+6*T1),] + Z[(i+5*T1),] + Z[(i+6*T1),])/2./two
    W[(j+14*T1),] <- W[(j+13*T1),]
    W[(j+15*T1),] <- (W[(i+7*T1),] + Z[(i+7*T1),])/two
    # gamma_i gamma_5 (AA) -> symmetric -cosh
    W[(j+16*T1),] <- -(W[(i+8*T1),] + Z[(i+8*T1),])/two
    W[(j+17*T1),] <- -(W[(i+9*T1),] + W[(i+10*T1),] + Z[(i+9*T1),] + Z[(i+10*T1),])/2./two
    W[(j+18*T1),] <- (W[(j+17*T1),])
    W[(j+19*T1),] <- -(W[(i+11*T1),] + Z[(i+11*T1),])/two
    # (4A + A4)/2 for LL, LF, FL, FF -> antisymmetric sinh
    W[(j+20*T1),] <- (W[(i+20*T1),] - Z[(i+20*T1),] + W[(i+24*T1),] - Z[(i+24*T1),])/2./two
    W[(j+21*T1),] <- (W[(i+21*T1),] - Z[(i+21*T1),] + W[(i+26*T1),] - Z[(i+26*T1),])/2./two
    W[(j+22*T1),] <- (W[(i+22*T1),] - Z[(i+22*T1),] + W[(i+25*T1),] - Z[(i+25*T1),])/2./two
    W[(j+23*T1),] <- (W[(i+23*T1),] - Z[(i+23*T1),] + W[(i+27*T1),] - Z[(i+27*T1),])/2./two
    # (4A + A4)/2 -> sinh
    W[(j+24*T1),] <- W[(j+20*T1),]
    W[(j+25*T1),] <- W[(j+21*T1),]
    W[(j+26*T1),] <- W[(j+22*T1),]
    W[(j+27*T1),] <- W[(j+23*T1),]
    # (VA + AV)/2 use AV -> cosh
    W[(j+28*T1),] <- (W[(i+28*T1),] + Z[(i+28*T1),] + W[(i+32*T1),] + Z[(i+32*T1),])/2./two
    W[(j+29*T1),] <- (W[(i+29*T1),] + Z[(i+29*T1),] + W[(i+34*T1),] + Z[(i+34*T1),])/2./two
    W[(j+30*T1),] <- (W[(i+30*T1),] + Z[(i+30*T1),] + W[(i+33*T1),] + Z[(i+33*T1),])/2./two
    W[(j+31*T1),] <- (W[(i+31*T1),] + Z[(i+31*T1),] + W[(i+35*T1),] + Z[(i+35*T1),])/2./two
    # (VA + AV)/2 -> cosh
    W[(j+32*T1),] <- W[(j+28*T1),]
    W[(j+33*T1),] <- W[(j+29*T1),]
    W[(j+34*T1),] <- W[(j+30*T1),]
    W[(j+35*T1),] <- W[(j+31*T1),]
  }
  return(invisible(W))
}


fitmasses.vector <- function(Cor, Err, t1, t2, Time, par=c(1.,0.1,0.12),
                         N=2, no.masses=1, no=1, kludge=TRUE) {
  Thalf <- Time/2
  T1 <- Thalf+1
  t1p1 <- (t1+1)
  t2p1 <- (t2+1)
  tr <- (t2-t1+1)

  if(no.masses == 1) {
#    fit <- optim(par, ChiSqr.1mass, method="BFGS", Thalf=Thalf,
#                     x=c((t1):(t2)), y=Cor, err=Err, tr=tr, N=N)
    fit <- gsl_fit_correlator_matrix(par, Thalf=Thalf,
                    x=c((t1):(t2)), y=Cor, err=Err, tr = tr, N=N)
    return(abs(fit$par[N+1]))
  }
  else if (no.masses == 2) {
#    fit <- optim(par, ChiSqr.2mass, method="BFGS", Thalf=Thalf,
#                     x=c((t1):(t2)), y=Cor, err=Err, tr=tr, N=N, kludge=kludge)
    fit <- gsl_fit_correlator_matrix(par, Thalf=Thalf,
                    x=c((t1):(t2)), y=Cor, err=Err, tr = tr, N=N, no_masses = 2)

    return(sort(abs(fit$par[c((N+1),(2*N+2))]))[no])
  }
  else if (no.masses == 3) {
    fit <- optim(par, ChiSqr.3mass, method="BFGS", Thalf=Thalf,
                     x=c((t1):(t2)), y=Cor, err=Err, tr=tr, N=N, kludge=kludge)
    return(sort(abs(fit$par[c((N+1),(2*N+2),(3*N+3))]))[no])
  }
}

fitmasses.vector.boot <- function(Z, d, Err, t1, t2, Time, par=c(1.,0.1,0.12),
                                  N=2, no.masses=1, kludge=TRUE) {
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
#    fit <- optim(par, ChiSqr.1mass, method="BFGS", Thalf=Thalf,
#                     x=c((t1):(t2)), y=Cor, err=Err, tr=tr, N=N)
    fit <- gsl_fit_correlator_matrix(par, Thalf=Thalf,
                    x=c((t1):(t2)), y=Cor, err=Err, tr = tr, N=N)
    return(c(fit$par[c(1:N)], abs(fit$par[N+1]), fit$value))
  }
  else if (no.masses == 2) {
#    fit <- optim(par, ChiSqr.2mass, method="BFGS", Thalf=Thalf,
#                     x=c((t1):(t2)), y=Cor, err=Err, tr=tr, N=N, kludge=kludge)
    fit <- gsl_fit_correlator_matrix(par, Thalf=Thalf,
                    x=c((t1):(t2)), y=Cor, err=Err, tr = tr, N=N, no_masses = 2)

    sort.ind <- order(fit$par[c((N+1),(2*N+2))])
    return(c(fit$par[c(((sort.ind[1]-1)*(N+1)+1):((sort.ind[1])*(N+1)-1))],
             abs(fit$par[sort.ind[1]*(N+1)]),
             fit$par[c(((sort.ind[2]-1)*(N+1)+1):((sort.ind[2])*(N+1)-1))],
             abs(fit$par[sort.ind[2]*(N+1)]),
             fit$value))
  }
  else if (no.masses == 3) {
    fit <- optim(par, ChiSqr.3mass, method="BFGS", Thalf=Thalf,
                     x=c((t1):(t2)), y=Cor, err=Err, tr=tr, N=N, kludge=kludge)
    sort.ind <- order(fit$par[c((N+1),(2*N+2),(3*N+3))])
    return(c(fit$par[c(((sort.ind[1]-1)*(N+1)+1):((sort.ind[1])*(N+1)-1))],
             abs(fit$par[sort.ind[1]*(N+1)]),
             fit$par[c(((sort.ind[2]-1)*(N+1)+1):((sort.ind[2])*(N+1)-1))],
             abs(fit$par[sort.ind[2]*(N+1)]),
             fit$par[c(((sort.ind[3]-1)*(N+1)+1):((sort.ind[3])*(N+1)-1))],
             abs(fit$par[sort.ind[3]*(N+1)]),
             fit$value))
  }

}

