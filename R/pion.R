pion <- function(cmicor, mu=0.1, kappa=0.156, t1, t2, S=1.5, pl=FALSE, skip=0,
                variational=list(ta=3, tb=4, N=6), ind.vec=c(1,3,4,5),
                no.masses=1, matrix.size=2, boot.R=99, boot.l=10, tsboot.sim="geom",
                method="uwerr", mass.guess, par.guess) {
  
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
    par.guess <- c(1.,0.8,0.1,0.1,0.1,0.1, 0.1,0.1,0.1,0.1,0.1, 0.1,0.1,0.1,0.1,0.1,) 
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

  Z <- array(cmicor[((Skip):Length),ind.vec[3]], 
             dim=c(nrObs*(T1)*4,(length(cmicor[((Skip):Length),ind.vec[3]])/(nrObs*(T1)*4))))
  # negative times
  W <- array(cmicor[((Skip):Length),ind.vec[4]], 
             dim=c(nrObs*(T1)*4,(length(cmicor[((Skip):Length),ind.vec[4]])/(nrObs*(T1)*4))))

  rm(cmicor)
  W <- arrangeCor.pion(T1=T1, W=W, Z=Z)
  rm(Z)

  options(show.error.messages = FALSE)
  pion.eff.ll <- effectivemass(from=(t1+1), to=(t2+1), Time, W[1:T1,] , pl=FALSE, S=1.5)
  pion.eff.lf <- effectivemass(from=(t1+1), to=(t2+1), Time, W[(T1+1):(2*T1),] , pl=FALSE, S=1.5)
  pion.eff.ff <- effectivemass(from=(t1+1), to=(t2+1), Time, W[(3*T1+1):(4*T1),] , pl=FALSE, S=1.5) 
  options(show.error.messages = TRUE)
  
  pion.eff <- data.frame(t=pion.eff.ll$t, mll=pion.eff.ll$mass, dmll=pion.eff.ll$dmass,
                         mlf=pion.eff.lf$mass, dmlf=pion.eff.lf$dmass,
                         mff=pion.eff.ff$mass, dmff=pion.eff.ff$dmass)
  
  Cor <- rep(0., times=9*4*T1)
  E <- rep(0., times=9*4*T1)
  
  for(i in 1:(9*4*T1)) {
    Cor[i] = mean(W[(i),])
    E[i] = uwerrprimary(W[(i),], pl=F)$dvalue
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
  if(FALSE) {

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
    print(par)
  }
  
  variational.masses <-  -log(abs(variational.solve$values[variational.sortindex]))/(tb-ta)
  rm(C1, C2, C3, ta, tb, N)

  # Index vector of data to be used in the analysis
  ii <- c((t1p1):(t2p1), (t1p1+T1):(t2p1+T1), (t1p1+3*T1):(t2p1+3*T1))
  if(matrix.size > 2) {
    ii <- c(ii, (t1p1+4*T1):(t2p1+4*T1), (t1p1+5*T1):(t2p1+5*T1),
            (t1p1+6*T1):(t2p1+6*T1), (t1p1+7*T1):(t2p1+7*T1),
            (t1p1+12*T1):(t2p1+12*T1), (t1p1+13*T1):(t2p1+13*T1),
            (t1p1+15*T1):(t2p1+15*T1))
  }
  if(matrix.size > 4) {
    ii <- c(ii, (t1p1+16*T1):(t2p1+16*T1), (t1p1+17*T1):(t2p1+17*T1),
            (t1p1+19*T1):(t2p1+19*T1),
            (t1p1+20*T1):(t2p1+20*T1), (t1p1+21*T1):(t2p1+21*T1),
            (t1p1+22*T1):(t2p1+22*T1), (t1p1+23*T1):(t2p1+23*T1),
            (t1p1+28*T1):(t2p1+28*T1), (t1p1+29*T1):(t2p1+29*T1),
            (t1p1+30*T1):(t2p1+30*T1), (t1p1+31*T1):(t2p1+31*T1))
  }

  #BFGS
  if(no.masses == 1) {
    pionfit <- optim(par, ChiSqr.1mass, method="BFGS", control=list(trace=0),Thalf=Thalf,
                    x=c((t1):(t2)), y=Cor[ii], err=E[ii], tr = (t2-t1+1), N=matrix.size)
    fit.mass <- abs(pionfit$par[matrix.size+1])
    fit.fpi <- 2*kappa*2*mu/sqrt(2)*abs(pionfit$par[1])/sqrt(fit.mass^3)
    if(matrix.size > 2) {
      fit.pcac <- abs(pionfit$par[5])*pionfit$par[3]/pionfit$par[1]/2.
    }

  }
  else if(no.masses == 2) {
    pionfit <- optim(par, ChiSqr.2mass, method="BFGS", control=list(trace=0),Thalf=Thalf,
                    x=c((t1):(t2)), y=Cor[ii], err=E[ii], tr = (t2-t1+1), N=matrix.size)
    fit.mass <- sort(abs(pionfit$par[c((matrix.size+1),(2*matrix.size+2))]))
  }
  else if(no.masses > 2) {
    pionfit <- optim(par, ChiSqr.3mass, method="BFGS", control=list(trace=0),Thalf=Thalf,
                    x=c((t1):(t2)), y=Cor[ii], err=E[ii], tr = (t2-t1+1), N=matrix.size)
    fit.mass <- sort(abs(pionfit$par[c((matrix.size+1),(2*matrix.size+2),(3*matrix.size+3))]))
  }
  if(pionfit$convergence != 0) {
    warning("optim did not converge for pionfit!", call.=F)
  }

  fit.dof <- (t2-t1+1)*3-length(pionfit$par)
  fit.chisqr <- pionfit$value

  if(pl) {
    plot.effmass(m=fit.mass, ll=pion.eff.ll, lf=pion.eff.lf, ff=pion.eff.ff)
  }

  fit.uwerrm <- NULL
  fit.uwerrf <- NULL
  fit.uwerrpcac <- NULL
  fit.uwerrm2 <- NULL
  fit.uwerrm3 <- NULL
  fit.boot <- NULL
  fit.tsboot <- NULL
  if(method == "uwerr" || method == "all") {
    fit.uwerrm <- uwerr(f=fitmasses.pion, data=W[ii,], S=S, pl=pl,
                        Time=Time, t1=t1, t2=t2, Err=E[ii], par=par, N=matrix.size, no.masses=no.masses)

    fit.uwerrf <- uwerr(f=fitf.pion, data=W[ii,], S=S, pl=pl,
                        Time=Time, t1=t1, t2=t2, Err=E[ii], par=par, N=matrix.size, no.masses=no.masses,
                        kappa=kappa, mu=mu)
    if(matrix.size > 2) {
      fit.uwerrpcac <- uwerr(f=fitmpcac.pion, data=W[ii,], S=S, pl=pl,
                             Time=Time, t1=t1, t2=t2, Err=E[ii], par=par, N=matrix.size, no.masses=no.masses,
                             kappa=kappa, mu=mu)
    }
    
    if(no.masses == 2) {
      fit.uwerrm2 <- uwerr(f=fitmasses.pion, data=W[ii,], S=S, pl=pl,
                           Time=Time, t1=t1, t2=t2, Err=E[ii], par=par, N=matrix.size, no.masses=no.masses, no=2)
    }
    if(no.masses > 2) {
      fit.uwerrm3 <- uwerr(f=fitmasses.pion, data=W[ii,], S=S, pl=pl,
                           Time=Time, t1=t1, t2=t2, Err=E[ii], par=par, N=matrix.size, no.masses=no.masses, no=3)
    }
  }
  if(method == "boot" || method == "all") {
    fit.boot <- boot(data=t(W[ii,]), statistic=fit.pion.boot, R=boot.R, stype="i",
                     Time=Time, t1=t1, t2=t2, Err=E[ii], par=par, N=matrix.size, no.masses=no.masses,
                     kappa=kappa, mu=mu)

    fit.tsboot <- tsboot(tseries=t(W[ii,]), statistic=fit.pion.boot, R=boot.R, l=boot.l, sim=tsboot.sim,
                         Time=Time, t1=t1, t2=t2, Err=E[ii], par=par, N=matrix.size, no.masses=no.masses,
                         kappa=kappa, mu=mu)
  }

  
  Chi <- rep(0., times=9*4*T1)
  Fit <- rep(0., times=9*4*T1)

  jj <-  c(t1p1:t2p1)
  Fit[jj] <- pionfit$par[1]^2*CExp(m=fit.mass, Time=2*Thalf, x=jj-1)
  Fit[jj+T1] <- pionfit$par[1]*pionfit$par[2]*CExp(m=fit.mass, Time=2*Thalf, x=jj-1)
  Fit[jj+2*T1] <- pionfit$par[1]*pionfit$par[2]*CExp(m=fit.mass, Time=2*Thalf, x=jj-1)
  Fit[jj+3*T1] <- pionfit$par[2]*pionfit$par[2]*CExp(m=fit.mass, Time=2*Thalf, x=jj-1)
  if(matrix.size > 2) {
    Fit[jj+20*T1] <- pionfit$par[1]*pionfit$par[3]*CExp(m=fit.mass, Time=2*Thalf, x=jj-1, sign=-1.)
    Fit[jj+21*T1] <- pionfit$par[1]*pionfit$par[4]*CExp(m=fit.mass, Time=2*Thalf, x=jj-1, sign=-1.)
    Fit[jj+22*T1] <- pionfit$par[2]*pionfit$par[3]*CExp(m=fit.mass, Time=2*Thalf, x=jj-1, sign=-1.)
    Fit[jj+23*T1] <- pionfit$par[2]*pionfit$par[4]*CExp(m=fit.mass, Time=2*Thalf, x=jj-1, sign=-1.)

    Fit[jj+16*T1] <- pionfit$par[3]*pionfit$par[3]*CExp(m=fit.mass, Time=2*Thalf, x=jj-1)
    Fit[jj+17*T1] <- pionfit$par[3]*pionfit$par[4]*CExp(m=fit.mass, Time=2*Thalf, x=jj-1)
    Fit[jj+18*T1] <- pionfit$par[3]*pionfit$par[4]*CExp(m=fit.mass, Time=2*Thalf, x=jj-1)
    Fit[jj+19*T1] <- pionfit$par[4]*pionfit$par[4]*CExp(m=fit.mass, Time=2*Thalf, x=jj-1)
  }
  if(matrix.size > 4) {
    Fit[jj+12*T1] <- pionfit$par[5]*pionfit$par[5]*CExp(m=fit.mass, Time=2*Thalf, x=jj-1)
    Fit[jj+13*T1] <- pionfit$par[5]*pionfit$par[6]*CExp(m=fit.mass, Time=2*Thalf, x=jj-1)
    Fit[jj+14*T1] <- pionfit$par[5]*pionfit$par[6]*CExp(m=fit.mass, Time=2*Thalf, x=jj-1)
    Fit[jj+15*T1] <- pionfit$par[6]*pionfit$par[6]*CExp(m=fit.mass, Time=2*Thalf, x=jj-1)

    Fit[jj+4*T1] <- pionfit$par[1]*pionfit$par[5]*CExp(m=fit.mass, Time=2*Thalf, x=jj-1, sign=-1.)
    Fit[jj+5*T1] <- pionfit$par[1]*pionfit$par[6]*CExp(m=fit.mass, Time=2*Thalf, x=jj-1, sign=-1.)
    Fit[jj+6*T1] <- pionfit$par[2]*pionfit$par[5]*CExp(m=fit.mass, Time=2*Thalf, x=jj-1, sign=-1.)
    Fit[jj+7*T1] <- pionfit$par[2]*pionfit$par[6]*CExp(m=fit.mass, Time=2*Thalf, x=jj-1, sign=-1.)
    
    Fit[jj+28*T1] <- pionfit$par[3]*pionfit$par[5]*CExp(m=fit.mass, Time=2*Thalf, x=jj-1)
    Fit[jj+29*T1] <- pionfit$par[3]*pionfit$par[6]*CExp(m=fit.mass, Time=2*Thalf, x=jj-1)
    Fit[jj+30*T1] <- pionfit$par[4]*pionfit$par[5]*CExp(m=fit.mass, Time=2*Thalf, x=jj-1)
    Fit[jj+31*T1] <- pionfit$par[4]*pionfit$par[6]*CExp(m=fit.mass, Time=2*Thalf, x=jj-1)
  }
  
  Chi[ii] <- (Fit[ii]-Cor[ii])/E[ii]
  
  res <- list(fitresult=pionfit, t1=t1, t2=t2, N=length(W[1,]), Time=Time,
              fitdata=data.frame(t=(jj-1), Fit=Fit[ii], Cor=Cor[ii], Err=E[ii], Chi=Chi[ii]),
              uwerrresultmps=fit.uwerrm, uwerrresultmps2=fit.uwerrm2, uwerrresultmps3=fit.uwerrm3,
              uwerrresultfps=fit.uwerrf, uwerrresultmpcac=fit.uwerrpcac,
              boot=fit.boot, tsboot=fit.tsboot,
              effmass=pion.eff, kappa=kappa, mu=mu,
              variational.masses=variational.masses, no.masses=no.masses,
              matrix.size = matrix.size)
  attr(res, "class") <- c("pionfit", "list")  
  return(invisible(res))
}

arrangeCor.pion <- function(T1, W, Z) {

  for(i in 1:(T1)) {
    two <- 2.
    if(i==1 || i==(T1)) {
      # Take care of zeros in the correlators when summing t and T-t+1
      two <- 1.
    }

    # PP for LL, (LF + FL)/2, FF -> symmetric cosh
    W[i,] <- (W[i,] + Z[i,])/two
    W[(i+T1),] <- (W[(i+T1),] + W[(i+2*T1),] + Z[(i+T1),] + Z[(i+2*T1),])/two/2.
    W[(i+2*T1),] <- W[(i+T1),]
    W[(i+3*T1),] <- (W[(i+3*T1),] + Z[(i+3*T1),])/two
    # PA for LL, LF, FL, FF -> antisymmetric -sinh
    W[(i+4*T1),] <- (W[(i+4*T1),] - Z[(i+4*T1),])/two
    W[(i+5*T1),] <- (W[(i+5*T1),] - Z[(i+5*T1),])/two
    W[(i+6*T1),] <- (W[(i+6*T1),] - Z[(i+6*T1),])/two
    W[(i+7*T1),] <- (W[(i+7*T1),] - Z[(i+7*T1),])/two
    # use again PA -> -sinh
    W[(i+8*T1),] <- W[(i+4*T1),]
    W[(i+9*T1),] <- W[(i+5*T1),]
    W[(i+10*T1),] <- W[(i+6*T1),]
    W[(i+11*T1),] <- W[(i+7*T1),]
    # AA for LL, (LF + FL)/2, FF -> symmetric cosh
    W[(i+12*T1),] <- (W[(i+12*T1),] + Z[(i+12*T1),])/two
    W[(i+13*T1),] <- (W[(i+13*T1),] + W[(i+14*T1),] + Z[(i+13*T1),] + Z[(i+14*T1),])/two/2.
    W[(i+14*T1),] <-  W[(i+13*T1),]
    W[(i+15*T1),] <- (W[(i+15*T1),] + Z[(i+15*T1),])/two
    # 44 for LL, LF + FL/2, FF -> -cosh
    W[(i+16*T1),] <- -(W[(i+16*T1),] + Z[(i+16*T1),])/two
    W[(i+17*T1),] <- -(W[(i+17*T1),] + W[(i+17*T1),] + Z[(i+18*T1),] + Z[(i+18*T1),])/two/2.
    W[(i+18*T1),] <-   W[(i+17*T1),]
    W[(i+19*T1),] <- -(W[(i+19*T1),] + Z[(i+19*T1),])/two
    # P4 for LL, LF, FL, FF -> antisymmetric -sinh
    W[(i+20*T1),] <- (W[(i+20*T1),] - Z[(i+20*T1),])/two
    W[(i+21*T1),] <- (W[(i+21*T1),] - Z[(i+21*T1),])/two
    W[(i+22*T1),] <- (W[(i+22*T1),] - Z[(i+22*T1),])/two
    W[(i+23*T1),] <- (W[(i+23*T1),] - Z[(i+23*T1),])/two
    # 4P use again P4 -> -sinh
    W[(i+24*T1),] <- W[(i+20*T1),]
    W[(i+25*T1),] <- W[(i+21*T1),]
    W[(i+26*T1),] <- W[(i+22*T1),]
    W[(i+27*T1),] <- W[(i+23*T1),]
    # A4 use 4A -> cosh
    W[(i+28*T1),] <- (W[(i+32*T1),] + Z[(i+32*T1),])/two
    W[(i+29*T1),] <- (W[(i+34*T1),] + Z[(i+34*T1),])/two
    W[(i+30*T1),] <- (W[(i+33*T1),] + Z[(i+33*T1),])/two
    W[(i+31*T1),] <- (W[(i+35*T1),] + Z[(i+35*T1),])/two
    # 4A for LL, LF, FL, FF -> cosh
    W[(i+32*T1),] <- W[(i+28*T1),]
    W[(i+33*T1),] <- W[(i+29*T1),]
    W[(i+34*T1),] <- W[(i+30*T1),]
    W[(i+35*T1),] <- W[(i+31*T1),]
  }
  return(invisible(W))
}

fitmasses.pion <- function(Cor, Err, t1, t2, Time, par=c(1.,0.1,0.12),
                         N=2, no.masses=1, no=1, kludge=FALSE) {
  Thalf <- Time/2
  T1 <- Thalf+1
  t1p1 <- (t1+1)
  t2p1 <- (t2+1)
  tr <- (t2-t1+1)

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

fitf.pion <- function(Cor, Err, t1, t2, Time, par=c(1.,0.1,0.12),
                         N=2, no.masses=1, no=1, kappa, mu, kludge=FALSE) {
  Thalf <- Time/2
  T1 <- Thalf+1
  t1p1 <- (t1+1)
  t2p1 <- (t2+1)
  tr <- (t2-t1+1)

  if(no.masses == 1) {
    fit <- optim(par, ChiSqr.1mass, method="BFGS", Thalf=Thalf,
                     x=c((t1):(t2)), y=Cor, err=Err, tr=tr, N=N)
    sort.ind <- c(1)
  }
  else if (no.masses == 2) {
    fit <- optim(par, ChiSqr.2mass, method="BFGS", Thalf=Thalf,
                     x=c((t1):(t2)), y=Cor, err=Err, tr=tr, N=N, kludge=kludge)
    sort.ind <- order(fit$par[c((N+1),(2*N+2))])
  }
  else if (no.masses == 3) {
    fit <- optim(par, ChiSqr.3mass, method="BFGS", Thalf=Thalf,
                     x=c((t1):(t2)), y=Cor, err=Err, tr=tr, N=N, kludge=kludge)
    sort.ind <- order(fit$par[c((N+1),(2*N+2),(3*N+3))])
  }
#  fit.fpi <- 2*kappa*2*mu/sqrt(2)*abs(fit$par[(sort.ind[1]-1)*(N+1)+1])/sqrt(abs(fit$par[sort.ind[1]*(N+1)])^3)
  fit.fpi <- abs(fit$par[(sort.ind[1]-1)*(N+1)+1])/sqrt(abs(fit$par[sort.ind[1]*(N+1)])^3)
  return(fit.fpi)
}

fitmpcac.pion <- function(Cor, Err, t1, t2, Time, par=c(1.,0.1,0.12),
                          N=2, no.masses=1, no=1, kappa, mu, kludge=FALSE) {
  Thalf <- Time/2
  T1 <- Thalf+1
  t1p1 <- (t1+1)
  t2p1 <- (t2+1)
  tr <- (t2-t1+1)

  if(no.masses == 1) {
    fit <- optim(par, ChiSqr.1mass, method="BFGS", Thalf=Thalf,
                     x=c((t1):(t2)), y=Cor, err=Err, tr=tr, N=N)
    sort.ind <- c(1)
  }
  else if (no.masses == 2) {
    fit <- optim(par, ChiSqr.2mass, method="BFGS", Thalf=Thalf,
                     x=c((t1):(t2)), y=Cor, err=Err, tr=tr, N=N, kludge=kludge)
    sort.ind <- order(fit$par[c((N+1),(2*N+2))])
  }
  else if (no.masses == 3) {
    fit <- optim(par, ChiSqr.3mass, method="BFGS", Thalf=Thalf,
                     x=c((t1):(t2)), y=Cor, err=Err, tr=tr, N=N, kludge=kludge)
    sort.ind <- order(fit$par[c((N+1),(2*N+2),(3*N+3))])
  }
  fit.pcac <- abs(fit$par[sort.ind[1]*(N+1)])*fit$par[(sort.ind[1]-1)*(N+1)+3]/fit$par[(sort.ind[1]-1)*(N+1)+1]/2.
  return(fit.pcac)
}


fit.pion.boot <- function(Z, d, Err, t1, t2, Time, par=c(1.,0.1,0.12),
                          N=2, no.masses=1, kludge=FALSE, kappa, mu) {
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
    fit <- optim(par, ChiSqr.1mass, method="BFGS", Thalf=Thalf,
                     x=c((t1):(t2)), y=Cor, err=Err, tr=tr, N=N)
    sort.ind <- c(1)
    fit.fpi <- 2*kappa*2*mu/sqrt(2)*abs(fit$par[(sort.ind[1]-1)*(N+1)+1])/sqrt(abs(fit$par[sort.ind[1]*(N+1)])^3)
    if(N > 2) {
      fit.pcac <- abs(fit$par[sort.ind[1]*(N+1)])*fit$par[(sort.ind[1]-1)*(N+1)+3]/fit$par[(sort.ind[1]-1)*(N+1)+1]/2.
      return(c(abs(fit$par[N+1]), fit.fpi, fit$par[c(1:N)],
               fit.pcac,
               fit$value))
    }
    else{
      return(c(abs(fit$par[N+1]), fit.fpi, fit$par[c(1:N)],
               fit$value))
    }
  }
  else if (no.masses == 2) {
    fit <- optim(par, ChiSqr.2mass, method="BFGS", Thalf=Thalf,
                     x=c((t1):(t2)), y=Cor, err=Err, tr=tr, N=N, kludge=kludge)
    sort.ind <- order(fit$par[c((N+1),(2*N+2))])
    fit.fpi <- 2*kappa*2*mu/sqrt(2)*abs(fit$par[(sort.ind[1]-1)*(N+1)+1])/sqrt(abs(fit$par[sort.ind[1]*(N+1)])^3)
    if(N > 2) {
      fit.pcac <- abs(fit$par[sort.ind[1]*(N+1)])*fit$par[(sort.ind[1]-1)*(N+1)+3]/fit$par[(sort.ind[1]-1)*(N+1)+1]/2.
      
      return(c(abs(fit$par[sort.ind[1]*(N+1)]), fit.fpi,
               fit$par[c(((sort.ind[1]-1)*(N+1)+1):((sort.ind[1])*(N+1)-1))],
               fit.pcac,
               abs(fit$par[sort.ind[2]*(N+1)]),
               fit$par[c(((sort.ind[2]-1)*(N+1)+1):((sort.ind[2])*(N+1)-1))],
               fit$value))
    }
    else {
      return(c(abs(fit$par[sort.ind[1]*(N+1)]), fit.fpi,
               fit$par[c(((sort.ind[1]-1)*(N+1)+1):((sort.ind[1])*(N+1)-1))],
               abs(fit$par[sort.ind[2]*(N+1)]),
               fit$par[c(((sort.ind[2]-1)*(N+1)+1):((sort.ind[2])*(N+1)-1))],
               fit$value))
    }
  }
  else if (no.masses == 3) {
    fit <- optim(par, ChiSqr.3mass, method="BFGS", Thalf=Thalf,
                     x=c((t1):(t2)), y=Cor, err=Err, tr=tr, N=N, kludge=kludge)
    sort.ind <- order(fit$par[c((N+1),(2*N+2),(3*N+3))])
    fit.fpi <- 2*kappa*2*mu/sqrt(2)*abs(fit$par[(sort.ind[1]-1)*(N+1)+1])/sqrt(abs(fit$par[sort.ind[1]*(N+1)])^3)
    if(N > 2) {
      fit.pcac <- abs(fit$par[sort.ind[1]*(N+1)])*fit$par[(sort.ind[1]-1)*(N+1)+3]/fit$par[(sort.ind[1]-1)*(N+1)+1]/2.
      
      return(c(abs(fit$par[sort.ind[1]*(N+1)]), fit.fpi,
               fit$par[c(((sort.ind[1]-1)*(N+1)+1):((sort.ind[1])*(N+1)-1))],
               fit.pcac,
               abs(fit$par[sort.ind[2]*(N+1)]),
               fit$par[c(((sort.ind[2]-1)*(N+1)+1):((sort.ind[2])*(N+1)-1))],
               abs(fit$par[sort.ind[3]*(N+1)]),
               fit$par[c(((sort.ind[3]-1)*(N+1)+1):((sort.ind[3])*(N+1)-1))],
               fit$value))
    }
    else {
      return(c(abs(fit$par[sort.ind[1]*(N+1)]), fit.fpi,
               fit$par[c(((sort.ind[1]-1)*(N+1)+1):((sort.ind[1])*(N+1)-1))],
               abs(fit$par[sort.ind[2]*(N+1)]),
               fit$par[c(((sort.ind[2]-1)*(N+1)+1):((sort.ind[2])*(N+1)-1))],
               abs(fit$par[sort.ind[3]*(N+1)]),
               fit$par[c(((sort.ind[3]-1)*(N+1)+1):((sort.ind[3])*(N+1)-1))],
               fit$value))
    }
  }
}


