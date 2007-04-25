# This is a 2x2 Fit to the gamma_i gamma_4 correlator only

rho <- function(cmicor, mu=0.1, kappa=0.156, t1, t2, S=1.5, pl=FALSE, skip=0,
                variational=list(ta=4, tb=5, N=6), ind.vec=c(1,3,4,5),
                no.masses=1, matrix.size = 2) {
  
  if(missing(cmicor)) {
    stop("Error! Data is missing!")
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
  W <- arrangeCor.vector(T1=T1, W=W, Z=Z)
  rm(Z)

  rho.eff.ll <- effectivemass(from=(t1+1), to=(t2+1), Time, W[1:T1,] , pl=FALSE, S=1.5)
  rho.eff.lf <- effectivemass(from=(t1+1), to=(t2+1), Time, W[(T1+1):(2*T1),] , pl=FALSE, S=1.5)
  rho.eff.ff <- effectivemass(from=(t1+1), to=(t2+1), Time, W[(3*T1+1):(4*T1),] , pl=FALSE, S=1.5) 

  rho.eff <- data.frame(t=rho.eff.ll$t, mll=rho.eff.ll$mass, dmll=rho.eff.ll$dmass,
                         mlf=rho.eff.lf$mass, dmlf=rho.eff.lf$dmass,
                         mff=rho.eff.ff$mass, dmff=rho.eff.ff$dmass)
  
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
               left.vectors[(1:matrix.size),i],
               -log(abs(variational.solve$values[variational.sortindex[2]]))/(tb-ta)) 
    }
  }
  print(par)
  
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
    rhofit <- optim(par, ChiSqr.1mass, method="BFGS", control=list(trace=0),Thalf=Thalf,
                    x=c((t1):(t2)), y=Cor[ii], err=E[ii], tr = (t2-t1+1), N=matrix.size)
  }
  else if(no.masses == 2) {
    rhofit <- optim(par, ChiSqr.2mass, method="BFGS", control=list(trace=0),Thalf=Thalf,
                    x=c((t1):(t2)), y=Cor[ii], err=E[ii], tr = (t2-t1+1), N=matrix.size)
  }
  if(rhofit$convergence != 0) {
    warning("optim did not converge for rhofit!", call.=F)
  }
  print(rhofit)
  fit.mass <- abs(rhofit$par[matrix.size+1])

  fit.dof <- (t2-t1+1)*3-length(rhofit$par)
  fit.chisqr <- rhofit$value

  if(pl) {
    plot.effmass(m=fit.mass, ll=rho.eff.ll, lf=rho.eff.lf, ff=rho.eff.ff)
  }

  fit.uwerrm <- uwerr(f=fit2by2mass2, data=W[ii,], S=S, pl=pl,
                      Time=Time, t1=t1, t2=t2, Err=E[ii], par=par, N=matrix.size, no.masses=no.masses)
  fit.uwerrm2 <- NULL
  if(no.masses == 2) {
    fit.uwerrm2 <- uwerr(f=fit2by2mass2, data=W[ii,], S=S, pl=pl,
                        Time=Time, t1=t1, t2=t2, Err=E[ii], par=par, N=matrix.size, no.masses=no.masses, no=2)
  }
  Chi <- rep(0., times=9*4*T1)
  Fit <- rep(0., times=9*4*T1)

  jj <-  c(t1p1:t2p1)
  Fit[jj] <- rhofit$par[1]^2*CExp(m=fit.uwerrm$value, Time=2*Thalf, x=jj-1)
  Fit[jj+T1] <- rhofit$par[1]*rhofit$par[2]*CExp(m=fit.uwerrm$value, Time=2*Thalf, x=jj-1)
  Fit[jj+2*T1] <- rhofit$par[1]*rhofit$par[2]*CExp(m=fit.uwerrm$value, Time=2*Thalf, x=jj-1)
  Fit[jj+3*T1] <- rhofit$par[2]*rhofit$par[2]*CExp(m=fit.uwerrm$value, Time=2*Thalf, x=jj-1)
  if(matrix.size > 2) {
    Fit[jj+20*T1] <- rhofit$par[1]*rhofit$par[3]*CExp(m=fit.uwerrm$value, Time=2*Thalf, x=jj-1, sign=-1.)
    Fit[jj+21*T1] <- rhofit$par[1]*rhofit$par[4]*CExp(m=fit.uwerrm$value, Time=2*Thalf, x=jj-1, sign=-1.)
    Fit[jj+22*T1] <- rhofit$par[2]*rhofit$par[3]*CExp(m=fit.uwerrm$value, Time=2*Thalf, x=jj-1, sign=-1.)
    Fit[jj+23*T1] <- rhofit$par[2]*rhofit$par[4]*CExp(m=fit.uwerrm$value, Time=2*Thalf, x=jj-1, sign=-1.)

    Fit[jj+16*T1] <- rhofit$par[3]*rhofit$par[3]*CExp(m=fit.uwerrm$value, Time=2*Thalf, x=jj-1)
    Fit[jj+17*T1] <- rhofit$par[3]*rhofit$par[4]*CExp(m=fit.uwerrm$value, Time=2*Thalf, x=jj-1)
    Fit[jj+18*T1] <- rhofit$par[3]*rhofit$par[4]*CExp(m=fit.uwerrm$value, Time=2*Thalf, x=jj-1)
    Fit[jj+19*T1] <- rhofit$par[4]*rhofit$par[4]*CExp(m=fit.uwerrm$value, Time=2*Thalf, x=jj-1)
  }
  Chi[ii] <- (Fit[ii]-Cor[ii])/E[ii]
  
  res <- list(fitresult=rhofit, t1=t1, t2=t2, N=length(W[1,]), Time=Time,
              fitdata=data.frame(t=(jj-1), Fit=Fit[ii], Cor=Cor[ii], Err=E[ii], Chi=Chi[ii]),
              uwerrresultmv=fit.uwerrm, uwerrresultmv2=fit.uwerrm2,
              effmass=rho.eff, kappa=kappa, mu=mu,
              variational.masses=variational.masses, no.masses=no.masses,
              matrix.size = matrix.size)
  attr(res, "class") <- c("cfit", "list")  
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

getNxNmatrix <- function(Cor, T1, t, N=2) {
  C1 <- array(0., dim=c(N,N))
  C1[1,1] = Cor[t]
  C1[1,2] = Cor[(t+T1)]
  C1[2,1] = Cor[(t+2*T1)]
  C1[2,2] = Cor[(t+3*T1)]
  if(N > 2) {
    C1[1,3] = Cor[(t+20*T1)]
    C1[3,1] = C1[1,3]
    C1[1,4] = Cor[(t+21*T1)]
    C1[4,1] = C1[1,4]
    C1[2,3] = Cor[(t+22*T1)]
    C1[3,2] = C1[2,3]
    C1[2,4] = Cor[(t+23*T1)]
    C1[4,2] = C1[2,4]
    C1[3,3] = Cor[(t+16*T1)]
    C1[3,4] = Cor[(t+17*T1)]
    C1[4,3] = C1[3,4]
    C1[4,4] = Cor[(t+19*T1)]
  }
  if(N > 4) {
    C1[5,5] = Cor[(t+12*T1)]
    C1[5,6] = Cor[(t+13*T1)]
    C1[6,5] = C1[5,6]
    C1[6,6] = Cor[(t+15*T1)]
    C1[1,5] = Cor[(t+4*T1)]
    C1[5,1] = C1[1,5]
    C1[1,6] = Cor[(t+5*T1)]
    C1[6,1] = C1[1,6]
    C1[2,5] = Cor[(t+6*T1)]
    C1[5,2] = C1[2,5]
    C1[2,6] = Cor[(t+7*T1)]
    C1[6,2] = C1[2,6]
    C1[3,5] = Cor[(t+28*T1)]
    C1[5,3] = C1[3,5]
    C1[3,6] = Cor[(t+29*T1)]
    C1[6,3] = C1[3,6]
    C1[4,5] = Cor[(t+30*T1)]
    C1[5,4] = C1[4,5]
    C1[4,6] = Cor[(t+31*T1)]
    C1[6,4] = C1[4,6]
  }
  return(invisible(C1))
}

fit2by2mass2 <- function(Cor, Err, t1, t2, Time, par=c(1.,0.1,0.12), N=2, no.masses=1, no=1) {
  Thalf <- Time/2
  T1 <- Thalf+1
  t1p1 <- (t1+1)
  t2p1 <- (t2+1)

  if(no.masses == 1) {
    fit <- optim(par, ChiSqr.1mass, method="BFGS", Thalf=Thalf,
                     x=c((t1):(t2)), y=Cor, err=Err, tr = (t2-t1+1), N=N)
    return(abs(fit$par[N+1]))
  }
  else if (no.masses == 2) {
    fit <- optim(par, ChiSqr.2mass, method="BFGS", Thalf=Thalf,
                     x=c((t1):(t2)), y=Cor, err=Err, tr = (t2-t1+1), N=N)
    return(sort(abs(fit$par[c((N+1),(2*N+2))]))[no])
  }
  else if (no.masses == 3) {
    fit <- optim(par, ChiSqr.3mass, method="BFGS", Thalf=Thalf,
                     x=c((t1):(t2)), y=Cor, err=Err, tr = (t2-t1+1), N=N)
    return(sort(abs(fit$par[c((N+1),(2*N+2),(3*N+3))]))[no])
  }
}

CExp <- function(m, Time, x, sign=1.) {
  return(0.5*(exp(-abs(m)*x) + sign*exp(-abs(m)*(Time-x))))
}

ChiSqr.1mass <- function(par, Thalf, x, y, err, tr, N=2) {
  # index of mass
  # l <- length(par)/no.masses
  ii <- c(1:tr)
  Sumall <- 0.
  if(N > 1) {
    # 44
    Sumall = Sumall + sum(((y[ii]
      - par[1]*par[1]*(CExp(m=abs(par[N+1]), Time=2*Thalf, x=x)))/err[ii])^2)
    
    Sumall = Sumall + sum(((y[ii+tr]
      - par[1]*par[2]*(CExp(m=abs(par[N+1]), Time=2*Thalf, x=x)))/err[ii+tr])^2)

    Sumall = Sumall + sum(((y[ii+2*tr]
      - par[2]*par[2]*(CExp(m=abs(par[N+1]), Time=2*Thalf, x=x)))/err[ii+2*tr])^2)
    
  }
  if(N > 2) {
    # 4A (sinh!)
    Sumall = Sumall + sum(((y[ii+3*tr]
      - par[1]*par[3]*(CExp(m=abs(par[N+1]), Time=2*Thalf, x=x, sign=-1.)))/err[ii+3*tr])^2)
    Sumall = Sumall + sum(((y[ii+4*tr]
      - par[1]*par[4]*(CExp(m=abs(par[N+1]), Time=2*Thalf, x=x, sign=-1.)))/err[ii+4*tr])^2)
    Sumall = Sumall + sum(((y[ii+5*tr]
      - par[2]*par[3]*(CExp(m=abs(par[N+1]), Time=2*Thalf, x=x, sign=-1.)))/err[ii+5*tr])^2)
    Sumall = Sumall + sum(((y[ii+6*tr]
      - par[2]*par[4]*(CExp(m=abs(par[N+1]), Time=2*Thalf, x=x, sign=-1.)))/err[ii+6*tr])^2)
    
    # AA
    Sumall = Sumall + sum(((y[ii+7*tr]
      - par[3]*par[3]*(CExp(m=abs(par[N+1]), Time=2*Thalf, x=x)))/err[ii+7*tr])^2)
    Sumall = Sumall + sum(((y[ii+8*tr]
      - par[3]*par[4]*(CExp(m=abs(par[N+1]), Time=2*Thalf, x=x)))/err[ii+8*tr])^2)
    Sumall = Sumall + sum(((y[ii+9*tr]
      - par[4]*par[4]*(CExp(m=abs(par[N+1]), Time=2*Thalf, x=x)))/err[ii+9*tr])^2)

  }
  if(N > 4) {
    # VV
    Sumall = Sumall + sum(((y[ii+10*tr]
      - par[5]*par[5]*(CExp(m=abs(par[N+1]), Time=2*Thalf, x=x)))/err[ii+10*tr])^2)
    Sumall = Sumall + sum(((y[ii+11*tr]
      - par[5]*par[6]*(CExp(m=abs(par[N+1]), Time=2*Thalf, x=x)))/err[ii+11*tr])^2)
    Sumall = Sumall + sum(((y[ii+12*tr]
      - par[6]*par[6]*(CExp(m=abs(par[N+1]), Time=2*Thalf, x=x)))/err[ii+12*tr])^2)    

    # 4V (sinh!)
    Sumall = Sumall + sum(((y[ii+13*tr]
      - par[1]*par[5]*(CExp(m=abs(par[N+1]), Time=2*Thalf, x=x, sign=-1.)))/err[ii+13*tr])^2)
    Sumall = Sumall + sum(((y[ii+14*tr]
      - par[1]*par[6]*(CExp(m=abs(par[N+1]), Time=2*Thalf, x=x, sign=-1.)))/err[ii+14*tr])^2)
    Sumall = Sumall + sum(((y[ii+15*tr]
      - par[2]*par[5]*(CExp(m=abs(par[N+1]), Time=2*Thalf, x=x, sign=-1.)))/err[ii+15*tr])^2)
    Sumall = Sumall + sum(((y[ii+16*tr]
      - par[2]*par[6]*(CExp(m=abs(par[N+1]), Time=2*Thalf, x=x, sign=-1.)))/err[ii+16*tr])^2)

    # AV cosh
    Sumall = Sumall + sum(((y[ii+17*tr]
      - par[3]*par[5]*(CExp(m=abs(par[N+1]), Time=2*Thalf, x=x)))/err[ii+17*tr])^2)
    Sumall = Sumall + sum(((y[ii+18*tr]
      - par[3]*par[6]*(CExp(m=abs(par[N+1]), Time=2*Thalf, x=x)))/err[ii+18*tr])^2)
    Sumall = Sumall + sum(((y[ii+19*tr]
      - par[4]*par[5]*(CExp(m=abs(par[N+1]), Time=2*Thalf, x=x)))/err[ii+19*tr])^2)
    Sumall = Sumall + sum(((y[ii+20*tr]
      - par[4]*par[6]*(CExp(m=abs(par[N+1]), Time=2*Thalf, x=x)))/err[ii+20*tr])^2)
  }
  return(Sumall)
}


ChiSqr.2mass <- function(par, Thalf, x, y, err, tr, N=2) {
  # index of mass
  l <- length(par)/2
  ii <- c(1:tr)
  m1 <- N+1
  m2 <- 2*N+2
  Sumall <- 0.
  if(N > 1) {
    Sumall = Sumall + sum(((y[ii]
      - par[1]*par[1]*(CExp(m=par[m1], Time=2*Thalf, x=x))
      - par[N+2]*par[N+2]*(CExp(m=par[m2], Time=2*Thalf, x=x)))/err[ii])^2)
    
    Sumall = Sumall + sum(((y[ii+tr]
      - par[1]*par[2]*(CExp(m=par[m1], Time=2*Thalf, x=x))
      - par[N+2]*par[N+3]*(CExp(m=par[m2], Time=2*Thalf, x=x)))/err[ii+tr])^2)

    Sumall = Sumall + sum(((y[ii+2*tr]
      - par[2]*par[2]*(CExp(m=par[m1], Time=2*Thalf, x=x))
      - par[N+3]*par[N+3]*(CExp(m=par[m2], Time=2*Thalf, x=x)))/err[ii+2*tr])^2)
    
  }
  if(N > 2) {
    # 4A (sinh!)
    Sumall = Sumall + sum(((y[ii+3*tr]
      - par[1]*par[3]*(CExp(m=abs(par[m1]), Time=2*Thalf, x=x, sign=-1.))
      - par[N+2]*par[N+4]*(CExp(m=abs(par[m2]), Time=2*Thalf, x=x, sign=-1.)))/err[ii+3*tr])^2)
    Sumall = Sumall + sum(((y[ii+4*tr]
      - par[1]*par[4]*(CExp(m=abs(par[m1]), Time=2*Thalf, x=x, sign=-1.))
      - par[N+2]*par[N+5]*(CExp(m=abs(par[m2]), Time=2*Thalf, x=x, sign=-1.)))/err[ii+4*tr])^2)
    Sumall = Sumall + sum(((y[ii+5*tr]
      - par[2]*par[3]*(CExp(m=abs(par[m1]), Time=2*Thalf, x=x, sign=-1.))
      - par[N+3]*par[N+4]*(CExp(m=abs(par[m2]), Time=2*Thalf, x=x, sign=-1.)))/err[ii+5*tr])^2)
    Sumall = Sumall + sum(((y[ii+6*tr]
      - par[2]*par[4]*(CExp(m=abs(par[m1]), Time=2*Thalf, x=x, sign=-1.))
      - par[N+3]*par[N+5]*(CExp(m=abs(par[m2]), Time=2*Thalf, x=x, sign=-1.)))/err[ii+6*tr])^2)
    
    # AA (no a1 in A)
    Sumall = Sumall + sum(((y[ii+7*tr]
      - par[3]*par[3]*(CExp(m=abs(par[m1]), Time=2*Thalf, x=x))
      - par[N+4]*par[N+4]*(CExp(m=abs(par[m2]), Time=2*Thalf, x=x)))/err[ii+7*tr])^2)
    Sumall = Sumall + sum(((y[ii+8*tr]
      - par[3]*par[4]*(CExp(m=abs(par[m1]), Time=2*Thalf, x=x))
      - par[N+4]*par[N+5]*(CExp(m=abs(par[m2]), Time=2*Thalf, x=x)))/err[ii+8*tr])^2)
    Sumall = Sumall + sum(((y[ii+9*tr]
      - par[4]*par[4]*(CExp(m=abs(par[m1]), Time=2*Thalf, x=x))
      - par[N+5]*par[N+5]*(CExp(m=abs(par[m2]), Time=2*Thalf, x=x)))/err[ii+9*tr])^2)

  }
  if(N > 4) {
    # VV
    Sumall = Sumall + sum(((y[ii+10*tr]
      - par[5]*par[5]*(CExp(m=abs(par[m1]), Time=2*Thalf, x=x))
      - par[N+6]*par[N+6]*(CExp(m=abs(par[m2]), Time=2*Thalf, x=x)))/err[ii+10*tr])^2)
    Sumall = Sumall + sum(((y[ii+11*tr]
      - par[5]*par[6]*(CExp(m=abs(par[m1]), Time=2*Thalf, x=x))
      - par[N+6]*par[N+7]*(CExp(m=abs(par[m2]), Time=2*Thalf, x=x)))/err[ii+11*tr])^2)
    Sumall = Sumall + sum(((y[ii+12*tr]
      - par[6]*par[6]*(CExp(m=abs(par[m1]), Time=2*Thalf, x=x))
      - par[N+7]*par[N+7]*(CExp(m=abs(par[m2]), Time=2*Thalf, x=x)))/err[ii+12*tr])^2)    

    # 4V (sinh!)
    Sumall = Sumall + sum(((y[ii+13*tr]
      - par[1]*par[5]*(CExp(m=abs(par[m1]), Time=2*Thalf, x=x, sign=-1.))
      - par[N+2]*par[N+6]*(CExp(m=abs(par[m2]), Time=2*Thalf, x=x, sign=-1.)))/err[ii+13*tr])^2)
    Sumall = Sumall + sum(((y[ii+14*tr]
      - par[1]*par[6]*(CExp(m=abs(par[m1]), Time=2*Thalf, x=x, sign=-1.))
      - par[N+2]*par[N+7]*(CExp(m=abs(par[m2]), Time=2*Thalf, x=x, sign=-1.)))/err[ii+14*tr])^2)
    Sumall = Sumall + sum(((y[ii+15*tr]
      - par[2]*par[5]*(CExp(m=abs(par[m1]), Time=2*Thalf, x=x, sign=-1.))
      - par[N+3]*par[N+6]*(CExp(m=abs(par[m2]), Time=2*Thalf, x=x, sign=-1.)))/err[ii+15*tr])^2)
    Sumall = Sumall + sum(((y[ii+16*tr]
      - par[2]*par[6]*(CExp(m=abs(par[m1]), Time=2*Thalf, x=x, sign=-1.))
      - par[N+3]*par[N+7]*(CExp(m=abs(par[m2]), Time=2*Thalf, x=x, sign=-1.)))/err[ii+16*tr])^2)

    # AV cosh
    Sumall = Sumall + sum(((y[ii+17*tr]
      - par[3]*par[5]*(CExp(m=abs(par[m1]), Time=2*Thalf, x=x))
      - par[N+4]*par[N+6]*(CExp(m=abs(par[m2]), Time=2*Thalf, x=x)))/err[ii+17*tr])^2)
    Sumall = Sumall + sum(((y[ii+18*tr]
      - par[3]*par[6]*(CExp(m=abs(par[m1]), Time=2*Thalf, x=x))
      - par[N+4]*par[N+7]*(CExp(m=abs(par[m2]), Time=2*Thalf, x=x)))/err[ii+18*tr])^2)
    Sumall = Sumall + sum(((y[ii+19*tr]
      - par[4]*par[5]*(CExp(m=abs(par[m1]), Time=2*Thalf, x=x))
      - par[N+5]*par[N+6]*(CExp(m=abs(par[m2]), Time=2*Thalf, x=x)))/err[ii+19*tr])^2)
    Sumall = Sumall + sum(((y[ii+20*tr]
      - par[4]*par[6]*(CExp(m=abs(par[m1]), Time=2*Thalf, x=x))
      - par[N+5]*par[N+7]*(CExp(m=abs(par[m2]), Time=2*Thalf, x=x)))/err[ii+20*tr])^2)
  }
  return(Sumall)
}
