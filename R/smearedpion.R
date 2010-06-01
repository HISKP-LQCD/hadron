smearedpion <- function(cmicor, mu1=0.0035, mu2=0.0035, kappa=0.161240, t1, t2, S=1.5,
                        skip=0, ind.vec=c(1,3,4,5), boot.R=99, boot.l=10, tsboot.sim="geom",
                        nrep=1, method="uwerr", debug=FALSE, par, fit.routine="optim", seed=123456) {
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
  nrOp <- 2
  if(missing(nrep)) {
    nrep <- c(length(cmicor[((Skip):Length),ind.vec[3]])/(nrObs*(T1)*nrOp))
  }
  else {
    skip <- 0
    if(sum(nrep) != length(cmicor[((Skip):Length),ind.vec[3]])/(nrObs*(T1)*nrOp)) {
      stop("sum of replica differs from total no of measurements!")
    }
  }
  Z <- array(cmicor[((Skip):Length),ind.vec[3]], 
             dim=c(nrObs*(T1)*nrOp,(length(cmicor[((Skip):Length),ind.vec[3]])/(nrObs*(T1)*nrOp))))
  ## negative times
  W <- array(cmicor[((Skip):Length),ind.vec[4]], 
             dim=c(nrObs*(T1)*nrOp,(length(cmicor[((Skip):Length),ind.vec[4]])/(nrObs*(T1)*nrOp))))

  rm(cmicor)
  ## averaged positivie and negative times and store in W
  W <- getCor(T1=T1, W=W, Z=Z, type=c("cosh", "cosh"))
  rm(Z)

  ## compute effective masses
  eff.sl <- effectivemass(from=(t1+1), to=(t2+1), Time, W[1:T1,] , pl=FALSE, S=1.5, nrep=nrep)
  eff.ss <- effectivemass(from=(t1+1), to=(t2+1), Time, W[(T1+1):(2*T1),] , pl=FALSE, S=1.5, nrep=nrep)
  options(show.error.messages = TRUE)

  mass.eff <- data.frame(t=eff.sl$t, mll=eff.sl$mass, dmll=eff.sl$dmass,
                         mlf=eff.ss$mass, dmlf=eff.ss$dmass)

  ## generate correlators plus errors
  Cor <- rep(0., times=nrOp*T1)
  E <- rep(0., times=nrOp*T1)
  
  for(i in 1:(nrOp*T1)) {
    Cor[i] <- mean(W[(i),])
    tmpe <- try(uwerrprimary(W[(i),], pl=F, nrep=nrep)$dvalue, TRUE)
    if(!inherits(tmpe, "try-error")) E[i] = tmpe
    else {
      warning("error of correlator replaced by naive estimate!\n", call.=F)
      E[i] = sd(W[(i),])/sqrt(length(W[(i),]))
    }
  }
  ## initial fit parameters from effective masses
  if(missing(par)) {
    par <- numeric(3)
    par[3] <- mass.eff$mll[1]
    par[2] <- sqrt(abs(Cor[(t1+1+T1)]*exp(par[3]*(t1+1))))
    par[1] <- abs(Cor[(t1+1)]*exp(mass.eff$mll[1]*(t1+1)))/par[2]
  }
  if(length(par) != 3) {
    par <- par[1:3]
  }
  ## indices in Cor and E for fitting
  ii <- c((t1p1):(t2p1), (t1p1+T1):(t2p1+T1))

  pionfit <- optim(par, ChiSqr.smeared, method="BFGS", control=list(trace=0, parscale=c(1.,1.,1.)), Thalf=Thalf,
                   x=c((t1):(t2)), y=Cor[ii], err=E[ii], tr = (t2-t1+1))
  if(fit.routine != "gsl") {
    pionfit <- optim(pionfit$par, ChiSqr.smeared, method="BFGS", control=list(trace=0, parscale=1/pionfit$par), Thalf=Thalf,
                     x=c((t1):(t2)), y=Cor[ii], err=E[ii], tr = (t2-t1+1))
  }
  else {
    pionfit <- gsl_fit_smeared_correlator(pionfit$par, Thalf=Thalf,
                                          x=c((t1):(t2)), y=Cor[ii], err=E[ii], tr = (t2-t1+1))
  }
  
  fit.mass <- abs(pionfit$par[3])
  fit.f <- 4.*kappa*(mu1+mu2)/2.*pionfit$par[1]/sqrt(fit.mass)^3/sqrt(2.)
  fit.sinhf <- 4.*kappa*(mu1+mu2)/2.*pionfit$par[1]/sqrt(fit.mass)/sinh(fit.mass)/sqrt(2.)
  fit.dof <- (t2-t1+1)*2-length(pionfit$par)
  fit.chisqr <- pionfit$value
  cat("mass =", abs(pionfit$par[3]), "fps =", fit.f, "fps(sinh) =", fit.sinhf, "chisqr/dof =", fit.chisqr,"/",fit.dof, "=", fit.chisqr/fit.dof, "\n")
  if(debug) {
    plot.effmass(m=fit.mass, ll=eff.ss, lf=eff.sl)
  }

  ## now do error analysis
  fit.uwerrm <- NULL
  fit.uwerrf <- NULL
  fit.boot <- NULL
  fit.tsboot <- NULL

  ## with UWERR method
  if(method == "uwerr" || method == "all") {
    fit.uwerrm <- uwerr(f=fitmasses.smeared, data=W[ii,], S=S, pl=debug, nrep=nrep,
                        Time=Time, t1=t1, t2=t2, Err=E[ii], par=pionfit$par,
                        fit.routine=fit.routine)

    fit.uwerrf <- uwerr(f=fitf.smeared, data=W[ii,], S=S, pl=debug, nrep=nrep,
                        Time=Time, t1=t1, t2=t2, Err=E[ii], par=pionfit$par,
                        fit.routine=fit.routine)
  }
  ## or bootstrap
  if(method == "boot" || method == "all") {
    set.seed(seed)
    fit.boot <- boot(data=t(W[ii,]), statistic=fit.smeared.boot, R=boot.R, stype="i",
                     Time=Time, t1=t1, t2=t2, Err=E[ii], par=pionfit$par,
                     kappa=kappa, mu1=mu1, mu2=mu2, fit.routine=fit.routine)
    
    fit.tsboot <- tsboot(tseries=t(W[ii,]), statistic=fit.smeared.boot, R=boot.R, l=boot.l, sim=tsboot.sim,
                         Time=Time, t1=t1, t2=t2, Err=E[ii], par=pionfit$par,
                         kappa=kappa, mu1=mu1, mu2=mu2, fit.routine=fit.routine)
  }

  ## compute Chi for each t value
  Chi <- rep(0., times=nrOp*T1)
  Fit <- rep(0., times=nrOp*T1)

  jj <-  c(t1p1:t2p1)
  Fit[jj] <- pionfit$par[1]*pionfit$par[2]*CExp(m=fit.mass[1], Time=2*Thalf, x=jj-1)
  Fit[jj+T1] <- pionfit$par[2]*pionfit$par[2]*CExp(m=fit.mass[1], Time=2*Thalf, x=jj-1)

  Chi[ii] <- (Fit[ii]-Cor[ii])/E[ii]
  ## return results to caller.
  res <- list(fitresult=pionfit, t1=t1, t2=t2, N=length(W[1,]), Time=Time, dof=fit.dof,
              fitdata=data.frame(t=(jj-1), Fit=Fit[ii], Cor=Cor[ii], Err=E[ii], Chi=Chi[ii]),
              uwerrresultmps=fit.uwerrm, uwerrresultfps=fit.uwerrf, 
              boot=fit.boot, tsboot=fit.tsboot, method=method,
              boot.R=boot.R, boot.l=boot.l, tsboot.sim=tsboot.sim,
              effmass=mass.eff, kappa=kappa, mu1=mu1, mu2=mu2, seed=seed,
              ##fit.routine=fit.routine, variational.masses=variational.masses, no.masses=no.masses, res.var=res.var
              matrix.size = 1, nrep=nrep)
  attr(res, "class") <- c("smearedfit", "list")  
  return(invisible(res))
}


fitmasses.smeared <- function(Cor, Err, t1, t2, Time, par=c(1.,0.1,0.12),
                              fit.routine="gsl") {
  Thalf <- Time/2
  T1 <- Thalf+1
  t1p1 <- (t1+1)
  t2p1 <- (t2+1)
  tr <- (t2-t1+1)

  fit <- optim(par*rnorm(length(par), mean=1., sd=0.001), ChiSqr.smeared, method="BFGS", Thalf=Thalf,
               control=list(trace=0, parscale=c(1,1,1)),
               x=c((t1):(t2)), y=Cor, err=Err, tr=tr)
  if(fit.routine != "gsl") {
    fit <- optim(fit$par, ChiSqr.smeared, method="BFGS", Thalf=Thalf,
                 control=list(trace=0, parscale=1./fit$par),
                 x=c((t1):(t2)), y=Cor, err=Err, tr=tr)
  }
  else {
    fit <- gsl_fit_smeared_correlator(fit$par, Thalf=Thalf,
                                      x=c((t1):(t2)), y=Cor, err=Err, tr = tr)
  }

  return(abs(fit$par[3]))
}

fitf.smeared <- function(Cor, Err, t1, t2, Time, par=c(1.,0.1,0.12),
                         fit.routine="gsl") {
  Thalf <- Time/2
  T1 <- Thalf+1
  t1p1 <- (t1+1)
  t2p1 <- (t2+1)
  tr <- (t2-t1+1)

  fit <- optim(par*rnorm(length(par), mean=1., sd=0.001), ChiSqr.smeared, method="BFGS", Thalf=Thalf,
               control=list(trace=0, parscale=c(1,1,1)),
               x=c((t1):(t2)), y=Cor, err=Err, tr=tr)
  if(fit.routine != "gsl") {
    fit <- optim(fit$par, ChiSqr.smeared, method="BFGS", Thalf=Thalf,
                 control=list(trace=0, parscale=1./fit$par),
                 x=c((t1):(t2)), y=Cor, err=Err, tr=tr)
  }
  else {
    fit <- gsl_fit_smeared_correlator(fit$par, Thalf=Thalf,
                                      x=c((t1):(t2)), y=Cor, err=Err, tr = tr)
  }

  return(abs(fit$par[1]/sqrt(abs(fit$par[3]))^3))
}

fit.smeared.boot <- function(Z, d, Err, t1, t2, Time, par=c(1.,0.1,0.12),
                             kappa, mu1, mu2,
                             fit.routine="gsl") {
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
  fit <- optim(par*rnorm(length(par), mean=1., sd=0.001), ChiSqr.smeared, method="BFGS", Thalf=Thalf,
               control=list(trace=0, parscale=c(1,1,1)),
               x=c((t1):(t2)), y=Cor, err=Err, tr=tr)
  if(fit.routine != "gsl") {
    fit <- optim(fit$par, ChiSqr.smeared, method="BFGS", Thalf=Thalf,
                 control=list(trace=0, parscale=1./fit$par),
                 x=c((t1):(t2)), y=Cor, err=Err, tr=tr)
    
  }
  else {
    fit <- gsl_fit_smeared_correlator(fit$par, Thalf=Thalf,
                                      x=c((t1):(t2)), y=Cor, err=Err, tr = tr)
  }

  fit.fpi <- 2*kappa*2*(mu1+mu2)/2./sqrt(2)*abs(fit$par[1])/sqrt(abs(fit$par[3])^3)
  return(c(abs(fit$par[3]), fit.fpi, fit$par[c(1:2)],
           fit$value))
}

summary.smearedfit <- function(fit) {

  kappa <- fit$kappa
  mu1 <- fit$mu1
  mu2 <- fit$mu2
  t1 <- fit$t1
  t2 <- fit$t2
  fit.mass <- abs(fit$fitresult$par[3])
  fit.fpi <- 2*kappa*2*(mu1+mu2)/2/sqrt(2)*abs(fit$fitresult$par[1])/sqrt(fit.mass^3)
  fit.sinhfpi <- 2*kappa*2*(mu1+mu2)/2/sqrt(2)*abs(fit$fitresult$par[1])/sqrt(fit.mass)/sinh(fit.mass)
  fit.chisqr <- fit$fitresult$value
  fit.dof <- length(fit$fitdata$t)-length(fit$fitresult$par)
  
  cat("mu1    = ", mu1, "\n")
  cat("mu2    = ", mu2, "\n")
  cat("kappa  = ", kappa, "\n")
  cat("Nr of measurements = ", fit$N, "\n")
  cat("No of replica = ", length(fit$nrep), "\n")
  cat("no or measurements per replicum: ", fit$nrep, "\n")
  cat("fitrange = ", t1, "-", t2, "\n")
  cat("chi^2    = ", fit.chisqr, "\n")
  cat("dof    = ", fit.dof, "\n")
  cat("chi^2/dof = ", fit.chisqr/fit.dof, "\n")
  
  cat("\nmps    = ", fit.mass, "\n", sep="\t")
  cat("fps    = ", fit.fpi, "\n", sep="\t")
  cat("fps    = ", fit.sinhfpi, " (from sinh definition)\n", sep="\t")
  cat("P_L    =", fit$fitresult$par[1], "\n", sep="\t")
  cat("P_S    =", fit$fitresult$par[2], "\n", sep="\t")

  if(!is.null(fit$uwerrresultmps)) {
    cat("\n--- Autocorrelation analysis for m_ps ---\n")
    cat("\nS        = ", fit$uwerrresultmps$S, "\n")
    cat("mps      = ", fit$uwerrresultmps$value, "\n")
    cat("dmps     = ", fit$uwerrresultmps$dvalue, "\n")
    cat("ddmps    = ", fit$uwerrresultmps$ddvalue, "\n")
    cat("tauint   = ", fit$uwerrresultmps$tauint, "\n")
    cat("dtauint  = ", fit$uwerrresultmps$dtauint, "\n")
    cat("Wopt     = ", fit$uwerrresultmps$Wopt, "\n")
    if(fit$uwerrresultmps$R>1) {
      cat("Qval     =", fit$uwerrresultmps$Qval, "\n")
    }
  }
  if(!is.null(fit$uwerrresultfps)) {
    cat("\n--- Autocorrelation analysis for f_ps ---\n")    
    cat("\nS        = ", fit$uwerrresultfps$S, "\n")
    cat("fps      = ", fit$uwerrresultfps$value*2*kappa*2*(mu1+mu2)/2./sqrt(2), "\n")
    cat("dfps     = ", fit$uwerrresultfps$dvalue*2*kappa*2*(mu1+mu2)/2./sqrt(2), "\n")
    cat("ddfps    = ", fit$uwerrresultfps$ddvalue*2*kappa*2*(mu1+mu2)/2./sqrt(2), "\n")
    cat("tauint   = ", fit$uwerrresultfps$tauint, "\n")
    cat("dtauint  = ", fit$uwerrresultfps$dtauint, "\n")
    cat("Wopt     = ", fit$uwerrresultfps$Wopt, "\n")
    if(fit$uwerrresultfps$R>1) {
      cat("Qval     =", fit$uwerrresultfps$Qval, "\n")
    }
  }

  if(!is.null(fit$boot)) {
    cat("--- Bootstrap analysis  ---\n")
    cat("---", fit$boot$R, "samples  ---\n")
    cat("          mean        -err           +err            stderr        bias\n")
    fit$boot.ci <- boot.ci(fit$boot, type = c("norm"), index=1)
    cat("mps    = ", fit$boot$t0[1], "(", (fit$boot.ci$normal[1,2]-fit$boot$t0[1])/1.96
	, ",", -(fit$boot$t0[1]-fit$boot.ci$normal[1,3])/1.96, ")", sd(fit$boot$t[,1]),
	mean(fit$boot$t[,1])-fit$boot$t0[1],"\n")
    
    fit$boot.ci <- boot.ci(fit$boot, type = c("norm"), index=2)
    cat("fps    = ", fit$boot$t0[2], "(", (fit$boot.ci$normal[1,2]-fit$boot$t0[2])/1.96
	, ",", -(fit$boot$t0[2]-fit$boot.ci$normal[1,3])/1.96, ")", sd(fit$boot$t[,2]),
	mean(fit$boot$t[,2])-fit$boot$t0[2], "\n")

    fit$boot.ci <- boot.ci(fit$boot, type = c("norm"), index=3)
    cat("P_L    = ", fit$boot$t0[3], "(", (fit$boot.ci$normal[1,2]-fit$boot$t0[3])/1.96
	, ",", -(fit$boot$t0[3]-fit$boot.ci$normal[1,3])/1.96, ")", sd(fit$boot$t[,3]),
	mean(fit$boot$t[,3])-fit$boot$t0[3], "\n")
    
    fit$boot.ci <- boot.ci(fit$boot, type = c("norm"), index=4)
    cat("P_S    = ", fit$boot$t0[4], "(", (fit$boot.ci$normal[1,2]-fit$boot$t0[4])/1.96
	, ",", -(fit$boot$t0[4]-fit$boot.ci$normal[1,3])/1.96, ")", sd(fit$boot$t[,4]),
	mean(fit$boot$t[,4])-fit$boot$t0[4],"\n")
  }
  if(!is.null(fit$tsboot)) {
    cat("\n--- Bootstrap analysis with blocking ---\n")
    cat("---", fit$boot$R, "samples  ---\n")
    cat("--- block size", fit$tsboot$l, "---\n")
    fit$tsboot.ci <- boot.ci(fit$tsboot, type = c("norm"), index=1)
    cat("mps    = ", fit$tsboot$t0[1], "(", (fit$tsboot.ci$normal[1,2]-fit$tsboot$t0[1])/1.96
	, ",", -(fit$tsboot$t0[1]-fit$tsboot.ci$normal[1,3])/1.96, ")", sd(fit$tsboot$t[,1]),
	mean(fit$tsboot$t[,1])-fit$tsboot$t0[1], "\n")
    
    fit$tsboot.ci <- boot.ci(fit$tsboot, type = c("norm"), index=2)
    cat("fps    = ", fit$tsboot$t0[2], "(", (fit$tsboot.ci$normal[1,2]-fit$tsboot$t0[2])/1.96
	, ",", -(fit$tsboot$t0[2]-fit$tsboot.ci$normal[1,3])/1.96, ")", sd(fit$tsboot$t[,2]),
	mean(fit$tsboot$t[,2])-fit$tsboot$t0[2], "\n")

    fit$tsboot.ci <- boot.ci(fit$tsboot, type = c("norm"), index=3)
    cat("P_L    = ", fit$tsboot$t0[3], "(", (fit$tsboot.ci$normal[1,2]-fit$tsboot$t0[3])/1.96
	, ",", -(fit$tsboot$t0[3]-fit$tsboot.ci$normal[1,3])/1.96, ")", sd(fit$tsboot$t[,3]),
	mean(fit$tsboot$t[,3])-fit$tsboot$t0[3], "\n")
    
    fit$tsboot.ci <- boot.ci(fit$tsboot, type = c("norm"), index=4)
    cat("P_F    = ", fit$tsboot$t0[4], "(", (fit$tsboot.ci$normal[1,2]-fit$tsboot$t0[4])/1.96
	, ",", -(fit$tsboot$t0[4]-fit$tsboot.ci$normal[1,3])/1.96, ")", sd(fit$tsboot$t[,4]),
	mean(fit$tsboot$t[,4])-fit$tsboot$t0[4], "\n")
  }
}
