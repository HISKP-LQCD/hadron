
compRpipi <- function(c4, c2, Thalf) {
  tt <- c(1:Thalf)
  return((c4[tt] - c4[tt+1])/(c2[tt]^2 - c2[tt+1]^2))
}

## fit formel: R(t+1/2) = A*(cosh(dE*t') +sinh(dE*t')*coth(2*Mpi*t'))
## t' = t+1/2-T/2
Rfn <- function(par, tp, m) {
  return(par[1]*(cosh(par[2]*tp) + sinh(par[2]*tp)/tanh(2*m*tp)))
}


read.pipidata <- function(path="./", T, Rformat=FALSE, Tformat=FALSE) {
  pion.cor <- cf()
  pipi.cor <- cf()
  if(Rformat) {
    load(paste(path, "pion.cor.Rdata", sep="") )
    load(paste(path, "pipi.cor.Rdata", sep="") )
  }
  else if(Tformat) {
    ind.vector <- c(2,2)
    pipi.cor <- readtextcf("pipi_pipi_A1_corr_TP0_00.dat", T=T, check.t=1, path=path, ind.vector=ind.vector)
    pion.cor <- readtextcf("pi_corr_p0.dat", T=T, check.t=1, path=path, ind.vector=ind.vector)
  }
  else {
    files <- Sys.glob( paste(path, "C2_pi+-*", sep="") )
    pion.cor <- readbinarycf(files, obs=0, T=T)
    files <- Sys.glob( paste(path, "C4_1*", sep="") )
    pipi1 <- readbinarycf(files, obs=0, T=T)
    files <- Sys.glob( paste(path, "C4_2*", sep="") )
    pipi2 <- readbinarycf(files, obs=0, T=T)
    files <- Sys.glob( paste(path, "C4_3*", sep="") )
    pipi3 <- readbinarycf(files, obs=0, T=T)
    pipi3 <- mul.cf(pipi3, a=2.)
    rm(files)
    save(pion.cor, file=paste(path, "pion.cor.Rdata", sep="") )
    save(pipi1, file=paste(path, "pipi1.Rdata", sep="") )
    save(pipi2, file=paste(path, "pipi2.Rdata", sep="") )
    save(pipi3, file=paste(path, "pipi3.Rdata", sep="") )
    pipi.cor <- pipi1 + pipi2 - pipi3
    save(pipi.cor, file=paste(path, "pipi.cor.Rdata", sep="") )
  }
  return(invisible(list(pion.cor=pion.cor, pipi.cor=pipi.cor)))
}

## Luescher Formula to Order 1/L^5
## a0 from deltaE and L
## using
## deltaE = - 4 pi a0 / m/L^3 (1 -2.837297 a0/L + 6.375183 a0^2/L^2)

deltaEvL <- function(a0, deltaE, L, m, debug=FALSE) {
  if(debug) cat("1/L^4:", -2.837297*a0/L, " 1/L^5:", 6.375183*a0^2/L^2, "\n")
  return(deltaE + 4*pi*a0/(m*L^3)*(1 - 2.837297*a0/L + 6.375183*a0^2/L^2))
  ##return(deltaE + 4*pi*a0/(m*L^3)*(1 - 2.837297*a0/L + 0*6.375183*a0^2/L^2))
  ##return(deltaE + 4*pi*a0/(m*L^3))
}


## extrac scattering length Luescher formula using
## the ratio to determine energy shift values
run.pipi.analysis.ratio <- function(data, boot.R=999, boot.l=1, useCov=TRUE,
                              t1=13, t2=31, tr1, tr2, L, ens="") {
  
  if(missing(t1) || missing(t2)) {
    stop("t1 and t2 must be specified, aborting...\n")
  }
  if(missing(tr1)) {
    tr1 <- t1
  }
  if(missing(tr2)) {
    tr2 <- t2
  }
  if(missing(data)) {
    stop("data is missing, aborting...\n")
  }
  pion.cor <- data$pion.cor
  pipi.cor <- data$pipi.cor

  Thalf <- pion.cor$Time/2
  if(missing(L)) {
    L <- Thalf
  }
  if(t2 > Thalf || tr2 > Thalf) {
    stop("t2 or tr2 not in range, aborting...\n")
  }

  pipi.cor <- bootstrap.cf(pipi.cor, boot.R=boot.R, boot.l=boot.l)

  if(!inherits(pion.cor, "effectivemassfit")) {
    pion.cor <- bootstrap.cf(pion.cor, boot.R=boot.R, boot.l=boot.l)
    pion.effmass <- bootstrap.effectivemass(pion.cor, boot.R=boot.R, boot.l=boot.l, type="solve")
    pion.effmass <- fit.effectivemass(pion.effmass, t1=t1, t2=t2, useCov=useCov)
  }
  
  ## this is R(t+1/2)
  Rpipi <- compRpipi(c4=pipi.cor$cf0, c2=pion.cor$cf0, Thalf=Thalf)

  Rpipi.tsboot <- array(NA, dim=c(boot.R, Thalf))
  for(i in c(1:boot.R)) {
    Rpipi.tsboot[i,] <- compRpipi(pipi.cor$cf.tsboot$t[i,], pion.cor$cf.tsboot$t[i,], Thalf=Thalf)
  }
  ## error of the ratio
  dRpipi <- apply(Rpipi.tsboot, 2, sd)

  ## fit function includes covariance matrix M
  fitfn <- function(par, y, t, m, Thalf, M) {
    tp <- t - Thalf 
    z <- Rfn(par, tp, m)
    return((z-y) %*% M %*% (z-y))
  }
  fitfn2 <- function(par, y, t, m, Thalf, L) {
    z <- Rfn(par, t-Thalf, m)
    return(L %*% (z-y))
  }
  ## fit range
  tt <- c(tr1:tr2)
  ## ratio is displaced by 1/2
  tphys <- tt-0.5
  ## build the covariance matrix from bootstrap samples, if wanted
  M <- diag(1./dRpipi[tt]^2)
  if(useCov) M <- invertCovMatrix(Rpipi.tsboot[,tt], boot.samples=TRUE)

  par <- c(1.6,0.01)
  ## fit on the original data
  lm.avail <- require(minpack.lm)

  if(lm.avail) {
    ML <- chol(M)
    opt.res <- nls.lm(par, fn=fitfn2,  L=ML, Thalf=Thalf, t=tphys, y=Rpipi[tt], m=pion.effmass$opt.res$par[1])
  }
  else {
    opt.res <- optim(par, fn = fitfn, method="BFGS", M=M, Thalf=Thalf, t=tphys, y=Rpipi[tt], m=pion.effmass$opt.res$par[1])
    opt.res <- optim(opt.res$par, fn = fitfn, method="BFGS", M=M, Thalf=Thalf, t=tphys, y=Rpipi[tt],
                     m=pion.effmass$opt.res$par[1], control=list(parscale=1./opt.res$par,ndeps=rep(1.e-6, times=length(par))))
  }
  data$chisq <- fitfn(opt.res$par, M=M, Thalf=Thalf, t=tphys, y=Rpipi[tt], m=pion.effmass$opt.res$par[1])

  a0 <- uniroot(deltaEvL, c(0.,-6.), deltaE=opt.res$par[2], L=L, m=pion.effmass$opt.res$par[1])$root
  mpia0 <- pion.effmass$opt.res$par[1]*a0

  par <- opt.res$par
  opt.tsboot <- array(NA, dim=c(boot.R, 5))
  for(i in 1:boot.R) {
    if(lm.avail) {
      opt <- nls.lm(par, fn=fitfn2,  L=ML, Thalf=Thalf, t=tphys, y=Rpipi.tsboot[i,tt], m=pion.effmass$massfit.tsboot[i,1])
    }
    else {
      opt <- optim(par, fn = fitfn,
                   t=tphys, Thalf=Thalf, control=list(ndeps=rep(1.e-6, times=length(par)), parscale=1/par),
                   method="BFGS", M=M, y = Rpipi.tsboot[i,tt], m=pion.effmass$massfit.tsboot[i,1])
    }
    opt.tsboot[i, c(1,2)] <- opt$par
    opt.tsboot[i, 3] <- fitfn(opt$par, t=tphys, Thalf=Thalf, M=M, y = Rpipi.tsboot[i,tt], m=pion.effmass$massfit.tsboot[i,1])
               
    if(opt$par[2] < 0) opt.tsboot[i, 5] <- NA
    else opt.tsboot[i, 5] <- uniroot(deltaEvL, c(0.,-6.), deltaE=opt$par[2], L=L, m=pion.effmass$massfit.tsboot[i,1])$root
    opt.tsboot[i, 4] <- pion.effmass$massfit.tsboot[i,1]*opt.tsboot[i, 5]
  }

  ## this is deltaE as a function of a0, m and L
  lfn2 <- function(a0, L, m, debug=FALSE) {
    if(debug) cat("1/L^4:", -2.837297*a0/L, " 1/L^5:", 6.375183*a0^2/L^2, "\n")
    return(- 4*pi*a0/(m*L^3)*(1 - 2.837297*a0/L + 6.375183*a0^2/L^2))
    ##return(deltaE + 4*pi*a0/(m*L^3)*(1 - 2.837297*a0/L + 0*6.375183*a0^2/L^2))
    ##return(deltaE + 4*pi*a0/(m*L^3))
  }
  

  data$pion.cor <- pion.cor
  data$pipi.cor <- pipi.cor
  data$Rpipi <- Rpipi
  data$dRpipi <- dRpipi
  data$opt.res <- opt.res
  data$opt.tsboot <- opt.tsboot
  data$t1 <- t1
  data$t2 <- t2
  data$useCov <- useCov
  data$tr1 <- tr1
  data$tr2 <- tr2
  data$boot.R <- boot.R
  data$boot.l <- boot.l
  data$L <- L
  data$T <- pion.cor$Time
  data$ens <- ens
#  data$chisq <- opt.res$value
  data$dof <- length(tt)-2
  data$Qval <- 1-pchisq(data$chisq, data$dof)
  data$pionQval <- pion.effmass$Qval
  data$pion.effmass <- pion.effmass
  data$mpia0 <- mpia0
  data$a0 <- a0
  data$mpi <- pion.effmass$opt.res$par[1]
  data$lfn <- lfn2
  data$tphys <- tphys
  attr(data, "class") <- c("pipi", class(data))
  return(invisible(data))
}

summary.pipi <- function(pipi) {
  summary(pipi$pion.effmass)

  cat("\n ** Luescher Analysis **\n\n")
  cat("fitrange from", pipi$tr1-0.5, "to", pipi$tr2-0.5, "\n")
  cat("correlated fit = ", pipi$useCov, "\n")
  chisq <- pipi$opt.res$value
  dof <- pipi$dof
  cat("chi^2 = ", pipi$chisq, "\n")
  cat("dof = ", pipi$dof, "\n")
  cat("chi^2/dof = ", pipi$chisq/pipi$dof, "\n")
  cat("Qval = ", pipi$Qval, "\n")
  cat("Qval of pion fit = ", pipi$pionQval, "\n\n")

#  cat("contributions by the orders in L relativ to L cubed\n")
#  data$lfn(a0=pipi$mpia0/pipi$mpi, m=pipi$mpi, L=pipi$L, debug=TRUE)
#  cat("\n")
  
  cat("deltaE = ", pipi$opt.res$par[2], "+-", sd(pipi$opt.tsboot[,2]), "(", 100*sd(pipi$opt.tsboot[,2])/pipi$opt.res$par[2], "%)", "\n")
  cat("a_0    = ", pipi$a0, "+-", sd(pipi$opt.tsboot[, 5], na.rm=TRUE), "(", 100*sd(pipi$opt.tsboot[, 5], na.rm=TRUE)/abs(pipi$a0), "%)\n")
  cat("mpi*a_0 = ", pipi$mpia0, "+-", sd(pipi$opt.tsboot[, 4], na.rm=TRUE), "(", 100*sd(pipi$opt.tsboot[, 4], na.rm=TRUE)/abs(pipi$mpia0), "%)\n")
}

plot.pipi <- function(pipi, ylim, ...) {

  plot(pipi$pion.effmass, ylim=c(0.9*pipi$pion.effmass$opt.res$par[1], 1.1*pipi$pion.effmass$opt.res$par[1]), xlab=c("t/a"), ylab=c("aMeff"), main=paste("pion effectivemass ", pipi$t1, "-", pipi$t2, " p=", pipi$pionQval, sep=""))
  if(interactive()) X11()
  qqnorm(pipi$opt.tsboot[,1], main=paste("qqnorm amplitude ", pipi$tr1, "-", pipi$tr2, " p=", pipi$Qval, sep=""))
  if(interactive()) X11()
  qqnorm(pipi$opt.tsboot[,2], main=paste("qqnorm deltaE ", pipi$tr1, "-", pipi$tr2, " p=", pipi$Qval, sep=""))
  if(interactive()) X11()
  qqplot(pipi$opt.tsboot[,3], rchisq(n=pipi$boot.R, df=pipi$dof), main=paste("qqplot chisq", pipi$tr1, "-", pipi$tr2, " p=", pipi$Qval, sep=""))
  Thalf <- pipi$T/2
  if(missing(ylim)) {
    rg <- range(pipi$Rpipi[c(5:Thalf)])
#    ylim <- c(rg[1]-0.1*(rg[2]-rg[1]), rg[2]+0.1*(rg[2]-rg[1]))
    ylim <- c(0.96*rg[1], 1.04*rg[2])
  }
  if(interactive()) X11()
  plotwitherror(c(1:(Thalf))-0.5, pipi$Rpipi, pipi$dRpipi, ylim=ylim, xlim=c(5,Thalf), ylab=c("R(t)"), xlab=c("t"), main=paste("Ratio", pipi$tr1, "-", pipi$tr2, " p=", pipi$Qval, sep=""), ...)
  
  x <- seq(pipi$tphys[1], pipi$tphys[length(pipi$tphys)]+1, 0.01)
  rat <- Rfn(pipi$opt.res$par, tp=x-Thalf, m=pipi$mpi)
  lines(x=x, y=rat, col="red")
  
  write.table(data.frame(x=x, rat=rat), file=paste(pipi$ens, ".", pipi$tr1, "-", pipi$tr2,".ratiofit.dat", sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE)
  write.table(data.frame(t=c(1:Thalf)-0.5, pipi$Rpipi, pipi$dRpipi), file=paste(pipi$ens, ".", pipi$tr1, "-", pipi$tr2, ".ratio.dat", sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE)
  
  
  x <- seq(0, 0.06, 0.001)
  dE <- pipi$lfn(pipi$mpia0/pipi$mpi, 1/x, pipi$mpi)
  
  write.table(data.frame(x=x, dE=dE), file=paste(pipi$ens, ".", pipi$tr1, "-", pipi$tr2, ".dEovL.dat", sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE)
}
