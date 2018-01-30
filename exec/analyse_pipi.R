compRpipi <- function(c4, c2, Thalf) {
  tt <- c(1:Thalf)
  return((c4[tt] - c4[tt+1])/(c2[tt]^2 - c2[tt+1]^2))
}

## fit formel: R(t+1/2) = A*(cosh(dE*t') +sinh(dE*t')*coth(2*Mpi*t'))
## t' = t+1/2-T/2
Rfn <- function(par, tp, m) {
  return(par[1]*(cosh(par[2]*tp) + sinh(par[2]*tp)/tanh(2*m*tp)))
}


read.pipidata <- function(path="./", T, Rformat=FALSE) {
  if(!Rformat) {
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
  else {
    load(paste(path, "pion.cor.Rdata", sep="") )
    load(paste(path, "pipi.cor.Rdata", sep="") )
  }
  return(invisible(list(pion.cor=pion.cor, pipi.cor=pipi.cor)))
}

run.pipi.analysis <- function(data, boot.R=999, boot.l=1, useCov=TRUE,
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

  pion.cor <- bootstrap.cf(pion.cor, boot.R=boot.R, boot.l=boot.l)


  pion.effmass <- bootstrap.effectivemass(pion.cor, boot.R=boot.R, boot.l=boot.l, type="acosh")
  pion.effmass <- fit.effectivemass(pion.effmass, t1=t1, t2=t2, useCov=useCov)

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
    return(sum((z-y) %*% M %*% (z-y)))
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
  opt.res <- optim(par, fn = fitfn, method="BFGS", M=M, Thalf=Thalf, t=tphys, y=Rpipi[tt], m=pion.effmass$opt.res$par[1])
  opt.res <- optim(opt.res$par, fn = fitfn, method="BFGS", M=M, Thalf=Thalf, t=tphys, y=Rpipi[tt],
                   m=pion.effmass$opt.res$par[1], control=list(parscale=1./opt.res$par))

  ## a0 from deltaE and L
  ## using
  ## deltaE = - 4 pi a0 / m/L^3 (1 -2.837297 a0/L + 6.375183 a0^2/L^2)
  lfn <- function(a0, deltaE, L, m, debug=FALSE) {
    if(debug) cat("1/L^4:", -2.837297*a0/L, " 1/L^5:", 6.375183*a0^2/L^2, "\n")
    return(deltaE + 4*pi*a0/(m*L^3)*(1 - 2.837297*a0/L + 6.375183*a0^2/L^2))
    ##return(deltaE + 4*pi*a0/(m*L^3)*(1 - 2.837297*a0/L + 0*6.375183*a0^2/L^2))
    ##return(deltaE + 4*pi*a0/(m*L^3))
  }

  mpia0 <- pion.effmass$opt.res$par[1]*uniroot(lfn, c(-.1,-4.), deltaE=opt.res$par[2], L=L, m=pion.effmass$opt.res$par[1])$root
  par <- opt.res$par
  opt.tsboot <- array(NA, dim=c(boot.R,4))
  for(i in 1:boot.R) {
    opt <- optim(par, fn = fitfn,
                 t=tphys, Thalf=Thalf,
                 method="BFGS", M=M, y = Rpipi.tsboot[i,tt], m=pion.effmass$massfit.tsboot[i,1])
    opt <- optim(opt$par, fn = fitfn,
                 control=list(parscale=1/opt$par), t=tphys, Thalf=Thalf,
                 method="BFGS", M=M, y = Rpipi.tsboot[i,tt], m=pion.effmass$massfit.tsboot[i,1])
    opt.tsboot[i, c(1,2)] <- opt$par
    opt.tsboot[i, 3] <- opt$value
    opt.tsboot[i, 4] <- pion.effmass$massfit.tsboot[i,1]*uniroot(lfn, c(0.,-6.), deltaE=opt$par[2], L=L, m=pion.effmass$massfit.tsboot[i,1])$root
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
  data$chisq <- opt.res$value
  data$dof <- length(tt)-2
  data$qval <- 1-pchisq(data$chisq, data$dof)
  data$pion.effmass <- pion.effmass
  data$mpia0 <- mpia0
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
  qval <- pipi$qval
  cat("chi^2 = ", pipi$chisq, "\n")
  cat("dof = ", pipi$dof, "\n")
  cat("chi^2/dof = ", pipi$chisq/pipi$dof, "\n")
  cat("Qval = ", pipi$qval, "\n\n")

  cat("contributions by the orders in L relativ to L cubed\n")
  data$lfn(a0=pipi$mpia0/pipi$mpi, m=pipi$mpi, L=pipi$L, debug=TRUE)
  cat("\n")
  
  cat("deltaE = ", pipi$opt.res$par[2], "+-", sd(pipi$opt.tsboot[,2]), "(", 100*sd(pipi$opt.tsboot[,2])/pipi$opt.res$par[2], "%)", "\n")
  cat("mpi*a_0 = ", pipi$mpia0, "+-", sd(pipi$opt.tsboot[, 4]), "(", 100*sd(pipi$opt.tsboot[, 4])/abs(pipi$mpia0), "%)\n")


}

plot.pipi <- function(pipi, ...) {

  qqnorm(pipi$opt.tsboot[,1], main=c("qqnorm amplitude"))
  X11()
  qqnorm(pipi$opt.tsboot[,2], main=c("qqnorm deltaE"))
  X11()
  qqplot(pipi$opt.tsboot[,3], rchisq(n=pipi$boot.R, df=pipi$dof), main=c("qqplot chisq"))
  Thalf <- pipi$T/2
  ylim <- c(1.6,1.9)
  plotwitherror(c(1:(Thalf))-0.5, pipi$Rpipi, pipi$dRpipi, ylim=ylim, xlim=c(5,Thalf), ylab=c("R(t)"), xlab=c("t"), ...)
  
  x <- seq(pipi$tphys[1], pipi$tphys[length(pipi$tphys)]+1, 0.01)
  rat <- Rfn(pipi$opt.res$par, tp=x-Thalf, m=pipi$mpi)
  lines(x=x, y=rat, col="red")

  write.table(data.frame(x=x, rat=rat), file=paste(pipi$ens, ".ratiofit.dat", sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE)
  write.table(data.frame(t=c(1:Thalf)-0.5, pipi$Rpipi, pipi$dRpipi), file=paste(pipi$ens, ".ratio.dat", sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE)
  
  
  x <- seq(0, 0.06, 0.001)
  dE <- pipi$lfn(pipi$mpia0/pipi$mpi, 1/x, pipi$mpi)
  
  write.table(data.frame(x=x, dE=dE), file=paste(pipi$ens, ".dEovL.dat", sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE)


}
