cfeffectivemass <- function(cf, Thalf, type="solve") {
  if(length(dim(cf)) == 2) {
    Cor <- apply(cf, 2, mean)
  }
  else Cor <- cf
  effMass <- rep(0., (Thalf-1))
  if(type == "acosh") {
    for(t in c(2:(Thalf-1))) {
      effMass[t] <- acosh((Cor[t+1] + Cor[t-1])/2./Cor[t])
    }
  }
  else {
    ii <- c(1:(Thalf-1))
    iip1 <- ii+1
    Ratio <- Cor[ii]/Cor[iip1]
    if(type == "log") {
      effMass <- log(Ratio)
    }
    else {
      for(t in ii) {
        effMass[t] <- invcosh(Ratio[t], timeextent=2*Thalf, t=t)
      }
    }
  }
  return(invisible(effMass))
}

bootstrap.effectivemass <- function(cf, boot.R=400, boot.l=20, seed=12345, type="solve") {

  N <- length(cf$cf[,1])
  effMass <- cfeffectivemass(cf$cf, cf$Time/2, type=type)
  ## we set the seed for reproducability and correlation
  set.seed(seed)
  ## and bootstrap the fit
  effMass.tsboot <- tsboot(tseries=cf$cf, statistic=cfeffectivemass, R=boot.R, l=boot.l, sim="geom",
                           Thalf=cf$Time/2, type=type)
  deffMass=apply(effMass.tsboot$t, 2, sd, na.rm=TRUE)
  ret <- list(t=c(1:(cf$T/2-1)), effMass=effMass, deffMass=deffMass, effMass.tsboot=effMass.tsboot,
              opt.res=NULL, t1=NULL, t2=NULL, type=type, useCov=NULL, boot.R=boot.R, boot.l=boot.l,
              massfit.tsboot=NULL, Time=cf$Time, N=N)
  attr(ret, "class") <- c("effectivemass", class(ret))
  return(ret)
}

fit.effectivemass <- function(cf, t1, t2, useCov=FALSE) {
  if(missing(cf) || !("effectivemass" %in% class(cf))) {
    stop("cf is missing or must be of class \"effectivemass\"! Aborting...!\n")
  }
  if(missing(t1) || missing(t2)) {
    stop("t1 and t2 must be specified! Aborting...!\n")
  }
  cf$t1 <- t1
  cf$t2 <- t2
  cf$useCov <- useCov
  ii <- c((t1):(t2))
  attach(cf) 

  ## here we generate the inverse covariance matrix, if required
  ## otherwise take inverse errors squared
  M <- diag(1/deffMass[ii]^2)
  if(useCov) {
    ## compute correlation matrix and compute the correctly normalised inverse
    CovMatrix <- cov(effMass.tsboot$t[,ii])
    cov.svd <- svd(CovMatrix)
    ## replace smallest eigenvalues by their mean, if needed
    if(floor(sqrt(length(cov.svd$d))) < length(effMass.tsboot$t[,1])) {
      cov.svd$d[floor(sqrt(length(cov.svd$d))):length(cov.svd$d)] <-
        mean(cov.svd$d[floor(sqrt(length(cov.svd$d))):length(cov.svd$d)])
    }
    D <- diag(1/cov.svd$d)
    M <- cov.svd$v %*% D %*% t(cov.svd$u)
  }

  par <- c(effMass[t1])
  opt.res <- optim(par, fn = function(par, y, M) { (y-par[1]) %*% M %*% (y-par[1])},
                   method="BFGS", M=M, y = effMass[ii])
  par <- opt.res$par
  ## now we bootstrap the fit
  massfit.tsboot <- array(0, dim=c(boot.R, 2))
  for(i in 1:boot.R) {
    opt <- optim(par, fn = function(par, y, M) { (y-par[1]) %*% M %*% (y-par[1])},
                 method="BFGS", M=M, y = effMass.tsboot$t[i,ii])
    massfit.tsboot[i, 1] <- opt$par[1]
    massfit.tsboot[i, 2] <- opt$value
  }

  detach(cf)
  cf$massfit.tsboot <- massfit.tsboot
  cf$opt.res <- opt.res
  attr(cf, "class") <- c("effectivemassfit", class(cf))
  return(invisible(cf))
}

summary.effectivemass <- function(effMass) {
  attach(effMass)
  cat("\n ** effective mass values **\n\n")
  cat("no. measurements\t=\t", N, "\n")
  cat("boot.R\t=\t", boot.R, "\n")
  cat("boot.l\t=\t", boot.l, "\n")
  cat("Time extend\t=\t", Time, "\n")
  cat("values with errors:\n\n")
  print(data.frame(t= t, m = effMass, dm = deffMass))
  detach(effMass)
}

summary.effectivemassfit <- function(effMass, verbose=FALSE) {
  attach(effMass)
  cat("\n ** Result of effective mass analysis **\n\n")
  cat("no. measurements\t=\t", N, "\n")
  cat("type\t=\t", type, "\n")
  cat("boot.R\t=\t", boot.R, "\n")
  cat("boot.l\t=\t", boot.l, "\n")
  cat("Time extend\t=\t", Time, "\n")
  cat("time range from", t1, " to ", t2, "\n")
  if(verbose) {
    cat("values with errors:\n\n")
    print(data.frame(t= t, m = effMass, dm = deffMass))
  }
  cat("correlated fit\t=\t", useCov, "\n")
  cat("m\t=\t", opt.res$par[1], "\n")
  cat("dm\t=\t", sd(massfit.tsboot[,1]), "\n")
  cat("chisqr\t=\t", opt.res$value, "\n")
  cat("dof\t=\t", t2-t1, "\n")
  cat("chisqr/dof=\t",
      opt.res$value/(t2-t1), "\n")
  cat("Quality of the fit (p-value):",   1-pchisq(opt.res$value, t2-t1), "\n")

  detach(effMass)
}

print.effectivemassfit <- function(effMass, verbose=FALSE) {
  summary(effMass, verbose=verbose)
}

plot.effectivemass <- function(effMass, ref.value, ...) {
  op <- options()
  options(warn=-1)
  plotwitherror(effMass$t, effMass$effMass, effMass$deffMass, ...)
  options(op)
  if(!missing(ref.value)) {
    abline(h=ref.value, col=c("darkgreen"), lwd=c(3))
  }
  if(!is.null(effMass$opt.res)) {
    arrows(x0=effMass$t1, y0=effMass$opt.res$par[1], x1=effMass$t2, y1=effMass$opt.res$par[1], col=c("red"), length=0)
    arrows(x0=effMass$t1, y0=effMass$opt.res$par[1]+sd(effMass$massfit.tsboot[,1]),
           x1=effMass$t2, y1=effMass$opt.res$par[1]+sd(effMass$massfit.tsboot[,1]), col=c("red"), length=0, lwd=c(1))
    arrows(x0=effMass$t1, y0=effMass$opt.res$par[1]-sd(effMass$massfit.tsboot[,1]),
           x1=effMass$t2, y1=effMass$opt.res$par[1]-sd(effMass$massfit.tsboot[,1]), col=c("red"), length=0, lwd=c(1))
  }
}
