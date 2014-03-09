cfeffectivemass <- function(cf, Thalf, type="solve", nrObs=1) {
  ## cf is supposed to have nrObs*(Thalf+1) elements in time direction
  if(length(dim(cf)) == 2) {
    Cor <- apply(cf, 2, mean)
  }
  else {
    Cor <- cf
  }

  if(length(Cor) != nrObs*(Thalf+1)) {
    stop("cf does not have the correct time extend in cfeffectivemass! Aborting...!\n")
  }
  tt <- c(1:(nrObs*(Thalf+1)))
  cutii <- c()
  cutii2 <- c()
  for(i in 1:nrObs) {
    cutii <- c(cutii, (i-1)*(Thalf+1)+1, i*(Thalf+1))
    cutii2 <- c(cutii2, i*(Thalf+1))
  }
  t2 <- tt[-cutii2]
  
  effMass <- rep(NA, nrObs*(Thalf+1))
  if(type == "acosh") {
    t <- tt[-cutii]
    effMass[t] <- acosh((Cor[t+1] + Cor[t-1])/2./Cor[t])
  }
  else {
    t <- tt[c(1:(length(tt)-1))]
    Ratio <- Cor[t]/Cor[t+1]
    if(type == "log") {
      effMass[t2] <- log(Ratio[t2])
    }
    else {
      for(t in t2) {
        effMass[t] <- invcosh(Ratio[t], timeextent=2*Thalf, t=(t %% (Thalf+1)))
      }
    }
  }
  return(invisible(effMass[t2]))
}

bootstrap.effectivemass <- function(cf, boot.R=400, boot.l=20, seed=12345, type="solve") {

  if(!any(class(cf) == "cf")) {
    stop("bootstrap.effectivemass requires an object of class cf as input! Aborting!\n")
  }

  if(!cf$boot.samples) {
    cf <- bootstrap.cf(cf, boot.R=boot.R, boot.l=boot.l, seed=seed)
  }
  ## number of measurements
  N <- length(cf$cf[,1])
  ## number of time slices (hopefully in units of T/2+1)
  Nt <- length(cf$cf0)
  nrObs <- floor(Nt/(cf$Time/2+1))
  ## we run on the original data first
  effMass <- cfeffectivemass(cf$cf0, cf$Time/2, type=type, nrObs=nrObs)
  ## now we do the same on all samples
  effMass.tsboot <- t(apply(cf$cf.tsboot$t, 1, cfeffectivemass, cf$Time/2, type=type, nrObs=nrObs))

  deffMass=apply(effMass.tsboot, 2, sd, na.rm=TRUE)
  ret <- list(t=c(1:(cf$Time/2)),
              effMass=effMass, deffMass=deffMass, effMass.tsboot=effMass.tsboot,
              opt.res=NULL, t1=NULL, t2=NULL, type=type, useCov=NULL, invCovMatrix=NULL,
              boot.R=boot.R, boot.l=boot.l,
              massfit.tsboot=NULL, Time=cf$Time, N=N, nrObs=nrObs, dof=NULL)
  attr(ret, "class") <- c("effectivemass", class(ret))
  return(ret)
}

fit.effectivemass <- function(cf, t1, t2, useCov=FALSE, replace.na=TRUE) {
  if(missing(cf) || !any(class(cf) == "effectivemass" )) {
    stop("cf is missing or must be of class \"effectivemass\"! Aborting...!\n")
  }
  if(missing(t1) || missing(t2)) {
    stop("t1 and t2 must be specified! Aborting...!\n")
  }
  cf$t1 <- t1
  cf$t2 <- t2
  cf$useCov <- useCov
  cf$replace.na <- replace.na
  attach(cf)
  
  ## create an index array for the fit range
  ## the '+1' for Fortran index convention
  ## t1 and t2 can be in range 0-T/2
  ii <- c()
  for(i in 1:nrObs) {
    ii <- c(ii, ((i-1)*Time/2+t1+1):((i-1)*Time/2+t2+1))
  }
  ## get rid of the NAs for the fit
  ii.na <- which(is.na(effMass[ii]))
  ii <- ii[-ii.na]
  cf$ii <- ii
  cf$dof <-  length(ii)-1

  ## here we generate the inverse covariance matrix, if required
  ## otherwise take inverse errors squared
  M <- diag(1/deffMass[ii]^2)
  
  if(useCov) {
    ## compute correlation matrix and compute the correctly normalised inverse
    M <- invertCovMatrix(effMass.tsboot[,ii], boot.samples=TRUE)
  }

  cf$invCovMatrix <- M
  par <- c(effMass[t1])
  opt.res <- optim(par, fn = function(par, y, M) { sum((y-par[1]) %*% M %*% (y-par[1]))},
                   method="BFGS", M=M, y = effMass[ii])
  opt.res <- optim(opt.res$par, fn = function(par, y, M) { sum((y-par[1]) %*% M %*% (y-par[1]))},
                   control=list(parscale=1/opt.res$par),
                   method="BFGS", M=M, y = effMass[ii])
  par <- opt.res$par
  ## now we bootstrap the fit
  massfit.tsboot <- array(0, dim=c(boot.R, 2))

  ## now we replace all NAs by randomly chosen other bootstrap values
  tb.save <-  effMass.tsboot
  if(replace.na) {
    for(k in ii) {
      jj <- which(is.na(effMass.tsboot[,k]))
      rj <- sample.int(n=length(effMass.tsboot[jj, k]), size=length(jj), replace=FALSE)
      effMass.tsboot[jj, k] <- effMass.tsboot[rj, k]
    }
  }
  for(i in 1:boot.R) {
    opt <- optim(par, fn = function(par, y, M) { return(sum((y-par[1]) %*% M %*% (y-par[1])))},
                 control=list(parscale=1/par),
                 method="BFGS", M=M, y = effMass.tsboot[i,ii])
    massfit.tsboot[i, 1] <- opt$par[1]
    massfit.tsboot[i, 2] <- opt$value
  }
  effMass.tsboot <- tb.save
  rm(tb.save)
  detach(cf)
  cf$massfit.tsboot <- massfit.tsboot
  cf$opt.res <- opt.res
  attr(cf, "class") <- c("effectivemassfit", class(cf))
  return(invisible(cf))
}

summary.effectivemass <- function(effMass) {
  cat("\n ** effective mass values **\n\n")
  cat("no. measurements\t=\t", effMass$N, "\n")
  cat("boot.R\t=\t", effMass$boot.R, "\n")
  cat("boot.l\t=\t", effMass$boot.l, "\n")
  cat("Time extend\t=\t", effMass$Time, "\n")
  cat("total NA count in bootstrap samples:\t", length(which(is.na(effMass$effMass.tsboot))), "\n")
  cat("values with errors:\n\n")
  print(data.frame(t= effMass$t, m = effMass$effMass, dm = effMass$deffMass))
}

summary.effectivemassfit <- function(effMass, verbose=FALSE) {
  cat("\n ** Result of effective mass analysis **\n\n")
  cat("no. measurements\t=\t", effMass$N, "\n")
  cat("type\t=\t", effMass$type, "\n")
  cat("boot.R\t=\t", effMass$boot.R, "\n")
  cat("boot.l\t=\t", effMass$boot.l, "\n")
  cat("Time extend\t=\t", effMass$Time, "\n")
  cat("NA count in fitted bootstrap samples:\t", length(which(is.na(effMass$effMass.tsboot[,effMass$ii]))),
      "(",length(which(is.na(effMass$effMass.tsboot[,effMass$ii])))/ length(effMass$effMass.tsboot[,effMass$ii]), "%)\n")
  cat("NAs replaced in fit:", effMass$replace.na, "\n")
  cat("time range from", effMass$t1, " to ", effMass$t2, "\n")
  cat("No correlation functions", effMass$nrObs, "\n")
  if(verbose) {
    cat("values with errors:\n\n")
    print(data.frame(t= effMass$t, m = effMass$effMass, dm = effMass$deffMass))
  }
  cat("correlated fit\t=\t", effMass$useCov, "\n")
  cat("m\t=\t", effMass$opt.res$par[1], "\n")
  cat("dm\t=\t", sd(effMass$massfit.tsboot[,1]), "\n")
  cat("chisqr\t=\t", effMass$opt.res$value, "\n")
  cat("dof\t=\t", effMass$dof, "\n")
  cat("chisqr/dof=\t",
      effMass$opt.res$value/effMass$dof, "\n")
  cat("Quality of the fit (p-value):",   1-pchisq(effMass$opt.res$value, effMass$dof), "\n")

}

print.effectivemassfit <- function(effMass, verbose=FALSE) {
  summary(effMass, verbose=verbose)
}

plot.effectivemass <- function(effMass, ref.value, col,...) {
  if(missing(col)) {
    col <- c("black", palette(rainbow(max(effMass$nrObs,4))))
  }
  op <- options()
  options(warn=-1)
  t <- effMass$t
  plotwitherror(t, effMass$effMass[t], effMass$deffMass[t], ...)
  if(effMass$nrObs > 1) {
    for(i in 1:(effMass$nrObs-1)) {
      plotwitherror(t, effMass$effMass[t+i*effMass$Time/2], effMass$deffMass[t+i*effMass$Time/2], rep=TRUE, col=col[i+1], ...)
    }
  }
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
