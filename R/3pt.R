convert2cf <- function(data, symmetric=TRUE) {
  sign <- +1.
  if(!symmetric) {
    sign <- -1.
  }
  Time <- max(data[[1]])+1
  Thalf <- Time/2
  Thalfp1 <- Thalf+1
  i1 <- seq(2,Time/2)
  i2 <- seq(Time, Time/2+2)
  X <- array(data[[2]], dim=c(Time, length(data[[2]])/Time))
  X[i1,] <- 0.5*(X[i1,] + sign*X[i2,])
  X <- X[c(1:Thalfp1),]
  ret <- list(cf=t(X), icf=NULL, Time=Time, nrStypes=1, nrObs=1, boot.samples=FALSE)
  attr(ret, "class") <- c("cf", class(ret))
  return(invisible(ret))
}

averx <- function(data3pt, data2pt, 
                  boot.R=400, boot.l=2, piont1, piont2, useCov=FALSE,
                  t1, t2) {

  if(missing(data3pt) || missing(data2pt)) {
    stop("Error! Data is missing!")
  }
  if(missing(t1) || missing(t2)) {
    stop("Error! t1 and t2 must be specified!")
  }
  if(missing(piont1) || missing(piont2)) {
    stop("Error! piont1 and piont2 must be specified!")
  }

  ## convert to modern format
  Cf2pt <- convert2cf(data2pt)
  Cf3pt <- convert2cf(data3pt)
  Cf3pt <- mul.cf(Cf3pt, a=-1.)
  
  ## bootstrap the data
  Cf2pt <- bootstrap.cf(Cf2pt, boot.R=boot.R, boot.l=boot.l)
  Cf3pt <- bootstrap.cf(Cf3pt, boot.R=boot.R, boot.l=boot.l)
  
  
  ## Determine the pion mass
  effmass <- bootstrap.effectivemass(Cf2pt, boot.R=boot.R, boot.l=boot.l, type="acosh")
  effmass <- fit.effectivemass(effmass, t1=piont1, t2=piont2, useCov=useCov)


  ## fit interval
  ii <- c((t1+1):(t2+1))
  ## error weights
  w <- 1/apply(Cf3pt$cf.tsboot$t[,ii], 2, sd)
  plateau <- weighted.mean(x=Cf3pt$cf0[ii], w=w)
  plateau.tsboot <- apply(Cf3pt$cf.tsboot$t[,ii], 1, weighted.mean, w=w)
  Thalfp1 <- Cf2pt$Time/2+1
  averx <- plateau/effmass$opt.res$par[1]/Cf2pt$cf0[Thalfp1]
  daverx <- sd(plateau.tsboot/effmass$massfit.tsboot[,1]/Cf2pt$cf.tsboot$t[,Thalfp1])

  res <- list(averx=averx, daverx=daverx, plateau=plateau, plateau.tsboot=plateau.tsboot,
              effmass=effmass, Cf2pt=Cf2pt, Cf3pt=Cf3pt,
              t1=t1, t2=t2, piont1=piont1, piont2=piont2,
              boot.R=boot.R, boot.l=boot.l, ii=ii)
  attr(res, "class") <- c("averx", "list")  
  return(invisible(res))
}

summary.averx <- function(averx) {
  summary(averx$effmass)
  cat("\nAnalysis for <x>\n\n")
  cat("<x>      =", averx$averx, "\n")
  cat("based on", length(averx$Cf3pt$cf[,1]), "measurements\n")
  cat("error    =", averx$daverx, "\n")
}

print.averx <- function(averx) {
  summary.averx(averx)
}

plot.averx <- function(averx) {
  Thalfp1 <- averx$Cf2pt$time/2+1
  plot(averx$effmass, ylim=c(averx$effmass$opt.res$par[1]/2, 3/2*averx$effmass$opt.res$par[1]))
  X11()
  plot(averx$Cf3pt, xlab=c("t/a"), ylab=c("C3pt"), ylim=c(0, 2*averx$plateau))
  arrows(x0=averx$t1, y0=averx$plateau,
         x1=averx$t2, y1=averx$plateau, col=c("red"), length=0)
  arrows(x0=averx$t1, y0=averx$plateau+sd(averx$plateau.tsboot),
         x1=averx$t2, y1=averx$plateau+sd(averx$plateau.tsboot),
         col=c("red"), length=0, lwd=c(1))
  arrows(x0=averx$t1, y0=averx$plateau-sd(averx$plateau.tsboot),
         x1=averx$t2, y1=averx$plateau-sd(averx$plateau.tsboot),
         col=c("red"), length=0, lwd=c(1))
}

