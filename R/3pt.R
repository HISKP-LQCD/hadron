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

  ## now a constrained fit to the matrix
  matrixfit <- matrixfit(Cf2pt, t1=piont1-1, t2=piont2+1, symmetrise=TRUE, useCov=useCov,
                         matrix.size=1, parlist=array(c(1,1), dim=c(2,1)))

  ## which we can use to determine the 2pt correlator at T/2
  Thalfp1 <- Cf2pt$Time/2+1
  Amp <- matrixfit$opt.res$par[2]
  Mass <- matrixfit$opt.res$par[1]
  Cf2ptThalf <- 0.5*Amp^2*(exp(-Mass*(Cf2pt$Time-Thalfp1)) + exp(-Mass*Thalfp1))
  
  ## fit interval
  ii <- c((t1+1):(t2+1))
  ## error weights
  w <- 1/apply(Cf3pt$cf.tsboot$t[,ii], 2, sd)


  ## here we generate the inverse covariance matrix, if required
  ## otherwise take inverse errors squared
  M <- diag(w^2)

  if(useCov) {
    ## compute correlation matrix and compute the correctly normalised inverse
    ##M <- invertCovMatrix(Cf3pt$cf.tsboot$t[,ii], boot.samples=TRUE)
    M <- invertCovMatrix(Cf3pt$cf[,ii], boot.samples=FALSE, boot.l=boot.l)
  }
  fn <- function(par, y, M) { sum((y-par[1]) %*% M %*% (y-par[1]))}

  par <- Cf3pt$cf0[Cf2pt$Time/4]
  opt.res <- optim(par, fn = fn,
                   method="BFGS", M=M, y = Cf3pt$cf0[ii])
  opt.res <- optim(opt.res$par, fn = fn,
                   control=list(parscale=1/opt.res$par),
                   method="BFGS", M=M, y = Cf3pt$cf0[ii])
  par <- opt.res$par
  plateau <- par[1]
  plateau.tsboot <- array(NA, dim=c(boot.R,2))
  for(i in 1:boot.R) {
    opt <- optim(par, fn = fn,
                 control=list(parscale=1/par),
                 method="BFGS", M=M, y = Cf3pt$cf.tsboot$t[i,ii])
    plateau.tsboot[i,1] <- opt$par[1]
    plateau.tsboot[i,2] <- opt$value
  }
  ##  plateau <- weighted.mean(x=Cf3pt$cf0[ii], w=w)
  ##  plateau.tsboot <- apply(Cf3pt$cf.tsboot$t[,ii], 1, weighted.mean, w=w)


  averx <- plateau/effmass$opt.res$par[1]/Cf2pt$cf0[Thalfp1]
  averxfit <- plateau/matrixfit$opt.res$par[1]/Cf2ptThalf
  daverx <- sd(plateau.tsboot[,1]/effmass$massfit.tsboot[,1]/Cf2pt$cf.tsboot$t[,Thalfp1])
  daverxfit <- sd(plateau.tsboot[,1]/matrixfit$opt.tsboot[1,]/(0.5*matrixfit$opt.tsboot[2,]^2*(exp(-matrixfit$opt.tsboot[1,]*(Cf2pt$Time-Thalfp1)) + exp(-matrixfit$opt.tsboot[1,]*Thalfp1))))

  res <- list(averx=averx, daverx=daverx, plateau=plateau, plateau.tsboot=plateau.tsboot,
              effmass=effmass, Cf2pt=Cf2pt, Cf3pt=Cf3pt, matrixfit=matrixfit,
              t1=t1, t2=t2, piont1=piont1, piont2=piont2, chisqr=opt.res$value, dof=length(ii)-1,
              boot.R=boot.R, boot.l=boot.l, ii=ii, useCov=useCov,
              averxfit=averxfit, daverxfit=daverxfit, invCovMatrix=M)
  attr(res, "class") <- c("averx", "list")  
  return(invisible(res))
}

summary.averx <- function(averx) {
  summary(averx$effmass)
  cat("\n")
  summary(averx$matrixfit)
  cat("\nAnalysis for <x>\n\n")
  cat("based on", length(averx$Cf3pt$cf[,1]), "measurements\n")
  cat("correlated fit\t=\t", averx$useCov, "\n")
  cat("fitrange from", averx$t1, " to ", averx$t2, "\n")
  cat("chisqr\t=\t", averx$chisqr, "\n")
  cat("dof\t=\t", averx$dof, "\n")
  cat("chisqr/dof=\t",
      averx$chisqr/averx$dof, "\n")
  cat("Quality of the fit (p-value):",   1-pchisq(averx$chisqr, averx$dof), "\n\n")

  cat("<x>      =", averx$averx, "\n")
  cat("error    =", averx$daverx, "\n")
  cat("Alternative (using fitted Cf2pt(t/2) ):\n")
  cat("<x>      =", averx$averxfit, "\n")
  cat("error    =", averx$daverxfit, "\n")  
}

print.averx <- function(averx) {
  summary.averx(averx)
}

