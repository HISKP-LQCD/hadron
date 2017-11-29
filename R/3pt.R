convert2cf <- function(data, symmetric=TRUE, symmetrise=TRUE) {
  sign <- +1.
  if(!symmetric) {
    sign <- -1.
  }
  Time <- max(data[[1]])+1
  X <- array(data[[2]], dim=c(Time, length(data[[2]])/Time))
  if(symmetrise){
    Thalf <- Time/2
    Thalfp1 <- Thalf+1
    i1 <- seq(2,Time/2)
    i2 <- seq(Time, Time/2+2)
    X[i1,] <- 0.5*(X[i1,] + sign*X[i2,])
    X <- X[c(1:Thalfp1),]
  }
  ret <- list(cf=t(X), icf=NULL, Time=Time, nrStypes=1, nrObs=1, boot.samples=FALSE, symmetrised=symmetrise)
  attr(ret, "class") <- c("cf", class(ret))
  return(invisible(ret))
}

pion_ff <- function(data3ptp0, data3ptp, data2ptp0, data2ptp,
                    boot.R=400, boot.l=2, t1.2pt, t2.2pt,
                    t1, t2, useCov=FALSE) {
    
  if(missing(data3ptp0) || missing(data3ptp) || missing(data2ptp0) || missing(data2ptp)) {
    stop("Error! Data is missing!")
  }
  if(missing(t1) || missing(t2)) {
    stop("Error! t1 and t2 must be specified!")
  }
  if(missing(t1.2pt)) {
    t1.2pt <- t1
  }
  if(missing(t2.2pt)) {
    t2.2pt <- t2
  }
  
  ## convert to cf
  Cf2ptp0 <- mul.cf(convert2cf(data2ptp0), 0.5)
  Cf2ptp <- mul.cf(convert2cf(data2ptp), 0.5)
  Cf3ptp0 <- convert2cf(data3ptp0, symmetric=FALSE)
  Cf3ptp <- convert2cf(data3ptp, symmetric=FALSE)

  ## we generate the appropriate ratios
  Cf3pt <- Cf3ptp/Cf3ptp0
  Cf2pt <- Cf2ptp0/Cf2ptp
  
  ## bootstrap the data
  Cf2ptp0 <- bootstrap.cf(Cf2ptp0, boot.R=boot.R, boot.l=boot.l)
  Cf3ptp0 <- bootstrap.cf(Cf3ptp0, boot.R=boot.R, boot.l=boot.l)
  Cf2ptp <- bootstrap.cf(Cf2ptp, boot.R=boot.R, boot.l=boot.l)
  Cf3ptp <- bootstrap.cf(Cf3ptp, boot.R=boot.R, boot.l=boot.l)

  ## and the ratios as well
  Cf2pt <- bootstrap.cf(Cf2pt, boot.R=boot.R, boot.l=boot.l)
  Cf3pt <- bootstrap.cf(Cf3pt, boot.R=boot.R, boot.l=boot.l)

  plateaufitZV <- fit.plateau2cf(Cf3ptp0, t1=t1, t2=t2, boot.samples=FALSE, boot.l=boot.l, boot.R=boot.R, useCov=useCov)
  plateaufitFF <- fit.plateau2cf(Cf3pt, t1=t1, t2=t2, boot.sample=FALSE, boot.l=boot.l, boot.R=boot.R, useCov=useCov)

  res <- list(Cf2ptratio=Cf2pt, Cf3ptratio=Cf3pt, Cf2ptp0=Cf2ptp0,
              Cf2ptp=Cf2ptp, Cf3ptp0=Cf3ptp0, Cf3ptp=Cf3ptp,
              plateaufitZV=plateaufitZV, plateaufitFF=plateaufitFF,
              boot.R=boot.R, boot.l=boot.l, t1=t1, t2=t2, useCov=useCov
              )

  attr(res, "class") <- c("pionff", "list")  
  return(invisible(res))
}

summary.pionff <- function(ff) {
  T <- ff$Cf2ptp0$Time
  Thalfp1 <- T/2+1

  cat("F(q^2) ", ff$plateaufitFF$plateau*ff$Cf2ptratio$cf0[Thalfp1], " ", sd(ff$plateaufitFF$plateau.tsboot[,1]*ff$Cf2ptratio$cf.tsboot$t[,Thalfp1]), "\n")
  cat("ZV     ", ff$Cf2ptp0$cf0[Thalfp1]/ff$plateaufitZV$plateau, " ", sd(ff$Cf2ptp0$cf.tsboot$t[,Thalfp1]/ff$plateaufitZV$plateau.tsboot[,1]), "\n")
}


averx <- function(data3pt, data2pt, pionfit, 
                  boot.R=400, boot.l=2, piont1, piont2, useCov=FALSE,
                  t1, t2, seed=123456, type="solve", force.optim=FALSE) {

  if(missing(data3pt) || missing(data2pt)) {
    stop("Error! Data is missing!")
  }
  if(missing(t1) || missing(t2)) {
    stop("Error! t1 and t2 must be specified!")
  }
  if(missing(piont1) || missing(piont2)) {
    stop("Error! piont1 and piont2 must be specified!")
  }

  Cf2pt <- data2pt
  ## convert to modern format if needed
  if(!inherits(Cf2pt, "cf")) {
    Cf2pt <- convert2cf(data2pt)
    ## bootstrap the data
    Cf2pt <- bootstrap.cf(Cf2pt, boot.R=boot.R, boot.l=boot.l, seed=seed)
  }

  Cf3pt <- data3pt
  if(!inherits(Cf3pt, "cf")) {
    Cf3pt <- convert2cf(data3pt)
    Cf3pt <- mul.cf(Cf3pt, a=-1.)
    
    Cf3pt <- bootstrap.cf(Cf3pt, boot.R=boot.R, boot.l=boot.l, seed=seed)
  }
  
  ## Determine the pion mass using effective masses
  ## effmass <- fit.effectivemass( bootstrap.effectivemass(Cf2pt, boot.R=boot.R, boot.l=boot.l, type=type), t1=piont1, t2=piont2-1, useCov=useCov)

  ## now a cosh fit to the 2pt correlator
  if(missing(pionfit)) {
    pionfit <- matrixfit(Cf2pt, t1=piont1, t2=piont2, symmetrise=TRUE, useCov=useCov,
                         matrix.size=1, parlist=array(c(1,1), dim=c(2,1)))
  }

  ## which we can use to determine the 2pt correlator at T/2
  Thalfp1 <- Cf2pt$Time/2+1
  Amp <- pionfit$opt.res$par[2]
  Mass <- pionfit$opt.res$par[1]
  Cf2ptThalf <- 0.5*Amp^2*(exp(-Mass*(Cf2pt$Time-Thalfp1)) + exp(-Mass*Thalfp1))
  
  ## fit interval
  ii <- c((t1+1):(t2+1))
  ## error weights
  w <- 1/apply(Cf3pt$cf.tsboot$t[,ii], 2, sd)

  ## here we generate the inverse covariance matrix, if required
  ## otherwise take inverse errors squared
  M <- diag(w^2)
  lm.avail <- require("minpack.lm")
  
  if(useCov) {
    ## compute correlation matrix and compute the correctly normalised inverse
    M <- invertCovMatrix(Cf3pt$cf[,ii], boot.samples=FALSE, boot.l=boot.l)
  }
  LM <- chol(M)

  fn <- function(par, y, M) { (y-par[1]) %*% M %*% (y-par[1])}
  fn.lm <- function(par, y, LM) {LM %*% (y-par[1])}
  
  par <- Cf3pt$cf0[Cf2pt$Time/4]
  if(lm.avail && !force.optim) {
    opt.res <- nls.lm(par, fn=fn.lm, LM=LM, y=Cf3pt$cf0[ii])
    chisq <- opt.res$rsstrace[length(opt.res$rsstrace)]
  }
  else {
    opt.res <- optim(par, fn = fn,
                     method="BFGS", M=M, y = Cf3pt$cf0[ii])
    opt.res <- optim(opt.res$par, fn = fn,
                     control=list(parscale=1/opt.res$par),
                     method="BFGS", M=M, y = Cf3pt$cf0[ii])
    chisq <- opt.res$value
  }
  par <- opt.res$par
  plateau <- par[1]
  plateau.tsboot <- array(NA, dim=c(boot.R,2))
  for(i in 1:boot.R) {
    if(lm.avail && !force.optim) {
      opt <- nls.lm(par, fn=fn.lm, LM=LM, y=Cf3pt$cf.tsboot$t[i,ii])
      plateau.tsboot[i,2] <- opt$rsstrace[length(opt$rsstrace)]
    }
    else {
      opt <- optim(par, fn = fn,
                   control=list(parscale=1/par),
                   method="BFGS", M=M, y = Cf3pt$cf.tsboot$t[i,ii])
      plateau.tsboot[i,2] <- opt$value
    }
    plateau.tsboot[i,1] <- opt$par[1]

  }
  plateau.wm <- weighted.mean(x=Cf3pt$cf0[ii], w=w)
  plateau.wm.tsboot <- apply(Cf3pt$cf.tsboot$t[,ii], 1, weighted.mean, w=w)


  averx <- plateau/pionfit$opt.res$par[1]/Cf2pt$cf0[Thalfp1]
  averxfit <- plateau/pionfit$opt.res$par[1]/Cf2ptThalf
  daverx <- sd(plateau.tsboot[,1]/pionfit$opt.tsboot[1,]/Cf2pt$cf.tsboot$t[,Thalfp1])
  daverxfit <- sd(plateau.tsboot[,1]/pionfit$opt.tsboot[1,]/(0.5*pionfit$opt.tsboot[2,]^2*(exp(-pionfit$opt.tsboot[1,]*(Cf2pt$Time-Thalfp1)) + exp(-pionfit$opt.tsboot[1,]*Thalfp1))))
  averx.wm <- plateau.wm/pionfit$opt.res$par[1]/Cf2pt$cf0[Thalfp1]
  daverx.wm <- sd(plateau.wm.tsboot/pionfit$opt.tsboot[1,]/Cf2pt$cf.tsboot$t[,Thalfp1])
  
  res <- list(averx=averx, daverx=daverx, plateau=plateau, plateau.tsboot=plateau.tsboot,
              Cf2pt=Cf2pt, Cf3pt=Cf3pt, matrixfit=pionfit,
              t1=t1, t2=t2, piont1=piont1, piont2=piont2, chisqr=chisq, dof=length(ii)-1,
              boot.R=boot.R, boot.l=boot.l, ii=ii, useCov=useCov,
              averxfit=averxfit, daverxfit=daverxfit, invCovMatrix=M,
              plateau.wm=plateau.wm, plateau.wm.tsboot=plateau.wm.tsboot,
              averx.wm=averx.wm, daverx.wm=daverx.wm, Qval=1-pchisq(chisq, length(ii)-1),
              pionQval = pionfit$Qval)
  attr(res, "class") <- c("averx", "list")  
  return(invisible(res))
}

summary.averx <- function(averx) {
  cat("\n")
  summary(averx$matrixfit)
  cat("\nAnalysis for <x>\n\n")
  cat("based on", length(averx$Cf3pt$cf[,1]), "measurements\n")
  cat("correlated fit\t=\t", averx$useCov, "\n")
  cat("fitrange from", averx$t1, " to ", averx$t2, "\n")
  cat("boot.R =", averx$boot.R, "\n")
  cat("boot.l =", averx$boot.l, "\n")
  cat("chisqr\t=\t", averx$chisqr, "\n")
  cat("dof\t=\t", averx$dof, "\n")
  cat("chisqr/dof=\t",
      averx$chisqr/averx$dof, "\n")
  cat("Quality of the plateau fit (p-value):", averx$Qval, "\n")
  cat("Quality of the pion fit (p-value):", averx$pionQval, "\n\n")

  cat("Using Cf2pt(T/2) data\n")
  cat("<x>      =", averx$averx, "\n")
  cat("error    =", averx$daverx, "\n")
  cat("Alternative (using fitted Cf2pt(T/2) ):\n")
  cat("<x>      =", averx$averxfit, "\n")
  cat("error    =", averx$daverxfit, "\n")  
  cat("Alternative (using weighted average over plateau)\n")
  cat("<x>      =", averx$averx.wm, "\n")
  cat("error    =", averx$daverx.wm, "\n")  
}

print.averx <- function(averx) {
  summary.averx(averx)
}

