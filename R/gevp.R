## t0 starts counting at 0!
##                           ( 1 , 2 )
## order c(1,2,3,4) goes to 
##                           ( 3 , 4 )

gevp <- function(cf, Time, t0=1, matrix.size=2, element.order=c(1,2,3,4),
                 for.tsboot=TRUE, sort.type="vectors") {
  
  if(matrix.size^2 != length(element.order)) {
    stop("matrix.size^2 must be Aborting...\n")
  }
  if(t0 < 0 || t0 > (Time/2-2)) {
    stop("t0 must be in between 0 and T/2-2. Aborting ...\n")
  }
  if(!any(c("values", "vectors") == sort.type)) {
    stop("possible values for sort.ype are values or vectors. Aborting\n")
  }
  Thalf <- Time/2
  if(length(dim(cf)) == 2) {
    Cor <- apply(cf, 2, mean)
  }
  else {
    Cor <- cf
  }

  ## need to check consistency of cf here!
  ## can only operate on a square matrix

  ## number of correlators in cf
  Ncor <- length(Cor)/(Thalf+1)
  if(length(element.order) != matrix.size^2) {
    stop("gevp can only operate on square matrices, please adjust element.order! Aborting!\n")
  }
  if(max(element.order) > Ncor) {
    stop("element.order tries to index beyond the available correlators in cf! Aborting...\n")
  }
  ## index array for indexing the linear data
  ii <- c()
  for(i in 1:Ncor) {
    ii <- c(ii, (i-1)*(Thalf+1)+1)
  }
  ## re-order as to match the input order
  ii <- ii[element.order]
  
  evalues <-  array(NA, dim=c(Thalf+1, matrix.size))
  evectors <- array(NA, dim=c(Thalf+1, matrix.size, matrix.size))
  amplitudes <- array(NA, dim=c(Thalf+1, matrix.size, matrix.size))
  ## matrix at t=t0 (ii takes care of the indices starting at 1 and not 0)
  ## and symmetrise
  cM <- 0.5*matrix(Cor[ii+t0], nrow=matrix.size, ncol=matrix.size)
  cM <- cM + t(cM)
  ## check for positive definiteness
  ev.cM <- eigen(cM, symmetric=TRUE, only.values = TRUE)
  if(any(ev.cM$values < 0)) {
    stop("gevp: matrix at t0 is not positive definite. Aborting...\n")
  }
  ## compute Cholesky factorisation
  L <- chol(cM)
  invL <- solve(L)
  
  ## now the time dependence for t > t0
  ## we need to multiply from the left with t(invL) and from the right with invL
  for(t in (t0+1):(Thalf)) {
    ## matrix at t and symmetrise
    cM <- 0.5*matrix(Cor[ii+t], nrow=matrix.size, ncol=matrix.size)
    cM <- cM + t(cM)
    ## determine eigenvalues and vectors
    variational.solve <- eigen(t(invL) %*% cM %*% invL,
                               symmetric=TRUE, only.values = FALSE, EISPACK=FALSE)
    ## sort depending on input by values or vectors
    sortindex <- c()
    if(sort.type == "values" || (t < t0+2)) {
      sortindex <- order(variational.solve$values, decreasing=TRUE)
    }
    else {
      ## compute the scalar product of eigenvectors with those at t0+1
      idx <- apply(abs( t(variational.solve$vectors) %*% evectors[t0+2,,] ), 1, order, decreasing=TRUE)
      sortindex <- idx[,1]
      if(!all(order(sortindex) == c(1:matrix.size))) {
        sortindex <- order(variational.solve$values, decreasing=TRUE)
      }
    }
    evalues[t+1,] <- variational.solve$values[sortindex]
    evectors[t+1,,] <- invL %*% variational.solve$vectors[, sortindex]
    tmp <- matrix(Cor[ii+t], nrow=matrix.size, ncol=matrix.size) %*% evectors[t+1,,]
    ## t(evectors[t+1,,]) %*% tmp must be proportional to delta_ij
    ## these are the amplitudes up to a factor sqrt(exp(-mt) \pm exp(-m(T-t)))
    ## diag(t(evectors[t+1,,]) %*% tmp) might get negative due to fluctuations
    ## we set them to NA first
    d <- diag(t(evectors[t+1,,]) %*% tmp)
    d[d < 0] <- NA
    amplitudes[t+1,,] <- t(t(tmp)/sqrt(d))
    rm(tmp)
  }
  evalues[t0+1,] <- 1.
  ## in case of bootstrapping everything (eigenvalues and eigenvectors)
  ## is concatenated into a single vector
  if(for.tsboot) {
    return(c(as.vector(evalues), as.vector(amplitudes), as.vector(evectors)))
  }
  
  return(invisible(list(evalues=evalues, evectors=evectors, amplitudes=amplitudes)))
}


bootstrap.gevp <- function(cf, t0=1, boot.R=400, boot.l=2, matrix.size=2,
                           element.order=c(1,2,3,4), seed=1234, sort.type="vectors") {
  ## number of measurements
  if(!any(class(cf) == "cf")) {
    stop("bootstrap.gevp requires an object of class cf as input! Aborting!\n")
  }
  if(!any(c("values", "vectors") == sort.type)) {
    stop("possible values for sort.ype are values or vectors. Aborting\n")
  }
  N <- length(cf$cf[,1])
  if(!cf$boot.samples) {
    cf <- bootstrap.cf(cf, boot.R=boot.R, boot.l=boot.l, seed=seed)
  }
  res <- gevp(cf$cf0, Time=cf$Time, t0, matrix.size, element.order, for.tsboot=FALSE, sort.type=sort.type)


  gevp.tsboot <- t(apply(cf$cf.tsboot$t, 1, gevp, Time=cf$Time, t0=t0,
                         matrix.size=matrix.size, element.order=element.order,
                         for.tsboot=TRUE, sort.type=sort.type))

  ## gevp.tsboot contains first the N*(Thalf+1) eigenvalues
  ## and the the N*N*(Thalf+1) eigenvectors
  
  ret <- list(cf=cf, res.gevp=res, gevp.tsboot=gevp.tsboot, boot.R=cf$boot.R, boot.l=cf$boot.l, matrix.size=matrix.size)
  class(ret) <- c("gevp", class(ret))
  return(invisible(ret))
}

gevp2cf <- function(gevp, id=1) {
  if(!inherits(gevp, "gevp")) {
    stop("gevp2cf requires an element of class gevp as input. Aborting...\n")
  }
  if(id > gevp$matrix.size || id < 1) {
    stop("gevp2cf: id must be <= matrix.size and > 0. Aborting...\n")
  }
  cf <- list()
  cf$N <- length(gevp$cf$cf[,1])
  cf$cf0 <- gevp$res.gevp$evalues[,id]
  cf$boot.samples <- TRUE
  cf$nrStypes <- 1
  cf$nrObs <- 1
  cf$Time <- gevp$cf$Time
  tt <- (id-1)*(cf$Time/2+1)+seq(1, cf$Time/2+1)
  cf$cf.tsboot <- list()
  cf$cf.tsboot$t <- gevp$gevp.tsboot[,tt]
  cf$cf.tsboot$t0 <- gevp$res.gevp$evalues[,id]
  cf$id <- id
  attr(cf, "class") <- c("cf", class(cf))
  return(invisible(cf))
}

gevp2amplitude <- function(gevp, mass, id=1, op.id=1, type="cosh", t1, t2, useCov=TRUE) {
  if(id > gevp$matrix.size || id < 1 || op.id > gevp$matrix.size || op.id < 1) {
    stop("gevp2cf: id and op.id must be <= matrix.size and > 0. Aborting...\n")
  }
  if(missing(t1) || missing(t2)) {
    if(inherits(mass, "effectivemassfit") || inherits(mass, "matrixfit")) {
      t1 <- mass$t1
      t2 <- mass$t2
    }
    else {
      stop("gevp2amplitude: t1 and/or t2 missing... Aborting\n")
    }
  }
  if((t2 <= t1) || (t1 < 0) || (t2 > (gevp$cf$Time/2-1))) {
    stop("gevp2amplitude: t1 < t2 and both in 0...T/2-1 is required. Aborting...\n")
  }
  
  sign <- +1
  if(type != "cosh") sign <- -1.
  if(!inherits(gevp, "gevp")) {
    stop("gevp2amplitude requires an element of class gevp as input. Aborting...\n")
  }
  if(is.numeric(mass)){
    if(length(mass) != gevp$boot.R) {
      stop(paste("gevp2amplitude: gevp and mass differ in number of bootstrap samples:", gevp$boot.R, "and", length(mass), ". Aborting...\n"))
    }
    m0 <- mean(mass)
    m <-  mass
  }
  else if(inherits(mass, "effectivemassfit") || inherits(mass, "matrixfit")) {
    if(gevp$boot.R != mass$boot.R) {
      stop(paste("gevp2amplitude: gevp and mass differ in number of bootstrap samples:", gevp$boot.R, "and", mass$boot.R, ". Aborting...\n"))
    }
    m0 <- mass$opt.res$par[1]
    m  <- mass$massfit.tsboot[,1]
  }
  else {
    stop("gevp2amplitude requires a numeric vector or an object either of type effectivemassfit or matrixfit as input. Abortgin...\n")
  }
  T <- gevp$cf$Time
  t <- c(0:(T/2))
  amplitude <- abs(gevp$res.gevp$amplitudes[,id,op.id])/sqrt(.5*(exp(-m0*t)+ sign*exp(-m0*(T-t))))
  tt <- gevp$matrix.size*(T/2+1) + ((id-1)*gevp$matrix.size+(op.id-1))*(T/2+1) + seq(1, T/2+1)
  amplitude.tsboot <- array(NA, dim=c(gevp$boot.R, T/2+1))
  for(i in c(1:gevp$boot.R)) {
    amplitude.tsboot[i,] <- abs(gevp$gevp.tsboot[i,tt])/sqrt(.5*(exp(-m[i]*t)+ sign*exp(-m[i]*(T-t))))
  }
  damplitude <- apply(amplitude.tsboot, 2, sd)
  
  ## now we perform a constant fit
  ii <- c((t1+1):(t2+1))

  M <- diag(1/damplitude[ii]^2)
  if(useCov) {
    ## compute correlation matrix and compute the correctly normalised inverse
    M <- invertCovMatrix(amplitude.tsboot[,ii], boot.samples=TRUE)
  }
  ## the chisqr function
  fn <- function(par, y, M) { sum((y-par[1]) %*% M %*% (y-par[1]))}
  
  par <- c(amplitude[t1+1])
  opt.res <- optim(par, fn = fn,
                   method="BFGS", M=M, y = amplitude[ii])
  opt.res <- optim(opt.res$par, fn = fn,
                   control=list(parscale=1/opt.res$par),
                   method="BFGS", M=M, y = amplitude[ii])
  meanAmplitude <- par[1]
  par <- opt.res$par

  meanAmplitude.tsboot <- array(0, dim=c(gevp$boot.R, 2))
  for(i in 1:gevp$boot.R) {
    opt <- optim(par, fn = fn,
                 control=list(parscale=1/par),
                 method="BFGS", M=M, y = amplitude.tsboot[i,ii])
    meanAmplitude.tsboot[i, 1] <- opt$par[1]
    meanAmplitude.tsboot[i, 2] <- opt$value
  }

  res <- list(amplitude=amplitude,
              amplitude.tsboot=amplitude.tsboot,
              damplitude=damplitude,
              meanAmplitude=meanAmplitude,
              meanAmplitude.tsboot=meanAmplitude.tsboot,
              chisqr = opt.res$value,
              dof=t2-t1, t1=t1, t2=t2,
              mass=mass, gevp=gevp,
              id=id, op.id=op.id,
              Time=T, m0=m0, m0.tsboot=m, useCov=useCov,
              Qval=1-pchisq(opt.res$value, t2-t1)
              )
  attr(res, "class") <- c("gevp.amplitude", class(res))
  return(invisible(res))
}

summary.gevp.amplitude <- function(amp) {
  cat("\n ** Result of a GEVP analysis for the amplitude **\n\n")
  cat("time range from", amp$t1, " to ", amp$t2, "\n")
  cat("mass:\n")
  cat("m \t=\t", amp$m0, "\n")
  cat("dm\t=\t", sd(amp$m0.tsboot), "\n")
  cat("\nAmplitude:\n")
  cat("operator id:", amp$op.id, "\n")
  cat("state id   :", amp$id, "\n")
  cat(" P[", amp$id, ",", amp$op.id, "] = ", amp$meanAmplitude, "\n")
  cat("dP[", amp$id, ",", amp$op.id, "] = ", sd(amp$meanAmplitude.tsboot[,1]), "\n")
  cat("\n")
  cat("boot.R\t=\t", amp$gevp$boot.R, "\n")
  cat("boot.l\t=\t", amp$gevp$boot.l, "\n")
  cat("useCov\t=\t", amp$useCov, "\n")
  cat("chisqr\t=\t", amp$chisqr, "\ndof\t=\t", amp$dof, "\nchisqr/dof=\t",
      amp$chisqr/amp$dof, "\n")
  cat("Quality of the fit (p-value):", amp$Qval, "\n")
  if(any(names(bla) == "fps")) {
    cat("\nDecay Constant (derived quantity):\n")
    cat("mu1 \t=\t", amp$mu1, "\n")
    cat("mu2 \t=\t", amp$mu2, "\n")
    if(amp$normalisation == "cmi") cat("kappa\t=\t", amp$kappa,"\n")
    cat("fps \t=\t", amp$fps, "\n")
    cat("dfps\t=\t", sd(amp$fps.tsboot), "\n")
  }
}

plot.gevp.amplitude <- function(amp, ...) {
  plotwitherror(c(0:(amp$Time/2)), amp$amplitude, amp$damplitude, ...)
  arrows(x0=amp$t1, y0=amp$meanAmplitude,
         x1=amp$t2, y1=amp$meanAmplitude, col=c("red"), length=0)
  arrows(x0=amp$t1, y0=amp$meanAmplitude+sd(amp$meanAmplitude.tsboot[,1]),
         x1=amp$t2, y1=amp$meanAmplitude+sd(amp$meanAmplitude.tsboot[,1]),
         col=c("red"), length=0, lwd=c(1))
  arrows(x0=amp$t1, y0=amp$meanAmplitude-sd(amp$meanAmplitude.tsboot[,1]),
         x1=amp$t2, y1=amp$meanAmplitude-sd(amp$meanAmplitude.tsboot[,1]),
         col=c("red"), length=0, lwd=c(1))
}
