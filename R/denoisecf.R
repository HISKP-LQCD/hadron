projectHankel <- function(X, N = dim(X)[1], sN = 1, verbose=FALSE) {
  stopifnot(sN > 0 & N %% sN != 0)
  n <- N/sN
  Y <- X
  if(sN == 1) {
    for(j in c(2:(2*N-1))) {
      Y[row(X) == j-col(X)] = mean(X[row(X) == j-col(X)])
    }
  }
  else {
    ii <- seq(from=1, to=N, by=sN)
    for(i in c(1:sN)) {
      for(j in c(1:sN)) {
        Z <- X[ii+i-1, ii+j-1]
        for(k in c(2:(2*n-1))) {
          Z[row(Z) == k-col(Z)] = mean(Z[row(Z) == k-col(Z)])
        }
        Y[ii+i-1, ii+j-1] <- Z
      }
    }
    Y <- 0.5(Y + t(Y))
  }
  if(verbose) {
    cat("Hankel projection diff:", sum(abs(X-Y)), "\n")
  }
  return(Y)
}

projectPSD <- function(X, N = dim(X)[1], verbose=FALSE) {
  ## symmetrise
  X <- 0.5*(X + t(X))
  X.eigen <- eigen(X, symmetric=TRUE)
  ii <- which(X.eigen$values < 0)
  nevsum <- sum(X.eigen$values[ii])
  X.eigen$values[ii] <- 0
  Z <- X.eigen$vectors %*% diag(X.eigen$values) %*% t(X.eigen$vectors)
  if(verbose) {
    cat("PSD neg. ev sum:", nevsum, "\n")
    cat("PSD diff norm:", sqrt(sum((Z-X)^2)), "\n")
  }
  return(Z)
}

projectTopPSD <- function(X, N=dim(X)[1], verbose=FALSE) {
  Z <- X
  Xtop <- X[c(1:(N-1)), c(2:N)]
  Y <- projectPSD(Xtop, verbose=verbose)
  Z[c(1:(N-1)), c(2:N)] <- Y
  Z[c(2:N),1] <- Y[1,]
  Z[N, c(1:(N-1))] <- Y[,N-1]
  if(verbose) {
    cat("PSD top diff norm:", sqrt(sum((Z-X)^2)), "\n")
  }
  return(Z)
}

projectDensity <- function(X, N=dim(X)[1], density, verbose=FALSE) {
  if(verbose) {
    cat("Density projection diff:", sum(abs(X[1,1]-1+density) + abs(X[N,N]-density)), "\n")
  }
  X[1,1] <- density[1]
  if(length(density)>1) {
    X[N,N] <- density[2]
  }
  return(X)
}

Hankel2cf <- function(H, N=dim(H)[1]) {
  return( c(H[1,], H[N, c(2:N)]))
}

dykstraIteration <- function(cf, N, sN=1, verbose=FALSE, tol=1.e-15, niter=10,
                             densityIndex=c(1), element.order=c(1,2,3,4), Lcf=NULL) {

  dens <- cf[densityIndex]
  n <- length(cf)
  stopifnot(n >= 2*N-1)
  H <- array(NA, dim=c(N, N))
  neff <- 2*N-1
  if(sN == 1) {
    H <- hadron:::hankel.matrix(n=N, z=cf[1:neff])
  }
  else {
    for(i in c(1:sN)) {
      for(j in c(1:sN)) {
        cor.id <- element.order[(i-1)*sN + j]
        H[ii+i-1,ii+j-1] <- hadron:::hankel.matrix(n=N/sN, z=cf[1 + (cor.id-1)*Lcf])
      }
    }
    ## symmetrise
    H <- 0.5*(H + t(H))
  }
  X <- H
  Xtmp <- X
  Y <- array(0, dim=c(4, N, N))
  normsqr <- 10000
  m <- 0
  while(normsqr > tol & m < niter) {
    normsqr <- 0
    for(p in c(1:4)) {
      Xtmp <- X
      if(p == 1) {
        X <- projectDensity(X=Xtmp-Y[p,,], density=dens, verbose=verbose)
      }
      if(p == 2) {
        X <- projectHankel(X=Xtmp-Y[p,,], verbose=verbose, sN=sN)
      }
      if(p == 3) {
        X <- projectPSD(X=Xtmp-Y[p,,], verbose=verbose)
      }
      if(p == 4) {
        X <- projectTopPSD(X=Xtmp-Y[p,,], verbose=verbose)
      }
      Ytmp <- Y[p,,]
      Y[p,,] <- X - (Xtmp - Ytmp)
      normsqr <- normsqr + sum((Ytmp-Y[p,,])^2)
    }
    m <- m+1
    if(verbose) cat(m, "tol", normsqr, "correction", sum((Xtmp-X)^2), "\n\n")
  }
  if(n > neff) {
    return(invisible(c(Hankel2cf(X), rep(NA, times=n-neff))))
  }
  return(invisible(Hankel2cf(X)))
}


#' @title denoise a correlation function
#'
#' @description
#'   t.b.w.
#' 
#' @param cf object of type \link{cf}, which needs to be bootstrapped before
#' @param n Integer. dimension of Hankel matrix to be build from
#'   2n-1 elements of cf. cf must have more than 2n-2 elements.
#' 
#' @param tol numeric. tolerance threshold for denoising
#' @param niter integer. maximal number of Dykstra denoising iterations
#' @param densityIndex integer vector. the (time) indices of the elements in cf
#'   to fix during Dykstra iteration
#'
#' @references "Denoising of imaginary time response functions with Hankel projections"
#'       Yang Yu, Alexander F. Kemper, Chao Yang, Emanuel Gull,
#'       https://doi.org/10.1103/PhysRevResearch.6.L032042
#' @author Carsten Urbach, \email{curbach@@gmx.de}
#' 
#' @details
#'   t.b.w.
#' 
#' @return
#' object of type \link{cf} with original data and bootstrap data for the correlation
#' function denoised
#' 
#' @family hankel
#' @export
denoise.cf <- function(cf, n, tol=1.e-15, niter=10, densityIndex=c(1), errortype="dbboot") {
  stopifnot(inherits(cf, 'cf_meta'))
  stopifnot(inherits(cf, 'cf_boot'))
  stopifnot(errortype %in% c("outlier-removal", "normal", "dbboot"))
  if(errortype == "dbboot" && !inherits(cf, 'cf_dbboot')) {
    stop("errortype dbboot requires double bootstrapped cf\n")
  }
  
  Nmax <- length(cf$cf0)
  Neff <- 2*n-1
  stopifnot(Nmax>Neff)

  cf$cf0 <- dykstraIteration(cf$cf0, N=n, tol=tol, niter=niter, densityIndex=densityIndex)
  if(errortype=="dbboot" ) {
    cf$doubleboot$cf <- aperm(apply(X=cf$doubleboot$cf, MARGIN=c(1L, 2L), FUN=dykstraIteration, N=n, tol=tol, niter=niter, densityIndex=densityIndex),
                              perm=c(2,3,1))
    cf$cf.tsboot$t <- apply(cf$doubleboot$cf, MARGIN=c(1L,3L), FUN=median, na.rm=TRUE)
  }
  else {
    cf$cf.tsboot$t <- t(apply(X=cf$cf.tsboot$t, MARGIN=1L, FUN=dykstraIteration, N=n, tol=tol, niter=niter, densityIndex=densityIndex))
  }
  
  if(errortype == "outlier-removal") {
    remove_outliers <- function(x, probs=c(0.25,0.75)) {
      Q <- quantile(x, probs=probs, na.rm=TRUE)
      iqr <- Q[2]-Q[1]
      x[x<(Q[1]-1.5*iqr) | x > (Q[2] + 1.5*iqr)] <- NA
      return(invisible(x))
    }
    cf$cf.tsboot$t  <- apply(cf$cf.tsboot$t, MARGIN=2L, FUN=remove_outliers)
  }

  cf$tsboot.se <- apply(cf$cf.tsboot$t, MARGIN = 2L, FUN = cf$error_fn, na.rm=TRUE)
  cf$denoised <- TRUE
  
  return(invisible(cf))
}
