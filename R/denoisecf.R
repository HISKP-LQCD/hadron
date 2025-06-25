hankelise <- function(X, N = dim(X)[1], sN = 1, verbose=FALSE) {
  stopifnot(sN > 0 & N %% sN == 0)
  n <- N/sN
  Y <- X
  if(sN == 1) {
    for(j in c(2:(2*N-1))) {
      Y[row(X) == j-col(X)] = mean(X[row(X) == j-col(X)])
    }
  }
  else {
##    for(i in c(2:(2*N-1))) {
##      for(k in c(0:(2*sN-1))) {
##        Y[(row(X) == i-col(X)) & ((((row(X)-1) %% sN) + ((col(X) -1) %% sN)) == k )] = mean(X[(row(X) == i-col(X)) & ((((row(X)-1) %% sN) + ((col(X) -1) %% sN)) == k )])
##      }
    ##    }
    ii <- seq(1, N, sN)-1
    for(i in c(1:sN)) {
      for(j in c(1:sN)) {
        Y[ii+i, ii+j] <- hankelise(X=Y[ii+i, ii+j], sN=1)
      }
    }
    for(i in c(1:sN)) {
      j <- i+1
      while(j <= sN) {
        Y[ii+i, ii+j] <- 0.5*(Y[ii+i, ii+j] + Y[ii+j, ii+i])
        Y[ii+j, ii+i] <- Y[ii+i, ii+j]
        j <- j+1
      }
    }
  }
  if(verbose) {
    cat("Hankel projection diff:", sum(abs(X-Y)), "\n")
  }
  return(Y)
}

projectPSD <- function(X, N = dim(X)[1], verbose=FALSE, cutNoise=TRUE) {
  ## symmetrise
  X <- 0.5*(X + t(X))
  X.eigen <- eigen(X, symmetric=TRUE)
  ii <- which(X.eigen$values < 0)
  if(cutNoise) ii <- which(X.eigen$values < abs(min(X.eigen$values)))
  nevsum <- sum(abs(X.eigen$values[ii]))
  if(verbose) {
    cat("PSD no removed evs:", length(ii), " min ", min(X.eigen$values), "max", max(X.eigen$values[ii]), "\n")
    cat("PSD neg. ev sum:", nevsum, "\n")
  }
  X.eigen$values[ii] <- 0
  Z <- X.eigen$vectors %*% diag(X.eigen$values) %*% t(X.eigen$vectors)
  if(verbose) {
    cat("PSD diff norm:", sqrt(sum((Z-X)^2)), "\n")
  }
  return(Z)
}

projectTopPSD <- function(X, N=dim(X)[1], sN=1, verbose=FALSE) {
  Z <- X
  Y <- projectPSD(X=X[c(1:(N-sN)), c((sN+1):N)], verbose=verbose)
  Z[c(1:(N-sN)), c((sN+1):N)] <- Y
  Z[c((sN+1):N),c(1:sN)] <- t(Y[c(1:sN),, drop=FALSE])
  Z[c((N-sN+1):N), c(1:(N-sN))] <- t(Y[,c((N-2*sN+1):(N-sN)), drop=FALSE])
  if(verbose) {
    cat("PSD top diff norm:", sqrt(sum((Z-X)^2)), "\n")
  }
  return(Z)
}

projectDensity <- function(X, N=dim(X)[1], D, sN=1, verbose=FALSE) {
  if(verbose) {
    cat("Density projection diff:", sum(abs(X[c(1:sN),c(1:sN)]-D)), "\n")
  }
  X[c(1:sN),c(1:sN)] <- D
  return(X)
}

Hankel2cf <- function(H, N=dim(H)[1], sN=1, element.order=c(1,2,3,4), cf.orig, Lcf, t0p1=1) {

  neff <- 2*N/sN-1
  ii <- seq(from=1, to=N, by=sN)-1
  cfii1 <- seq(from=t0p1, to=t0p1-1+N/sN, by=1)
  cfii2 <- seq(from=t0p1-1+N/sN, to=t0p1-1+2*N/sN-1, by=1)
  for(i in c(1:sN)) {
    for(j in c(1:sN)) {
      cor.id <- element.order[(i-1)*sN + j]
      ## take first and last block row of H
      cf.orig[cfii1 + (cor.id-1)*Lcf] <- H[i, ii+j]
      cf.orig[cfii2 + (cor.id-1)*Lcf] <- H[N-(sN-i), ii+j]
    }
  }
  return(cf.orig)
}

cf2Hankel <- function(cf, N, sN, t0p1=1, Lcf, element.order, symmetrise=TRUE) {
  H <- array(NA, dim=c(N, N))
  neff <- 2*N/sN-1
  cfii <- seq(from=t0p1, to=neff, by=1)
  if(sN == 1) {
    return(hankel.matrix(n=N, z=cf[cfii]))
  }
  ii <- seq(from=1, to=N, by=sN)
  for(i in c(1:sN)) {
    for(j in c(1:sN)) {
      cor.id <- element.order[(i-1)*sN + j]
      H[ii+i-1,ii+j-1] <- hankel.matrix(n=N/sN, z=cf[cfii + (cor.id-1)*Lcf])
    }
  }
  ## symmetrise
  if(symmetrise) H <- 0.5*(H + t(H))
  return(H)
}

isBlockHankel <- function(H, sN=1) {
  blockHankel <- TRUE
  N <- dim(H)[1]
  neff <- N/sN
  for(i in c(1:sN)) {
    for(j in c(i+1, N)) {
      if(j > N) break
      blockHankel <- abs(H)
    }
  }
  return(blockHankel)
}

dykstraIteration <- function(cf, N, sN=1, verbose=FALSE, tol=1.e-15, niter=10,
                             element.order=c(1,2,3,4), Lcf=length(cf), pmax=3,
                             cutNoise=TRUE, finalHankelisation=FALSE,
                             t0=0) {

  if(sN > 1 & is.null(Lcf)) {
    stop("In dykstraIteration: for sN>1, Lcf needs to be an integer\n")
  }
  stopifnot(length(element.order) >= sN^2)
  stopifnot(Lcf-t0*sN^2 >= 2*N/sN-1)
  t0p1 <- t0+1
  H <- cf2Hankel(cf=cf, N=N, sN=sN, t0p1=t0p1, Lcf=Lcf, element.order=element.order)
  if(verbose > 1) {
    cat("Testing Hankel2cf and cf2Hankel being inverse operations\n")
    tmp <- Hankel2cf(H=H, Lcf=Lcf, sN=sN, t0p1=t0p1, element.order=element.order, cf.orig=cf)
    cat(" -> sum of |cf - Hankel2cf(cf2Hankel(cf))| =", sum(abs(tmp-cf)), "\n")
  }
  ## build the zero timeslice matrix to fix
  D <- H[c(1:sN), c(1:sN)]

  X <- H
  Xtmp <- X
  Y <- array(0, dim=c(4, N, N))
  normsqr <- 10000
  m <- 0
  while(normsqr > tol & m < niter) {
    normsqrold <- normsqr
    normsqr <- 0
    for(p in c(1:pmax)) {
      Xtmp <- X
      if(p == 1) {
        X <- projectDensity(X=Xtmp-Y[p,,], D=D, sN=sN, verbose=(verbose > 1))
      }
      if(p == 2) {
        X <- hankelise(X=Xtmp-Y[p,,], verbose=(verbose > 1), sN=sN)
      }
      if(p == 3) {
        X <- projectPSD(X=Xtmp-Y[p,,], verbose=(verbose > 1), cutNoise=cutNoise)
      }
      if(p == 4) {
        X <- projectTopPSD(X=Xtmp-Y[p,,], sN=sN, verbose=(verbose > 1))
      }
      Ytmp <- Y[p,,]
      Y[p,,] <- X - (Xtmp - Ytmp)
      normsqr <- normsqr + sum((Ytmp-Y[p,,])^2)
      if(verbose > 1) cat(p, " ", normsqr, "\n")
    }
    m <- m+1
    if(verbose) cat(m, "tol", normsqr, "correction", sum((Xtmp-X)^2), "\n")
    if(m > niter & normsqr > 5*normsqrold) {
      X <- Xtmp
      break
    }
  }
  if(finalHankelisation) {
    X <- hankelise(X=X, verbose=(verbose > 1), sN=sN)
  }
  return(invisible(Hankel2cf(H=X, Lcf=Lcf, sN=sN, t0p1=t0p1, element.order=element.order, cf.orig=cf)))
}


#' @title denoise a correlation function
#'
#' @description
#'   t.b.w.
#' 
#' @param cf object of type \link{cf}, which needs to be bootstrapped before
#' @param n Integer. dimension of Hankel matrix to be build from
#'   2n-1 elements of cf. cf must have more than 2n-2 elements.
#' @param N Integer. Maximal time index in correlation function to be used in
#'                   Hankel matrix
#' 
#' @param tol numeric. tolerance threshold for denoising
#' @param niter integer. maximal number of Dykstra denoising iterations
#' @param submatrix.size Integer. Submatrix size to be used in build
#'   of Hankel matrices.
#' @param errortype string. Determines the treatment of the bootstrap
#'   histograms to determine the statistical error on eigenvalues. Can
#'   be: 1. 'outlier-removal' for which outliers are removed according to
#'   the 0.25 and 0.75 quantiles and the inter-quantile-range,
#'   i.e. only values are kept which are in the interval
#'   \eqn{[Q_25-1.5IQR, Q_75+1.5IQR]}
#'   and the error is computed from the standard deviation of the bootstrap distribution.
#'   2. 'quantiles' for which the error is estimated from the difference
#'   between the 0.32 and 0.68 quantile of the original bootstrap distribution
#'   3. 'dbboot' which works only, if the 'cf' is double bootstrapped. It will
#'   estimate the error from the true error of the median
#' @param element.order Integer vector. specifies how to fit the \code{n} linearly ordered single
#'    correlators into the correlator
#'    matrix for submatrix.size > 1. \code{element.order=c(1,2,3,4)} leads to a matrix
#'    \code{matrix(cf[element.order], nrow=2)}.
#'    Matrix elements can occur multiple times, such as \code{c(1,2,2,3)} for the symmetric case,
#'    for example.
#' @param verbose bool or integer. triggers verbose output in iteration on original data. 'TRUE' or '1'
#'    is the lowest level, 'verbose=2' triggers more verbose output.
#' @param cutNoise bool. if set to 'TRUE' also the smallest positive eigenvalues
#'   will be removed from the spectrum. In detail, all eigenvalues with magnitude
#'   smaller than 'abs(min(eigenvalues))' will be removed. If 'FALSE', only negative
#'   eigenvalues will be removed.
#' @param finalHankelisation bool. Apply one more 'hankelisation' as the final step
#'   before conversion back to cf. Per default the last step would be the projection to a
#'   positive semidefinite matrix.
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
denoise.cf <- function(cf, n, N = (cf$Time/2+1), tol=1.e-15, niter=10, errortype="dbboot", 
                       submatrix.size=1, element.order=c(1,2,3,4), verbose=FALSE,
                       cutNoise=TRUE, finalHankelisation=FALSE) {
  stopifnot(inherits(cf, 'cf_meta'))
  stopifnot(inherits(cf, 'cf_boot'))
  stopifnot(errortype %in% c("outlier-removal", "normal", "dbboot"))
  if(errortype == "dbboot" && !inherits(cf, 'cf_dbboot')) {
    stop("errortype dbboot requires double bootstrapped cf\n")
  }
  
  Nmax <- length(cf$cf0)
  Neff <- 2*n-1
  stopifnot(Nmax>Neff)

  cf$cf0 <- dykstraIteration(cf$cf0, N=n, sN=submatrix.size, element.order=element.order, tol=tol, niter=niter, verbose=verbose, Lcf=N, cutNoise=cutNoise, finalHankelisation=finalHankelisation)
  
  if(errortype=="dbboot" ) {
    cf$doubleboot$cf <- aperm(apply(X=cf$doubleboot$cf, MARGIN=c(1L, 2L), FUN=dykstraIteration, N=n, sN=submatrix.size, element.order=element.order, tol=tol, niter=niter, Lcf=N, cutNoise=cutNoise, finalHankelisation=finalHankelisation),
                              perm=c(2,3,1))
    cf$cf.tsboot$t <- apply(cf$doubleboot$cf, MARGIN=c(1L,3L), FUN=median, na.rm=TRUE)
  }
  else {
    cf$cf.tsboot$t <- t(apply(X=cf$cf.tsboot$t, MARGIN=1L, FUN=dykstraIteration, N=n, sN=submatrix.size, element.order=element.order, tol=tol, niter=niter, Lcf=N, cutNoise=cutNoise, finalHankelisation=finalHankelisation))
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
