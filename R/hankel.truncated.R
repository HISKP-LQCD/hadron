## Get the eigensystem of the full Hankel matrix and
## return the GEVP solution given a truncation dimension
spectrum.truncated.gevp <- function(ev.cM, n, deltat, submatrix.size,
                                    truncation.dim, error.weights, symmetric) {
  n.full <- n + deltat*submatrix.size
  ii0 <- 1:n
  ii.shift <- ii0 + deltat*submatrix.size
  ii1 <- rev((1:n.full)[order(abs(ev.cM$values))])[1:truncation.dim]

  M.0 <- ev.cM$vectors[ii0,ii1]
  M.t <- ev.cM$vectors[ii.shift,ii1]
  if(symmetric){
    chi <- sqrt(error.weights[ii0]^2 + error.weights[ii.shift]^2)
    M.bar <- M.0 + M.t
  } else {
    chi <- error.weights[ii.shift]
    M.bar <- M.0
  }
  M.bar <- chi * M.bar
  M.0 <- t(M.bar) %*% (chi * M.0)
  M.t <- t(M.bar) %*% (chi * M.t)
  M <- try(solve(M.0, M.t), TRUE)

  if(!inherits(M, "try-error")) {
    M.eigen <- try(eigen(M, symmetric=FALSE, only.values=TRUE), TRUE)
    if(!inherits(M.eigen, "try-error")) {
      return(invisible(M.eigen$values))
    } else {
      warning("eigen failed in gevp.truncated.hankel\n")
    }
  } else {
    warning("inversion failed in gevp.truncated.hankel\n")
  }

  return(invisible(rep(NA, truncation.dim)))
}

## Calculate the coefficients for reconstructing the correlators
## given the decay eigenvalues lambda = exp(-E deltat)
coeffs.truncated.gevp <- function(cf.mat, t0, deltat, Delta, lambda, submatrix.size,
                                  truncation.dim, error.weights) {
  times <- (t0 + (0:(nrow(cf.mat)-1))*Delta)/deltat
  t.min <- min(times)
  t.max <- max(times)
  vandermonde <- outer(times, lambda, function(t, a) ifelse(abs(a) < 1, a^(t-t.min), a^(t-t.max)))
  scale <- ifelse(abs(lambda) < 1, lambda^(-t.min), lambda^(-t.max))

  mat.coeffs <- sapply(1:submatrix.size^2, function(i) {
                         w <- error.weights[, i]^2
                         M.i <- Conj(t(vandermonde)) %*% (w * cf.mat[, i])
                         M.v <- Conj(t(vandermonde)) %*% (w * vandermonde)
                         M <- try(solve(M.v, M.i), TRUE)
                         if(!inherits(M, "try-error")) {
                           return(M)
                         } else {
                           warning("inversion failed in coeffs.truncated.gevp\n")
                         }
                         return(invisible(rep(NA, length(lambda))))
                                  })
  if(truncation.dim == 1) mat.coeffs <- t(mat.coeffs)

  return(invisible(mat.coeffs * scale))
}

## Approximate the coefficient matrix by a rank-1 matrix and return the vector
vec.coeffs.truncated.gevp <- function(mat.coeffs, truncation.dim, submatrix.size){
  if(truncation.dim == 1) mat.coeffs <- t(mat.coeffs)
  vec.coeffs <- apply(as.matrix(mat.coeffs), 1, function(M) {
                       M.mat <- matrix(M, nrow=submatrix.size)
                       M.mat <- 0.5 * (M.mat + Conj(t(M.mat)))
                       vec <- sqrt(diag(M.mat))
                       vec <- vec * M.mat[,1] / abs(M.mat[,1])
                       return(vec)
                                  })
  vec.coeffs <- t(vec.coeffs)

  return(invisible(vec.coeffs))
}

## Given the decay eigenvalues lambda = exp(-E deltat) and the coefficients
## reconstruct the full correlator matrix
reconstruct.correlators <- function(lambda, times, coeffs, lambda0=lambda){
  stopifnot(length(lambda) == length(lambda0))
  t.min <- min(times)
  t.max <- max(times)
  vandermonde <- outer(times, lambda, function(t, a) ifelse(abs(a) < 1, a^(t-t.min), a^(t-t.max)))
  scale <- ifelse(abs(lambda0) < 1, lambda0^(t.min), lambda0^(t.max))

  mat.coeffs <- coeffs * scale
  cor <- vandermonde %*% mat.coeffs

  return(invisible(cor))
}

#' @title GEVP method based on truncated Hankel matrices.
#' 
#' @description
#' Alternative method to determine energy levels from correlation
#'   matrices. A so-called Hankel matrix is generated from an input
#'   real numeric vector, truncated via SVD
#'   and a generalised eigenvalue problem is solved then.
#'
#' @param cf Numeric vector (this will generally be the time slices of a correlation function).
#' @param t0 Integer. Initial time value of the GEVP, must be in between 0 and
#'    \code{Time/2-2}. Default is 1.
#' @param deltat Integer. Time shift to be used to build the Hankel matrix.
#' @param n Integer. Size of the Hankel matrices to generate. This needs to include the factor of
#'   'submatrix.size'.
#' @param N Integer. Maximal time index in correlation function to be used in
#'                   Hankel matrix.
#' @param max.truncation Integer. Maximal truncation dimension to be used. Default is
#'   \code{n*submatrix.size}, the maximal possible value.
#' @param submatrix.size Integer. Submatrix size to be used in build
#'   of Hankel matrices. Submatrix.size > 1 is experimental.
#' @param element.order Integer vector. specifies how to fit the \code{n} linearly ordered single
#'    correlators into the correlator
#'    matrix for submatrix.size > 1. \code{element.order=c(1,2,3,4)} leads to a matrix
#'    \code{matrix(cf[element.order], nrow=2)}.
#'    Matrix elements can occur multiple times, such as \code{c(1,2,2,3)} for the symmetric case,
#'    for example.
#' @param Delta integer. Delta is the time shift used in the Hankel matrix.
#' @param get.coeffs boolean. If 'TRUE', correlator coefficients are also calculated,
#'    not only the decay spectrum.
#' @param effTime integer. Per default it is set to 'N'. It is only
#'   relevant for 'submatrix.size>1', and must contain the effective
#'   time extent of a single correlator, i.e. the spacing
#'   separating the different single correlator sequences in 'cf'.
#' @param error.weights boolean or numeric vector. If 'FALSE', no error weighting
#'   is applied. If 'TRUE', the inverse standard error of the correlator is used as weights.
#'   If a numeric vector is given, it must be of the same length as \code{cf$cf0} and
#'   contains the weights to be used.
#' @param symmetric boolean. If 'TRUE', the energy spectrum is guaranteed to be symmetric about 0.
#'   Default is \code{cf$symmetrised}.
#' @return
#' List object containing the following entries:
#' \item{cfii}{Integer vector. The time indices of the correlator used in the Hankel matrix.}
#' \item{spectrum}{Numeric matrix. The decay eigenvalues for truncation dimensions
#'   1 to \code{max.truncation}.}
#' \item{singular.values}{Numeric vector. The singular values of the full Hankel matrix,
#'   sorted by absolute value.}
#' \item{coefficients}{(if \code{get.coeffs=TRUE}) Numeric array. The correlator coefficients
#'   for truncation dimensions 1 to \code{max.truncation}.}
#' 
#' @family hankel
#' @export
gevp.truncated.hankel <- function(cf, t0=1, deltat=1, n, N, max.truncation=n,
                                  submatrix.size=1, element.order=c(1,2,3,4),
                                  Delta=1, get.coeffs=FALSE,
                                  effTime=N, error.weights=FALSE, symmetric=TRUE) {
  stopifnot((t0 >= 0) && (n > 0) && (N > 0) && (Delta > 0) && (submatrix.size > 0) && (max.truncation <= n))
  stopifnot((t0 + 1 + 2*(n/submatrix.size-1)*Delta + deltat) <= N)
  stopifnot(length(element.order) >= submatrix.size^2)
  stopifnot(length(error.weights) == 1 || length(error.weights) == length(cf))

  n.full <- n + deltat*submatrix.size
  hankel.dim <- n.full/submatrix.size

  cM0 <- array(NA, dim=c(n.full, n.full))
  ii <- seq(from=1, to=n.full, by=submatrix.size)

  t0p1 <- t0+1
  cfii <- t0p1 + (0:(2*hankel.dim - 2))*Delta

  for(i in c(1:submatrix.size)) {
    for(j in c(1:submatrix.size)) {
      cor.id <- element.order[(i-1)*submatrix.size+j]
      cM0[ii+(i-1), ii+(j-1)] <- hankel.matrix(n=hankel.dim, z=cf[cfii + (cor.id-1)*effTime])
    }
  }
  if(submatrix.size > 1) {
    ## symmetrise
    cM0 <- 0.5*(cM0 + t(cM0))
  }

  if(length(error.weights) == 1){
    error.weights <- matrix(1, nrow=length(cfii), ncol=submatrix.size^2)
    outer.weights <- rep(1, n.full)
    inner.weights <- outer.weights
  } else{
    error.mat <- matrix(error.weights, nrow=effTime)
    error.weights <- error.mat[cfii, element.order, drop=FALSE]
    error.mat.diag <- error.mat[, element.order[seq(1, submatrix.size^2, by=submatrix.size+1)], drop=FALSE]
    outer.weights <- c(t(error.mat.diag[cfii[(1:hankel.dim)*2 - 1],]))
    inner.weights <- outer.weights / c(sapply(1:hankel.dim, function(k){ rep((hankel.dim - abs(hankel.dim+1-2*k))^(1/4), submatrix.size) }))
  }

  cM0 <- t(inner.weights * t(inner.weights * cM0))
  ev.cM <- eigen(cM0, symmetric=TRUE, only.values = FALSE)
  ev.cM$vectors <- 1/inner.weights * ev.cM$vectors

  spectrum <- array(NA, dim=c(max.truncation, max.truncation))
  for(truncation.dim in 1:max.truncation){
    spectrum[truncation.dim, 1:truncation.dim] <-
      spectrum.truncated.gevp(ev.cM=ev.cM, n=n, deltat=deltat,
                              submatrix.size=submatrix.size,
                              truncation.dim=truncation.dim,
                              error.weights=outer.weights,
                              symmetric=symmetric)
  }

  res <- list(spectrum=spectrum, singular.values=rev(ev.cM$values[order(abs(ev.cM$values))]), cfii=cfii)

  if(get.coeffs) {
    coefficients <- array(NA, dim=c(max.truncation, max.truncation, submatrix.size^2))
    vec.coeffs <- array(NA, dim=c(max.truncation, max.truncation, submatrix.size))
    chi2 <- c()
    cf.mat <- matrix(cf, nrow=effTime)[cfii, element.order, drop=FALSE]
    for(truncation.dim in 1:max.truncation){
      coefficients[truncation.dim, 1:truncation.dim, ] <-
        coeffs.truncated.gevp(cf.mat, t0=t0, deltat=deltat, Delta=Delta,
                              lambda=spectrum[truncation.dim, 1:truncation.dim],
                              submatrix.size=submatrix.size,
                              truncation.dim=truncation.dim,
                              error.weights=error.weights)
      vec.coeffs[truncation.dim, 1:truncation.dim, ] <-
        vec.coeffs.truncated.gevp(coefficients[truncation.dim, 1:truncation.dim, ],
                                  truncation.dim=truncation.dim, submatrix.size=submatrix.size)
      cor.reconstructed <- reconstruct.correlators(lambda=spectrum[truncation.dim, 1:truncation.dim],
                                                   times=(cfii-1)/deltat,
                                                   coeffs=coefficients[truncation.dim, 1:truncation.dim,])
      chi <- (cf.mat - cor.reconstructed)*error.weights
      chi2[truncation.dim] <- sum(abs(chi)^2)
    }
    res$coefficients <- coefficients
    res$vec.coeffs <- vec.coeffs
    res$chi2 <- chi2
  }

  return(res)
}

#' @title Truncated PGEVM 
#' 
#' @description
#' Alternative method to determine energy levels from correlation
#'   matrices. A so-called Hankel matrix is generated from an input
#'   \link{cf} object, truncated via SVD and a generalised eigenvalue
#'   problem is solved then. This is the function to call.
#'   It will perform a bootstrap analysis. 
#'
#' @param cf object of type \link{cf}
#' @param deltat Integer. value of deltat used in the hankel GEVP. Default is 1. Used
#'   \code{t0fixed=FALSE}
#' @param Delta integer. Delta is the time shift used in the Hankel matrix.
#' @param N Integer. Maximal time index in correlation function to be used in
#'                   Hankel matrix.
#' @param t0 Integer. Initial time value of the GEVP, must be in between 0 and
#'    \code{Time/2-n}. Default is 1. Used when \code{t0fixed=TRUE}.
#' @param n Integer. Maximal size of the Hankel matrices to generate.
#'   Total Hankel matrix dimension will be \code{n*submatrix.size}.
#' @param submatrix.size Integer. Submatrix size to be used in build
#'   of Hankel matrices.
#' @param element.order Integer vector. specifies how to fit the \code{n} linearly ordered single
#'    correlators into the correlator
#'    matrix for \code{submatrix.size > 1}. E.g. \code{element.order=c(1,2,3,4)} leads to a matrix
#'    \code{matrix(cf[element.order], nrow=2)}.
#'    Matrix elements can occur multiple times, such as \code{c(1,2,2,3)} for the symmetric case.
#' @param max.truncation Integer. Maximal truncation dimension to be used. Default is
#'   \code{n*submatrix.size}, the maximal possible value.
#' @param error.weights boolean or numeric vector. If 'FALSE', no error weighting
#'   is applied. If 'TRUE', the inverse standard error of the correlator is used as weights.
#'   If a numeric vector is given, it must be of the same length as \code{cf$cf0} and
#'   contains the weights to be used.
#' @param symmetric boolean. If 'TRUE', the energy spectrum is guaranteed to be symmetric about 0.
#'   Default is \code{cf$symmetrised}.
#' @param bootstrap.coeffs boolean. If 'TRUE', the correlator coefficients are also
#'   calculated for each bootstrap sample, not only the original data.
#' @param eps numeric. Threshold for the singular value in the SVD to be considered for the
#'   proposed truncation dimension returned as \code{opt.idx}. Default is 1e-15.
#' 
#' @details
#' tbw
#'
#' @return
#' List object of classes "PGEVM" and "truncated.pgevm". The eigenvalues are stored in a
#' numeric vector \code{evs}, the corresponding samples in \code{evs.tsboot}.
#'
#' @family hankel
#' @export
bootstrap.truncated.pgevm <- function(cf, deltat=1, Delta=1, N = (cf$Time/2+1), t0 = 1,
                                      n = floor(((N - 1 - t0 - deltat)/Delta)/2 + 1),
                                      submatrix.size=1, element.order=1,
                                      max.truncation = n*submatrix.size, error.weights=FALSE, symmetric=cf$symmetrised,
                                      bootstrap.coeffs=FALSE, eps=1e-15) {
  stopifnot(inherits(cf, 'cf_meta'))
  stopifnot(inherits(cf, 'cf_boot'))
  dbboot <- inherits(cf, 'cf_dbboot')
  max.truncation <- min(max.truncation, n*submatrix.size)

  if(length(error.weights) == 1 & all(error.weights)) {
    error.weights <- 1/cf$tsboot.se
  }

  ## we need the inter-correlator spacing in 'cf'
  ## for gevp.truncated.hankel, because 'N' can be different now
  effTime <- cf$Time/2+1
  if(!cf$symmetrised) {
    effTime <- cf$Time
  }

  t0p1 <- t0 + 1
  boot.R <- cf$boot.R

  ## the last correlator element entering H(t+delta t) is
  ## C(t0+delta t + (2n-1)Delta)
  ## the last available element is N-1
  ## thus see n in argument list
  if(n < 1) n <- 1

  evs.dbboot <- array()
  coefficients.tsboot <- array()
  coefficients.dbboot <- array()
  vec.coeffs.tsboot <- array()
  vec.coeffs.dbboot <- array()
  dbboot.R <- c()
  if(dbboot) {
    dbboot.R <- cf$doubleboot$dbboot.R
  }

  evs <- gevp.truncated.hankel(cf$cf0, t0=t0, deltat=deltat, Delta=Delta, get.coeffs=TRUE,
                               n=n*submatrix.size, N=N, max.truncation=max.truncation,
                               submatrix.size=submatrix.size, element.order=element.order,
                               effTime=effTime, error.weights=error.weights, symmetric=symmetric)
  res.tsboot <- apply(cf$cf.tsboot$t, MARGIN=1L, FUN=function(cf0, ...) gevp.truncated.hankel(cf0, ...),
                      t0=t0, deltat=deltat, Delta=Delta, get.coeffs=bootstrap.coeffs,
                      n=n*submatrix.size, N=N, max.truncation=max.truncation,
                      submatrix.size=submatrix.size, element.order=element.order,
                      effTime=effTime, error.weights=error.weights, symmetric=symmetric)
  evs.tsboot <- array(t(sapply(res.tsboot, FUN=function(x) x$spectrum)),
                      dim=c(boot.R, max.truncation, max.truncation))
  if(bootstrap.coeffs){
    coefficients.tsboot <- array(t(sapply(res.tsboot, FUN=function(x) c(x$coefficients, x$vec.coeffs))),
                                 dim=c(boot.R, max.truncation, max.truncation, submatrix.size^2+submatrix.size))
    vec.coeffs.tsboot <- coefficients.tsboot[,,,(submatrix.size^2+1):(submatrix.size^2+submatrix.size), drop=FALSE]
    coefficients.tsboot <- coefficients.tsboot[,,,1:(submatrix.size^2), drop=FALSE]
  }
  if(dbboot) {
    res.dbboot <- apply(cf$doubleboot$cf, MARGIN=c(1L,2L), FUN=function(cf0, ...) gevp.truncated.hankel(cf0, ...),
                        t0=t0, deltat=deltat, Delta=Delta, get.coeffs=bootstrap.coeffs,
                        n=n*submatrix.size, N=N, max.truncation=max.truncation,
                        submatrix.size=submatrix.size, element.order=element.order,
                        effTime=effTime, error.weights=error.weights, symmetric=symmetric)
    evs.dbboot <- array(t(sapply(res.dbboot, FUN=function(x) x$spectrum)),
                        dim=c(boot.R, dbboot.R, max.truncation, max.truncation))
    if(bootstrap.coeffs){
      coefficients.dbboot <- array(t(sapply(res.dbboot, FUN=function(x) c(x$coefficients, x$vec.coeffs))),
                                   dim=c(boot.R, dbboot.R, max.truncation, max.truncation, submatrix.size^2+submatrix.size))
      vec.coeffs.dbboot <- coefficients.dbboot[,,,,(submatrix.size^2+1):(submatrix.size^2+submatrix.size), drop=FALSE]
      coefficients.dbboot <- coefficients.dbboot[,,,,1:(submatrix.size^2), drop=FALSE]
    }
  }

  truncation.error <- abs(evs$singular.values[-1] / evs$singular.values[-max.truncation])
  opt.idx <- which.min(truncation.error)
  opt.idx[2] <- min(which(evs$singular.values < eps) - 1, max.truncation)
  dof <- (2*(n+deltat)-1)*submatrix.size^2 - (submatrix.size+1)*(1:max.truncation)

  ret <- list(cf=cf,
              evs=evs$spectrum,
              evs.tsboot=evs.tsboot,
              evs.dbboot=evs.dbboot,
              singular.values=evs$singular.values,
              coefficients=evs$coefficients,
              vec.coeffs=evs$vec.coeffs,
              coefficients.tsboot=coefficients.tsboot,
              coefficients.dbboot=coefficients.dbboot,
              vec.coeffs.tsboot=vec.coeffs.tsboot,
              vec.coeffs.dbboot=vec.coeffs.dbboot,
              chi2=evs$chi2,
              cfii=evs$cfii,
              opt.idx=opt.idx,
              truncation.error=truncation.error,
              dof=dof,
              eps=eps,
              boot.R=boot.R,
              boot.l=cf$boot.l,
              seed=cf$seed,
              t0=t0,
              max.truncation=max.truncation,
              submatrix.size=submatrix.size,
              element.order=element.order,
              error.weights=error.weights,
              Delta=Delta,
              deltat=deltat,
              effTime=effTime,
              n=c(1:max.truncation),
              N=N)
  class(ret) <- c("PGEVM", "truncated.pgevm", class(ret))
  return(invisible(ret))
}

#' @title pgevm2bootstrapfit
#'
#' @param pgevm an object of class 'PGEVM' generated by 'bootstrap.truncated.pgevm'
#' @param truncation.dim integer. The truncation dimension to be used.
#'   Default is the most likely optimal truncation dimension `pgevm$opt.idx`.
#' @param errortype string. Determines the treatment of the bootstrap
#'   histograms to determine the statistical error on fit result. Can
#'   be: 1. 'outlier-removal' for which outliers are removed according to
#'   the 0.25 and 0.75 quantiles and the inter-quantile-range,
#'   i.e. only values are kept which are in the interval
#'   \eqn{[Q_25-1.5IQR, Q_75+1.5IQR]}
#'   and the error is computed from the standard deviation of the bootstrap distribution.
#'   2. 'std-dev' for which the error is estimated from the standard deviation.
#' @family hankel
#' @seealso input is generated via \link{bootstrap.truncated.pgevm}
#' See also \link{bootstrap.nlsfit}
#'
#' @return
#' Returns an object of S3 class `bootstrapfit`.
#' 
#' @export
pgevm2bootstrapfit <- function(pgevm, truncation.dim=pgevm$opt.idx[1], errortype="outlier-removal") {
  stopifnot(inherits(pgevm, "PGEVM"))

  if(errortype == "outlier-removal"){
    error.function <- function(x, probs=c(0.25,0.75), na.rm=TRUE) {
      Q <- quantile(x, probs=probs, na.rm=na.rm)
      iqr <- Q[2]-Q[1]
      x[x<(Q[1]-1.5*iqr) | x > (Q[2] + 1.5*iqr)] <- NA
      return(invisible(sd(x, na.rm=na.rm)))
    }
  }else{
    error.function <- sd
  }

  range <- 1:truncation.dim
  lambda <- pgevm$evs[truncation.dim, range]
  basic.res <- list(x=(0:(pgevm$N-1))/pgevm$deltat, boot.R=pgevm$boot.R, errormodel="yerrors",
                    par.guess=range, t0=lambda,
                    t=as.matrix(pgevm$evs.tsboot[, truncation.dim, range]),
                    useCov=FALSE, chisqr=pgevm$chi2[truncation.dim], dof=pgevm$dof[truncation.dim],
                    error.function=error.function, mask=pgevm$cfii,
                    tofn=list(coeffs=pgevm$coefficients[truncation.dim, range, ]))
  attr(basic.res, "class") <- c("bootstrapfit", "PGEVM", class(basic.res))

  bootstrap.coeffs <- any(!is.na(pgevm$coefficients.tsboot))
  if(bootstrap.coeffs){
    basic.res$tofn$coeffs.tsboot <- pgevm$coefficients.tsboot[, truncation.dim, range, , drop=FALSE]
  }else{
    basic.res$tofn$lambda0 <- lambda
  }

  res <- lapply(seq(pgevm$element.order), function(i) {
                  res <- basic.res
                  res$y <- pgevm$cf$cf0[1:pgevm$N + (pgevm$element.order[i]-1)*pgevm$effTime]
                  res$bsamples <- pgevm$cf$cf.tsboot$t[,1:pgevm$N + (pgevm$element.order[i]-1)*pgevm$effTime]
                  res$dy <- pgevm$cf$tsboot.se[1:pgevm$N + (pgevm$element.order[i]-1)*pgevm$effTime]
                  if(bootstrap.coeffs){
                    res$fn <- function(par, x, boot.r, coeffs, coeffs.tsboot, ...){
                      if(boot.r > 0) coeffs <- coeffs.tsboot[boot.r,,,]
                      return(Re(reconstruct.correlators(lambda=c(par), times=x, coeffs=coeffs)[,i]))
                    }
                  }else{
                    res$fn <- function(par, x, boot.r, ...) Re(reconstruct.correlators(lambda=c(par), times=x, ...)[,i])
                  }
                  return(invisible(res))
                    })
  return(invisible(res))
}

#' plot.truncated.pgevm
#'
#' @param x Object of type `truncated.pgevm`.
#' @param ... Graphical parameters to be passed on.
#'
#' @return
#' No return value.
#' 
#' @export
plot.truncated.pgevm <- function(x, ...) {
  par.old <- par(mfrow = c(1,2))
  plot(abs(x$singular.values), log="y", ylab="|s|")
  points(-Re(x$singular.values), col="red")
  legend("topright", legend=c("positive",  "negative", "opt. index 1", "opt. index 2"), col=c("black", "red", "black", "black"), pch=c(1, 1, NA, NA), lty=c(NA, NA, 2,3), bty="n")
  abline(v=x$opt.idx[1], lty=2)
  abline(v=x$opt.idx[2], lty=3)
  rcvalues <- cumsum(rev(abs(x$singular.values)))
  N <- length(rcvalues)
  totalsum <- rcvalues[N]
  rcvalues <- rcvalues/totalsum
  plot(rcvalues, log="y", ylab=c("cumsum(rev(|s|))/sum(|s|)"), xlab=c("rev(Index)"))
  abline(h=rcvalues[N-x$opt.idx[1]], lty=2)
  abline(h=rcvalues[N-x$opt.idx[2]], lty=3)
  legend("topleft", legend=c("opt. cut 1", "opt. cut 2"), lty=c(2,3), bty="n")
  par(par.old)
}
