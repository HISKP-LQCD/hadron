## creates a square Hankel matrix of size nxn
hankel.matrix <- function(n, z){
  outer(1:n,1:n,
        function(x,y) z[x+y-1])
}

#' @title GEVP method based on Hankel matrices.
#' 
#' @description
#' Alternative method to determine energy levels from correlation
#'   matrices. A so-called Hankel matrix is generated from an input
#'   real numeric vector and a generalised eigenvalue problem is solved
#'   then.
#'
#' @param cf Numeric vector (this will generally be the time slices of a correlation function).
#' @param t0values Integer vector. The t0 values to sum over.
#' @param deltat Integer. The value of the time shift to use to build the Hankel matrices.
#' @param n Integer. Size of the Hankel matrices to generate
#' @param N Integer. Maximal time index in correlation function to be used in
#'                   Hankel matrix
#'
#' @family hankel
#' @return
#' A complex vector of length \code{n + n^2} which contains the eigenvalues in the first
#' \code{n} elements and the eigenvectors in the remaining \code{n^2} elements.
#' 
#' A vector of NAs of \code{n + n^2} is returend in case the QR decomposition fails.
gevp.hankel_summed <- function(cf, t0values=c(1), deltat = 1, n, N) {

  cM1 <- array(0., dim=c(n, n))
  cM2 <- cM1

  ## build the two Hankel matrices
  for(i in c(1:length(t0values))) {
    cM1 <- cM1 + hankel.matrix(n=n, z=cf[(t0values[i]):(N)])
    cM2 <- cM2 + hankel.matrix(n=n, z=cf[(t0values[i]+deltat):(N)])
  }
  qr.cM1 <- qr(cM1)

  M <- try(qr.coef(qr.cM1, cM2), TRUE)
  if(inherits(M, "try-error")) {
    return(rep(NA, times=(n+n^2)))
  }

  M.eigen <- try(eigen(M, symmetric=FALSE, only.values=FALSE), TRUE)
  if(inherits(M.eigen, "try-error")) {
    return(rep(NA, times=(n+n^2)))
  }
  return(invisible(c(M.eigen$values, as.vector(M.eigen$vectors))))
}

#' @title GEVP method based on Hankel matrices. 
#' 
#' @description
#' Alternative method to determine energy levels from correlation
#'   matrices. A so-called Hankel matrix is generated from an input
#'   \link{cf} object and a generalised eigenvalue problem is solved
#'   then. This is the function to call. It will perform a bootstrap
#'   analysis. 
#'
#' @param cf object of type \link{cf}
#' @param t0values Integer vector. The t0 values to sum over. Default is \code{c(1:max)}.
#' All elements must be larger or equal to zero and smaller or equal than
#' \code{max=N-2*n-deltat}
#' @param deltat Integer. value of deltat used in the hankel GEVP. Default is 1.
#' 
#' @param n Integer. Size of the Hankel matrices to generate, default is 2.
#' @param N Integer. Maximal time index in correlation function to be used in
#'                   Hankel matrix
#'
#' @details
#' See `vignette(name="hankel", package="hadron")`
#'
#' @return
#' List object of class "hankel.summed". The eigenvalues are stored in a
#' numeric vector \code{t0}, the corresonding samples in \code{t}. The reference input
#' times \code{t0values} is stored as \code{t0values} in the returned list. In addition,
#' \code{deltat} is stored in the returned list.
#'
#' @family hankel
#' @export
#' @examples
#'
#' data(correlatormatrix)
#' correlatormatrix <- bootstrap.cf(correlatormatrix, boot.R=99, boot.l=1, seed=132435)
#' t0 <- 4
#' correlatormatrix.gevp <- bootstrap.gevp(cf=correlatormatrix, t0=t0, element.order=c(1,2,3,4))
#' pc1 <- gevp2cf(gevp=correlatormatrix.gevp, id=1)
#' pc1.hankel <- bootstrap.hankel_summed(cf=pc1, t0=c(1:15), n=2)
bootstrap.hankel_summed <- function(cf, t0values=c(1:(N-2*n-deltat)), deltat=1,
                                    n=2, N = cf$Time/2+1) {
  stopifnot(inherits(cf, 'cf_meta'))
  stopifnot(inherits(cf, 'cf_boot'))
  stopifnot(all(t0values>=0))
  stopifnot(max(t0values) <= N-2*n-deltat)
  
  boot.R <- cf$boot.R
  evs <- rep(NA, times=n + n^2)
  evs.tsboot <- array(NA, dim=c(boot.R, n + n^2))

  evs <- gevp.hankel_summed(cf$cf0, t0values=t0values,
                            deltat=deltat, n=n, N=N)
  evs.tsboot <- t(apply(cf$cf.tsboot$t, 1, gevp.hankel_summed,
                        n=n, N=N, deltat=deltat,
                        t0values=t0values))
                                
  ret <- list(cf=cf,
              t0=evs[c(1:n)],
              t=evs.tsboot[ ,c(1:n)],
              vectors=evs[c((n+1):(n+n^2))],
              vectors.tsboot=evs.tsboot[ ,c((n+1):(n+n^2))],
              boot.R=boot.R,
              boot.l=cf$boot.l,
              seed=cf$seed,
              t0values=t0values,
              deltat=deltat,
              n=n,
              N=N)
  class(ret) <- c("hankel_summed", class(ret))
  return(invisible(ret))
}

#' summary.hankel_summed
#'
#' @param object Object of type "hankel_summed" generated
#' by \link{bootstrap.hankel_summed}
#' @param ... Generic parameters to pass on.
#'
#' @return
#' Returns invisibly a data frame with columns `E`, `dE` (energies and their standard errors)
#' `q16` and `q84` the 16 and 84 percent quantiles.
#' 
#' @export
summary.hankel_summed <- function(object, ...) {
  X <- data.frame(E   =-log(object$t0)/object$deltat,
                  dE  =apply(-log(object$t)/object$deltat, 2, FUN=sd, na.rm=TRUE),
                  q16 =apply(-log(object$t)/object$deltat, 2, FUN=quantile, na.rm=TRUE, probs=0.159),
                  q84 =apply(-log(object$t)/object$deltat, 2, FUN=quantile, na.rm=TRUE, probs=0.841))
  cat("Summed Hankel analysis\n")
  cat("t0 values", object$t0values, "\n")
  print(X)
  return(invisible(X))
}

#' @title GEVP method based on Hankel matrices.
#' 
#' @description
#' Alternative method to determine energy levels from correlation
#'   matrices. A so-called Hankel matrix is generated from an input
#'   real numeric vector and a generalised eigenvalue problem is solved
#'   then.
#'
#' @param cf Numeric vector (this will generally be the time slices of a correlation function).
#' @param t0 Integer. Initial time value of the GEVP, must be in between 0 and
#'    \code{Time/2-2}. Default is 1.
#' @param n Integer. Size of the Hankel matrices to generate
#' @param N Integer. Maximal time index in correlation function to be used in
#'                   Hankel matrix
#' @param deltat Integer. Time shift to be used to build the Hankel matrix
#' @param submatrix.size Integer. Submatrix size to be used in build
#'   of Hankel matrices. Submatrix.size > 1 is experimental.
#' @param element.order Integer vector. specifies how to fit the \code{n} linearly ordered single
#'    correlators into the correlator
#'    matrix for submatrix.size > 1. \code{element.order=c(1,2,3,4)} leads to a matrix
#'    \code{matrix(cf[element.order], nrow=2)}.
#'    Matrix elements can occur multiple times, such as \code{c(1,2,2,3)} for the symmetric case,
#'    for example.
#' @param Delta integer. Delta is the time shift used in the Hankel matrix.
#' 
#' @return
#' A complex vector of length \code{n + n^2} which contains the eigenvalues in the first
#' \code{n} elements and the eigenvectors in the remaining \code{n^2} elements.
#' 
#' A vector of NAs of \code{n + n^2} is returend in case the QR decomposition fails.
#' 
#' @family hankel
gevp.hankel <- function(cf, t0=1, deltat=1, n, N,
                        submatrix.size=1, element.order=c(1,2,3,4),
                        Delta=1) {
  stopifnot(t0 >= 0 && n > 0 && N > 0 && Delta > 0)
  stopifnot(t0 + 1 + 2*(n/submatrix.size-1)*Delta + deltat <= N)
  t0p1 <- t0+1
  cM1 <- array(NA, dim=c(n, n))
  cM2 <- cM1
  ii <- seq(from=1, to=n, by=submatrix.size)

  for(i in c(1:submatrix.size)) {
    for(j in c(1:submatrix.size)) {
      cor.id <- element.order[(i-1)*submatrix.size+j]
      cM1[ii+i-1,ii+j-1] <- hankel.matrix(n=n/submatrix.size, z=cf[seq(from=(cor.id-1)*N+t0p1, to=(cor.id)*N, by=Delta)])
      cM2[ii+i-1,ii+j-1] <- hankel.matrix(n=n/submatrix.size, z=cf[seq(from=(cor.id-1)*N+t0p1+deltat, to=(cor.id)*N, by=Delta)])
    }
  }
  if(submatrix.size > 1) {
    ## symmetrise
    cM1 <- 0.5*(cM1 + t(cM1))
    cM2 <- 0.5*(cM2 + t(cM2))
  }
  ev.cM <- eigen(cM1, symmetric=TRUE, only.values = TRUE)
  positive <- TRUE
  if(any(ev.cM$values <= 0)) positive <- FALSE
  M <- matrix()
  if(positive) {
    ## compute Cholesky factorisation
    invL <- solve(chol(cM1))
    M <- t(invL) %*% cM2 %*% invL
  }
  if(!positive) {
    ## QR decomposition
    qr.cM1 <- qr(cM1)
    M <- try(qr.coef(qr.cM1, cM2), TRUE)
  }
  if(inherits(M, "try-error")) {
    return(rep(NA, times=(n+n^2)))
  }

  M.eigen <- try(eigen(M, symmetric=positive, only.values=FALSE), TRUE)
  if(inherits(M.eigen, "try-error")) {
    return(rep(NA, times=(n+n^2)))
  }
  return(invisible(c(M.eigen$values, as.vector(M.eigen$vectors))))
}

#' @title GEVP method based on Hankel matrices. 
#' 
#' @description
#' Alternative method to determine energy levels from correlation
#'   matrices. A so-called Hankel matrix is generated from an input
#'   \link{cf} object and a generalised eigenvalue problem is solved
#'   then. This is the function to call. It will perform a bootstrap
#'   analysis. 
#'
#' @param cf object of type \link{cf}
#' @param t0 Integer. Initial time value of the GEVP, must be in between 0 and
#'    \code{Time/2-n}. Default is 1. Used when \code{t0fixed=TRUE}.
#' @param deltat Integer. value of deltat used in the hankel GEVP. Default is 1. Used
#'   \code{t0fixed=FALSE}
#' @param submatrix.size Integer. Submatrix size to be used in build
#'   of Hankel matrices. Submatrix.size > 1 is experimental.
#' @param n Integer. Size of the Hankel matrices to generate
#' @param N Integer. Maximal time index in correlation function to be used in
#'                   Hankel matrix
#' @param t0fixed Integer. If set to \code{TRUE}, keep t0 fixed and vary deltat, otherwise
#'    keep deltat fixed and vary t0.
#' @param element.order Integer vector. specifies how to fit the \code{n} linearly ordered single
#'    correlators into the correlator
#'    matrix for submatrix.size > 1. \code{element.order=c(1,2,3,4)} leads to a matrix
#'    \code{matrix(cf[element.order], nrow=2)}.
#'    Matrix elements can occur multiple times, such as \code{c(1,2,2,3)} for the symmetric case,
#'    for example.
#' @param Delta integer. Delta is the time shift used in the Hankel matrix.
#' 
#' @details
#' See `vignette(name="hankel", package="hadron")`
#'
#' @return
#' List object of class "hankel". The eigenvalues are stored in a
#' numeric vector \code{t0}, the corresonding samples in \code{t}. The reference input
#' time \code{t0} is stored as \code{reference_time} in the returned list.
#'
#' @family hankel
#' @export
#' @examples
#'
#' data(correlatormatrix)
#' correlatormatrix <- bootstrap.cf(correlatormatrix, boot.R=99, boot.l=1, seed=132435)
#' t0 <- 4
#' correlatormatrix.gevp <- bootstrap.gevp(cf=correlatormatrix, t0=t0, element.order=c(1,2,3,4))
#' pc1 <- gevp2cf(gevp=correlatormatrix.gevp, id=1)
#' pc1.hankel <- bootstrap.hankel(cf=pc1, t0=1, n=2)
#' hpc1 <- hankel2cf(hankel=pc1.hankel, id=1)
#' plot(hpc1, log="y")
#' heffectivemass1 <- hankel2effectivemass(hankel=pc1.hankel, id=1)
bootstrap.hankel <- function(cf, t0=1, n=2, N = (cf$Time/2+1),
                             t0fixed=TRUE, deltat=1, Delta=1,
                             submatrix.size=1, element.order=1) {
  stopifnot(inherits(cf, 'cf_meta'))
  stopifnot(inherits(cf, 'cf_boot'))
  stopifnot((n %% submatrix.size) == 0)
  if(!t0fixed) {
    stopifnot(deltat > 0)
  }
  else {
    stopifnot(t0 > 0)
  }
  ## R/Fortran index convention
  t0p1 <- t0 + 1
  boot.R <- cf$boot.R
  evs <- array(NA, dim=c(N, n + n^2))
  evs.tsboot <- array(NA, dim=c(boot.R, N, n + n^2))

  if(t0fixed) {
    for(deltat in c(1:(N-1-t0-2*(n/submatrix.size-1)*Delta))) {
      evs[deltat+t0p1, ] <- gevp.hankel(cf$cf0, t0=t0,
                                        n=n, N=N, deltat=deltat,
                                        submatrix.size=submatrix.size, element.order=element.order,
                                        Delta=Delta)
      
      evs.tsboot[, deltat+t0p1, ] <- t(apply(cf$cf.tsboot$t, 1, gevp.hankel, t0=t0,
                                             n=n, N=N, deltat=deltat,
                                             submatrix.size=submatrix.size, element.order=element.order,
                                             Delta=Delta))
    }
  }
  else {
    for(t02 in c(1:(N-1-deltat-2*(n/submatrix.size-1)*Delta))) {
      evs[t02, ] <- gevp.hankel(cf$cf0, t0=t02,
                                n=n, N=N, deltat=deltat,
                                submatrix.size=submatrix.size, element.order=element.order,
                                Delta=Delta)
      
      evs.tsboot[, t02, ] <- t(apply(cf$cf.tsboot$t, 1, gevp.hankel, t0=t02,
                                     n=n, N=N, deltat=deltat,
                                     submatrix.size=submatrix.size, element.order=element.order,
                                     Delta=Delta))
    }
  }
  ret <- list(cf=cf,
              t0=evs[ ,c(1:n), drop=FALSE],
              t=evs.tsboot[ ,, c(1:n), drop=FALSE],
              vectors=evs[ ,c((n+1):(n+n^2)), drop=FALSE],
              vectors.tsboot=evs.tsboot[ ,, c((n+1):(n+n^2)), drop=FALSE],
              boot.R=boot.R,
              boot.l=cf$boot.l,
              seed=cf$seed,
              reference_time=t0,
              t0fixed=t0fixed,
              submatrix.size=submatrix.size,
              element.order=element.order,
              Delta=Delta,
              n=n,
              N=N)
  if(!t0fixed) {
    ret$reference_time=deltat
  }
  class(ret) <- c("hankel", class(ret))
  return(invisible(ret))
}

#' @title plot_hankel_spectrum
#'
#' @description
#' produces a scatter plot of the complex \eqn{-\log}{-log}
#' of the eigenvalues produced by the
#' \link{bootstrap.hankel} method. In addition, produces a
#' histogramm of all real and positive eigenvalues after computing
#' \eqn{-\log(ev)/\delta t}{-log(ev)/delta t} in the range
#' (0,1) and determines its mode.
#' 
#' @param hankel object as returned from \link{bootstrap.hankel}
#' @param deltat Integer. Time shift at which to plot 
#' @param id Integer vector. Indices of eigenvalues to be plotted.
#' Must be part of \code{c(1:hankel$n)}.
#' 
#' @family hankel
#'
#' @return
#' No return value.
plot_hankel_spectrum <- function(hankel, deltat=1, id=c(1:hankel$n)) {
  n <- hankel$n
  stopifnot(max(id) <= n & min(id) >= 1)
  col <- c("black", rainbow(n=n-1))
  xlim <- range(c(Re(-log(hankel$t)/deltat), Re(-log(hankel$t0)/deltat)), na.rm=TRUE)
  ylim <- range(c((Im(-log(hankel$t)) %% pi)/deltat,
                  (Im(-log(hankel$t0)) %% pi)/deltat),
                na.rm=TRUE)
  if(max(abs(ylim)) < 0.001) ylim=c(-0.01, 0.01)
  plot(NA, xlim=xlim, ylim=ylim,
       xlab="Re", ylab="Im")

  tt <- hankel$reference_time+deltat
  for(i in id) {
    points(x=Re(-log(hankel$t[, tt, i])/deltat),
           y=(Im(-log(hankel$t[, tt, i])) %% pi)/deltat,
           col=col[i], bg="white")
    points(x=Re(-log(hankel$t0[tt, i])/deltat),
           y=(Im(-log(hankel$t0[tt, i])) %% pi)/deltat,
           col=col[i], bg=col[i])
  }
  tmp <- hankel$t
  tmp[Im(tmp) != 0] <- NA
  tmp[Re(tmp) < 0] <- NA
  tmp[Re(tmp) > 1] <- NA
  tmp <- Re(-log(tmp[, tt, ])/deltat)
  tmp[tmp > 1]  <- NA
  new_window_if_appropriate()
  hist(tmp, xlim=c(0, 1), main="Histogram of Samples",
       xlab="E", breaks=seq(0, 1, 0.02))
  mode<-density(tmp, na.rm=TRUE)$x[which.max(density(tmp, na.rm=TRUE)$y)]
  new_window_if_appropriate()
  plot(density(tmp, na.rm=TRUE))
  message("Mode of this density:", mode, "deltat:", deltat, "\n")
}

#' @title hankel2cf
#'
#' @param hankel object as returned from \link{bootstrap.hankel}
#' @param id Integer. ID of the principal correlator to extract
#' @param eps Numeric. Cut-off: if the imaginary part of the generalised
#' eigenvalues is larger than eps, the eigenvalue is discarded.
#' @param range Numeric vector. Value-range for the real part of the eigenvalues
#' (not the energies). If outside this range, the eigenvalue will be discarded
#' @param id Integer. Index of eigenvalue to consider,
#' \eqn{1\leq id\leq n}{1 <= id <= n}.
#' @param sort.type the sort algorithm to be used to sort the eigenvalues. This can be
#' either simply "values", or the eigenvector information is used in addition with
#' "vectors"
#' @param sort.t0 Boolean. Whether to use the eigenvector at t0 or the one at deltat-1
#' for sorting
#' 
#' @family hankel
#' @seealso input is generated via \link{bootstrap.hankel}
#' alternatively use \link{hankel2effectivemass}. For the `cf` class see
#' \link{cf}
#'
#' @return
#' Returns an object of S3 class `cf`.
#' @export
hankel2cf <- function(hankel, id=c(1), range=c(0,1), eps=1.e-16,
                      sort.type="values", sort.t0=TRUE) {
  stopifnot(inherits(hankel, "hankel"))
  stopifnot((id <= hankel$n && id >= 1))
  stopifnot(length(id) == 1)
  stopifnot(sort.type %in% c("values", "vectors", "mindist"))
  N <- hankel$N
  n <- hankel$n
  submatrix.size <- hankel$submatrix.size
  Delta <- hankel$Delta
  
  range <- range(range)
  
  reftime <- hankel$reference_time
  ## Base `cf` properties.
  mycf <- cf_meta(nrObs = 1,
                  Time = hankel$cf$Time,
                  nrStypes = 1,
                  symmetrised = hankel$cf$symmetrised)

  ## Add the `cf_boot` mixin.
  cf.tsboot <- list(
      t0 = array(NA, dim=c(N)),
      t = array(NA, dim=c(hankel$boot.R, N))
  )

  cf.tsboot$t0[reftime+1] <- 1
  cf.tsboot$t[, reftime+1] <- 1
  if(sort.type == "values" || sort.type=="mindist") {
    idpass <- id
    if(hankel$t0fixed) {
      idpass <- NA
    }
    .fn <- function(evs, range, eps, id) {
      if(is.na(id)) {
        ret <- Re(evs[which.min(abs(Re(evs)-mean(range)))])
        if(length(ret) == 0) return(NA)
        return(ret)
      }
      ii <- which(abs(Im(evs)) <= eps & Re(evs) > range[1]
                  & Re(evs) < range[2])
      return(Re(evs[ii[id]]))
    }
    .fn2 <- function(evs, refvalue, eps) {
      evs2 <- evs[which(abs(Im(evs)) <= eps)]
      ret <- Re(evs2[which.min(abs(Re(evs2)-refvalue))])
      if(length(ret) == 0) return(NA)
      return(ret)
    }
    for(deltat in c(1:(N-1-reftime-2*(n/submatrix.size - 1)*Delta))) {
      cf.tsboot$t0[deltat+reftime] <- .fn(evs=hankel$t0[deltat+reftime, , drop = FALSE],
                                          range=range, eps=eps, id=idpass)
      if(sort.type=="values") {
        cf.tsboot$t[, deltat+reftime] <- apply(X=hankel$t[, deltat+reftime, , drop = FALSE],
                                               MARGIN=1, FUN=.fn,
                                               range=range, eps=eps, id=idpass)
      }
      else { ## mindist
        cf.tsboot$t[, deltat+reftime] <- apply(X=hankel$t[, deltat+reftime, , drop = FALSE],
                                               MARGIN=1, FUN=.fn2, eps=eps,
                                               refvalue=cf.tsboot$t0[deltat+reftime])
      }
    }
  }
  if(sort.type == "vectors") {
    ## obtain indices of "real" eigenvalues in the appropriate range
    ## chosing the one with maximal overlap in the eigenvectors
    .fnii <- function(evs, range, eps) {
      return(which(abs(Im(evs)) <= eps & Re(evs) > range[1]
                & Re(evs) < range[2]))
    }
    .fn3 <- function(evs, ii, id, n,
                    vectors, v0) {
      if(length(ii) == 0) {
        return(NA)
      }
      if(is.na(v0[1])) {
        return(Re(evs[ii[id]]))
      }
      sp <- c()
      ## chose the one with maximal overlap in the eigenvectors
      for(i in c(1:length(ii))) {
        sp[i] <- abs(sum(Conj(v0) * vectors[c(((ii[i]-1)*n+1) : (ii[i]*n))]))
      }
      return(Re(evs[ii[order(sp, decreasing=TRUE, na.last=TRUE)[id]]]))
    }

    ## precompute reference eigenvector at deltat == 1
    ## and use method "values" at deltat == 1
    ii <- .fnii(evs=hankel$t0[reftime+1, , drop = FALSE],
                range=range, eps=eps)
    v0 <- NA
    ## orthogonality is for <v_i, H(t0) v0_j>
    if(length(ii) != 0) v0 <- hankel.matrix(n=n, z=hankel$cf$cf0[reftime+1]) %*%
                          hankel$vectors[reftime+1, c(((ii[id]-1)*n+1) : (ii[id]*n))]
    cf.tsboot$t0[reftime+1] <- .fn3(evs=hankel$t0[1+reftime, , drop = FALSE],
                                   ii=ii, id=id, n=n,
                                   vectors=NA, v0=NA)
    v0.tsboot <- array(NA, dim=c(hankel$boot.R, n))
    for(rr in c(1:hankel$boot.R)) {
      ii <- .fnii(evs=hankel$t[rr, reftime+1, , drop = FALSE],
                  range=range, eps=eps)
      v0.tsboot[rr,] <- NA
      if(length(ii) != 0) {
        v0.tsboot[rr,] <- hankel.matrix(n=n, z=hankel$cf$cf.tsboot$t[rr, reftime+1]) %*%
          hankel$vectors.tsboot[rr, reftime+1, c(((ii[id]-1)*n+1) : (ii[id]*n))]
      }
      cf.tsboot$t[rr, reftime+1] <- .fn3(evs=hankel$t[rr, 1+reftime, , drop = FALSE],
                                        ii=ii, id=id, n=n,
                                        vectors=NA, v0=NA)
    }
    ## now use pre-computed vectors
    for(deltat in c(2:(N-1-reftime-2*(n/submatrix.size-1)*Delta))) {
      ii <- .fnii(evs=hankel$t0[reftime+deltat, , drop = FALSE],
                  range=range, eps=eps)
      cf.tsboot$t0[reftime+deltat] <- .fn3(evs=hankel$t0[reftime+deltat, , drop = FALSE],
                                          ii=ii, id=id, n=n,
                                          vectors=hankel$vectors[reftime+deltat,],
                                          v0=v0)
      if(!sort.t0 && !is.na(cf.tsboot$t0[reftime+deltat])) {
        k <- which(Re(hankel$t0[reftime+deltat, ]) == cf.tsboot$t0[reftime+deltat])
        v0 <- hankel.matrix(n=n, z=hankel$cf$cf0[reftime+deltat]) %*%
          hankel$vectors[reftime+deltat, c(((k-1)*n+1) : (k*n))]
      }
      for(rr in c(1:hankel$boot.R)) {
        ii <- .fnii(evs=hankel$t[rr, reftime+deltat, , drop = FALSE],
                    range=range, eps=eps)
        cf.tsboot$t[rr, reftime+deltat] <- .fn3(evs=hankel$t[rr, reftime+deltat, , drop = FALSE],
                                               ii=ii, id=id, n=n,
                                               vectors=hankel$vectors.tsboot[rr, reftime+deltat,],
                                               v0=v0.tsboot[rr,])
        if(!sort.t0 && !is.na(cf.tsboot$t[rr, reftime+deltat])) {
          k <- which(Re(hankel$t[rr, reftime+deltat,]) == cf.tsboot$t[rr, reftime+deltat])
          v0.tsboot[rr,] <- hankel.matrix(n=n, z=hankel$cf$cf.tsboot$t[rr, reftime+deltat]) %*%
            hankel$vectors.tsboot[rr, reftime+deltat, c(((k-1)*n+1) : (k*n))]
        }
      }
    }
  }

  mycf <- cf_boot(.cf = mycf,
                  boot.R = hankel$boot.R,
                  boot.l = hankel$boot.l,
                  seed = hankel$seed,
                  sim = hankel$cf$sim,
                  endcorr = hankel$cf$endcorr,
                  cf.tsboot = cf.tsboot,
                  resampling_method = hankel$cf$resampling_method)
  
  mycf <- cf_principal_correlator(.cf = mycf,
                                  id = id,
                                  gevp_reference_time = hankel$reference_time)
  mycf$hankel_eps <- eps
  mycf$hankel_range <- range
  return(invisible(mycf))
}

#' @title hankel2effectivemass
#'
#' @inheritParams hankel2cf
#' @param type Character vector. Type of effective mass to use.
#' Must be in \code{c("log", "acosh")}
#' 
#' @family hankel
#' @seealso input is generated via \link{bootstrap.hankel}
#' alternatively use \link{hankel2cf}. See also
#' \link{bootstrap.effectivemass}
#'
#' @return
#' Returns an object of S3 class `effectivemass`.
#' 
#' @export
hankel2effectivemass  <- function(hankel, id=c(1), type="log",
                                  range=c(0,1), eps=1.e-16,
                                  sort.type="values", sort.t0=TRUE) {
  stopifnot(inherits(hankel, "hankel"))
  stopifnot(length(id) == 1)
  stopifnot((id <= hankel$n && id >= 1))
  range <- range(range)
  
  tmax <- hankel$cf$Time/2
  if(!hankel$cf$symmetrised){
    tmax <- hankel$cf$Time-1
  }

  pc <- hankel2cf(hankel=hankel, id=id, range=range, eps=eps,
                  sort.type=sort.type, sort.t0=sort.t0)

  effMass <- c()
  effMass.tsboot <- c()
  if(hankel$t0fixed) {
    deltat <- c(1:(tmax+1))-hankel$reference_time
    if(type == "log") {
      effMass <- -log(pc$cf.tsboot$t0)/deltat
      effMass.tsboot <- -log(pc$cf.tsboot$t)/t(array(deltat, dim=rev(dim(pc$cf.tsboot$t))))
    }
    if(type == "acosh") {
      effMass <- -acosh(pc$cf.tsboot$t0)/deltat
      effMass.tsboot <- -acosh(pc$cf.tsboot$t)/t(array(deltat, dim=rev(dim(pc$cf.tsboot$t))))
    }
  }
  else {
    deltat <- hankel$reference_time
    effMass <- -log(pc$cf.tsboot$t0)/deltat
    effMass.tsboot <- -log(pc$cf.tsboot$t)/deltat
  }
  deffMass <- apply(effMass.tsboot, 2, hankel$cf$error_fn, na.rm=TRUE)

  ret <- list(t.idx=c(1:(tmax+1)),
              effMass=effMass, deffMass=deffMass, effMass.tsboot=effMass.tsboot,
              opt.res=NULL, t1=NULL, t2=NULL, type=type, useCov=NULL, CovMatrix=NULL, invCovMatrix=NULL,
              boot.R = hankel$boot.R, boot.l = hankel$boot.l, seed = hankel$seed,
              massfit.tsboot=NULL, Time=hankel$cf$Time, N=1, nrObs=1, dof=NULL,
              chisqr=NULL, Qval=NULL
             )
  ret$cf <- hankel$cf
  ret$t0 <- effMass
  ret$t <- effMass.tsboot
  ret$se <- apply(ret$t, MARGIN=2L, FUN=hankel$cf$error_fn, na.rm=TRUE)
  attr(ret, "class") <- c("effectivemass", class(ret))
  return(ret)
}

#' @title hankeldensity2effectivemass
#'
#' @param hankel object as returned from \link{bootstrap.hankel}
#' @param range Numeric vector. Value-range for the real part of the eigenvalues.
#' If outside this range, the eigenvalue will be discarded
#' @param method Character vector. Method to be used to determine the central value
#' of the effective mass. Must be "median" (default) or "density"
#' @description
#'
#' computes the density of all bootstrap replicates of effective masses
#'
#' @seealso \link{bootstrap.effectivemass}, \link{hankel2effectivemass}
#' @return
#' Returns an object of S3 class `effectivemass`.
#' #' 
hankeldensity2effectivemass <- function(hankel, range=c(0,1),
                                        method="median") {

  stopifnot(inherits(hankel, "hankel"))
  stopifnot(method %in% c("median", "density"))
  
  N <- hankel$N
  n <- hankel$n
  reftime <- hankel$reference_time
  
  
  tmax <- hankel$cf$Time/2
  if(!hankel$cf$symmetrised){
    tmax <- hankel$cf$Time-1
  }

  effMass <- rep(NA, times=N)
  effMass.tsboot <- array(NA, dim=c(hankel$boot.R, N))
  se <- rep(NA, times=N)
  for(deltat in c(1:(N-2-reftime-2*n))) {
    tmp <- as.vector(hankel$t[, deltat+reftime,])
  
    ii <- which(Im(tmp) == 0 &
                Re(tmp) > 0 &
                Re(tmp) < 1)
    tmp <- -log(Re(tmp[ii]))/deltat
    ii <- which(tmp <= range[2] & tmp >= range[1])
    tmp <- tmp[ii]
    
    if(method=="median") {
      effMass[reftime+deltat] <- median(tmp)
    }
    else {
      hdens <- density(tmp, na.rm=TRUE, from=range[1], to=range[2])  
      effMass[reftime+deltat] <- hdens$x[which.max(hdens$y)]
    }
    se[reftime+deltat] <- hankel$cf$error_fn(tmp)
  }
  ret <- list(t.idx=c(1:(tmax+1)),
              effMass=effMass, deffMass=se, effMass.tsboot=effMass.tsboot,
              opt.res=NULL, t1=NULL, t2=NULL, type="log", useCov=NULL, CovMatrix=NULL, invCovMatrix=NULL,
              boot.R = hankel$boot.R, boot.l = hankel$boot.l, seed = hankel$seed,
              massfit.tsboot=NULL, Time=hankel$cf$Time, N=1, nrObs=1, dof=NULL,
              chisqr=NULL, Qval=NULL
             )
  ret$cf <- hankel$cf
  ret$t0 <- effMass
  ret$t <- effMass.tsboot
  ret$se <- se
  attr(ret, "class") <- c("effectivemass", class(ret))
  return(ret)
}

#' Resample bootstrap samples in Hankel effmass
#'
#' The bootstrap distribution in the Hankel effective mass can be quite broad
#' due to outliers and long tails. These screw with proper error estimation.
#' Therefore it can be useful to trim these tails. Just trimming a bootstrap
#' distribution would lead to less samples, therefore we do a parametric
#' resampling.
#'
#' The central values are also inferred from the distribution because they
#' often are outliers themselves. The new central value is the middle between
#' the upper and lower quantile, making the resulting distribution symmetric.
#'
#' Half the distance between the quantiles is taken to be the error, therefore
#' the quantiles are chosen at 16 and 84 percent to match the standard
#' deviation. All points that are more than “distance” errors away from the new
#' central value are taken to be outliers.
#'
#' @param hankel_effmass Hankel effective mass from `hankel2effmass`.
#' @param distance Numeric, threshold for marking outliers.
#'
#' @return
#' The Hankel effmass object is returned with the same fields, the numbers
#' have been changed.
#'
#' Additionally there are the followi1ng fields:
#'
#' - `cov_full` contains the full covariance matrix as determined from all the
#'   data. This will be skewed by the outliers.
#' - `finite_count` gives the number of non-outliers per time slice.
#' - `complete_count` gives the numbers of complete cases if all outliers are
#'   taken out. This number is often zero because the late time slices contain
#'   lots of outliers due to the noise.
#'- `cov_3sigma_pairwise` is the covariance matrix using only the non-outliers
#'  and removing NAs in a pairwise fashion, using the maximum of the data. This
#'  is the covariance matrix that is used for the resampling.
#'  
#' In case that no time slices had a finite error estimate, this function
#' returns just `NA`.
#'
#' @export
resample_hankel <- function (hankel_effmass, distance = 5.0) {
    boot_R <- nrow(hankel_effmass$effMass.tsboot)
    lower <- apply(hankel_effmass$effMass.tsboot, 2, quantile, 0.16, na.rm = TRUE)
    upper <- apply(hankel_effmass$effMass.tsboot, 2, quantile, 0.84, na.rm = TRUE)
    val <- (upper + lower) / 2
    err <- (upper - lower) / 2
    
    use <- !is.na(err)
    
    if (sum(use[2:length(use)]) == 0) {
        return (NA)
    }
    
    hankel_effmass$cov_full <- cov(hankel_effmass$effMass.tsboot)
    all_outlier <- matrix(NA, nrow = nrow(hankel_effmass$effMass.tsboot), ncol = ncol(hankel_effmass$effMass.tsboot))
    
    for (t in 1:ncol(hankel_effmass$effMass.tsboot)) {
        outlier <- hankel_effmass$effMass.tsboot[, t] > val[t] + distance * err[t] |
            hankel_effmass$effMass.tsboot[, t] < val[t] - distance * err[t]
        all_outlier[, t] <- outlier
        hankel_effmass$effMass.tsboot[outlier, t] <- NA
    }
    
    hankel_effmass$finite_count <- apply(!all_outlier, 2, sum, na.rm = TRUE)
    hankel_effmass$complete_count <- sum(complete.cases(all_outlier))
    hankel_effmass$cov_3sigma_pairwise <- cov(hankel_effmass$effMass.tsboot, use = 'pairwise.complete.obs')
    
    hankel_effmass$effMass <- val
    hankel_effmass$deffMass <- err
    hankel_effmass$effMass.tsboot[, 1:ncol(hankel_effmass$effMass.tsboot)] <- NA
    hankel_effmass$effMass.tsboot[, use] <- parametric.bootstrap.cov(boot_R, val[use], hankel_effmass$cov_3sigma_pairwise[use, use, drop = FALSE], 1)
    
    # Don't even ask …
    hankel_effmass$t0 <- hankel_effmass$effMass
    hankel_effmass$t <- hankel_effmass$effMass.tsboot
    hankel_effmass$se <- hankel_effmass$deffMass

    return (hankel_effmass)
}
