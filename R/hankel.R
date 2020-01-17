## creates a square Hankel matrix of size nxn
hankel.matrix <- function(n, z){
  outer(1:n,1:n,
        function(x,y) z[x+y-1])
}

gevp.hankel.old <- function(cf, t0, deltat = 1, n, N, id=c(1),
                            range=c(0.,1.), eps=0.001, debug=FALSE, cor.id=1) {
  if(t0+2*n+deltat+1 > N) {
    stop("t0+n+deltat > N\n")
  }
  ## build the two Hankel matrices
  cM1 <- hankel.matrix(n=n, z=cf[((cor.id-1)*N+t0):((cor.id)*N)])
  cM2 <- hankel.matrix(n=n, z=cf[(((cor.id-1)*N+t0)+deltat):((cor.id)*N)])

  qr.cM1 <- qr(cM1)
  M <- qr.coef(qr.cM1, cM2)
  ## M = cM1^-1 * cM2

  M.eigen <- eigen(M, symmetric=FALSE, only.values=TRUE)

  ## reduce to real eigenvalues with real part larger than 1
  ii <- which(abs(Im(M.eigen$values)) < eps*abs(Re(M.eigen$values)) & Re(M.eigen$values) > range[1] & Re(M.eigen$values) < range[2] )
  if(debug) print(M.eigen$values[which(abs(Im(M.eigen$values)) < eps*abs(Re(M.eigen$values)))])
  return(sort(Re(M.eigen$values[ii]), decreasing=TRUE)[id])
  
}

#' @title GEVP method based on Hankel matrices.
#' 
#' @description
#' Alternative method to determine energy levels from correlation
#'   matrices. A so-called Hankel matrix is generated from an input
#'   real numeric vector and a generalised eigenvalue problem is solved
#'   then.
#'
#' @param cf Numeric vector.
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
#'    Double indexing is allowed.
#'
#'
#' @return
#' A complex vector of length \code{n + n^2} with firsth the eigenvalues and
#' then the eigenvectors is returned
#' 
#' A vector of NAs of \code{n + n^2} is returend in case the QR decomposition fails.
#' 
#' @family hankel
gevp.hankel <- function(cf, t0=1, deltat=1, n, N,
                        submatrix.size=1, element.order=c(1,2,3,4)) {
  if(t0+2*n+deltat > N) {
    stop("t0+n+deltat > N\n")
  }
  t0 <- t0+1
  cM1 <- array(NA, dim=c(submatrix.size*n, submatrix.size*n))
  cM2 <- cM1
  ii <- seq(from=1, to=submatrix.size*n, by=submatrix.size)

  for(i in c(1:submatrix.size)) {
    for(j in c(1:submatrix.size)) {
      cor.id <- element.order[(i-1)*submatrix.size+j]
      cM1[ii+i-1,ii+j-1] <- hankel.matrix(n=n, z=cf[((cor.id-1)*N+t0):((cor.id)*N)])
      cM2[ii+i-1,ii+j-1] <- hankel.matrix(n=n, z=cf[(((cor.id-1)*N+t0)+deltat):((cor.id)*N)])
    }
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
#' @param t0     initial time value of the GEVP, must be in between 0 and
#'    \code{Time/2-2}. Default is 1.
#' @param n Integer. Size of the Hankel matrices to generate
#' @param N Integer. Maximal time index in correlation function to be used in
#'                   Hankel matrix
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
#' correlatormatrix <- bootstrap.cf(correlatormatrix, boot.R=400, boot.l=1, seed=132435)
#' t0 <- 4
#' correlatormatrix.gevp <- bootstrap.gevp(cf=correlatormatrix, t0=t0, element.order=c(1,2,3,4))
#' pc1 <- gevp2cf(gevp=correlatormatrix.gevp, id=1)
#' pc1.hankel <- bootstrap.hankel(cf=pc1, t0=1, n=2)
#' hpc1 <- hankel2cf(hankel=pc1.hankel, id=1)
#' plot(hpc1, log="y")
#' heffectivemass1 <- hankel2effectivemass(hankel=pc1.hankel, id=1)
bootstrap.hankel <- function(cf, t0, n=2, N = cf$Time/2+1, id=c(1)) {
  stopifnot(inherits(cf, 'cf_meta'))
  stopifnot(inherits(cf, 'cf_boot'))

  boot.R <- cf$boot.R
  evs <- array(NA, dim=c(N, n + n^2))
  evs.tsboot <- array(NA, dim=c(boot.R, N, n + n^2))

  
  for(deltat in c(1:(N-2-t0-2*n))) {
    evs[deltat+t0, ] <- gevp.hankel(cf$cf0, t0=t0,
                                    n=n, N=N, deltat=deltat,
                                    submatrix.size=1, element.order=c(1))

    evs.tsboot[, deltat+t0, ] <- t(apply(cf$cf.tsboot$t, 1, gevp.hankel, t0=t0,
                                         n=n, N=N, deltat=deltat,
                                         submatrix.size=1, element.order=1))
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
              n=n,
              N=N)
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
#' @export
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
  hadron:::new_window_if_appropriate()
  hist(tmp, xlim=c(0, 1), main="Histogram of Samples",
       xlab="E", breaks=seq(0, 1, 0.02))
  mode<-density(tmp, na.rm=TRUE)$x[which.max(density(tmp, na.rm=TRUE)$y)]
  hadron:::new_window_if_appropriate()
  plot(density(tmp, na.rm=TRUE))
  cat("Mode of this density:", mode, "deltat:", deltat, "\n")
}

#' @title hankel2cf
#'
#' @param hankel object as returned from \link{bootstrap.hankel}
#' @param id Integer. ID of the principal correlator to extract
#' @param eps Numeric. Cut-off: if the imaginary part of the generalised
#' eigenvalues is larger than eps, the eigenvalue is discarded.
#' @param range Numeric vector. Value-range for the real part of the eigenvalues.
#' If outside this range, the eigenvalue will be discarded
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
#' alternatively use \link{hankel2effectivemass}.
#'
#' @export
hankel2cf <- function(hankel, id=c(1), range=c(0,1), eps=1.e-16,
                      sort.type="values", sort.t0=TRUE) {
  stopifnot(inherits(hankel, "hankel"))
  stopifnot((id <= hankel$n && id >= 1))
  stopifnot(length(id) == 1)
  stopifnot(sort.type %in% c("values", "vectors"))
  N <- hankel$N
  n <- hankel$n
  
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

  cf.tsboot$t0[reftime] <- 1
  cf.tsboot$t[, reftime] <- 1
  if(sort.type == "values") {
    .fn <- function(evs, range, eps, id) {
      ii <- which(abs(Im(evs)) <= eps & Re(evs) > range[1]
                  & Re(evs) < range[2])
      return(Re(evs[ii[id]]))
    }
    
    for(deltat in c(1:(N-2-reftime-2*n))) {
      cf.tsboot$t0[deltat+reftime] <- .fn(evs=hankel$t0[deltat+reftime, , drop = FALSE],
                                          range=range, eps=eps, id=id)
      cf.tsboot$t[, deltat+reftime] <- apply(X=hankel$t[, deltat+reftime, , drop = FALSE],
                                             MARGIN=1, FUN=.fn,
                                             range=range, eps=eps, id=id)
    }
  }
  if(sort.type == "vectors") {
    ## obtain indices of "real" eigenvalues in the appropriate range
    ## chosing the one with maximal overlap in the eigenvectors
    .fnii <- function(evs, range, eps) {
      return(which(abs(Im(evs)) <= eps & Re(evs) > range[1]
                & Re(evs) < range[2]))
    }
    .fn <- function(evs, ii, id, n,
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
    ## and use method "values" at detlat == 1
    ii <- .fnii(evs=hankel$t0[reftime+1, , drop = FALSE],
                range=range, eps=eps)
    v0 <- NA
    if(length(ii) != 0) v0 <- hankel$vectors[reftime+1, c(((ii[id]-1)*n+1) : (ii[id]*n))]
    cf.tsboot$t0[reftime+1] <- .fn(evs=hankel$t0[1+reftime, , drop = FALSE],
                                   ii=ii, id=id, n=n,
                                   vectors=NA, v0=NA)
    v0.tsboot <- array(NA, dim=c(hankel$boot.R, n))
    for(rr in c(1:hankel$boot.R)) {
      ii <- .fnii(evs=hankel$t[rr, reftime+1, , drop = FALSE],
                  range=range, eps=eps)
      v0.tsboot[rr,] <- NA
      if(length(ii) != 0) {
        v0.tsboot[rr,] <- hankel$vectors.tsboot[rr, reftime+1, c(((ii[id]-1)*n+1) : (ii[id]*n))]
      }
      cf.tsboot$t[rr, reftime+1] <- .fn(evs=hankel$t[rr, 1+reftime, , drop = FALSE],
                                        ii=ii, id=id, n=n,
                                        vectors=NA, v0=NA)
    }
    ## now use pre-computed vectors
    for(deltat in c(2:(N-2-reftime-2*n))) {
      ii <- .fnii(evs=hankel$t0[reftime+deltat, , drop = FALSE],
                  range=range, eps=eps)
      cf.tsboot$t0[reftime+deltat] <- .fn(evs=hankel$t0[reftime+deltat, , drop = FALSE],
                                          ii=ii, id=id, n=n,
                                          vectors=hankel$vectors[reftime+deltat,],
                                          v0=v0)
      if(!sort.t0 && !is.na(cf.tsboot$t0[reftime+deltat])) {
        k <- which(Re(hankel$t0[reftime+deltat, ]) == cf.tsboot$t0[reftime+deltat])
        v0 <- hankel$vectors[reftime+deltat, c(((k-1)*n+1) : (k*n))]
      }
      for(rr in c(1:hankel$boot.R)) {
        ii <- .fnii(evs=hankel$t[rr, reftime+deltat, , drop = FALSE],
                    range=range, eps=eps)
        cf.tsboot$t[rr, reftime+deltat] <- .fn(evs=hankel$t[rr, reftime+deltat, , drop = FALSE],
                                               ii=ii, id=id, n=n,
                                               vectors=hankel$vectors.tsboot[rr, reftime+deltat,],
                                               v0=v0.tsboot[rr,])
        if(!sort.t0 && !is.na(cf.tsboot$t[rr, reftime+deltat])) {
          k <- which(Re(hankel$t[rr, reftime+deltat,]) == cf.tsboot$t[rr, reftime+deltat])
          v0.tsboot[rr,] <- hankel$vectors.tsboot[rr, reftime+deltat, c(((k-1)*n+1) : (k*n))]
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
#' alternatively use \link{hankel2cf}.
#'
#' @export
hankel2effectivemass  <- function(hankel, id=c(1), type="log",
                                  range=c(0,1), eps=1.e-16,
                                  sort.type="values", sort.t0=TRUE) {
  stopifnot(inherits(hankel, "hankel"))
  stopifnot(length(id) == 1)
  stopifnot((id <= hankel$n && id >= 1))
  stopifnot(type %in% c("log", "acosh"))
  stopifnot(sort.type %in% c("values", "vectors"))
  
  tmax <- hankel$cf$Time/2
  if(!hankel$cf$symmetrised){
    tmax <- hankel$cf$Time-1
  }

  pc <- hankel2cf(hankel=hankel, id=id, range=range, eps=eps,
                  sort.type=sort.type, sort.t0=sort.t0)

  deltat <- c(1:(tmax+1))-hankel$reference_time
  if(type == "log") {
    effMass <- -log(pc$cf.tsboot$t0)/deltat
    effMass.tsboot <- -log(pc$cf.tsboot$t)/t(array(deltat, dim=rev(dim(pc$cf.tsboot$t))))
  }
  if(type == "acosh") {
    effMass <- -acosh(pc$cf.tsboot$t0)/deltat
    effMass.tsboot <- -acosh(pc$cf.tsboot$t)/t(array(deltat, dim=rev(dim(pc$cf.tsboot$t))))
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

#' hankeldensity2effectivemass
#'
#' @param hankel object as returned from \link{bootstrap.hankel}
#' @param range Numeric vector. Value-range for the real part of the eigenvalues.
#' If outside this range, the eigenvalue will be discarded
#' @param method Character vector. Method to be used to determine the central value
#' of the effective mass. Must be "median" (default) or "density"
#' @description
#' 
#' @export
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
