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
#'   \link{cf} object and a generalised eigenvalue problem is solved
#'   then.
#'
#' @inheritParams bootstrap.hankel
#' @param debug Boolean. Enable debug output.
#' @param deltat Integer. Time shift to be used to build the Hankel matrix
#' @param submatrix.size Integer. Submatrix size to be used in build
#'   of Hankel matrices. Submatrix.size > 1 is experimental.
#' @param element.order Integer vector. specifies how to fit the \code{n} linearly ordered single
#'    correlators into the correlator
#'    matrix for submatrix.size > 1. \code{element.order=c(1,2,3,4)} leads to a matrix
#'    \code{matrix(cf[element.order], nrow=2)}.
#'    Double indexing is allowed.
#' 
#' @family hankel
gevp.hankel <- function(cf, t0=1, deltat=1, n, N, eps=0.0001, range=c(0,1),
                        submatrix.size=1, element.order=c(1,2,3,4), id=c(1),
                        debug=FALSE) {
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

  ## M = cM1^-1 * cM2
  M <- try(qr.coef(qr.cM1, cM2), TRUE)
  if(inherits(M, "try-error")) {
    return(c(NA))
  }
  M[abs(M) < 1.e-12] <- 0
  M.eigen <- try(eigen(M, symmetric=FALSE, only.values=TRUE), TRUE)
  if(inherits(M.eigen, "try-error")) {
    return(NA)
  }
  ii <- which(abs(Im(M.eigen$values)) < eps & Re(M.eigen$values) > range[1]
              & Re(M.eigen$values) < range[2])
  if(id == "all") {
    return(invisible(Re(M.eigen$values[ii])))
  }
  return(invisible(Re(M.eigen$values[ii[id]])))
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
#' @param n Integer. Size of the submatrix Hankel matrices to generate
#' @param N Integer. Maximal time index in correlation function to be used in
#'                   Hankel matrix
#' @param eps Numeric. Cut-off: if the imaginary part of the generalised
#' eigenvalues is larger than eps, the eigenvalue is discarded.
#' @param range Numeric vector. Value-range of eigenvalues to be considered
#' @param id Integer. Vector of indices of eigenvalues to consider.
#'
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
#' 
#' @family hankel
bootstrap.hankel <- function(cf, t0, n=2, N = cf$Time/2+1, id=c(1), range=c(0,1), eps=0.001) {

  stopifnot(inherits(cf, 'cf_meta'))
  stopifnot(inherits(cf, 'cf_boot'))

  boot.R <- cf$boot.R
  evs <- array(NA, dim=c(N, length(id)))
  evs.tsboot <- array(NA, dim=c(boot.R, N, length(id)))
  
  for(deltat in c(1:(N-2-t0-2*n))) {
    evs[deltat+t0, ] <- gevp.hankel(cf$cf0, t0=t0, n=n, N=N, id=id, range=range, eps=eps, debug=FALSE, deltat=deltat,
                       submatrix.size=1, element.order=c(1))
    
    evs.tsboot[, deltat+t0, ] <- t(apply(cf$cf.tsboot$t, 1, gevp.hankel, t0=t0,
                                         n=n, N=N, id=id, range=range, eps=eps, deltat=deltat,
                                         submatrix.size=1, element.order=1))
  }
  ret <- list(cf=cf,
              res.hankel=evs,
              hankel.tsboot=evs.tsboot,
              boot.R=boot.R,
              boot.l=cf$boot.l,
              seed=cf$seed,
              t0=t0,
              n=n,
              N=N,
              eps=eps,
              range=range)
  class(ret) <- c("hankel", class(ret))
  return(invisible(ret))
}


#' @title hankel2cf
#'
#' @param hankel object of type \link{hankel}
#' @param id Integer. ID of the principal correlator to extract
#' 
#' @family hankel
hankel2cf <- function(hankel, id=1) {
  stopifnot(inherits(hankel, "hankel"))
  stopifnot((id <= hankel$n && id >= 1))
  
  ## Base `cf` properties.
  cf <- cf_meta(nrObs = 1,
                Time = hankel$cf$Time,
                nrStypes = 1,
                symmetrised = hankel$cf$symmetrised)

  ## Add the `cf_boot` mixin.
  cf.tsboot <- list(t = hankel$hankel.tsboot[, , id],
                    t0 = hankel$res.hankel[,id])
  
  cf <- cf_boot(cf,
                boot.R = hankel$boot.R,
                boot.l = hankel$boot.l,
                seed = hankel$seed,
                sim = hankel$cf$sim,
                endcorr = TRUE,
                cf.tsboot = cf.tsboot,
                resampling_method = hankel$cf$resampling_method)
  
  cf <- cf_principal_correlator(cf,
                                id = id,
                                gevp_reference_time = hankel$t0)
  
  return(invisible(cf))
}

#' @title hankel2effectivemass
#'
#' @param hankel object of type \link{hankel}
#' @param id Integer. ID of the principal correlator to extract
#' @param type Character vector. Type of effective mass to use.
#' 
#' @family hankel
hankel2effectivemass  <- function(hankel, id=c(1), type="log") {
  stopifnot(inherits(hankel, "hankel"))
  stopifnot(length(id) == 1)
  stopifnot((id <= hankel$n && id >= 1))
  stopifnot(type %in% c("log", "acosh"))
  
  tmax <- hankel$cf$Time/2
  if(!hankel$cf$symmetrised){
    tmax <- hankel$cf$Time-1
  }

  nrObs <- 1

  deltat <- c(1:(tmax+1))-hankel$t0
  if(type == "log") {
    effMass <- -log(hankel$res.hankel[,id])/deltat
    effMass.tsboot <- -log(hankel$hankel.tsboot[,, id])/t(array(deltat, dim=rev(dim(hankel$hankel.tsboot[,, id]))))
  }
  if(type == "acosh") {
    effMass <- -acosh(hankel$res.hankel[,id])/deltat
    effMass.tsboot <- -acosh(hankel$hankel.tsboot[,, id])/t(array(deltat, dim=rev(dim(hankel$hankel.tsboot[,, id]))))
  }
  deffMass <- apply(effMass.tsboot, 2, hankel$cf$error_fn, na.rm=TRUE)

  ret <- list(t.idx=c(1:(tmax)),
              effMass=effMass, deffMass=deffMass, effMass.tsboot=effMass.tsboot,
              opt.res=NULL, t1=NULL, t2=NULL, type=type, useCov=NULL, CovMatrix=NULL, invCovMatrix=NULL,
              boot.R = hankel$boot.R, boot.l = hankel$boot.l, seed = hankel$seed,
              massfit.tsboot=NULL, Time=hankel$cf$Time, N=1, nrObs=nrObs, dof=NULL,
              chisqr=NULL, Qval=NULL
             )
  ret$cf <- cf
  ret$t0 <- effMass
  ret$t <- effMass.tsboot
  ret$se <- apply(ret$t, MARGIN=2L, FUN=hankel$cf$error_fn, na.rm=TRUE)
  attr(ret, "class") <- c("effectivemass", class(ret))
  return(ret)
}
