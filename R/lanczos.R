

#' @title Lanczos method for LQCD correlators
#' 
#' @description
#'   Taking a single correlation function as input, the method
#'   determines the ground state energy plus its bootstrap uncertainty.
#'
#' @param cf object of type \link{cf}, optimally returned by
#'   \code{\link{bootstrap.cf}}
#' @param N Integer. Maximal time index in correlation function to be used in
#'                   Lanczos analysis
#' @param bias_correction boolean. If set to 'TRUE', the median of the bootstrap
#'   distribution is used as estimator for the energy values.
#' @param errortype string. Determines the treatment of the bootstrap
#'   histograms to determine the statistical error on eigenvalues. Can
#'   be: 1. 'outlier-removal' for which outliers are removed according to
#'   the 0.25 and 0.75 quantiles and the inter-quantile-range,
#'   i.e. only values are kept which are in the interval
#'   [Q_25-1.5IQR, Q_75+1.5IQR]
#'   and the error is computed from the standard deviation of the bootstrap distribution.
#'   2. 'quantiles' for which the error is estimated from the difference
#'   between the 0.32 and 0.68 quantile of the original bootstrap distribution
#' @param pivot boolean. If set to 'TRUE', the eigenvalues on the original data are used
#'   to find the "correct" eigenvalue on the bootstrap sample by the
#'   smallest distance.
#' @param probs numeric. Vector of probabilities for the error estimation method
#'   'quantiles'.
#' @seealso \code{\link{plot.effectivemass}}, \code{\link{bootstrap.effectivemass}}
#' @references M. Wagman, 'Lanczos, the transfer matrix, and the signal-to-noise problem',
#'   arXiv:2406.20009 
#' @return
#'   Returns an object of S3 class `effectivemass`.
#' 
#' @family lanczos
#' @export
#' @examples
#' data(pscor.sample)
#' newcf <- cf_orig(cf=t(array(pscor.sample[,2], dim=c(48, 316))))
#' newcf <- cf_meta(newcf, nrObs=1, Time=48, symmetrised=FALSE)
#' newcf.boot <- bootstrap.cf(newcf)
#' ncf.boot <- symmetrise.cf(newcf.boot)
#' ncf.effmass <- bootstrap.effectivemass(ncf.boot)
#' plot(ncf.effmass, ylim=c(0.1,0.2))
#' res <- bootstrap.lanczos(newcf.boot, N=newcf$Time)
#' plot(res, rep=TRUE, col="red", pch=22, xshift=0.2)
bootstrap.lanczos <- function(cf, N = (cf$Time/2+1), bias_correction=FALSE, errortype="outlier-removal", pivot=FALSE, probs=c(0.32,0.68)) {
  ## wrapper function, not yet bootstrapping...
  stopifnot(inherits(cf, 'cf_meta'))
  stopifnot(inherits(cf, 'cf_boot'))
  stopifnot(inherits(cf, 'cf_orig'))
  stopifnot(errortype %in% c("outlier-removal", "quantiles"))
  
  seed <- cf$seed
  boot.R <- cf$boot.R
  boot.l <- cf$boot.l

  res <- lanczos.solve(cf=cf$cf0, N=N, pivot=FALSE)
  lanczos.tsboot.orig <- t(apply(cf$cf.tsboot$t, 1, lanczos.solve, N=N, pivot=pivot, pivot_elements=res))
  lanczos.tsboot <- lanczos.tsboot.orig
  effMass <- -log(res)
  deffMass <- rep(NA, length(effMass))
  effMassMedian <- deffMass
  if(errortype=="outlier-removal") {
    remove_outliers <- function(x, probs=c(0.25,0.75)) {
      Q <- quantile(x, probs=probs, na.rm=TRUE)
      iqr <- Q[2]-Q[1]
      x[x<(Q[1]-1.5*iqr) | x > (Q[2] + 1.5*iqr)] <- NA
      return(invisible(x))
    }
    lanczos.tsboot <- apply(lanczos.tsboot.orig, 2, remove_outliers)
    deffMass <- apply(-log(lanczos.tsboot), 2L, cf$error_fn, na.rm=TRUE)
  }
  else if(errortype == "quantiles") {
    error_fn <- function(x, probs=c(0.32, 0.68)) {
      Q <- quantile(x, probs=probs)
      return(Q[2]-Q[1])
    }
    deffMass <- apply(-log(lanczos.tsboot), 2L, error_fn, probs=probs)
  }
  bias <- effMass - apply(-log(lanczos.tsboot), 2L, median, na.rm=TRUE)
  if(bias_correction) {
    effMass <- effMass - bias
  }
  ret <- list(t.idx=c(1:(length(res))), cf=cf, res.lanczos=res, bias=bias,
              lanczos.tsboot.orig=lanczos.tsboot.orig, lanczos.tsboot=lanczos.tsboot,
              effMass=effMass, deffMass=deffMass, effMass.tsboot=-log(lanczos.tsboot),
              opt.res=NULL, t1=NULL, t2=NULL, type="log", useCov=NULL, CovMatrix=NULL, invCovMatrix=NULL,
              boot.R = boot.R, boot.l = boot.l, seed = seed,
              massfit.tsboot=NULL, Time=cf$Time, nrObs=1, dof=NULL,
              chisqr=NULL, Qval=NULL
             )
  ret$t0 <- effMass
  ret$t <- ret$effMass.tsboot
  ret$se <- deffMass
  attr(ret, "class") <- c("effectivemass", "lanczos", class(ret))
  return(invisible(ret))
}

#' @title Lanczos solver
#' 
#' @description
#'   blub ...
#'
#' @param cf Numeric vector (this will generally be a correlation function or a bootstrap sample thereof).
#' @param N Integer. Maximal time index in correlation function to be used in
#'                   Lanczos analysis
#' @return
#'   tbw
#' 
#' @family lanczos
lanczos.solve <- function(cf, N, pivot=FALSE, pivot_elements=NULL) {
  ## container for the eigenvalues per m
  evs <- rep(NA, times=N/2)
  for(m in c(1:(N/2))) {
    Aj <- cf[2:N]/cf[1]
    Bj <- rep(0, times=length(Aj))
    Bjp1 <- Bj
    Gj <- Bj
    Gjp1 <- Gj
    Ajm1 <- Bj
    Ajp1 <- Bj
    
    alpha <- rep(NA, times=m)
    beta <- alpha
    gamma <- alpha
    alpha[1] <- Aj[1]
    beta[1] <- 0.
    gamma[1] <- 0.

    if(m > 1) {
      for(j in c(1:(m-1))) {
        ## eq.(53) suppl. mat.
        sr <- Aj[2]-alpha[j]^2-beta[j]*gamma[j]
        ## eq.(41) suppl. mat.
        rhojp1 <- sqrt(abs(sr))
        taujp1 <- sr/rhojp1
        ## below eq.(28) suppl. mat.
        kmax <- 2*(m-j) + 1
        ## B_{j+1}^k and G_{j+1}^k
        for(k in c(1:kmax)) {
          ## eqs. (44) and (45) suppl. mat.
          Gjp1[k] <- Aj[k+1] - alpha[j]*Aj[k] - gamma[j]*Bj[k]
          Bjp1[k] <- Aj[k+1] - alpha[j]*Aj[k] - beta[j]*Gj[k]
        }
        Gjp1 <- Gjp1/taujp1
        Bjp1 <- Bjp1/rhojp1
        ## A_{j+1}^k eq. (46) suppl. mat.
        for(k in c(1:kmax)) {
          Ajp1[k] <-  Aj[k+2] - 2*alpha[j]*Aj[k+1] + alpha[j]^2*Aj[k] +
            alpha[j]*(beta[j]*Gj[k] + gamma[j]*Bj[k]) -
            (beta[j]*Gj[k+1] + gamma[j]*Bj[k+1]) +
            gamma[j]*beta[j]*Ajm1[k]
        }
        Ajp1 <- Ajp1/rhojp1/taujp1
        ## eq. (38) suppl. mat.
        alpha[j + 1] <- Ajp1[1]
        beta[j + 1] <- Bjp1[1]
        gamma[j + 1] <- Gjp1[1]
        ## don't do the copying if not needed
        if(j == m-1) break
        Ajm1 <- Aj
        Aj <- Ajp1
        Bjm1 <- Bj
        Bj <- Bjp1
        Gjm1 <- Gj
        Gj <- Gjp1
      }
    }
    if(m == 1) {
      ## 1x1 case
      evs[m] <- alpha[1]
    }
    else{
      ##cat("m= ", m, "\n\n")
      ## construct mxm tridiagonal matrix
      M <- diag(alpha[1:m])
      ii <- c(1:(m-1))
      M[row(M) - col(M) == -1] <- beta[2:m]
      M[row(M) - col(M) == +1] <- gamma[2:m]
      ##print(M)
      
      ## eigensolve and extract the lowest eigenvalue
      eigvalues <- eigen(M, symmetric=FALSE, only.values = TRUE, EISPACK = FALSE)$values
      ##cat("\n")
      eigvalues <- Re(eigvalues[Im(eigvalues) == 0])
      eigvalues <- eigvalues[eigvalues > 0 & eigvalues < 1]
      ##if(m==10)
      ##  print(eigvalues)
      if(length(eigvalues) == 0) eigvalues[1] <- 0.01
      if(pivot) {
        ## find the eigenvalue closest to pivot_elements[m]
        evs[m] <- eigvalues[which.min(abs(eigvalues-pivot_elements[m]))]
      }
      else{
        evs[m] <- max(eigvalues, na.rm=TRUE)
      }
    }
  }
  return(evs)
}

