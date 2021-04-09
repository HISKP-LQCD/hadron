#' fits a plateau to an object of class \code{cf}
#' 
#' where applicable, a plateau is fitted to the averaged data in \code{cf}
#' using a (correlated) chisquare fit.
#' 
#' 
#' @param cf input object of class \code{cf}
#' @param t1 starting t-value for the fit
#' @param t2 final t-value for the fit.
#' @param useCov perform a correlated chisquare fit or not.
#' @return Returns a list with elements \item{plateau}{ the fitted plateau
#' value } \item{dplateau}{ its error }
#' @author Carsten Urbach \email{curbach@@gmx.de}
#' @seealso \code{\link{cf}}
#' @keywords bootstrap fit
#' @examples
#'
#' data(correlatormatrix)
#' cfnew <- extractSingleCor.cf(correlatormatrix, id=1)
#' cfnew <- bootstrap.cf(cfnew, boot.R=99, boot.l=1)
#' X <- fit.plateau2cf(cfnew, t1=13, t2=20)
#' @export fit.plateau2cf
fit.plateau2cf <- function(cf, t1, t2, useCov=FALSE) {
  stopifnot(inherits(cf, 'cf_meta'))
  stopifnot(inherits(cf, 'cf_boot'))

  boot.R <- cf$boot.R
  boot.l <- cf$boot.l
  ## fit interval
  ii <- c((t1+1):(t2+1))
  ## error weights
  w <- 1/apply(cf$cf.tsboot$t[,ii], 2, sd)

  ## here we generate the inverse covariance matrix, if required
  ## otherwise take inverse errors squared
  M <- diag(w^2)

  if(useCov) {
    ## compute correlation matrix and compute the correctly normalised inverse
    M <- invertCovMatrix(cf$cf.tsboot$t[,ii], boot.samples = TRUE)
  }
  fn <- function(par, y, M) { sum((y-par[1]) %*% M %*% (y-par[1]))}

  par <- cf$cf0[cf$Time/4]
  opt.res <- optim(par, fn = fn, lower=cf$cf0[cf$Time/4]-4*cf$tsboot.se[cf$Time/4], upper=cf$cf0[cf$Time/4]+4*cf$tsboot.se[cf$Time/4],
                   method="Brent", M=M, y = cf$cf0[ii])
  opt.res <- optim(opt.res$par, fn = fn, lower=cf$cf0[cf$Time/4]-4*cf$tsboot.se[cf$Time/4], upper=cf$cf0[cf$Time/4]+4*cf$tsboot.se[cf$Time/4], 
                   control=list(parscale=1/opt.res$par),
                   method="Brent", M=M, y = cf$cf0[ii])
  par <- opt.res$par
  plateau <- par[1]
  chisqr <- opt.res$value
  dof <- length(ii)-1
  plateau.tsboot <- array(NA, dim=c(boot.R,2))
  for(i in 1:boot.R) {
    opt <- optim(par, fn = fn, lower=cf$cf0[cf$Time/4]-4*cf$tsboot.se[cf$Time/4], upper=cf$cf0[cf$Time/4]+4*cf$tsboot.se[cf$Time/4], 
                 control=list(parscale=1/par),
                 method="Brent", M=M, y = cf$cf.tsboot$t[i,ii])
    plateau.tsboot[i,1] <- opt$par[1]
    plateau.tsboot[i,2] <- opt$value
  }

  dplateau <- sd(plateau.tsboot[,1])
  res <- list(plateau=plateau, dplateau=dplateau, plateau.tsboot=plateau.tsboot,
              t0=plateau, se=dplateau, t=plateau.tsboot,
              chisqr=chisqr, dof=dof, 
              cf=cf, t1=t1, t2=t2, boot.R=boot.R, boot.l=boot.l, useCov=useCov,
              invCovMatrix=M)
  return(invisible(res))
}
