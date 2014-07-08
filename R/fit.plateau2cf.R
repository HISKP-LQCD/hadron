fit.plateau2cf <- function(cf, t1, t2,
                           boot.samples=FALSE, boot.R=400, boot.l=2,
                           useCov=FALSE) {

  if(!cf$boot.samples) {
    cf <- bootstrap.cf(cf, boot.R=boot.R, boot.l=boot.l)
  }
  else {
    boot.R <- cf$boot.R
    boot.l <- cf$boot.l
  }
  ## fit interval
  ii <- c((t1+1):(t2+1))
  ## error weights
  w <- 1/apply(cf$cf.tsboot$t[,ii], 2, sd)

  ## here we generate the inverse covariance matrix, if required
  ## otherwise take inverse errors squared
  M <- diag(w^2)

  if(useCov) {
    ## compute correlation matrix and compute the correctly normalised inverse
    ##M <- invertCovMatrix(cf$cf.tsboot$t[,ii], boot.samples=TRUE)
    M <- invertCovMatrix(cf$cf[,ii], boot.samples=boot.samples, boot.l=boot.l)
  }
  fn <- function(par, y, M) { sum((y-par[1]) %*% M %*% (y-par[1]))}

  par <- cf$cf0[Cf2pt$Time/4]
  opt.res <- optim(par, fn = fn,
                   method="BFGS", M=M, y = cf$cf0[ii])
  opt.res <- optim(opt.res$par, fn = fn,
                   control=list(parscale=1/opt.res$par),
                   method="BFGS", M=M, y = cf$cf0[ii])
  par <- opt.res$par
  plateau <- par[1]
  chisqr <- opt.res$value
  dof <- length(ii)-1
  plateau.tsboot <- array(NA, dim=c(boot.R,2))
  for(i in 1:boot.R) {
    opt <- optim(par, fn = fn,
                 control=list(parscale=1/par),
                 method="BFGS", M=M, y = cf$cf.tsboot$t[i,ii])
    plateau.tsboot[i,1] <- opt$par[1]
    plateau.tsboot[i,2] <- opt$value
  }

  res <- list(plateau=plateau, dplateau=sd(plateau.tsboot[,1]),
              chisqr=chisqr, plateau.tsboot=plateau.tsboot,
              dof=dof, cf=cf, t1=t1, t2=t2, boot.R=boot.R, boot.l=boot.l, useCov=useCov,
              invCovMatrix=M)
  return(invisible(res))
}
