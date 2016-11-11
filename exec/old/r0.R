## r0/a = (r0/a)_chiral * [1 + C1 r0 * mR + C2 (r0 * mR)^2 + (D / (r0/
## a)_chiral)^2 + D1 a \mu / (r0/a)_chiral + D2 (a \mu)^2]

## par[1:N]       r0^chiral
## par[N+1:2*N]   ZP


r0fit <- function(npar, mdata, ZPdata, ii, expr, boot.R=10, seed=123456) {

  parguess <- c(4.56, 5.36, 6.72, 8.37, 0.43, 0.433, 0.453, 0.462)
  N <- length(ii)
  debug <- TRUE
  mini <- optim(par=c(parguess, rep(0., times=(npar-length(parguess)))), fn=r0chisqr, method="BFGS", hessian=TRUE,
                control=list(maxit=150, trace=debug, REPORT=50),
                mdata=mdata, ii=ii, ZPdata=ZPdata, expr=expr) 
  mini <- optim(par=mini$par, fn=r0chisqr, method="BFGS", hessian=TRUE,
                control=list(maxit=500, trace=debug, parscale=mini$par+0.01, REPORT=50),
                mdata=mdata, ii=ii, ZPdata=ZPdata, expr=expr)

  dof <- 0

  for(i in 1:N) {
    dof <- dof + length( na.omit( mdata[[i]]$r0a ) )
  }
  dof <- dof + N - npar
  
  set.seed(seed)
  ZPbootsamples <- array(0., dim=c(boot.R,N))
  for(i in 1:N) {
    ZPbootsamples[,i] <- rnorm(boot.R, mean=ZPdata$ZP[i], sd=ZPdata$dZP[i])
  }

  boots <- array(0., dim=c(boot.R, length(mini$par)))
  for(s in 1:boot.R) {
    df <- list(data.frame(mu=mdata[[1]]$mu,
                          r0a=rnorm(length(mdata[[1]]$r0a), mean=mdata[[1]]$r0a, sd=mdata[[1]]$dr0a),
                          dr0a=mdata[[1]]$dr0a))
    ZPdf <- data.frame(ZP=ZPbootsamples[s,], dZP=ZPdata$dZP[1:N])
    for(i in 2:N) {
      df[[i]] <- data.frame(mu=mdata[[i]]$mu,
                            r0a=rnorm(length(mdata[[i]]$r0a), mean=mdata[[i]]$r0a, sd=mdata[[i]]$dr0a),
                            dr0a=mdata[[i]]$dr0a)
    }
    mini.boot <- optim(par=mini$par, fn=r0chisqr, method="BFGS", hessian=TRUE,
                       control=list(maxit=500, trace=0, parscale=mini$par+0.01, REPORT=50),
                       mdata=df, ii=ii, ZPdata=ZPdf, expr=expr)
    boots[s,] <- mini.boot$par
  }
  boot.result <- data.frame(res=mini$par, dres=apply(boots, 2, sd))
  return(invisible(list(npar=npar, par=mini$par, mini=mini, boots=boots,
                        mdata=mdata, ZPdata=ZPdata, expr=expr, ii=ii, boot.R=boot.R,
                        chisqr=mini$value, dof=dof, seed=seed, boot.result=boot.result)))
}


r0chisqr <- function(par, expr, mdata, ZPdata, ii) {

  N <- length(ii)
  env <- new.env()
  assign("par", par, env)
  assign("N", N, env)
  chisum <- 0.
  for( i in 1:N) {
    assign("i", i, env)
    assign("mu", mdata[[i]]$mu[ii[[i]]], env)
    r0 <- eval(expr, envir=env)
    chisum <- chisum + sum(((mdata[[i]]$r0a[ii[[i]]]-r0)/mdata[[i]]$dr0[ii[[i]]])^2, na.rm=TRUE) +
      ((ZPdata$ZP[i]-par[N+i])/ZPdata$dZP[i])^2
  }
  rm(env)
  return(chisum)
}

average.r0fit <- function(list.fits, ii = c(1,2,3,4), av.weight=TRUE) {

  if(is.null(list.fits)) {
    list.fits <- dir(paste(".", sep=""), pattern="fit*", recursive=TRUE, full.names=TRUE)
  }
  
  e1 <- new.env()
  nl <- length(list.fits)
  
  ## first we copute the weights for each fit in list.fits
  ## we also build a list of all the fits -> fitlist
  weights <- rep(1., times=nl)
  i <- 1
  load(list.fits[i], e1)
  fit.cur <- get(ls(e1), e1)
  N <- length(fit.cur$mdata)
  boot.R <- fit.cur$boot.R
  if(av.weight) weights[i] <- 1-pgamma(fit.cur$chisqr/2, fit.cur$dof/2)
  fitlist <- list(fit.cur)
  rm(list=ls(e1), envir=e1)
  for(i in 2:nl) {
    load(list.fits[i], e1)
    fit.cur <- get(ls(e1), e1)
    if(boot.R > fit.cur$boot.R) {
      boot.R <- fit.cur$boot.R
    }
    cat(fit.cur$chisqr/2, fit.cur$dof/2, "\n")
    if(av.weight) weights[i] <- 1-pgamma(fit.cur$chisqr/2, fit.cur$dof/2)
                                        #    weights[i] <- 1./pgamma(fit.cur$result$chisqr/2, fit.cur$result$dof/2)
    fitlist[[i]] <- fit.cur
    rm(list=ls(e1), envir=e1)
  }
  cat(weights, "\n")
  ## bres will hold the results for each bootstrap sample
  bres <- array(0., dim=c(boot.R, length(ii)))
  tmp <- numeric(nl)
  
  for(j in 1:length(ii)) {
    for(i in 1:boot.R) {
      for(k in 1:nl) {
        tmp[k] <- fitlist[[k]]$boots[i, ii[j]]
      }
      bres[i, j] <- weighted.median(x = tmp, w=weights)
    }
  }
  nr <- length(ii)

  res <- numeric(nr)
  histres <- array(0., dim=c(nr, nl))
  for(i in 1:nr) {
    for(j in 1:nl) {
      histres[i,j] <- fitlist[[j]]$par[i]
    }
    res[i] <- weighted.median(histres[i,], w=weights, na.rm=TRUE)
  }
  for(i in 1:length(ii)) {
    cat("r0[", i, "]\t", signif(res[i], digits=4), "\t+-", signif(sd(bres[,i], na.rm=TRUE), digits=4), "\t+", signif(weighted.quantile(histres[i,], w=weights, prob=c(0.8427), na.rm=TRUE)-res[i], digits=4), "\t-", signif(-weighted.quantile(histres[i,], w=weights, prob=c(0.1573), na.rm=TRUE)+res[i], digits=4), "\t bias:", res[i]-mean(bres[,i], na.rm=TRUE), "\n")
  }
  for(i in 1:length(ii)) {
    printtab(res[i], sd(bres[,i], na.rm=TRUE))
  }
  cat("---\n")
  for(i in 1:length(ii)) {
    fac <- 1.
    printtab(res[i], sqrt(var(bres[,i], na.rm=TRUE) + (weighted.quantile(histres[i,], w=weights, prob=c(0.8427), na.rm=TRUE)-res[i])^2 + (-weighted.quantile(histres[i,], w=weights, prob=c(0.1573), na.rm=TRUE)+res[i])^2), c=fac)
  }
  
  return(invisible(list(bootres=bres, res=res, histres=histres, weights=weights)))

}
