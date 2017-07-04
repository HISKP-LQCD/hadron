getfit.boot <- function(Z, d, Err, t1, t2, Time, par=c(1.,0.12),
                     fit.routine="optim", sign) {
  Thalf <- Time/2
  T1 <- Thalf+1
  t1p1 <- (t1+1)
  t2p1 <- (t2+1)
  tr <- (t2-t1+1)
  Cor <- rep(0., times=length(Z[1,]))
  if(!missing(d)) {
    for(i in 1:length(Z[1,])) {
      Cor[i] = mean(Z[d,(i)])
    }
  }
  else {
    for(i in 1:length(Z[1,])) {
      Cor[i] = mean(Z[,(i)])
    }
  }

  fit <- optim(par, ChiSqr.singleCor, method="BFGS", Thalf=Thalf,
               x=c((t1):(t2)), y=Cor, err=Err, tr=tr, sign=sign)
  sort.ind <- c(1)
  return(c(abs(fit$par[2]), fit$par[1],
           fit$value))
}
