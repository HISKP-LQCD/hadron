fitmass <- function(Cor, Err, t1, t2, Time, par=c(1.,0.12), sign,
                    fit.routine="optim") {
  Thalf <- Time/2
  T1 <- Thalf+1
  t1p1 <- (t1+1)
  t2p1 <- (t2+1)
  tr <- (t2-t1+1)
  
  if(fit.routine != "gsl") {
    fit <- optim(par, ChiSqr.singleCor, method="BFGS", Thalf=Thalf,
                 x=c((t1):(t2)), y=Cor, err=Err, tr=tr, sign=sign)
  }
  else {
    fit <- gsl <- fit <- correlator <- matrix(par, Thalf=Thalf,
                                              x=c((t1):(t2)), y=Cor, err=Err, tr = tr, sign=sign)
  }
  
  return(abs(fit$par[2]))
}
