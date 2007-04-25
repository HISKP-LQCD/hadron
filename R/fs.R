fs <- function(data) {

  L <- data[,1]
  mps <- data[,2]
  err <- data[,3]
  fit.fs <- optim(par=c(0.15,1.,1.), ChiSqr.fs, method="BFGS", L=L, y=mps, err=err)
  plotwitherror(fit.fs$par[1]*L, mps, err, xlab="mps*L", ylab="FSrel")
  xfit <- seq(L[1], L[length(L)], 0.05)
  yfit <- fit.fs$par[1] + fit.fs$par[2]/xfit^(3/2)*exp(-fit.fs$par[3]*fit.fs$par[1]*xfit)
  lines(spline(xfit*fit.fs$par[1], yfit), col="blue", lty=1)
  abline(h=fit.fs$par[1], col="red")
  title("FS effects b=4.05, mu=0.006")
  return(invisible(list(data=data, res=fit.fs)))
}

ChiSqr.fs <- function(par, L, y, err) {
  return(
         sum(((y-par[1] - par[2]/L^(3/2)*exp(-par[3]*par[1]*L))/err)^2)
         )
}
