fs <- function(data) {

  L <- data[,1]
  mps <- data[,2]
  err <- data[,3]
  fps <- data[,4]
  dfps <-  data[,5]
  fit.mpi.fs <- optim(par=c(0.15,2.,1.), ChiSqr.fs, method="BFGS",
                      L=L, y=mps, err=err, control=list(maxit=1000))
  plotwitherror(fit.mpi.fs$par[1]*L, mps, err, xlab="mps*L", ylab="m_ps")
  xfit <- seq(L[1], L[length(L)], 0.05)
  yfit <- fit.mpi.fs$par[1] + fit.mpi.fs$par[2]/xfit^(3/2)*exp(-fit.mpi.fs$par[3]*fit.mpi.fs$par[1]*xfit)
  lines(spline(xfit*fit.mpi.fs$par[1], yfit), col="blue", lty=1)
  abline(h=fit.mpi.fs$par[1], col="red")
  title("FS effects in m_ps, single fit")

  fit.fpi.fs <- optim(par=c(0.15,-10.,1.), ChiSqr.fs.fps, method="BFGS",
                      L=L, y=fps, err=dfps, mpi=fit.mpi.fs$par[1], control=list(maxit=1000))
  new_window_if_appropriate()
  plotwitherror(fit.mpi.fs$par[1]*L, fps, dfps, xlab="mps*L", ylab="f_ps")
  yfit <- fit.fpi.fs$par[1] + fit.fpi.fs$par[2]/xfit^(3/2)*exp(-fit.fpi.fs$par[3]*fit.mpi.fs$par[1]*xfit)
  lines(spline(xfit*fit.mpi.fs$par[1], yfit), col="blue", lty=1)
  abline(h=fit.fpi.fs$par[1], col="red")
  title("FS effects in f_ps b=4.05, mu=0.006, single fit")

  fit.comb <- optim(par=c(fit.mpi.fs$par, fit.fpi.fs$par[1], fit.fpi.fs$par[2]), ChiSqr.fs.comb, method="BFGS",
                    L=L, mps=mps, fps=fps, dmps=err, dfps=dfps, control=list(maxit=1000))
  print(fit.comb)
  new_window_if_appropriate()
  plotwitherror(fit.mpi.fs$par[1]*L, mps, err, xlab="mps*L", ylab="m_ps")
  yfitm <- fit.comb$par[1] + fit.comb$par[2]/xfit^(3/2)*exp(-fit.comb$par[3]*fit.comb$par[1]*xfit)

  lines(spline(xfit*fit.comb$par[1], yfitm), col="blue", lty=1)
  abline(h=fit.comb$par[1], col="red")
  title("FS effects in m_ps, combined fit")

  new_window_if_appropriate()
  plotwitherror(fit.mpi.fs$par[1]*L, fps, dfps, xlab="mps*L", ylab="f_ps")
  yfitf <- fit.comb$par[4] + fit.comb$par[5]/xfit^(3/2)*exp(-fit.comb$par[3]*fit.comb$par[1]*xfit)
  lines(spline(xfit*fit.comb$par[1], yfitf), col="blue", lty=1)
  abline(h=fit.comb$par[4], col="red")
  title("FS effects in f_ps, combined fit")

  fit.comb2 <- optim(par=c(fit.mpi.fs$par[1:2]+0.2, fit.fpi.fs$par[1], fit.fpi.fs$par[2]), ChiSqr.fs.comb2, method="BFGS",
                    L=L, mps=mps, fps=fps, dmps=err, dfps=dfps)
  message("fit.com2\n")
  print(fit.comb2)
  new_window_if_appropriate()
  plotwitherror(fit.mpi.fs$par[1]*L, mps, err, xlab="mps*L", ylab="m_ps")
  yfitm <- fit.comb2$par[1] + fit.comb2$par[2]/xfit^(3/2)*exp(-fit.comb2$par[1]*xfit)

  lines(spline(xfit*fit.comb2$par[1], yfitm), col="blue", lty=1)
  abline(h=fit.comb2$par[1], col="red")
  title("FS effects in m_ps, combined fit (2)")

  new_window_if_appropriate()
  plotwitherror(fit.mpi.fs$par[1]*L, fps, dfps, xlab="mps*L", ylab="f_ps")
  yfitf <- fit.comb2$par[3] + fit.comb2$par[4]/xfit^(3/2)*exp(-fit.comb2$par[1]*xfit)
  lines(spline(xfit*fit.comb2$par[1], yfitf), col="blue", lty=1)
  abline(h=fit.comb2$par[3], col="red")
  title("FS effects in f_ps, combined fit (2)")


  fit.pow.mps <- optim(par=c(0.15,50.,3/2), ChiSqr.pow, method="BFGS",
                       L=L, y=mps, err=err)
  print(fit.pow.mps)
  new_window_if_appropriate()
  plotwitherror(fit.mpi.fs$par[1]*L, mps, err, xlab="mps*L", ylab="m_ps")
  yfit <- fit.pow.mps$par[1] + fit.pow.mps$par[2]/xfit^fit.pow.mps$par[3]
  lines(spline(xfit*fit.pow.mps$par[1], yfit), col="blue", lty=1)
  abline(h=fit.pow.mps$par[1], col="red")
  title("FS effects in m_ps, power law fit")

  fit.pow.fps <- optim(par=c(0.05,-5.,3/2), ChiSqr.pow, method="BFGS", control=list(maxit=10000),
                       L=L, y=fps, err=dfps)
  print(fit.pow.fps)
  new_window_if_appropriate()
  plotwitherror(fit.mpi.fs$par[1]*L, fps, dfps, xlab="mps*L", ylab="f_ps")
  yfit <- fit.pow.fps$par[1] + fit.pow.fps$par[2]/xfit^fit.pow.fps$par[3]
  lines(spline(xfit*fit.pow.fps$par[1], yfit), col="blue", lty=1)
  abline(h=fit.pow.fps$par[1], col="red")
  title("FS effects in f_ps, power law fit")

  return(invisible(list(data=data, res.mps=fit.mpi.fs, res.fps=fit.fpi.fs, res.comb=fit.comb, res.comb2=fit.comb2)))
}

ChiSqr.fs <- function(par, L, y, err) {
  return(
         sum(((y-par[1] - par[2]/L^(3/2)*exp(-par[3]*par[1]*L))/err)^2)
         )
}

ChiSqr.fs.fps <- function(par, L, y, err, mpi) {
  return(
         sum(((y-par[1] - par[2]/L^(3/2)*exp(-par[3]*mpi*L))/err)^2)
         )
}

ChiSqr.fs.comb <- function(par, L, mps, fps, dmps, dfps) {
  return(
         (sum(((mps-par[1] - par[2]/L^(3/2)*exp(-par[3]*par[1]*L))/dmps)^2)
          +
          sum(((fps-par[4] - par[5]/L^(3/2)*exp(-par[3]*par[1]*L))/dfps)^2))
         )
}

ChiSqr.fs.comb2 <- function(par, L, mps, fps, dmps, dfps) {
  return(
         (sum(((mps-par[1] - par[2]/L^(3/2)*exp(-par[1]*L))/dmps)^2)
          +
          sum(((fps-par[3] - par[4]/L^(3/2)*exp(-par[1]*L))/dfps)^2))
         )
}

ChiSqr.pow <- function(par, L, y, err) {
  return(
         sum(((y - par[1]-par[2]/L^par[3])/err)^2)
         )
}
