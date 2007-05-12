# $Id$
plotwitherror <- function(x, y, dy, ylim, ...) {
  if(missing(ylim)) {
    plot(x,y, ylim=c(min(y-2*dy),max(y+2*dy)), ...)
  }
  else {
    plot(x,y, ylim=ylim, ...)
  }
  arrows(x, y-dy, x, y+dy, length=0.01,angle=90,code=3)
}

plot.massfit <- function(data, xlab = "t", ylab = "m", ...) {
  plotwitherror(data$t,data$mass, data$dmass, xlab=xlab, ylab=ylab, ...)
}

plot.pionfit <- function(fit) {
  plot.cfit(fit)
}

plot.rhofit <- function(fit) {
  plot.cfit(fit)
}

plot.cfit <- function(fit) {
  fit.mass <- abs(fit$fitresult$par[fit$matrix.size+1])
  plot.effmass(m=fit.mass,
               ll=data.frame(t=fit$effmass$t, mass=fit$effmass$mll, dmass=fit$effmass$dmll),
               lf=data.frame(t=fit$effmass$t, mass=fit$effmass$mlf, dmass=fit$effmass$dmlf),
               ff=data.frame(t=fit$effmass$t, mass=fit$effmass$mff, dmass=fit$effmass$dmff))
  
  if(!is.null(fit$uwerrresultmpcac)) {
    plot(fit$uwerrresultmpcac)
  }
  if(!is.null(fit$uwerrresultmps)) {
    plot(fit$uwerrresultmps)
  }
  if(!is.null(fit$uwerrresultfps)) {
    plot(fit$uwerrresultfps)
  }
  if(!is.null(fit$uwerrresultmv)) {
    plot(fit$uwerrresultmv)
  }
  if(!is.null(fit$boot)) {
    X11()
    plot(fit$boot)
  }
  if(!is.null(fit$tsboot)) {
    X11()
    plot(fit$tsboot)
  }
  if(!is.null(fit$mv.boot)) {
    X11()
    plot(fit$mv.boot)
  }
  if(!is.null(fit$mv.tsboot)) {
    X11()
    plot(fit$mv.tsboot)
  }
}

plot.correlator <- function(data, xlab = "t", ylab = "C(t)", log="y", ...) {
  plotwitherror(data$t,data$corr, data$dcorr, xlab=xlab, ylab=ylab, log=log, ...)
}

plot.effmass <- function(m, ll, lf, ff, ...) {

  X11()
  plot.massfit(ff, ylab="m_eff", xlab="t")
  points((ll$t-0.2), ll$mass, pch=1, col="blue")
  arrows((ll$t-0.2), ll$mass-ll$dmass,
         (ll$t-0.2), ll$mass+ll$dmass, length=0.01,angle=90,code=3)
  points((lf$t+0.2), lf$mass, pch=2, col="red")
  arrows((lf$t+0.2), lf$mass-lf$dmass,
         (lf$t+0.2), lf$mass+lf$dmass, length=0.01,angle=90,code=3)
  lines(ll$t, rep(m, times=length(ll$t)))
}


plot.averx <- function(averx, ...) {
  X11()
  plotwitherror(c(0:(length(averx$data$Cor)-1)), averx$data$Cor/averx$mps, averx$data$Err/averx$mps, ...)
  abline(h=averx$averx, col="red")
  abline(h=averx$averx+averx$daverx, col="blue")
  abline(h=averx$averx-averx$daverx, col="blue")
  if(!is.null(averx$fit.uwerr)) {
    plot(averx$fit.uwerr)
  }
}
