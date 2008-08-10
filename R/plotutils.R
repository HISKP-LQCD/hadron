# $Id$
plotwitherror <- function(x, y, dy, ylim, rep=FALSE, col="black", ...) {
  if(missing(ylim)) {
    if(rep) {
      points(x, y, col=col, ...)
    }
    else {
      plot(x,y, ylim=c(min(y-2*dy, na.rm = TRUE),max(y+2*dy, na.rm = TRUE)), col=col, ...)
    }
  }
  else {
    plot(x,y, ylim=ylim, col=col, ...)
  }
  arrows(x, y-dy, x, y+dy, length=0.01,angle=90,code=3, col=col)
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

plot.b1fit <- function(fit) {
  plot.cfit(fit)
}

plot.ofit <- function(fit) {
  plot.cfit(fit)
}

plot.cfit <- function(fit) {
  fit.mass <- abs(fit$fitresult$par[fit$matrix.size+1])
  if(!is.null(fit$effmass$mll)) {
    plot.effmass(m=fit.mass,
                 ll=data.frame(t=fit$effmass$t, mass=fit$effmass$mll, dmass=fit$effmass$dmll),
                 lf=data.frame(t=fit$effmass$t, mass=fit$effmass$mlf, dmass=fit$effmass$dmlf),
                 ff=data.frame(t=fit$effmass$t, mass=fit$effmass$mff, dmass=fit$effmass$dmff))
  }
  if(!is.null(fit$effmass$m)) {
    plot.effmass(m=fit.mass, ll=data.frame(t=fit$effmass$t, mass=fit$effmass$m, dmass=fit$effmass$dm))
  }
  if(!is.null(fit$uwerrresultmpcac)) {
    plot(fit$uwerrresultmpcac, main=expression(m[PCAC]))
  }
  if(!is.null(fit$uwerrresultmps)) {
    plot(fit$uwerrresultmps, main=expression(m[PS]))
  }
  if(!is.null(fit$uwerrresultfps)) {
    plot(fit$uwerrresultfps, main=expression(f[PS]))
  }
  if(!is.null(fit$uwerrresultmv)) {
    plot(fit$uwerrresultmv, main=expression(m[V]))
  }
  if(!is.null(fit$boot)) {
    X11()
    plot(fit$boot, main="Bootstrap analysis for mps")
  }
  if(!is.null(fit$tsboot)) {
    X11()
    plot(fit$tsboot, main="TS Boostrap analysis for mps")
  }
  if(!is.null(fit$mv.boot)) {
    X11()
    plot(fit$mv.boot, main="Bootstrap analysis for mv")
  }
  if(!is.null(fit$mv.tsboot)) {
    X11()
    plot(fit$mv.tsboot, main="TS Bootstrap analysis for mv")
  }
  if(!is.null(fit$fitdata)) {
    X11()
    plot(fit$fitdata$Chi, main="Chi per data point", xlab="data point", ylab=expression(chi), type="h")
  }
  if(!is.null(fit$dpaopp)) {
    X11()
    plotwitherror(fit$dpaopp$t, fit$dpaopp$mass, fit$dpaopp$dmass,
                  main=expression(m[PCAC]), xlab="t", ylab=expression(m[PCAC]))
    abline(h=fit$fitresult$par[3]*fit$fitresult$par[2]/fit$fitresult$par[1]/2.)
  }
}

plot.correlator <- function(data, xlab = "t", ylab = "C(t)", log="y", ...) {
  plotwitherror(data$t,data$corr, data$dcorr, xlab=xlab, ylab=ylab, log=log, ...)
}

plot.effmass <- function(m, ll, lf, ff, ...) {

  X11()
  if(!missing(ff)) {
    plot.massfit(ff, ylab=expression(m[eff]), xlab="t")
    points((ll$t-0.2), ll$mass, pch=1, col="blue")
    arrows((ll$t-0.2), ll$mass-ll$dmass,
           (ll$t-0.2), ll$mass+ll$dmass, length=0.01,angle=90,code=3)
    points((lf$t+0.2), lf$mass, pch=2, col="red")
    arrows((lf$t+0.2), lf$mass-lf$dmass,
           (lf$t+0.2), lf$mass+lf$dmass, length=0.01,angle=90,code=3)
    lines(ll$t, rep(m, times=length(ll$t)))
  }
  else {
    plot.massfit(ll, ylab=expression(m[eff]), xlab="t")
    lines(ll$t, rep(m, times=length(ll$t)))
  }
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

#postscript(file = "test.eps", family="NimbusRom", paper="special", horizontal=FALSE, onefile=FALSE, width=6.996766, height=6.996766)
# linewidth, tick length and direction, margings, orientation
#opar <- par(lwd=0.5, tck=0.02, mgp=c(2,.5,0), las=1)
#plot(b, axes="F")
#main axis
#axis(1, lwd=0.5)
#axis(2, lwd=0.5)
#sub ticks
#axis(1, at=seq(0,100,10), tck=0.01, labels=F, lwd=0.5)
#axis(2, at=seq(0,100,10), tck=0.01, labels=F, lwd=0.5)
#box()
#dev.off()
#par(opar)
