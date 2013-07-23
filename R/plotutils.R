plotwitherror <- function(x, y, dy, ylim, dx, xlim, rep=FALSE, col="black", ...) {
  my.xlim <- NULL
  my.ylim <- NULL
  
  if(missing(xlim)) {
    if(!missing(dx)) {
      my.xlim <- c(min(x-2*dx, na.rm = TRUE),max(x+2*dx, na.rm = TRUE))
    }
  } else {
    my.xlim <- xlim
  }

  if(missing(ylim)) {
    if(!missing(dy)) {
      my.ylim <- c(min(y-2*dy, na.rm = TRUE),max(y+2*dy, na.rm = TRUE))
    }
  } else {
    my.ylim <- ylim
  }

  if(rep) {
    points(x, y, col=col, ...)
  }
  else if(!is.null(my.xlim) && !is.null(my.ylim)) {
    plot(x, y, ylim=my.ylim, xlim=my.xlim, col=col, ...)
  } 
  else if(is.null(my.xlim) && !is.null(my.ylim)) {
    plot(x, y, ylim=my.ylim, col=col, ...)
  }
  else if(!is.null(my.xlim) && is.null(my.ylim)) {
    plot(x, y, xlim=my.xlim, col=col, ...)
  }
  else {
    plot(x, y, col=col, ...)
  }

  options(show.error.messages = FALSE)
  if(!missing(dy)) {
    arrows(x, y-dy, x, y+dy, length=0.01,angle=90,code=3, col=col)
  }
  if(!missing(dx)) {
    arrows(x-dx,y,x+dx,y, length=0.01,angle=90,code=3, col=col)
  }
  options(show.error.messages = TRUE)
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
    abline(h=fit$fitresult$par[3]*fit$fitresult$par[2]/fit$fitresult$par[1]/2., col="red")
  }
  if(!is.null(fit$MChist.dpaopp)) {
    plot(seq(1, length(fit$MChist.dpaopp))+fit$skip, fit$MChist.dpaopp, type="l",
         main=expression(m[PCAC]), xlab=expression(t[HMC]), ylab=expression(m[PCAC]))
    abline(h=fit$fitresult$par[3]*fit$fitresult$par[2]/fit$fitresult$par[1]/2., col="red")
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
  else if(!missing(lf)) {
    plot.massfit(ll, ylab=expression(m[eff]), xlab="t")
    points((lf$t-0.2), lf$mass, pch=1, col="blue")
    arrows((lf$t-0.2), lf$mass-lf$dmass,
           (lf$t-0.2), lf$mass+lf$dmass, length=0.01,angle=90,code=3)
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

plot.outputdata <- function(data, skip=0, ...) {
  plaq.res <- uwerrprimary( data$V2[skip:length(data$V2)])
  dH.res <- uwerrprimary( exp(-data$V3[skip:length(data$V3)]))
  plot(data$V1, data$V2, type="l",
       main="Plaquette", xlab=expression(t[HMC]), ylab="P", col="red", ...)
  abline(h=plaq.res$value, col="blue")
  abline(h=plaq.res$value+plaq.res$dvalue)
  abline(h=plaq.res$value-plaq.res$dvalue)
  X11()
  plot(data$V1, data$V3, type="l",
       main=expression(paste(Delta, "H")), xlab=expression(t[HMC]), ylab=expression(paste(Delta, "H")), ...)
  return(invisible(list(data=data, plaq.res=plaq.res, dH.res = dH.res)))
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
