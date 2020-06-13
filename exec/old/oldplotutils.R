#' plot.pionfit
#'
#' @description
#' Generic function to plot an object of type `pionfit`
#'
#' @param x Object of type `pionfit`
#' @param ... Generic graphical parameter, ignored.
#'
#' @return
#' See \link{plot.cfit}
#' 
#' @export
plot.pionfit <- function(x, ...) {
  plot.cfit(x)
}

#' plot.rhofit
#'
#' @description
#' Generic function to plot an object of type `rhofit`
#'
#' @param x Object of type `rhofit`
#' @param ... Generic graphical parameter to be passed on to \link{plotwitherror}
#'
#' @return
#' See \link{plot.cfit}
#' 
#' @export
plot.rhofit <- function(x, ...) {
  plot.cfit(x)
}

#' plot.b1fit
#'
#' @description
#' Generic function to plot an object of type `b1fit`
#'
#' @param x Object of type `b1fit`
#' @param ... Generic graphical parameter to be passed on to \link{plotwitherror}
#'
#' @return
#' See \link{plot.cfit}
#' 
#' @export
plot.b1fit <- function(x, ...) {
  plot.cfit(x)
}

#' plot.c1fit
#'
#' @description
#' Generic function to plot an object of type `c1fit`
#'
#' @param x Object of type `c1fit`
#' @param ... Generic graphical parameter to be passed on to \link{plotwitherror}
#'
#' @return
#' No return value, only plots are generated.
#' 
#' @export
plot.cfit <- function(x, ...) {
  fit <- x
  fit.mass <- abs(fit$fitresult$par[fit$matrix.size+1])
  if(!is.null(fit$effmass$mll)) {
    new_window_if_appropriate()
    plot.effmass(m=fit.mass,
                 ll=data.frame(t=fit$effmass$t, mass=fit$effmass$mll, dmass=fit$effmass$dmll),
                 lf=data.frame(t=fit$effmass$t, mass=fit$effmass$mlf, dmass=fit$effmass$dmlf),
                 ff=data.frame(t=fit$effmass$t, mass=fit$effmass$mff, dmass=fit$effmass$dmff))
  }
  if(!is.null(fit$effmass$m)) {
    new_window_if_appropriate()
    plot.effmass(m=fit.mass, ll=data.frame(t=fit$effmass$t, mass=fit$effmass$m, dmass=fit$effmass$dm))
  }
  if(!is.null(fit$uwerrresultmpcac)) {
    new_window_if_appropriate()
    plot(fit$uwerrresultmpcac, main=expression(m[PCAC]))
  }
  if(!is.null(fit$uwerrresultmps)) {
    new_window_if_appropriate()
    plot(fit$uwerrresultmps, main=expression(m[PS]))
  }
  if(!is.null(fit$uwerrresultfps)) {
    new_window_if_appropriate()
    plot(fit$uwerrresultfps, main=expression(f[PS]))
  }
  if(!is.null(fit$uwerrresultmv)) {
    new_window_if_appropriate()
    plot(fit$uwerrresultmv, main=expression(m[V]))
  }
  if(!is.null(fit$boot)) {
    new_window_if_appropriate()
    plot(fit$boot, main="Bootstrap analysis for mps")
  }
  if(!is.null(fit$tsboot)) {
    new_window_if_appropriate()
    plot(fit$tsboot, main="TS Boostrap analysis for mps")
  }
  if(!is.null(fit$mv.boot)) {
    new_window_if_appropriate()
    plot(fit$mv.boot, main="Bootstrap analysis for mv")
  }
  if(!is.null(fit$mv.tsboot)) {
    new_window_if_appropriate()
    plot(fit$mv.tsboot, main="TS Bootstrap analysis for mv")
  }
  if(!is.null(fit$fitdata)) {
    new_window_if_appropriate()
    plot(fit$fitdata$Chi, main="Chi per data point", xlab="data point", ylab=expression(chi), type="h")
  }
  if(!is.null(fit$dpaopp)) {
    new_window_if_appropriate()
    plotwitherror(fit$dpaopp$t, fit$dpaopp$mass, fit$dpaopp$dmass,
                  main=expression(m[PCAC]), xlab="t", ylab=expression(m[PCAC]))
    abline(h=fit$fitresult$par[3]*fit$fitresult$par[2]/fit$fitresult$par[1]/2., col="red")
  }
  if(!is.null(fit$MChist.dpaopp)) {
    new_window_if_appropriate()
    plot(seq(1, length(fit$MChist.dpaopp))+fit$skip, fit$MChist.dpaopp, type="l",
         main=expression(m[PCAC]), xlab=expression(t[HMC]), ylab=expression(m[PCAC]))
    abline(h=fit$fitresult$par[3]*fit$fitresult$par[2]/fit$fitresult$par[1]/2., col="red")
  }
}
