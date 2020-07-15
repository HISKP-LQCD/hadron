#' compute.plotlims
#'
#' @description
#' Computes limits for plots
#'
#' @param val Numeric. Value.
#' @param logscale Boolean.
#' @param cumul.dval Numeric. Cumulative error.
#' @param cumul.mdval Numeric. Cumulative error.
#'
#' @return
#' The computed plot limits are returned as a two component
#' numeric vector.
compute.plotlims <- function(val, logscale, cumul.dval, cumul.mdval){
  tmp <- val - 0.1*abs(val)
  tmpp <- val + 0.1*abs(val)
  if(!is.null(cumul.dval)) {
    ## cumul.mdx is implicitly negative
    tmp <- val+2*apply(X=cumul.mdval,MARGIN=1,FUN=min,na.rm=TRUE)
    tmpp <- val+2*apply(X=cumul.dval,MARGIN=1,FUN=max,na.rm=TRUE)
  }
  if(logscale) {
    tmp <- tmp[ tmp > 0 ]
    tmpp <- tmpp[ tmpp > 0 ]
    if( ( all(is.na(tmp)) && all(is.na(tmpp)) ) | ( length(tmp) == 0 & length(tmpp) == 0 ) ){
      warning("compute.plotlims: log scale requested but there are no positive data, setting default range\n")
      tmp <- 10^(-6)
      tmpp <- 10^2
    }
  }
  range(c(as.vector(tmp),as.vector(tmpp)),na.rm=TRUE)
}

## there are two possibilities for one-dimensional vectors: the vector class or the array class
is.vectorial <- function(x) {
  ( is.vector(x) || length(dim(x))==1 )
}

errorpos <- function(dx,errsum.method="linear") {
  if( !any(errsum.method==c("linear","linear.quadrature","quadrature")) ){
    stop(sprintf("errorpos: called with unknown errsum.method %s",errsum.method))
  }
  norm <- 1

  if(errsum.method=="quadrature"){
    norm <- apply(X=dx,MARGIN=1,FUN=function(x){sqrt(sum(x^2))})
  }

  rval <- NULL
  if ( errsum.method=="linear.quadrature" ) {
    rval <- array(dim=c(nrow(dx),(ncol(dx)+2)))
  } else {
    rval <- array(dim=c(nrow(dx),(ncol(dx)+1)))
  }
  
  rval[,1] <- 0
  if (errsum.method=="linear.quadrature" || errsum.method=="linear") {
    rval[,2] <- dx[,1]
  } else if (errsum.method=="quadrature"){
    rval[,2] <- dx[,1]^2
  }
  
  ## compute error bar deviations from the columns of dx depending on the errsum.method
  ## so for the "quadrature" method we will have
  # norm[j] <- sqrt(sum(dx[j,]^2))
  # rval[j,] == c(0, dx[j,1]^2/norm[j], (dx[j,1]^2+dx[j,2]^2)/norm[j], (dx[j,1]^2+dx[j,2]^2+dx[j,3]^2)/norm[j],..,norm[j])
  if(ncol(dx)>1){ 
    for( i in 2:ncol(dx) ){
      ## when dx has only one row, dx[,1:i] is returned as a vector, making apply fail....
      if(nrow(dx)>1){
        rval[,(i+1)] <- apply( X=dx[,1:i],MARGIN=1,
                             FUN=function(x){
                                   if(errsum.method=="quadrature") return( sum(x^2) )
                                   else return( sum(x) )
                                 } )
      } else {
        if(errsum.method=="quadrature"){
          rval[,(i+1)] <- sum(dx[,1:i]^2)
        } else {
          rval[,(i+1)] <- sum(dx[,1:i])
        }
      }
    }
  }
  if(ncol(dx)>1 && errsum.method=="linear.quadrature"){
    ## unlike above, even if dx has only one row, this works fine
    rval[,(ncol(dx)+2)] <- apply(X=dx,MARGIN=1,FUN=function(x){ sqrt(sum(x^2)) })
  }
  rval <- rval/norm
  rval
}

#' Plot Command For XY Plots With Error Bars
#' 
#' Plot command for XY scatterplots based on plot and points which provides
#' support for multiple, non-symmetric error bars. Error bars are drawn as
#' vertical or horizontal lines originating from the point with narrow,
#' perpendicular lines at the end of the error bar (end caps). When multiple
#' errors are drawn, the width of the perpendicular line increases from the
#' innermost error bar to the outermost one. Different summation methods for
#' the individual errors are supported.
#' 
#' 
#' @param x vector of x coordinates
#' @param xlim limits for x-axis
#' @param y vector of y coordinates
#' @param ylim limits for y-axis
#' @param col colour of plotted data
#' @param dy one of: \itemize{ \item Vector of errors on y coordinates.
#' \item Array, matrix or data frame if multiple error bars are to be drawn,
#' such that each column refers to one error. The individual errors should be
#' provided as is, because they are summed internally to draw the final error
#' bars.  A given column can also be provided with 0 entries, in which case the
#' error bar will be drawn, but it will have zero length, such that only the
#' end caps for this error will be visible.  }
#' @param dx Same as \code{dy}, but for the x coordinate.
#' @param mdx Support for non-symmetric error bars. Same as \code{dx}, but for
#' errors in the negative x-direction. Errors should be provided as positive
#' numbers, the correct sign will be added internally. If not provided,
#' \code{dx} is used as a symmetric error.
#' @param mdy Same as \code{mdx} but for the y coordinate.
#' @param errsum.method Determines how the invidual errors should be summed for
#' display purposes. Valid argument values are: \itemize{ \item "linear"
#' \itemize{ \item Individual errors are summed linearly, such that the distance
#' from the point to the \eqn{i}'th error bar, \eqn{l_i}, is \deqn{ l_i =
#' \sum_{j=1}^i e_j } Hence, the third error bar, for example, would be located
#' at \deqn{ l_3 = e_1 + e_2 + e_3 } while the second error bar is at \deqn{
#' l_2 = e_1 + e_2 }
#' 
#' }
#' 
#' \item "quadrature" \itemize{ \item Individual errors are summed in quadrature
#' and error bars are drawn at the fractional position according to the
#' following formula: \deqn{ l_{max} = \sqrt{ \sum_{j=1}^{max} e_j^2 } } \deqn{
#' l_i = \sum_{j=1}^i e_j^2 / l_{max} }
#' 
#' }
#' 
#' \item "linear.quadrature" \itemize{ \item Errors are summed as for "linear",
#' but the total error summed in quadrature is also indicated as an end cap of
#' triple line width }
#' 
#' }
#' @param rep If set to \code{TRUE}, operate like "replot" in gnuplot. Allows
#' adding points with error bars to the current plot. Switches the underlying
#' plotting routine from \code{\link[base]{plot}} to
#' \code{\link[graphics]{points}}.
#' @param ...  any graphic options passed over to \code{\link[base]{plot}}
#' @return a plot with error bars is drawn on the current device
#' @author Carsten Urbach, \email{urbach@@hiskp.uni-bonn.de} \cr Bartosz
#' Kostrzewa, \email{bartosz.kostrzewa@@desy.de}
#' @seealso \code{\link[base]{plot}}, \code{\link[graphics]{points}}
#'
#' @return
#' Returns for convenience a list with elements `xlim` and `ylim` representing the
#' x- and y-limits chosen by the routine.
#' 
#' @examples
#' 
#' # Create some random data, set one error to zero.
#' x <- 1:50
#' y <- runif(50, 0, 1)
#' dy <- runif(50, 0.1, 0.2)
#' dy[4] <- 0
#' 
#' plotwitherror(x, y, dy)
#' 
#' @export plotwitherror
plotwitherror <- function(x, y, dy, ylim = NULL, dx, xlim = NULL, mdx, mdy, errsum.method="linear.quadrature", rep=FALSE, col="black", ...) {
  if(!missing(mdy) && missing(dy)){
    stop("plotwitherror: if 'mdy' is provided, 'dy' must be too (it can be 0)")
  }
  if(!missing(mdx) && missing(dx)){
    stop("plotwitherror: if 'mdx' is provided, 'dx' must be too (it can be 0)")
  }
  if( (!missing(dx) && is.null(ncol(dx)) && length(dx)>length(x) ) || ( !missing(mdx) && is.null(ncol(mdx)) && length(mdx)>length(x) ) ||
      (!missing(dy) && is.null(ncol(dy)) && length(dy)>length(y) ) || ( !missing(mdy) && is.null(ncol(dy)) && length(dy)>length(y) ) ) {
    stop("plotwitherror: if any dx, mdx, dy or mdy is a simple vector, it must be of the same length as the corresponing x/y, multiple errors must ALWAYS be specified in the form of arrays/data frames/matrices")
  }

  ## if no plotting limits are passed, we calculate them automatically
  ## for this to work properly, we need to know if any of the axes
  ## are requested in log scale
  dots <- list(...)
  ylog <- FALSE
  xlog <- FALSE
  if( 'log' %in% names(dots) ) {
    logidx <- which(names(dots) == 'log')
    if(any( grepl('y', dots[logidx]) ) ){ ylog <- TRUE }
    if(any( grepl('x', dots[logidx]) ) ){ xlog <- TRUE }
  }
  my.xlim <- c()
  my.ylim <- c()
 
  ## cumulative errors as computed with errsum.method 
  cumul.dx <- NULL
  cumul.mdx <- NULL
  cumul.dy <- NULL
  cumul.mdy <- NULL
  
  ## compute cumulative error bar positions, first convert vectorial dx, dy into arrays 
  ## see above for description of what errorpos does
  if(!missing(dx)){
    ## if dx is a simple vector, convert into an appropriately sized array
    if(is.vectorial(dx)){
      dx <- array(data=dx,dim=c(length(dx),1))
    }

    cumul.dx <- errorpos(dx,errsum.method)
    
    if(missing(mdx)){
      cumul.mdx <- -errorpos(dx,errsum.method)
    } else {
      if(is.vectorial(mdx)){
          mdx <- array(data=mdx,dim=c(length(mdx),1))
        }
      cumul.mdx <- -errorpos(mdx,errsum.method)
    }  
  }
  
  if(!missing(dy)){
    if(is.vectorial(dy)){
      dy <- array(data=dy,dim=c(length(dy),1))
    }
    cumul.dy <- errorpos(dy,errsum.method)
    
    if(missing(mdy)){
      cumul.mdy <- -errorpos(dy,errsum.method)
    } else {
      if(is.vectorial(mdy)){
        mdy <- array(data=mdy,dim=c(length(mdy),1))
      }
      cumul.mdy <- -errorpos(mdy,errsum.method)
    }  
  }

  if( is.null(xlim) ){
    my.xlim <- compute.plotlims(val=x, logscale=xlog, cumul.dval=cumul.dx, cumul.mdval=cumul.mdx)
  } else {
    my.xlim <- xlim
  }

  if( is.null(ylim) ){
    my.ylim <- compute.plotlims(val=y, logscale=ylog, cumul.dval=cumul.dy, cumul.mdval=cumul.mdy)
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

  ## We will need some langth in inches, so let's define the proportionality factor as
  ## size of the output in inches devided by the differnce of boundary coordinates.
  x.to.inches <- par("pin")[1]/diff(par("usr")[1:2])
  y.to.inches <- par("pin")[2]/diff(par("usr")[3:4])

  if(!is.null(cumul.dy)) {
    for(cumul.err in list(cumul.dy,cumul.mdy)){
      rng <- 2:ncol(cumul.err)
      if(ncol(cumul.err)>2 && errsum.method=="linear.quadrature") rng <- 2:(ncol(cumul.err)-1)
      ## this loop is necessary because the "length" parameter of "arrows" is not vectorial...
      ## so to accomodate the generalisations below, we need to draw the error for each point
      ## individually, it doesn't make it much slower
      for(rw in 1:length(y)){
        ## the length of the arrowhead lines will depend on the "level" (the more errors, the longer the arrowhead lines)
        arwhd.len <- 0.02
        clr <- col
        if(length(col)>1) clr <- col[rw]
        for(level in rng){
          start <- y[rw]+cumul.err[rw,(level-1)]
          end <- y[rw]+cumul.err[rw,level]
          if (!is.na(start) && !is.na(end)){
            if(par("ylog")){
              ## logarithmic scaling means distances are multiplicative
              if(start > 0 && end > 0){
                dy.in.inches <- abs(log10(end/start))*y.to.inches
                ## An arrow can only be plotted if it is longer than 1/1000 inches.
                if(dy.in.inches <= 1/1000) {
                  ## We want to plot an errorbar of minimum size anyway, so as to show the point is plotted with an error.
                  start <- y[rw] / 10^(1.01/2000/y.to.inches)
                  end <- y[rw] * 10^(1.01/2000/y.to.inches)
                }
              }else{
                dy.in.inches <- 0
              }
            }else{
              dy.in.inches <- abs(end-start)*y.to.inches
              if(dy.in.inches <= 1/1000) {
                start <- y[rw] - 1.01/2000/y.to.inches
                end <- y[rw] + 1.01/2000/y.to.inches
              }
            }

            if(dy.in.inches > 1/1000) {
              arrows(x[rw], start, x[rw], end, length=arwhd.len, angle=90, code=2, col=clr)
            }else if(dy.in.inches > 0){
              arrows(x[rw], start, x[rw], end, length=arwhd.len, angle=90, code=3, col=clr)
            }
            arwhd.len <- arwhd.len + 0.01
          }
        } 
        ## for the linear.quadrature method, show the total error as a line of triple thickness
        ## without drawing any "arrowstems"
        if(ncol(cumul.err)>2 && errsum.method=="linear.quadrature"){
          ## to be consistent, drawX/Ybars uses inches just like arrows
          arwhd.len <- arwhd.len + 0.02
          drawYbars(x=x[rw],y=y[rw],dy=cumul.err[rw,ncol(cumul.err)],length=arwhd.len,lwd=3,col=clr)
        }
      }
    }
  }
  if(!is.null(cumul.dx)) {
    for(cumul.err in list(cumul.dx,cumul.mdx)){
      rng <- 2:ncol(cumul.err)
      if(ncol(cumul.err)>2 && errsum.method=="linear.quadrature") rng <- 2:(ncol(cumul.err)-1)
      for(rw in 1:length(x)){
        arwhd.len <- 0.02
        clr <- col
        if(length(col)>1) clr <- col[rw]
        for(level in rng){
          start <- x[rw]+cumul.err[rw,(level-1)]
          end <- x[rw]+cumul.err[rw,level]
          if (!is.na(start) && !is.na(end)){
            if(par("xlog")){
              if(start > 0 && end > 0){
                dx.in.inches <- abs(log10(end/start))*x.to.inches
                if(dx.in.inches <= 1/1000) {
                  start <- x[rw] / 10^(1.01/2000/x.to.inches)
                  end <- x[rw] * 10^(1.01/2000/x.to.inches)
                }
              }else{
                dx.in.inches <- 0
              }
            }else{
              dx.in.inches <- abs(end-start)*x.to.inches
              if(dx.in.inches <= 1/1000) {
                start <- x[rw] - 1.01/2000/x.to.inches
                end <- x[rw] + 1.01/2000/x.to.inches
              }
            }

            if(dx.in.inches > 1/1000) {
              arrows(start, y[rw], end, y[rw], length=arwhd.len, angle=90, code=2, col=clr)
            }else if(dx.in.inches > 0){
              arrows(start, y[rw], end, y[rw], length=arwhd.len, angle=90, code=3, col=clr)
            }
            arwhd.len <- arwhd.len + 0.01
          }
        }
        if(ncol(cumul.err)>2 && errsum.method=="linear.quadrature"){
          arwhd.len <- arwhd.len + 0.02
          drawXbars(x=x[rw],y=y[rw],dx=cumul.err[rw,ncol(cumul.err)],length=arwhd.len,lwd=3,col=clr)
        }
      }
    } 
  }
  
  return(invisible(list(xlim=my.xlim, ylim=my.ylim)))
}

#' plothlinewitherror
#'
#' @description
#' plot a horizontal line with error band
#'
#' @param m Numeric. Mean value of the line to plot.
#' @param dp Numeric. Error up.
#' @param dm Numeric. Error down.
#' @param col String. Colour.
#' @param x0 Numeric. Left value of the range of the horizontal line.
#' @param x1 Numeric. Right value of the range of the horizontal line.
#'
#' @return
#' No return value, only graphics is generated.
#' 
#' @export
plothlinewitherror <- function(m, dp, dm, col=c("red"), x0, x1) {
  if(missing(dm)) {
    dm <- dp
  }
  arrows(x0=x0, y0=m, x1=x1, y1=m, col=col, length=0)
  arrows(x0=x0, y0=m+dp, x1=x1, y1=m+dp, col=col, length=0, lwd=c(1))
  arrows(x0=x0, y0=m-dm, x1=x1, y1=m-dm, col=col, length=0, lwd=c(1))
}

#' plot.massfit
#'
#' @description
#' Generic function to plot an object of type `massfit`
#'
#' @param x Object of type `massfit`
#' @param ... Generic graphical parameter to be passed on to \link{plotwitherror}
#' @param xlab String. Label for x-axis
#' @param ylab String. Lable for y-axis
#'
#' @return
#' See \link{plotwitherror}.
#'
#' @export
plot.massfit <- function(x, ..., xlab = "t", ylab = "m") {
  plotwitherror(x$t, x$mass, x$dmass, xlab=xlab, ylab=ylab, ...)
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


#' plot.ofit
#'
#' @description
#' Generic function to plot an object of type `ofit`
#'
#' @param x Object of type `ofit`
#' @param ... Generic graphical parameter to be passed on to \link{plotwitherror}
#'
#' @return
#' See \link{plot.cfit}
#' 
#' @export
plot.ofit <- function(x, ...) {
  plot.cfit(x)
}


#' plot.effmass
#'
#' @param x Object of class `effmass`
#' @param ... Graphical parameters to be passed on.
#' @param ll local-local effective mass object
#' @param lf local-fuzzed effective mass object
#' @param ff fuzzed-fuzzed effective mass object
#'
#' @return
#' No value returned, only plots are generated.
#' 
#' @export
plot.effmass <- function (x, ..., ll, lf, ff) {
  m <- x

  new_window_if_appropriate()
  if(!missing(ff)) {
    plot.massfit(ff, ylab=expression(m[eff]), xlab="t", ...)
    points((ll$t-0.2), ll$mass, pch=1, col="blue")
    arrows((ll$t-0.2), ll$mass-ll$dmass,
           (ll$t-0.2), ll$mass+ll$dmass, length=0.01,angle=90,code=3)
    points((lf$t+0.2), lf$mass, pch=2, col="red")
    arrows((lf$t+0.2), lf$mass-lf$dmass,
           (lf$t+0.2), lf$mass+lf$dmass, length=0.01,angle=90,code=3)
    lines(ll$t, rep(m, times=length(ll$t)))
  }
  else if(!missing(lf)) {
    plot.massfit(ll, ylab=expression(m[eff]), xlab="t", ...)
    points((lf$t-0.2), lf$mass, pch=1, col="blue")
    arrows((lf$t-0.2), lf$mass-lf$dmass,
           (lf$t-0.2), lf$mass+lf$dmass, length=0.01,angle=90,code=3)
    lines(ll$t, rep(m, times=length(ll$t)))
  }
  else {
    plot.massfit(ll, ylab=expression(m[eff]), xlab="t", ...)
    lines(ll$t, rep(m, times=length(ll$t)))
  }
}


#' Plots averx data
#'
#' @param x `averx` object
#' @param ... ignored
#'
#' @return
#' Returns the plotted data in from of a \link{data.frame} with named
#' columns `t` (the time index), `averx` the values of average x and
#' `daverx` the statistical error estimate.
#' 
#' @export
plot.averx <- function(x, ...) {
  averx <- x
  Thalfp1 <- averx$Cf2pt$Time/2+1
  #plot(averx$effmass, ylim=c(averx$effmass$opt.res$par[1]/2, 3/2*averx$effmass$opt.res$par[1]), main=c("Pion Effectivemass"), xlab=c("t/a"), ylab=c("a Meff"))
  
  new_window_if_appropriate()
  plot(averx$Cf3pt, xlab=c("t/a"), ylab=c("C3pt"), ylim=c(0, 2*averx$plateau), main=c("3pt ov 2pt"))
  plothlinewitherror(m=averx$plateau, dp=sd(averx$plateau.tsboot[,1]), dm=sd(averx$plateau.tsboot[,1]),
                    x0=averx$t1, x1=averx$t2)
  new_window_if_appropriate()
  plotwitherror(c(0:(averx$Cf2pt$Time/2)),
                averx$Cf3pt$cf0/averx$matrixfit$opt.res$par[1]/averx$Cf2pt$cf0[Thalfp1],
                apply(averx$Cf3pt$cf.tsboot$t/averx$matrixfit$opt.tsboot[1,]/averx$Cf2pt$cf.tsboot$t[,Thalfp1], 2, sd),
                ylim=c(0,0.6), xlab=c("t/a"), ylab=c("C3pt"), main=c("<x>")
                )
  new_window_if_appropriate()
  qqplot(averx$plateau.tsboot[,2], rchisq(n=averx$boot.R, df=averx$dof), main=paste("qqplot chisq"))
  return(invisible(data.frame(t=c(0:(averx$Cf2pt$Time/2)),
                              averx=averx$Cf3pt$cf0/averx$matrixfit$opt.res$par[1]/averx$Cf2pt$cf0[Thalfp1],
                              daverx=apply(averx$Cf3pt$cf.tsboot$t/averx$matrixfit$opt.tsboot[1,]/averx$Cf2pt$cf.tsboot$t[,Thalfp1], 2, sd)
                              )
                   )
         )
}

#' plot.pionff
#'
#' @description
#' Generic function to plot an object of type `pionff`
#'
#' @param x Object of type `pionff`
#' @param ... Generic graphical parameter to be passed on to \link{plotwitherror}
#'
#' @return
#' No return value, only plots are generated.
#' 
#' @export
plot.pionff <- function (x, ...) {
  ff <- x
  Time <- ff$Cf2ptp0$Time
  Thalfp1 <- Time/2+1
  plot(mul.cf(ff$Cf3ptp0, 1./ff$Cf2ptp0$cf0[Thalfp1]), main=c("1./Z_V"), xlab=c("t/a"), ylab=c("1/Z_V"))
  plothlinewitherror(m=1./ff$Cf2ptp0$cf0[Thalfp1]*ff$plateaufitZV$plateau, dp=sd(ff$plateaufitZV$plateau.tsboot[,1]/ff$Cf2ptp0$cf.tsboot$t[,Thalfp1]), dm=sd(ff$plateaufitZV$plateau.tsboot[,1]/ff$Cf2ptp0$cf.tsboot$t[,Thalfp1]),
                    x0=ff$plateaufitZV$t1, x1=ff$plateaufitZV$t2)

  new_window_if_appropriate()
  plot(mul.cf(ff$Cf3ptratio, ff$Cf2ptratio$cf0[Thalfp1]), ylim=c(0,1.3), main=c("F(q^2)"), xlab=c("t/a"), ylab=c("F(q^2)"))
  plothlinewitherror(m=ff$plateaufitFF$plateau*ff$Cf2ptratio$cf0[Thalfp1], dp=sd(ff$plateaufitFF$plateau.tsboot[,1]*ff$Cf2ptratio$cf.tsboot$t[,Thalfp1]), dm=sd(ff$plateaufitFF$plateau.tsboot[,1]*ff$Cf2ptratio$cf.tsboot$t[,Thalfp1]),
                    x0=ff$plateaufitFF$t1, x1=ff$plateaufitFF$t2)
}




#' Plot Command For Class Ouputdata
#' 
#' Generic plot routine for class \dQuote{ouputdata}. Currently it plots the
#' plaquette history and the history of \eqn{\Delta H}{Delta H}
#' 
#' 
#' @param x object of class \dQuote{outputdata} obtained from a read with
#' \code{readoutputdata}
#' @param skip number of trajectories to be skipped in analysis for plaquette
#' and \eqn{\exp(-\Delta H)}{exp(-Delta H)}.
#' @param ...  additional arguments passed to the generic plot function.
#' @return list containing the \dQuote{data}, an object of class \dQuote{uwerr}
#' called \dQuote{plaq.res} containing the statisical analysis for the
#' plaquette and a second object of type \dQuote{uwerr} called \dQuote{dH.res}
#' with the statisical analysis for \eqn{\exp(-\Delta }{exp(-Delta H)}\eqn{
#' H)}{exp(-Delta H)}.
#' @author Carsten Urbach, \email{curbach@@gmx.de}
#' @seealso \code{\link{readoutputdata}}, \code{\link{uwerr}}
#' @keywords methods hplot
#' @return
#' The plotted data is return in form of a \link{list} with named elements `data`
#' containing the input data, plaq.res an object returned by \link{uwerrprimary}
#' for the plaquette data dn `dH.res` an object returned by \link{uwerrprimary}
#' for \eqn{\Delta H}{Delta H}.
#' 
#' @export 
#' @examples
#' 
#' plaq <- readoutputdata(paste0(system.file(package="hadron"), "/extdata/output.data"))
#' plaq.plot <- plot(plaq, skip=100)
#' summary(plaq.plot$plaq.res)
#' 
plot.outputdata <- function (x, skip = 0, ...) {
  data <- x
  plaq.res <- uwerrprimary( data$V2[skip:length(data$V2)])
  dH.res <- uwerrprimary( exp(-data$V3[skip:length(data$V3)]))
  plot(data$V1, data$V2, type="l",
       main="Plaquette", xlab=expression(t[HMC]), ylab="P", col="red", ...)
  abline(h=plaq.res$value, col="blue")
  abline(h=plaq.res$value+plaq.res$dvalue)
  abline(h=plaq.res$value-plaq.res$dvalue)
  new_window_if_appropriate()
  plot(data$V1, data$V3, type="l",
       main=expression(paste(Delta, "H")), xlab=expression(t[HMC]), ylab=expression(paste(Delta, "H")), ...)
  return(invisible(list(data=data, plaq.res=plaq.res, dH.res = dH.res)))
}

## draw Y (X) error bars at coordinates y+dy (x+dx) (where dy (dx) can be negative)
## like for arrows, the bar width is specified in inches
## and additonal parameters (like lwd and col) can be passed
## to segments
drawYbars <- function(x,y,dy,length=0.01,...) {
  x.inch <- grconvertX(x=x,from="user",to="inches")
  xp.inch <- x.inch+length
  xm.inch <- x.inch-length
  xp <- grconvertX(x=xp.inch,from="inches",to="user")
  xm <- grconvertX(x=xm.inch,from="inches",to="user")
  segments(xp, y+dy,
           xm, y+dy, ...)
}

drawXbars <- function(x,y,dx,length=0.01,...) {
  y.inch <- grconvertY(y=y,from="user",to="inches")
  yp.inch <- y.inch+length
  ym.inch <- y.inch-length
  yp <- grconvertY(y=yp.inch,from="inches",to="user")
  ym <- grconvertY(y=ym.inch,from="inches",to="user")
  segments(x+dx, ym,
           x+dx, yp, ...)
}

new_window_if_appropriate <- function () {
  is_X11 <- grepl(pattern = "X11", x = names(dev.cur()), ignore.case = TRUE)
  is_null <- grepl(pattern = "null", x = names(dev.cur()), ignore.case = TRUE)

  if (interactive() && (is_X11 || is_null)) {
    message('Opening a new X11 window.\n')
    dev.new()
  }
}

#' pointswithslantederror
#'
#' plots data points with slanted error bars
#'
#' @description
#' This function plots points with x- and y-errors visualised as a
#' slanted errorbar. The length of the error bar represents x- and y-errors
#' added in quadrature. The slope of the error bar is positive of negative
#' depending on whether the correlation betwenn x and y is positive or
#' negative, respectively.
#'
#' @param x numeric vector. x-values
#' @param y numeric vector. y-values
#' @param dx numeric vector. x-standard errors
#' @param dy numeric vector. y-standard errors
#' @param cor numeric vector. Correlation coefficients between x- and y-
#' errors.
#' @param col the color of the points
#' @param bcol the color of the slanted error bars
#' @param ... further graphical parameters to be passed on to `points`
#'
#' @examples
#' x <- c(1:5)
#' y <- x^2
#' dx <- c(0.1, 0.2, 0.2, 0.1, 0.05)
#' dy <- c(0.05, 0.2, 0.1, 0.2, 0.1)
#' cor <- c(1, -1, -1, 1, 1)
#' plot(NA, xlim=range(x), ylim=range(y), xlab="y", ylab="y")
#' pointswithslantederror(x=x, y=y, dx=dx, dy=dy, cor=cor)
#' 
#' @export
pointswithslantederror <- function(x, y, dx, dy, cor, col="black", bcol="black", ...) {
  points(x=x, y=y, col=col, ...)
  l <- sqrt(dx^2+dy^2)
  fac <- rep(1, times=length(cor))
  fac[which(cor < 0)] <- -1
  arrows(code=0,
         x0=x-fac*dx, y0=y-dy,
         x1=x+fac*dx, y1=y+dy,
         col=bcol)
}
