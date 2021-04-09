#' plot_timeseries
#' 
#' @description
#' function to plot timeseries data, a corresponding histogram
#' and an error shading for an error analysis via uwerr
#'
#' @param dat Timeseries to analyse.
#' @param ylab Y-axis label.
#' @param xlab X-axis label.
#' @param plotsize Width and Height of plot.
#' @param titletext Text in the plot title.
#' @param hist.by Numeric. Stepping to compute the histogram breaks.
#' @param stat_range Optional integer vector of length 2. Start and end indices
#'        of the subset of `dat` to be plotted. If left empty, all of `dat` will be
#'        plotted.
#' @param pdf.filename String. PDF filename.
#' @param name String. Timeseries name.
#' @param hist.probs Optional numeric vector of length 2. Probability extrema to limit the width
#'        of the histogram or smoothed density plots. By default all data is used. Note: this
#'        has not effect on the analysis as a whole or other plots.
#' @param smooth_density Boolean. Instead of plotting a histogram, use a smoothed density.
#' @param errorband_color String. Colour of the error band.
#' @param type String. Plot type, see \link{plot} for details.
#' @param uwerr.S Numeric. `S` of the \link{uwerr} method to be used.
#' @param periodogram Boolean. Whether to show a periodogram.
#' @param debug Boolean. Generate debug output.
#' @param uw.summary Boolean. Generate an \link{uwerr} summary plot.
#' @param ... Generic graphical parameters to be passed on.
#'
#' @return
#' Returns a \link{data.frame} with named columns `val`, `dval`, `tauint`, `dtauint`, `Wopt`
#' and `stringsAsFactors`, see \link{uwerr}.
#' 
#' @export
plot_timeseries <- function(dat, 
                            ylab, plotsize, titletext, hist.by,
                            stat_range = c(1.0, length(dat$y)),
                            pdf.filename,
                            name="", xlab="$t_\\mathrm{MD}$", 
                            hist.probs=c(0.0,1.0), errorband_color=rgb(0.6,0.0,0.0,0.6),
                            type='l',
                            uwerr.S=2,
                            smooth_density=FALSE,
                            periodogram=FALSE,debug=FALSE,uw.summary=TRUE,...) {

  stopifnot(length(stat_range) == 2)
  stopifnot(length(hist.probs) == 2)

  yrange <- range(dat$y)

  stat_y <- dat$y[ seq(stat_range[1], stat_range[2]) ]
  
  uw.data <- uwerrprimary(stat_y, S=uwerr.S)
  if(debug) {
    print(paste("uw.",name,sep=""))
    print(summary(uw.data))
  }
  
  tikzfiles <- NULL
  if( !missing(pdf.filename) ){
    tikzfiles <- tikz.init(basename=pdf.filename,width=plotsize,height=plotsize)
  }

  op <- par(family="Palatino",cex.main=0.8,font.main=1)
  on.exit(par(op))
  par(mgp=c(2,1.0,0))

  # plot the timeseries
  plot(x=dat$t,
       xlim=range(dat$t),
       y=dat$y,
       ylab=ylab,
       type=type,
       xlab=xlab,
       main=titletext, ...)

  rect(xleft=dat$t[ stat_range[1] ],
       xright=dat$t[ stat_range[2] ],
       ytop=uw.data$value+uw.data$dvalue,
       ybottom=uw.data$value-uw.data$dvalue,
       border=FALSE, col=errorband_color)
  abline(h=uw.data$value,col="black",lwd=2)

  legend(x="topright",
         legend=sprintf("%s $=%s$",
                         ylab,
                         tex.catwitherror(x=uw.data$value,
                                          dx=uw.data$dvalue,
                                          digits=3,
                                          with.dollar=FALSE,
										  with.cdot=FALSE)
                         ),
         lty=1,
         pch=NA,
         col="red",
         bty='n')
  
  # plot the corresponding histogram
  hist.data <- NULL
  
  hist.breaks <- floor( ( max(stat_y)-min(stat_y) ) / uw.data$dvalue )
  if(!missing(hist.by)){
    hist.breaks <- floor( ( max(stat_y)-min(stat_y) ) / hist.by )
  } else {
    if(hist.breaks < 10 || hist.breaks > 150){
      hist.breaks <- 70
    }
  }
 
  if( smooth_density ){
    # determine an appropriate bandwidth for the density estimate
    # by finding the closest power of 2 smaller than half the
    # number of measurements
    n <- ifelse(length(stat_y)/2 > 512, 
                2^floor(log2(length(stat_y)/2)),
                512)
    d <- density(stat_y, bw = "SJ", n = n)
    plot(d, 
         xlim=quantile(stat_y,probs=hist.probs),
         main=titletext,
         xlab=ylab, 
         breaks=hist.breaks)
    ytop <- max(d$y)
    ybottom <- 0.0
  } else { 
    hist.data <- hist(stat_y,
                      xlim=quantile(stat_y,probs=hist.probs),
                      main=titletext,
                      xlab=ylab, 
                      breaks=hist.breaks)
    ytop <- max(hist.data$counts)
    ybottom <- 0.0
  }
  rect(ytop=ytop,
       ybottom=0,
       xright=uw.data$value+uw.data$dvalue,
       xleft=uw.data$value-uw.data$dvalue,
       border=FALSE, col=errorband_color)
  abline(v=uw.data$value, col="black", lwd=2) 

  # and a periodogram
  if(periodogram)
  {
    spec.pgram(x=stat_y,
               main=paste(ylab,paste("raw periodogram",titletext)))
  }

  # and the uwerr plots
  if(uw.summary){
    plot(uw.data, 
         main=paste(ylab,paste("UWErr analysis",titletext)),
         plot.hist=FALSE)
  }
   
  if(!missing(pdf.filename)){
    tikz.finalize(tikzfiles)
  }

  return(t(data.frame(val=uw.data$value, dval=uw.data$dvalue, tauint=uw.data$tauint, 
                      dtauint=uw.data$dtauint, Wopt=uw.data$Wopt, stringsAsFactors=FALSE)))
}

#' plot_eigenvalue_timeseries
#' 
#' @description
#' function to plot timeseries of eigenvlues, including minimum and maximum eigenvalue bands 
#'  as found in the monomial_0x.data files produced by tmLQCD
#' 
#' @param dat Timeseries to analyse.
#' @param ylab Y-axis label.
#' @param plotsize Width and Height of plot.
#' @param filelabel String. Label of the file.
#' @param titletext Text in the plot title.
#' @param stat_range range of statistics to use.
#' @param pdf.filename String. PDF filename.
#' @param errorband_color String. Colour of the error band.
#' @param debug Boolean. Generate debug output.
#'
#' @return
#' Returns a list with two named elements `mineval` and `maxeval` for the minimal
#' and the maximal eigenvalue, see \link{plot_timeseries}.
#' 
#' @export
plot_eigenvalue_timeseries <- function(dat,
                                       stat_range,
                                       ylab, plotsize, filelabel,titletext,
                                       pdf.filename,
                                       errorband_color=rgb(0.6,0.0,0.0,0.6),
                                       debug=FALSE) {
  if( missing(stat_range) ) { stat_range <- c(1,nrow(dat)) }
  yrange <- range(dat[,2:5])

  stat_min_ev <- dat[ seq(stat_range[1], stat_range[2]), "min_ev" ]
  stat_max_ev <- dat[ seq(stat_range[1], stat_range[2]), "max_ev" ]

  uw.min_ev <- uwerrprimary( stat_min_ev )
  uw.max_ev <- uwerrprimary( stat_max_ev )

  if(debug){
    print("uw.eval.min_ev")
    print(summary(uw.min_ev))
    print("uw.eval.max_ev")
    print(summary(uw.max_ev))
  }

  tikzfiles <- NULL
  if(!missing(pdf.filename)){
    tikzfiles <- tikz.init(basename=pdf.filename,width=plotsize,height=plotsize)
  }
  par_save <- par(mgp=c(2,1,0))
  on.exit(par(par_save))

  # plot the timeseries
  plot(x=dat$traj, xlim=range(dat$traj), 
       y=dat$min_ev, ylim=yrange, 
       t='l', 
       ylab=ylab, 
       xlab=expression(t[MD]), 
       main=titletext, log='y', tcl=0.02)

  lines(x=dat$traj, y=dat$max_ev)
 
  ## add the approximation interval
  lines(x=dat$traj, y=dat$ev_range_min, lty=2, col="darkgreen")
  lines(x=dat$traj, y=dat$ev_range_max, lty=2, col="darkgreen")
  
  # plot the corresponding histograms with error bands and mean values
  hist.min_ev <- hist(stat_min_ev, main=paste("min. eval",titletext),xlab="min. eval",tcl=0.02)
  rect(ytop=max(hist.min_ev$counts),
       ybottom=0,
       xright=uw.min_ev$value+uw.min_ev$dvalue,
       xleft=uw.min_ev$value-uw.min_ev$dvalue,border=FALSE,col=errorband_color)
  abline(v=uw.min_ev$value,col="black")                                                                                                   

  hist.max_ev <- hist(stat_max_ev, main=paste("max. eval",titletext), xlab="max. eval")
  rect(ytop=max(hist.max_ev$counts),
       ybottom=0,
       xright=uw.max_ev$value+uw.max_ev$dvalue,
       xleft=uw.max_ev$value-uw.max_ev$dvalue,border=FALSE,col=errorband_color)
  abline(v=uw.max_ev$value,col="black")                                                                                                   

  if(!missing(pdf.filename)){
    tikz.finalize(tikzfiles)
  }

  return(list(mineval=t(data.frame(val=uw.min_ev$value, dval=uw.min_ev$dvalue, tauint=uw.min_ev$tauint, 
                                   dtauint=uw.min_ev$dtauint, Wopt=uw.min_ev$Wopt, stringsAsFactors=FALSE)),
              maxeval=t(data.frame(val=uw.max_ev$value, dval=uw.max_ev$dvalue, tauint=uw.max_ev$tauint, 
                                   dtauint=uw.max_ev$dtauint, Wopt=uw.max_ev$Wopt, stringsAsFactors=FALSE)) ) )
              
}

