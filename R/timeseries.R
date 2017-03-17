# function to plot timeseries data, a corresponding histogram
# and an error shading for an error analysis via uwerr

plot_timeseries <- function(dat,xdat,pdf.filename,
                            ylab,plotsize,titletext,hist.by,
                            name="",xlab="$t_\\mathrm{MD}$",hist.probs=c(0.0,1.0),errorband_color=rgb(0.6,0.0,0.0,0.6),
                            stepsize=1,
                            uwerr.S=2,periodogram=FALSE,debug=FALSE,uw.summary=TRUE,...) {
  if(missing(xdat)) { xdat <- seq(1,length(dat),stepsize) }

  yrange <- range(dat)
  
  uw.data <- uwerrprimary(dat,S=uwerr.S)
  if(debug) {
    print(paste("uw.",name,sep=""))
    print(summary(uw.data))
  }
  
  tikzfiles <- tikz.init(basename=pdf.filename,width=plotsize,height=plotsize)
  op <- par(family="Palatino",cex.main=0.8,font.main=1)
  par(mgp=c(2,1.0,0))

  # plot the timeseries
  plot(x=xdat,xlim=range(xdat),y=dat,ylab=ylab,t='l',xlab=xlab,main=titletext,...)

  rect(xleft=range(xdat)[1],
       xright=range(xdat)[2],
       ytop=uw.data$value+uw.data$dvalue,
       ybottom=uw.data$value-uw.data$dvalue,border=FALSE,col=errorband_color)
  abline(h=uw.data$value,col="black",lwd=2)                                                                                                   
  # plot the corresponding histogram
  hist.data <- NULL
  
  hist.breaks <- floor( ( max(dat)-min(dat) ) / uw.data$dvalue )
  if(!missing(hist.by)){
    hist.breaks <- floor( ( max(dat)-min(dat) ) / hist.by )
  } else {
    if(hist.breaks < 10 || hist.breaks > 150){
      hist.breaks <- 70
    }
  }
  
  hist.data <- hist(dat,xlim=quantile(dat,probs=hist.probs),main=titletext,xlab=ylab, breaks=hist.breaks)
  rect(ytop=max(hist.data$counts),
       ybottom=0,
       xright=uw.data$value+uw.data$dvalue,
       xleft=uw.data$value-uw.data$dvalue,border=FALSE,col=errorband_color)
  abline(v=uw.data$value,col="black",lwd=2) 

  # and a periodogram
  if(periodogram)
  {
    spec.pgram(x=dat,main=paste(ylab,paste("raw periodogram",titletext)))
  }

  # and the uwerr plots
  if(uw.summary){
    plot(uw.data,main=paste(ylab,paste("UWErr analysis",titletext)),x11=FALSE,plot.hist=FALSE)
  }
   
  tikz.finalize(tikzfiles)

  return(t(data.frame(val=uw.data$value, dval=uw.data$dvalue, tauint=uw.data$tauint, 
                      dtauint=uw.data$dtauint, Wopt=uw.data$Wopt, stringsAsFactors=FALSE)))
}

# function to plot timeseries of eigenvlues, including minimum and maximum eigenvalue bands
# as found in the monomial_0x.data files produced by tmLQCD

plot_eigenvalue_timeseries <- function(dat,xdat,pdf.filename,
                            ylab,plotsize,filelabel,titletext,
                            stepsize=1,errorband_color=rgb(0.6,0.0,0.0,0.6),
                            debug=FALSE) {
  if(missing(xdat)) { xdat <- seq(1:nrow(dat),stepsize) }
  yrange <- range(dat[,2:5])
 
  uw.min_ev <- uwerrprimary(dat[,2])
  uw.max_ev <- uwerrprimary(dat[,3])

  if(debug){
    print("uw.eval.min_ev")
    print(summary(uw.min_ev))
    print("uw.eval.max_ev")
    print(summary(uw.max_ev))
  }

  tikzfiles <- tikz.init(basename=pdf.filename,width=plotsize,height=plotsize)
  par(mgp=c(2,1,0))

  # plot the timeseries
  plot(x=xdat,xlim=range(xdat),y=dat[,2],ylim=yrange,t='l',ylab=ylab,xlab=expression(t[MD]),main=titletext,log='y',tcl=0.02)
  lines(x=xdat,y=dat[,3])
 
  ## add the approximation interval
  lines(x=xdat,y=dat[,4],lty=2,col="darkgreen")
  lines(x=xdat,y=dat[,5],lty=2,col="darkgreen")
  
  # plot the corresponding histograms with error bands and mean values
  hist.min_ev <- hist(dat[,2],main=paste("min. eval",titletext),xlab="min. eval",tcl=0.02)
  rect(ytop=max(hist.min_ev$counts),
       ybottom=0,
       xright=uw.min_ev$value+uw.min_ev$dvalue,
       xleft=uw.min_ev$value-uw.min_ev$dvalue,border=FALSE,col=errorband_color)
  abline(v=uw.min_ev$value,col="black")                                                                                                   

  hist.max_ev <- hist(dat[,3],main=paste("max. eval",titletext),xlab="max. eval")
  rect(ytop=max(hist.max_ev$counts),
       ybottom=0,
       xright=uw.max_ev$value+uw.max_ev$dvalue,
       xleft=uw.max_ev$value-uw.max_ev$dvalue,border=FALSE,col=errorband_color)
  abline(v=uw.max_ev$value,col="black")                                                                                                   

  tikz.finalize(tikzfiles)

  return(list(mineval=t(data.frame(val=uw.min_ev$value, dval=uw.min_ev$dvalue, tauint=uw.min_ev$tauint, 
                                   dtauint=uw.min_ev$dtauint, Wopt=uw.min_ev$Wopt, stringsAsFactors=FALSE)),
              maxeval=t(data.frame(val=uw.max_ev$value, dval=uw.max_ev$dvalue, tauint=uw.max_ev$tauint, 
                                   dtauint=uw.max_ev$dtauint, Wopt=uw.max_ev$Wopt, stringsAsFactors=FALSE)) ) )
              
}

