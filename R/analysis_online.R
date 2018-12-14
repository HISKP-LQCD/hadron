# convenience function for analyzing online data from a tmLQCD run
### the various parameters are used to build a directory name into which R descends to read
### onlinemeas.xxxxxx and output.data
### but a directory can also be provided via 'rundir' with the understanding
### that L, T, kappa and mul must always be provided for the fits to work
### and the decay constant to be determined correctly
# 'addon' can be used to *append* arbitrary text to the directory name
# if it is construced on the fly. This is useful for analysing replicas
# which have identical parameters and thus need to be differentiated
# 'plaquette', 'dH', 'acc' and 'trajtime' control whether these are plotted
# cg_col indicates which column in output.data should be used to extract 
# a representative solver iteration counter

analysis_online <- function(L, T, t1, t2, beta, kappa, mul,
                            cg_col, evals, rundir, cg.ylim,
                            type="", csw=0, musigma=0, mudelta=0, muh=0, addon="",
                            skip=0, rectangle=TRUE,
                            plaquette=TRUE, dH=TRUE, acc=TRUE, trajtime=TRUE, omeas=TRUE,
                            plotsize=5, debug=FALSE, trajlabel=FALSE, title=FALSE,
                            pl=FALSE, method="uwerr", fit.routine="optim", oldnorm=FALSE, S=1.5,
                            stat_skip=0, omeas.samples=1, omeas.stride=1, omeas.avg=1,
                            omeas.start=0, omeas.stepsize=1, 
                            evals.start=0, evals.stepsize=1,
                            boot.R=1500, boot.l=2,
                            outname_suffix="", verbose=FALSE)
{
  # if we want to look at the online measurements in addition to output.data, we better provide
  # these parameters! 
  if( missing(L) | missing(T) | missing(beta) | missing(type) ){
    stop("L, T and beta must always be provided!\n");
  }
  if( omeas & (missing(t1) || missing(t2) || missing(kappa) || missing(mul)) ){
    stop("t1, t2, kappa and mul must be provided when 'omeas==TRUE'!");
  } else if(!omeas){
    kappa <- 0.0
    mul <- 0.0
    t1 <- 0
    t2 <- 0
  }

  # store analysis results in practical R format, replacing entries as new data is added
  resultsfile <- "omeas.summary.RData"
  
  resultsum <- list()
  if(file.exists(resultsfile)){
    cat("Loading analysis result database from ", resultsfile, "\n")
    load(resultsfile)
  }

  # vector with NAs to initialise result data frame
  navec <- t(data.frame(val=NA,dval=NA,tauint=NA,dtauint=NA,Wopt=NA,stringsAsFactors=FALSE))
  
  # set up data structure for analysis results 
  result <- list(params=data.frame(L=L,T=T,t1=t1,t2=t2,type=type,beta=beta,kappa=kappa,csw=csw,
                                   mul=mul,muh=muh,boot.l=boot.l,boot.R=boot.R,
                                   musigma=musigma,mudelta=mudelta,N.online=0,N.plaq=0,skip=skip,
                                   stat_skip=stat_skip,stringsAsFactors=FALSE),
                 obs=data.frame(mpcac_fit=navec, 
                                mpcac_mc=navec, 
                                mpi=navec, 
                                fpi=navec,
                                mpi_ov_fpi=navec, 
                                P=navec, 
                                dH=navec, 
                                expdH=navec, 
                                mineval=navec, 
                                maxeval=navec, 
                                CG.iter=navec, 
                                accrate=navec, 
                                trajtime=navec, stringsAsFactors=FALSE))
  
  errorband_color <- rgb(0.6,0.0,0.0,0.6)
  errorband_color2 <- rgb(0.0,0.0,0.6,0.6)
  
  if(missing(rundir)){
    rundir <- construct_rundir(type=type,beta=beta,L=L,T=T,kappa=kappa,mul=mul,
                             csw=csw,musigma=musigma,mudelta=mudelta,muh=muh,addon=addon,
                             debug=debug
                            )
  }
  
  filelabel <- rundir
  if(nchar(outname_suffix)>0){
    filelabel <- sprintf("%s_%s", filelabel, outname_suffix)
  }
      
  titletext <- NULL
  if(title) {
    titletext <- rundir
  } else {
    titletext <- ""
  }

  outfile <- sprintf("%s/output.data",rundir)
  plotcounter <- 0
  
  if(omeas){
    # read online measurements
    # get all omeas files that exist (hopefully in a consistent stepping)
    omeas.files <- getorderedfilelist(path=rundir, basename="onlinemeas", last.digits=6)
    # extract the trajectory numbers
    omeas.cnums <- getorderedconfignumbers(path=rundir, basename="onlinemeas", last.digits=6)
    # when the online measurements start later than traj 0 and we want to skip 'skip'
    # trajectories, the following should correspond to the correact indexing
    omeas.idx   <- c((1+as.integer(
                        (omeas.samples*ceiling((skip-omeas.start)/omeas.stepsize)))):length(omeas.files))
    omeas.files <- omeas.files[omeas.idx]
    omeas.cnums <- omeas.cnums[omeas.idx]
    pioncor <- readcmidatafiles( files=omeas.files, skip=0, 
                                 avg=omeas.avg, stride=omeas.stride, verbose=verbose )

    # when dealing with multi-sample online measurements, we need to thin out
    # the configuration numbers extracted above by stepping through them
    # with a stride of omeas.samples
    omeas.cnums <- omeas.cnums[seq(1, length(omeas.cnums), omeas.samples)]
    # add unique trajectory identifiers to the correlators
    pioncor <- cbind( pioncor, rep(omeas.cnums, each=3*(T/2+1) ) )

    result$params$N.online <- length(omeas.cnums)

    if(!any(class(pioncor)=='try-error')){
      # the correlation functions have been read externally, taking into account the measurement frequency
      # and possibly missing files. Therefore, skip=0! 
      onlineout <- onlinemeas(pioncor,
                              t1=t1,t2=t2,
                              kappa=kappa,
                              mu=mul,
                              skip=0,
                              method=method,
                              pl=pl,
                              fit.routine=fit.routine,
                              oldnorm=oldnorm,
                              S=S,
                              boot.R=boot.R,
                              boot.l=boot.l)

      if(debug){
        return(onlineout)
      }

      if(trajlabel){
        filelabel <- sprintf("%s_traj%06d-%06d",filelabel,min(omeas.cnums),max(omeas.cnums))
      }
      cat("Writing online measurements RData to ", sprintf("onlineout.%s.RData",filelabel), "\n")
      save(onlineout,file=sprintf("onlineout.%s.RData",filelabel))

      plotcounter <- plotcounter+1
      dpaopp_filename <- sprintf("%02d_dpaopp_%s",plotcounter,filelabel)
      result$obs$mpcac_mc <- plot_timeseries(dat=data.frame(y=onlineout$MChist.dpaopp,
                                                            t=omeas.cnums),
                                             stat_range=c( stat_skip+1, length(onlineout$MChist.dpaopp) ),
                                             pdf.filename=dpaopp_filename,
                                             ylab="$am_\\mathrm{PCAC}$",
                                             name="am_PCAC (MC history)",
                                             plotsize=plotsize,
                                             filelabel=filelabel,
                                             titletext=titletext,
                                             errorband_color=errorband_color)
                                             #ist.by=0.0002))
      
      # adjust autocorrelation times to be in terms of trajectories
      result$obs$mpcac_mc[3:5] <- result$obs$mpcac_mc[3:5]*omeas.stepsize

      lengthdpaopp <- length(onlineout$MChist.dpaopp)
      mindpaopp <- min(onlineout$MChist.dpaopp)
      maxdpaopp <- max(onlineout$MChist.dpaopp)
      
      mpcac_fit <- data.frame(val=(onlineout$fitresult$par[3]*onlineout$fitresult$par[2]/onlineout$fitresult$par[1]/2.),
                                  # note: the parameters are ordered differently between the result and the bootstrap samples
                                  dval=sd( onlineout$tsboot$t[,1]*onlineout$tsboot$t[,4]/(2*onlineout$tsboot$t[,3]) ),
                                  # no tauint from the fit
                                  tauint=NA, dtauint=NA, Wopt=NA, stringsAsFactors=FALSE)
      result$obs$mpcac_fit <- t(mpcac_fit)

      plotcounter <- plotcounter+1
      dpaopp_plateau_filename <- sprintf("%02d_dpaopp_plateau_%s",plotcounter,filelabel)
      tikzfiles <- tikz.init(basename=dpaopp_plateau_filename,width=plotsize,height=plotsize)
      op <- par(family="Palatino",cex.main=0.6,font.main=1)
      par(mgp=c(2,1,0))
      plotwitherror(x=onlineout$dpaopp$t,
                    y=onlineout$dpaopp$mass,dy=onlineout$dpaopp$dmass,t='p',
                    ylab="$am_\\mathrm{PCAC}$",
                    xlab="$t/a$",
                    main=titletext)
      rect(xleft=t1,
           xright=t2,
           ytop=mpcac_fit$val+mpcac_fit$dval,
           ybottom=mpcac_fit$val-mpcac_fit$dval,border=FALSE,col=errorband_color)
      abline(h=mpcac_fit$val,col="red",lwd=2)
      rect(xleft=t1,
           xright=t2,
           ytop=result$obs$mpcac_mc["val",]+result$obs$mpcac_mc["dval",],
           ybottom=result$obs$mpcac_mc["val",]-result$obs$mpcac_mc["dval",],border=FALSE,col=errorband_color2)
      abline(h=result$obs$mpcac_mc["val",],lwd=2,col="blue")
      legend(x="topright", bty='n', lty=1, lwd=4, col=c("red","blue"), 
            legend=c("$ M_\\mathrm{PS} |G_A| / 2|G_P| $ from 3-param ground-state fit",
                      sprintf("$ \\partial_0 \\langle A_0 P \\rangle / 2 \\langle P P \\rangle $ averaged from t=%d to t=%d",t1,t2) )
            )
      tikz.finalize(tikzfiles)
      
      plotcounter <- plotcounter+1
      mpi_plateau_filename <- sprintf("%02d_mpi_plateau_%s",plotcounter,filelabel)
      tikzfiles <- tikz.init(mpi_plateau_filename,width=plotsize,height=plotsize)
      op <- par(family="Palatino",cex.main=0.6,font.main=1)
      par(mgp=c(2,1,0))

      # sometimes the error computation for the effective mass fails
      # we deal with it here
      ploterror <- try(plotwitherror(x=onlineout$effmass$t,
                                     y=onlineout$effmass$m,dy=onlineout$effmass$dm,t='p',
                                     ylab="$aM_\\mathrm{PS}$",
                                     xlab="$t/a$",
                                     main=titletext),silent=FALSE)
      # and plot without errors if required
      if(inherits(ploterror,"try-error")) {
        plot(x=onlineout$effmass$t,y=onlineout$effmass$m)
      }
      rect(xleft=t1,
           xright=t2,
           ytop=onlineout$uwerrresultmps$value[1]+onlineout$uwerrresultmps$dvalue[1],
           ybottom=onlineout$uwerrresultmps$value[1]-onlineout$uwerrresultmps$dvalue[1],
           border=FALSE,
           col=errorband_color)
      abline(h=onlineout$uwerrresultmps$value[1],col="black")
      tikz.finalize(tikzfiles)

      result$obs$mpi <- t(data.frame(val=abs(onlineout$fitresultpp$par[2]),
                                    dval=onlineout$uwerrresultmps$dvalue[1],
                                    tauint=onlineout$uwerrresultmps$tauint[1]*omeas.stepsize,
                                    dtauint=onlineout$uwerrresultmps$dtauint[1]*omeas.stepsize,
                                    Wopt=onlineout$uwerrresultmps$Wopt[[1]]*omeas.stepsize, stringsAsFactors=FALSE) )

      result$obs$fpi <- t(data.frame(val=2*kappa*2*mul/sqrt(2)*abs(onlineout$fitresultpp$par[1])/
                                         (sqrt(onlineout$fitresultpp$par[2])*sinh(onlineout$fitresultpp$par[2])),
                                    dval=2*kappa*2*mul/sqrt(2)*onlineout$uwerrresultfps$dvalue[1],
                                    tauint=onlineout$uwerrresultfps$tauint[1]*omeas.stepsize,
                                    dtauint=onlineout$uwerrresultfps$dtauint[1]*omeas.stepsize,
                                    Wopt=onlineout$uwerrresultfps$Wopt[[1]]*omeas.stepsize, stringsAsFactors=FALSE) )

      if(method=="boot" | method=="all"){
        mpi_ov_fpi <- onlineout$tsboot$t[,1]/(2*kappa*2*mul/
                      sqrt(2)*(onlineout$tsboot$t[,3]/(sinh(onlineout$tsboot$t[,1])*sqrt(onlineout$tsboot$t[,1]))))
        result$obs$mpi_ov_fpi <- t(data.frame(val=mean(mpi_ov_fpi),
                                              dval=sd(mpi_ov_fpi),
                                              tauint=NA,
                                              dtauint=NA,
                                              Wopt=NA, stringsAsFactors=FALSE) )

      } else {
        # the storage format for the observables summary requires a bit of a transformation
        # to pass it to the naive error propagation function
        mpi_ov_fpi <- compute_ratio( as.list(result$obs$mpi[,1]), as.list(result$obs$fpi[,1]) )
        result$obs$mpi_ov_fpi <- t(data.frame(val=mpi_ov_fpi$val,
                                              dval=mpi_ov_fpi$dval,
                                              tauint=NA,
                                              dtauint=NA,
                                              Wopt=NA, stringsAsFactors=FALSE) )
      }

    } else { # (if error reading onlinemeas)
      stop("analysis_online: there was an error trying to read the output files (onlinemeas.xxxxxx/pionionline.dat)")
    }
  } # if(omeas)

  # something in the skip computation is odd, let's just solve it like this
  if(skip==0){
    shift <- 1
  } else {
    shift <- 0
  }

  outdat <- NULL
  tidx <- NULL
  if( plaquette || dH || trajtime || acc ) {
    # read output.data
    # in the first line of output.data which we want to consider 
    # (when the mass preconditioning is changed,
    # the number of columns may change so we need to be able to deal with that)
    no_columns <- max(count.fields(outfile,
                                   skip=skip))

    # if we have a rectange in the gauge action, the number of columns is different
    # in output.data
    tailvec <- c("acc","trajtime","rec")
    accshift <- length(tailvec)
    puregaugecols <- 7
    if(!rectangle){
      tailvec <- tailvec[1:2]
      accshift <- length(tailvec)
      puregaugecols <- 6
    }

    # the column names for the pure gauge case are given as follows 
    col.names <- c("traj","P","dH","expdH",tailvec)
    # but for the dynamical case, we have further information from the monomials
    if(no_columns > puregaugecols){
      col.names <- c("traj","P","dH","expdH",paste("V",5:(no_columns-accshift),sep=""),tailvec)
    }

    # read the entire file and then resize the number of columns to the target number
    # hopefully below, longer or short rows will be skipped because one generally skips
    # to the point where the trajectories become consistent 
    outdat <- read.table(outfile, fill=TRUE)
    outdat <- outdat[,1:no_columns]
    colnames(outdat) <- col.names

    tidx <- which(outdat$traj > skip)
    if(debug) print(outdat)
  }

  if(plaquette) {
    plotcounter <- plotcounter+1
    plaquette_filename <- sprintf("%02d_plaquette_%s",plotcounter,filelabel,title=filelabel)
    result$params$N.plaq <- length(tidx)-stat_skip
    result$obs$P <- plot_timeseries(dat=data.frame(y=outdat$P[tidx],
                                                   t=outdat$traj[tidx]),
                                    stat_range=c( 1+stat_skip, length( tidx ) ),
                                    pdf.filename=plaquette_filename,
                                    ylab="$ \\langle P \\rangle$" ,
                                    name="plaquette",
                                    plotsize=plotsize,
                                    filelabel=filelabel,
                                    titletext=titletext,
                                    errorband_color=errorband_color)
                                    #ist.by=0.00002))
  }
  if(dH) {
    plotcounter <- plotcounter+1
    dH_filename <- sprintf("%02d_dH_%s",plotcounter,filelabel)
    result$obs$dH <- plot_timeseries(dat=data.frame(y=outdat$dH[tidx],
                                                    t=outdat$traj[tidx]),
                                     stat_range=c( 1+stat_skip, length(tidx) ),
                                     pdf.filename=dH_filename,
                                     ylab="$ \\delta H $",
                                     name="dH",
                                     plotsize=plotsize,
                                     filelabel=filelabel,
                                     titletext=titletext,
                                     errorband_color=errorband_color,
                                     ylim=c(-2,3))
                                     #ist.by=0.2))

    plotcounter <- plotcounter+1
    expdH_filename <- sprintf("%02d_expdH_%s",plotcounter,filelabel)
    result$obs$expdH <- plot_timeseries(dat=data.frame(y=outdat$expdH[tidx],
                                                       t=outdat$traj[tidx]),
                                        stat_range=c( 1+stat_skip, length(tidx) ),
                                        pdf.filename=expdH_filename,
                                        ylab="$ \\exp(-\\delta H) $",
                                        name="exp(-dH)",
                                        plotsize=plotsize,
                                        filelabel=filelabel,
                                        titletext=titletext,
                                        errorband_color=errorband_color,
                                        hist.xlim=c(-2,4),
                                        ylim=c(-0,6))
                                        #ist.by=0.2))
  }
  if( !missing("cg_col") ) {
    plotcounter <- plotcounter+1
    cg_filename <- sprintf("%02d_cg_iter_%s", plotcounter, filelabel)
    result$obs$CG.iter <- plot_timeseries(dat=data.frame(y=outdat[tidx,cg_col],
                                                         t=outdat$traj[tidx]),
                                          stat_range=c( 1+stat_skip, length(tidx) ),
                                          pdf.filename=cg_filename,
                                          ylab="$N^\\mathrm{iter}_\\mathrm{CG}$",
                                          name="CG iterations",
                                          plotsize=plotsize,
                                          filelabel=filelabel,
                                          titletext=titletext,
                                          errorband_color=errorband_color)
                                          #ist.by=5))
  }
  if( !missing("evals") ) {
    plotcounter <- plotcounter+1
    ev_pdf_filename <- sprintf("%02d_evals_%02d_%s", plotcounter, evals, filelabel )
    ev_filename <- sprintf("%s/monomial-%02d.data", rundir, evals )
    
    evaldata <- tryCatch(read.table(ev_filename, stringsAsFactors=FALSE, 
                                    col.names=c("traj","min_ev","max_ev","ev_range_min","ev_range_max") ), 
                         error=function(e){ stop(sprintf("Reading of %s failed!",ev_filename)) } )

    eval.tidx <- which(evaldata$traj > skip)

    temp <- plot_eigenvalue_timeseries(dat=evaldata[eval.tidx,],
                                       stat_range=c( 1+stat_skip, length(eval.tidx) ),
                                       pdf.filename = ev_pdf_filename,
                                       ylab = "eigenvalue",
                                       plotsize=plotsize,
                                       filelabel=filelabel,
                                       titletext=titletext,
                                       errorband_color=errorband_color )
    temp$obs$mineval[3:5] <- temp$mineval[3:5]*evals.stepsize
    temp$obs$maxeval[3:5] <- temp$maxeval[3:5]*evals.stepsize
    
    result$obs$mineval <- temp$mineval
    result$obs$maxeval <- temp$maxeval
  }
  if( acc == TRUE ){
    # finally add acceptance rate
    plotcounter <- plotcounter+1
    accrate_filename <- sprintf("%02d_accrate_%s",plotcounter,filelabel,title=filelabel)
    result$obs$accrate <- plot_timeseries(dat=data.frame(y=outdat$acc[tidx],
                                                         t=outdat$traj[tidx]),
                                          stat_range=c( 1+stat_skip, length(tidx) ),
                                          pdf.filename=accrate_filename,
                                          ylab="Accept / Reject" ,
                                          name="accrate",
                                          plotsize=plotsize,
                                          filelabel=filelabel,
                                          titletext=titletext,
                                          errorband_color=errorband_color,
                                          hist.by=0.5)
  }
  if( trajtime == TRUE ){
    plotcounter <- plotcounter+1
    trajtime_filename <- sprintf("%02d_trajtime_%s",plotcounter,filelabel,title=filelabel)
    result$obs$trajtime <- plot_timeseries(data.frame(y=outdat$trajtime[tidx],
                                                      t=outdat$traj[tidx]),
                                           stat_range=c( 1+stat_skip, length(tidx) ),
                                           pdf.filename=trajtime_filename,
                                           ylab="Traj. time [s]" ,
                                           name="trajtime",
                                           plotsize=plotsize,
                                           filelabel=filelabel,
                                           titletext=titletext,
                                           errorband_color=errorband_color)
  }

  print(result$params)
  print(t(result$obs))

  # if the script "pdfcat" exists, concatenate the plots and remove the individual files
  if( Sys.which("pdfcat") != "" ) {
    commands <- c(sprintf("pdfcat analysis_%s.pdf ??_*_%s.pdf",filelabel,filelabel),
                  sprintf("rm -f ??_*_%s.pdf",filelabel) )
    for( command in commands ) {
      print(paste("calling",command))
      system(command=command)
    } 
  } else {
    print("pdfcat not found, not concatenating plots!")
  }
  if( Sys.which("pdfcrop") != "" ) {
    commands <- c(sprintf("pdfcrop analysis_%s.pdf analysis_%s.pdf",filelabel,filelabel))
    for( command in commands ) {
      print(paste("calling",command))
      system(command=command)
    } 
  } else {
    print("pdfcrop not found, not cropping plots!")
  }

  resultsum[[rundir]] <- result
  cat("Storing analysis result database in ", resultsfile, "\n")
  save(resultsum,file=resultsfile)

  return(invisible(result))
}

construct_rundir <- function(type,beta,L,T,kappa=0,mul=0,csw=0,musigma=0,mudelta=0,muh=0,addon="",debug=FALSE) {
  rundir <- NULL
  rundir <- sprintf("%s_b%s-L%dT%d",type,beta,L,T)

  if(csw!=0) {
    rundir <- sprintf("%s-csw%s",rundir,csw)
  }

  # mul=0 && kappa=0 means pure gauge, which is perhaps a bit misleading...
  if( mul != 0 && kappa != 0 ) {
    # silly sprintf prints numbers smaller than 0.001 in scientific notation unless g is used
    # on the other hand, for larger numbers, trailing zeroes are added...
    # so we use %g only in the former case and pass mu as a string otherwise!
    if(mul >= 1e-3) {
      rundir <- sprintf("%s-k%s-mul%s",rundir,kappa,mul)
    } else {
      rundir <- sprintf("%s-k%s-mul%g",rundir,kappa,mul)
    }
  }

  if(muh != 0) {
    rundir <- sprintf("%s-muh%s",rundir,muh)
  }

  if(musigma!=0){
    rundir <- sprintf("%s-musigma%s",rundir,musigma)
  }
  if(mudelta!=0){
    rundir <- sprintf("%s-mudelta%s",rundir,mudelta)
  }
  if(addon!=""){
    rundir <- sprintf("%s_%s",rundir,addon)
  }

  if(debug) {
    cat("Trying to read from directory:", rundir,"\n")
  }

  rundir
}
