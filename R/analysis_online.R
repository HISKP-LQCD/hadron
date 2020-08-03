append_pdf_filename <- function(basename, pdf_filenames){
  return( c(pdf_filenames, sprintf("%s.pdf", basename)) )
}

#' analysis_online
#'
#' @description
#' `analysis_online` is a function to analyse the online measurements
#'  and output files of the tmLQCD software, see references. The function
#'  operates on a subdirectory either passed via `rundir` or automatically
#'  constructed from the various function arguments. Depending on which
#'  parts of the analysis are requested, this subdirectory is expected
#'  to contain `onlinemeas.%06d` files with online correlator measurements,
#'  `output.data` containing the plaquette and energy violation, amongst others 
#'  and `monomial-%02d.data` with measurements of the extremal eigenvalues of the
#   non-degenerate twisted mass Dirac operator.
#' 
#' @param L integer. spatial lattice extent
#' @param Time integer. temporal lattice extent
#' @param t1 integer. initial time of fit range
#' @param t2 integer. end time of fit range
#' @param beta numeric. inverse squared gauge coupling
#' @param kappa numeric. hopping parameter
#' @param mul numeric. light sea twisted quark mass
#' @param cg_col integer. column of CG iteration counts from `output.data` to use
#' @param evals_id Integer. Monomial ID of the monomial for which eigenvalues are measured.
#'                          Function will attempt to open `monomial-%02d.data`.
#' @param rundir string. run directory. If not specified, run directory will
#'                       be constructed automatically. See \link{construct_onlinemeas_rundir} for
#'                       details.
#' @param cg.ylim numeric. y-limits for CG iteration counts
#' @param type string. Type specifier for the gauge action, this might be 'iwa' for Iwasaki,
#'                     for example.
#' @param csw numeric. clover coefficient
#' @param musigma numeric. average 1+1 sea twisted quark mass
#' @param mudelta numeric. splitting 1+1 sea twisted quark mass
#' @param muh numeric. "heavy" twisted mass in the case of a `n_f=2+2` run
#' @param addon string. addon to output filenames
#' @param skip integer. number of initial measurements to skip in analysis
#' @param rectangle boolean. If true, rectangle plaquettes are analysed
#' @param plaquette boolean. If true, square plaquettes are analysed
#' @param dH boolean. If true, delta H is analysed
#' @param acc boolean. If true, the acceptance rate is analysed
#' @param trajtime boolean. If true, the time per trajectory is analysed
#' @param omeas boolean. If true, online measurements are analysed (`onlinemeas.%06d`)
#' @param plotsize numeric. size of plots being generated
#' @param debug boolean. provide debug information
#' @param trajlabel boolean or string. If not `FALSE`, use as trajectory labels
#' @param title bolean or string. If not `FALSE`, use as main title of plots
#' @param pl boolean. If set to `TRUE` plots will be generated
#' @param method string. method to compute errors, can be "uwerr", "boot" or "all"
#' @param fit.routine string. minimisation routine for chisq, can be "optim"
#' @param oldnorm boolean. If `TRUE`, the function assumes that the `onlinemeas.%06d` are
#'                         in old tmLQCD normalisation.
#' @param S numeric. `S` parameter of \link{uwerr}
#' @param stat_skip integer. By passing this parameter, the various timeseries plots
#'                           will include `stat_skip` measurements,
#'                           but these will be skipped in the corresponding statistical analysis.
#'                           This maybe useful, for example, to visualise thermalisation.
#' @param omeas.samples integer. number of stochastic samples per online measurement
#' @param omeas.stride integer. stride length in the reading of online measurements
#' @param omeas.avg integer. Block average over this many subsequent measurements.
#' @param omeas.stepsize integer. Number of trajectories between online measurements. Autocorrelation
#'                                times of online measurement data will be scaled by this factor.
#' @param evals.stepsize integer. Numer of trajectories between (strange-charm Dirac opertoar) eigenvalue measurements.
#'                                Autocorrelation times of eigenvalues will be scaled by this factor.
#' @param boot.R integer. number of bootstrap samples to use in bootstrap-based parts of analysis.
#' @param boot.l integer. bootstrap block size
#' @param outname_suffix string. suffix for output files
#' @param verbose boolean. If `TRUE`, function produces verbose output.
#' #'
#' @return
#' a list is returned with all the accumulated results. Moreover, a
#'   PDF file with statistics and analytics is created and the results
#'   are written into .Rdata files. On the one hand, the result of the
#'   call to the \link{onlinemeas} function is written to
#'   `onlineout.%s.Rdata`, where `%s` is replaced with a label
#'   built from meta information based on the arguments above.
#'   On the other hand, summary data across many calls of this
#'   function is silently accumulated in the file
#'   `omeas.summary.Rdata` which contains the named list 'resultsum' with
#'   element names based on `rundir`.
#'
#' @references
#' K. Jansen and C. Urbach, Comput.Phys.Commun. 180 (2009) 2717-2738
#' 
#' @export
analysis_online <- function(L, Time, t1, t2, beta, kappa, mul,
                            cg_col, evals_id, rundir, cg.ylim,
                            type="", csw=0, musigma=0, mudelta=0, muh=0, addon="",
                            skip=0, rectangle=TRUE,
                            plaquette=TRUE, dH=TRUE, acc=TRUE, trajtime=TRUE, omeas=TRUE,
                            plotsize=5, debug=FALSE, trajlabel=FALSE, title=FALSE,
                            pl=FALSE, method="uwerr", fit.routine="optim", oldnorm=FALSE, S=1.5,
                            stat_skip=0, omeas.samples=1, omeas.stride=1, omeas.avg=1,
                            omeas.stepsize=1, 
                            evals.stepsize=1,
                            boot.R=1500, boot.l=2,
                            outname_suffix="", verbose=FALSE)
{
  staplr_avail <- requireNamespace("staplr")
  if( !staplr_avail ){
    stop("The 'staplr' package is required")
  }

  pdf_filenames <- c()

  # if we want to look at the online measurements in addition to output.data, we better provide
  # these parameters! 
  if( missing(L) | missing(Time) | missing(beta) | missing(type) ){
    stop("L, Time and beta must always be provided!\n");
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
    message("Loading analysis result database from ", resultsfile, "\n")
    load(resultsfile)
  }

  # vector with NAs to initialise result data frame
  navec <- t(data.frame(val=NA,dval=NA,tauint=NA,dtauint=NA,Wopt=NA,stringsAsFactors=FALSE))
  
  # set up data structure for analysis results 
  result <- list(params=data.frame(L=L,Time=Time,t1=t1,t2=t2,type=type,beta=beta,kappa=kappa,csw=csw,
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
    rundir <- construct_onlinemeas_rundir(type=type,beta=beta,L=L,Time=Time,kappa=kappa,mul=mul,
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
    omeas.idx <- which(omeas.cnums > skip)
    
    if( length(omeas.idx) < 1 ){
      stop(sprintf("After skipping %d trajectories, no online measurements are left!\n"))
    }


    omeas.files <- omeas.files[omeas.idx]
    omeas.cnums <- omeas.cnums[omeas.idx]
    pioncor <- readcmidatafiles( files=omeas.files, skip=0, 
                                 avg=omeas.avg, stride=omeas.stride, verbose=verbose )

    # when dealing with multi-sample online measurements, we need to thin out
    # the configuration numbers extracted above by stepping through them
    # with a stride of omeas.samples
    omeas.cnums <- omeas.cnums[seq(1, length(omeas.cnums), omeas.samples)]
    # add unique trajectory identifiers to the correlators
    pioncor <- cbind( pioncor, rep(omeas.cnums, each=3*(Time/2+1) ) )

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
      message("Writing online measurements RData to ", sprintf("onlineout.%s.RData",filelabel), "\n")
      save(onlineout,file=sprintf("onlineout.%s.RData",filelabel))

      plotcounter <- plotcounter+1
      dpaopp_filename <- sprintf("%02d_dpaopp_%s",plotcounter,filelabel)
      pdf_filenames <- append_pdf_filename(basename = dpaopp_filename,
                                           pdf_filenames = pdf_filenames)
      result$obs$mpcac_mc <- plot_timeseries(dat=data.frame(y=onlineout$MChist.dpaopp,
                                                            t=omeas.cnums),
                                             stat_range=c( stat_skip+1, length(onlineout$MChist.dpaopp) ),
                                             pdf.filename=dpaopp_filename,
                                             ylab="$am_\\mathrm{PCAC}$",
                                             name="am_PCAC (MC history)",
                                             plotsize=plotsize,
                                             filelabel=filelabel,
                                             titletext=titletext,
                                             errorband_color=errorband_color,
                                             smooth_density = TRUE)
      
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
      pdf_filenames <- append_pdf_filename(basename = dpaopp_plateau_filename,
                                           pdf_filenames = pdf_filenames)

      tikzfiles <- tikz.init(basename=dpaopp_plateau_filename,width=plotsize,height=plotsize)
      op <- par(family="Palatino",cex.main=0.6,font.main=1)
      on.exit(par(op))
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
      pdf_filenames <- append_pdf_filename(basename = mpi_plateau_filename,
                                           pdf_filenames = pdf_filenames)

      tikzfiles <- tikz.init(mpi_plateau_filename,width=plotsize,height=plotsize)
      op <- par(family="Palatino",cex.main=0.6,font.main=1)
      on.exit(par(op))
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
    
    # skip 'skip' lines. Hopefully the trajectory counter doesn't start at something much larger than 0
    # because otherwise we'll probably miss trajectories in the filtering below
    outdat <- read.table(outfile, skip=skip, fill=TRUE)
    no_rows <- nrow(outdat)

    # count the number of columns in the penultimate line of output.data
    # in the first line of output.data which we want to consider 
    # (when the mass preconditioning is changed,
    # the number of columns may change so we need to be able to deal with that)
    no_columns <- max(count.fields(outfile,
                                   skip=skip+no_rows-1))
    
    # restrict any outliers in outdat. Inconsistent lines will almost certainly
    # be filtered out below 
    outdat <- outdat[,1:no_columns]

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

    colnames(outdat) <- col.names

    # restrict to the lines which are definitely the trajectories we want to 
    # analyse
    tidx <- which(outdat$traj > skip)
    if(debug) print(outdat)
  }

  if(plaquette) {
    plotcounter <- plotcounter+1
    plaquette_filename <- sprintf("%02d_plaquette_%s",plotcounter,filelabel,title=filelabel)
    pdf_filenames <- append_pdf_filename(basename = plaquette_filename,
                                        pdf_filenames = pdf_filenames)

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
                                    errorband_color=errorband_color,
                                    smooth_density = TRUE)
  }
  if(dH) {
    plotcounter <- plotcounter+1
    dH_filename <- sprintf("%02d_dH_%s",plotcounter,filelabel)
    pdf_filenames <- append_pdf_filename(basename = dH_filename,
                                         pdf_filenames = pdf_filenames)

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
                                     smooth_density = TRUE,
                                     ylim=c(-2,3),
                                     hist.probs=c(0.01,0.99))

    plotcounter <- plotcounter+1
    expdH_filename <- sprintf("%02d_expdH_%s",plotcounter,filelabel)
    pdf_filenames <- append_pdf_filename(basename = expdH_filename,
                                         pdf_filenames = pdf_filenames)

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
                                        smooth_density = TRUE,
                                        ylim=c(0,6),
                                        hist.probs=c(0.01,0.99))
  }
  if( !missing("cg_col") ) {
    plotcounter <- plotcounter+1
    cg_filename <- sprintf("%02d_cg_iter_%s", plotcounter, filelabel)
    pdf_filenames <- append_pdf_filename(basename = cg_filename,
                                         pdf_filenames = pdf_filenames)

    result$obs$CG.iter <- plot_timeseries(dat=data.frame(y=outdat[tidx,cg_col],
                                                         t=outdat$traj[tidx]),
                                          stat_range=c( 1+stat_skip, length(tidx) ),
                                          pdf.filename=cg_filename,
                                          ylab="$N^\\mathrm{iter}_\\mathrm{CG}$",
                                          name="CG iterations",
                                          plotsize=plotsize,
                                          filelabel=filelabel,
                                          titletext=titletext,
                                          errorband_color=errorband_color,
                                          smooth_density = TRUE
                                          )
  }
  if( !missing("evals_id") ) {
    plotcounter <- plotcounter+1
    ev_pdf_filename <- sprintf("%02d_evals_%02d_%s", plotcounter, evals_id, filelabel )
    ev_filename <- sprintf("%s/monomial-%02d.data", rundir, evals_id )
    pdf_filenames <- append_pdf_filename(basename = ev_pdf_filename,
                                         pdf_filenames = pdf_filenames)
    
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
    accrate_filename <- sprintf("%02d_accrate_%s",plotcounter,filelabel)
    pdf_filenames <- append_pdf_filename(basename = accrate_filename,
                                         pdf_filenames = pdf_filenames)

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
                                          hist.by=0.5,
                                          smooth_density = FALSE)
  }
  if( trajtime == TRUE ){
    plotcounter <- plotcounter+1
    trajtime_filename <- sprintf("%02d_trajtime_%s",plotcounter,filelabel,title=filelabel)
    pdf_filenames <- append_pdf_filename(basename = trajtime_filename,
                                         pdf_filenames = pdf_filenames)


    result$obs$trajtime <- plot_timeseries(data.frame(y=outdat$trajtime[tidx],
                                                      t=outdat$traj[tidx]),
                                           stat_range=c( 1+stat_skip, length(tidx) ),
                                           pdf.filename=trajtime_filename,
                                           ylab="Traj. time [s]" ,
                                           name="trajtime",
                                           plotsize=plotsize,
                                           filelabel=filelabel,
                                           titletext=titletext,
                                           errorband_color=errorband_color,
                                           smooth_density = TRUE,
                                           hist.probs = c(0.01,0.99))
  }

  print(result$params)
  print(t(result$obs))

  # concatenate the individual files
  staplr::staple_pdf(input_files = pdf_filenames,
                     output_filepath = sprintf("analysis_%s.pdf", filelabel))
  # and remove them
  system(command=sprintf("rm -f ??_*_%s.pdf", filelabel))

  resultsum[[rundir]] <- result
  message("Storing analysis result database in ", resultsfile, "\n")
  save(resultsum,file=resultsfile)

  return(invisible(result))
}

#' Construct a run directory string for \link{analysis_online}
#'
#' @param type String. Short identifier for gauge action type. For example, `iwa` for Iwasaki gauge action.
#' @param beta Numeric. Inverse gauge coupling.
#' @param L Integer. Spatial lattice extent.
#' @param Time Integer. Temporal lattice extent.
#' @param kappa Numeric. Sea quark action hopping parameter.
#' @param mul Numeric. Sea light quark twisted mass.
#' @param csw Numeric. Sea quark action clover parameter.
#' @param musigma Numeric. Sea 1+1 "heavy" average twisted quark mass.
#' @param mudelta Numeric. Sea 1+1 "heavy" splitting twisted quark mass.
#' @param muh Numeric. In case of `n_f=2+2` run, "heavy" twisted quark mass.
#' @param addon String. Arbitratry string which will be suffixed to the constructed run directory.
#' @param debug Boolean. If `TRUE`, the constructed directory name is printed to screen.
#'
#' @return String. Directory name constructed out of the various function parameters. See source code for details.
#'
construct_onlinemeas_rundir <- function(type,beta,L,Time,kappa=0,mul=0,csw=0,musigma=0,mudelta=0,muh=0,addon="",debug=FALSE) {
  rundir <- NULL
  rundir <- sprintf("%s_b%s-L%dT%d",type,beta,L,Time)

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
    message("Trying to read from directory:", rundir,"\n")
  }

  return(rundir)
}
