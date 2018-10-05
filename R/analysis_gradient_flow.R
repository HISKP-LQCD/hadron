# these should also have references
ref_gf_scales <- list()
ref_gf_scales[["sqrt_t0_Eplaq"]] <- list(val=0.1416, dval=0.0008,
                                         label="t_0^{\\mathrm{plaq}}/a^2", 
                                         obs="tsqEplaq",
                                         obslabel="$\\langle t^2 E_{\\mathrm{plaq}}(t) \\rangle$",
                                         physobslabel="\\langle t_0^2 E_{\\mathrm{plaq}}(t_0) \\rangle")

ref_gf_scales[["sqrt_t0_Esym"]] <- list(val=0.1416, dval=0.0008,
                                        label="t_0^{\\mathrm{sym}}/a^2",
                                        obs="tsqEsym",
                                        obslabel="$\\langle t^2 E_{\\mathrm{sym}}(t) \\rangle$",
                                        physobslabel="\\langle t_0^2 E_{\\mathrm{sym}}(t_0) \\rangle")

ref_gf_scales[["w0_Wsym"]] <- list(val=0.1755, dval=0.0019,
                                   label="(w_0^{\\mathrm{sym}}/a)^2",
                                   obs="Wsym",
                                   obslabel="$\\langle W_{\\mathrm{sym}}(t) \\rangle$",
                                   physobslabel="\\langle W_{\\mathrm{sym}}(w_0^2) \\rangle")

analysis_gradient_flow <- function(path,outputbasename,basename="gradflow",read.data=TRUE,pl=FALSE,skip=0,start=0,scale=1,dbg=FALSE) {

  # much like for the analysis of online measurements, we store summary information
  # in the list "gradflow_resultsum" with elements named by "outputbasename"
  # such that when the analysis is rerun, the entries are replaced
  resultsfile <- "gradflow.summary.RData"
  gradflow_resultsum <- list()
  if(file.exists(resultsfile)){
    cat("Loading gradient flow results database from ", resultsfile, "\n")
    load(resultsfile)
  }

  if(read.data) {
    raw.gradflow <- readgradflow(path=path, skip=skip, basename=basename)
    save(raw.gradflow,file=sprintf("%s.raw.gradflow.Rdata",outputbasename),compress=FALSE)
  }else{
    cat(sprintf("Warning, reading data from %s.raw.gradflow.Rdata, if the number of samples changed, set read.data=TRUE to reread all output files\n",basename))
    load(sprintf("%s.raw.gradflow.Rdata",basename))
  }
  if(dbg==TRUE) print(raw.gradflow)

  t_vec <- unique(raw.gradflow$t)
  Ncol <- ncol(raw.gradflow)
  Nrow <- length(t_vec)
  # allocate some memory, for each observable, have space for the value, the error and the autocorrelation time
  # create a list with NULL rownams and sensible column names for this purpose
  subnames <- c("value","dvalue","ddvalue","tauint","dtauint")
  cnames <- colnames(raw.gradflow[,3:Ncol])
  outer.cnames <- t(outer(cnames,subnames,FUN=function(x,y) { paste(x,y,sep=".") }))
  grad.dimnames <- list()
  grad.dimnames[[1]] <- NULL
  grad.dimnames[[2]] <- c("t",as.vector(outer.cnames) )
  
  gradflow <- as.data.frame(matrix(data=NA,nrow=Nrow,ncol=(length(subnames)*(Ncol-2)+1),dimnames=grad.dimnames))
  for(i_row in 1:length(t_vec)){
    uwerr.gradflow <- apply(X=raw.gradflow[which(raw.gradflow$t==t_vec[i_row]),3:Ncol],MARGIN=2,FUN=uwerrprimary)
    summaryvec <- c(t_vec[i_row])
    for(i_col in 1:length(cnames)) {
      obs <- uwerr.gradflow[[cnames[i_col]]]
      summaryvec <- c(summaryvec, obs$value, obs$dvalue, obs$ddvalue, obs$tauint*scale, obs$dtauint*scale)
    }
    gradflow[i_row,] <- summaryvec
  }
  save(gradflow,file=sprintf("%s.result.gradflow.Rdata",outputbasename),compress=FALSE)
 
  gf_scales <- list()
  gf_latspacs <- list()
  gf_approx_scales <- list()
  for( i in 1:length(ref_gf_scales) ){
    val <- gradflow[,sprintf("%s.value", ref_gf_scales[[i]]$obs)]
    dval <- gradflow[,sprintf("%s.dvalue", ref_gf_scales[[i]]$obs)] 
    gf_scales[[i]] <- c(approx( x=val+dval,y=gradflow$t,xout=0.3)$y, 
                        approx( x=val, y=gradflow$t, xout=0.3 )$y, 
                        approx( x=val-dval, y=gradflow$t, xout=0.3)$y )

    ## determine which discrete value of t is closest to the scale in question
    for( tidx in 1:length(t_vec) ){
      if( t_vec[tidx] >= gf_scales[[i]][2] ){
        gf_approx_scales[[i]] <- t_vec[tidx]
        break
      }
    }

    # sqrt of our bounds to compute the lattice spacing
    sqrt_gf_scale <- sqrt(gf_scales[[i]])

    gf_latspacs[[i]] <- list()
    for( k in 1:length(sqrt_gf_scale) ){
      gf_latspacs[[i]][[k]] <- compute_ratio(dividend=ref_gf_scales[[i]],
                                             divisor=list(val=sqrt_gf_scale[k],
                                                          dval=0))
    }
  }
   
  if(pl) {
    tikzfiles <- tikz.init(basename=sprintf("%s.gradflow",outputbasename),width=4,height=4)
    for( i in 1:length(ref_gf_scales) ){
      scale_obs <- ref_gf_scales[[i]]$obs
      scale_obslabel <- ref_gf_scales[[i]]$obslabel
      scale_label <- ref_gf_scales[[i]]$label

      val <- gradflow[,sprintf("%s.value", scale_obs)]
      dval <- gradflow[,sprintf("%s.dvalue", scale_obs)]

      gf_scale <- gf_scales[[i]]
      gf_latspac <- gf_latspacs[[i]]

      # set up plot
      # xlim is set to 1.25 times the value of t corresponding to the upper
      # limit of the reference scale in lattice units
      plot(x=gradflow$t, 
           y=val,
           type='n',
           xlim=c(0,1.25*gf_scale[1]),
           ylim=c(0.0,0.4),
           xlab="$t/a^2$",
           ylab=scale_obslabel,
           las=1)
      # draw errorband
      poly.col <- rgb(red=1.0,green=0.0,blue=0.0,alpha=0.6)
      poly.x <- c(gradflow$t,rev(gradflow$t))
      poly.y <- c(val+dval,rev(val-dval))
      polygon(x=poly.x,y=poly.y,col=poly.col)
      abline(h=0.3)
      abline(v=gf_scale)
      lines(x=gradflow$t, y=val)
      legend(x="topleft",
             legend=c(sprintf("$%s = %s$",
                              scale_label,
                              tex.catwitherror(x=gf_scale[2],
                                               # for the error, we use the average of the plus and minus errors
                                               dx=0.5*(abs(gf_scale[3]-gf_scale[2]) + abs(gf_scale[2]-gf_scale[1])),
                                               digits=2,
                                               with.dollar=FALSE)
                              ),
                    sprintf("$a=%s$\\,fm", 
                            tex.catwitherror(x=gf_latspac[[2]]$val,
                                             dx=sqrt( (0.5* ( 
                                                             abs( gf_latspac[[3]]$val-gf_latspac[[2]]$val) +
                                                             abs( gf_latspac[[1]]$val-gf_latspac[[2]]$val) 
                                                             ) 
                                                       )^2 + gf_latspac[[2]]$dval^2 ),
                                             digits=2,
                                             with.dollar=FALSE)
                            )
                    ),
             bty='n')
      
      ### plot MD history of Wsym at t closest to scale
      approx_idx <- which(raw.gradflow$t==gf_approx_scales[[i]])
      tseries <- data.frame(y=raw.gradflow[approx_idx,scale_obs],
                            t=start + c( skip :(skip + 
                                                length(raw.gradflow[approx_idx,scale_obs]) - 1 ) )*scale )
      plot_timeseries(dat=tseries,
                      ylab=sprintf("%s$|_{t/a^2 = %.2f}$", 
                                   scale_obslabel,
                                   gf_approx_scales[[i]]),
                      titletext="")
    }
   
    if( any(cnames == "Qsym") ){
      ################ TOPOLOGICAL CHARGE ####################
      # set up plot
      plot(x=gradflow$t, 
           y=gradflow$Qsym.value,
           type='n',
           ylim=range(c(gradflow$Qsym.value+gradflow$Qsym.dvalue, gradflow$Qsym.value-gradflow$Qsym.dvalue)),
           xlim=c(0.01,max(gradflow$t)),
           xlab="$t/a^2$",
           ylab="$Q(t)$",
           las=1,
           log='x')
      # draw errorband
      poly.col <- rgb(red=1.0,green=0.0,blue=0.0,alpha=0.6)
      poly.x <- c(gradflow$t,rev(gradflow$t))
      poly.y <- c(gradflow$Qsym.value+gradflow$Qsym.dvalue,rev(gradflow$Qsym.value-gradflow$Qsym.dvalue))
      polygon(x=poly.x,y=poly.y,col=poly.col)
      lines(x=gradflow$t, y=gradflow$Qsym.value)
      
      # plot MD history of Q at maximal flow time
      tmax_idx <- which(raw.gradflow$t==max(raw.gradflow$t))
      tseries <- data.frame(y=raw.gradflow[tmax_idx,"Qsym"],
                            t=start + c( skip :( skip + length(raw.gradflow[tmax_idx,"Qsym"]) - 1 ) )*scale)
      plot_timeseries(dat=tseries,
                      ylab=sprintf("$Q\\left( t/a^2 = %.2f \\right)$",max(raw.gradflow$t)),
                      titletext="")
    }
    
    tikz.finalize(tikzfiles)
  }
}

