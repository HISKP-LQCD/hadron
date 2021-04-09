#' computes a disconnected correlation function from loops
#' 
#' The dimension of \code{cf$cf} and \code{cf$icf} must be \code{dim(Time, S, N)},
#' where \code{Time} is the time extent, \code{S} is the number of samples and
#' \code{N} the number of measurements (gauges). \code{cf2} is the same, but
#' needed only for cross-correlators.
#' 
#' If \code{subtract.vev=TRUE} the vev is estimated as the mean over all
#' gauges, samples and times available and subtracted from the original loop
#' data. (Same for \code{subtrac.vev2}.
#' 
#' The correlation is computed such as to avoid correlation between equal
#' samples, unless \code{nrSamples} is equal to 1.
#' 
#' \code{cf} and \code{cf2} must agree in \code{Time}, number of gauges and number
#' of samples. Matching of gauges is assumed. If this is not the case results
#' are wrong.
#' 
#' @param cf loop data as produced by \code{readcmidisc} or
#' \code{readbinarydisc}.
#' @param cf2 second set of loop data as produced by \code{readcmidisc} or
#' \code{readbinarydisc}. This is needed for cross-correlators
#' @param real use the real part \code{cf$cf}, if set to \code{TRUE}, otherwise
#' the imaginary part \code{cf$icf}.
#' @param real2 use the real part \code{cf2$cf}, if set to \code{TRUE},
#' otherwise the imaginary part \code{cf2$icf}.
#' @param subtract.vev subtract a vacuum expectation value. It will be
#' estimated as mean over all samples, gauges and times available.
#' @param subtract.vev2 subtract a vacuum expectation value for the second set
#' of loops. It will be estimated as mean over all samples, gauges and times
#' available.
#' @param subtract.equal subtract contributions of products computed on
#' identical samples. This will introduce a bias, if set to FALSE for missing
#' cf2 or if cf and cf2 are computed on the same set of random sources.
#' @param use.samples If set to an integer, only the specified number of
#' samples will be used for \code{cf}, instead of all samples.
#' @param use.samples2 Same like \code{use.samples}, but for \code{cf2}.
#' @param smeared use the loops instead of the local ones for \code{cf}.
#' @param smeared2 use the loops instead of the local ones for \code{cf2}.
#' @param type The correlation function can either be symmetric or
#' anti-symmetric in time. Anti-symmetric is of course only possible for
#' cross-correlators. In this case with \code{type="cosh"} it is assumed to be
#' symmetric, anti-symmetric otherwise.
#' @param verbose Print some debug output, like the VEVs of the loops.
#' @return Returns an object of type \code{cf} derived from a \code{list} with
#' elements \code{cf}, an array of dimension \code{dim(N, Time)}, where \code{N}
#' is the number of samples and \code{Time} the time extent, integers \code{Time}
#' for the time extent, \code{nrStypes} and \code{nrObs} for the available
#' smearing types and operators, and finally \code{nrSamples}, the number of
#' samples used to generate the correlation function \code{cf}.
#' @author Carsten Urbach, \email{curbach@@gmx.de}
#' @seealso \code{\link{readcmidisc}}, \code{\link{readbinarydisc}},
#' \code{\link{bootstrap.cf}}, \code{\link{add.cf}}, \code{\link{c.cf}}
#' @keywords correlator
#' @examples
#' 
#' data(loopdata)
#' Cpi0v4 <- computeDisc(cf=loopdata, real=TRUE, subtract.vev=TRUE)
#' Cpi0v4 <- bootstrap.cf(Cpi0v4, boot.R=99, boot.l=1, seed=14556)
#' 
#' @export computeDisc
computeDisc <- function(cf, cf2,
                        real=TRUE, real2 = TRUE,
                        smeared=FALSE, smeared2=FALSE,
                        subtract.vev=TRUE, subtract.vev2=TRUE,
                        subtract.equal = TRUE,
                        use.samples, use.samples2,
                        type="cosh", verbose=FALSE) {
  stopifnot(inherits(cf, 'cf_meta'))
  stopifnot(inherits(cf, 'cf_orig'))

  Time <- cf$Time
  ## extract the corresponding part of the correlation matrix
  tcf <- cf$cf
  if(smeared) {
    tcf <- cf$scf
  }
  if(!real) tcf <- cf$icf
  if(!real && smeared) tcf <- cf$sicf
  
  nrSamples <- cf$nrSamples
  nrSamples2 <- nrSamples
  if(!missing(use.samples) && !(use.samples > nrSamples) && (use.samples > 0) ) {
    nrSamples <- use.samples
  }
  sindex <- c(1:nrSamples)
  obs2 <- cf$obs
  conf.index <- cf$conf.index
  
  ## number of gauges
  N <- dim(tcf)[3]
  ## index array for t
  i <- c(1:Time)
  ## index array for t'
  i2 <- i
  ## space for the correlator
  Cf <- array(NA, dim=c(N, Time/2+1))

  vev <- 0.
  ## compute vev first
  ## mean over all gauges and times
  if(nrSamples == 1) vev <- mean(tcf)
  else vev <- mean(tcf[,sindex,])
  if(verbose) message("vev1 = ", vev, "\n")

  if(!subtract.vev) vev <- 0.

  if(missing(cf2)) {
    ## here we compute the actual correlation
    if(nrSamples != 1) {
      ## re-order data
      mtcf <- tcf - vev
      ## average over samples, tcf has dim(Time,N)
      tcf <- apply(mtcf[,sindex,], c(1,3), sum)
    }
    else{
      subtract.equal <- FALSE
      tcf <- tcf[,1,] - vev
    }
    ## need to run only to Time/2 because source and sink are equal
    ## only possible type is cosh
    for(dt in c(0:(Time/2))) {
      Cf[,1+dt] <- apply(tcf[i,]*tcf[i2,], 2, mean)
      ## subtract product of equal samples
      if(subtract.equal) Cf[,1+dt] <- Cf[,1+dt] - apply(apply(mtcf[i,sindex,]*mtcf[i2,sindex,], c(2,3), mean), 2, sum)
      ## shift the index array by 1 to the left
      i2 <- (i2) %% Time + 1
    }
    if(subtract.equal) Cf <- Cf/nrSamples/(nrSamples-1)
    else Cf <- Cf/nrSamples/(nrSamples)
  }
  ## now the more general case case of cross-correlators
  else {
    sign <- +1
    if(type != "cosh") sign <- -1
    
    nrSamples2 <- cf2$nrSamples
    if(!missing(use.samples2) && !(use.samples2 > nrSamples2) && (use.samples2 > 0) ) {
      nrSamples2 <- use.samples2
    }

    sindex2 <- c(1:nrSamples2)
    obs2 <- cf2$obs
    ## sanity checks
    if(nrSamples != nrSamples2 && subtract.equal) {
      warning("samples numbers are not equal for cf and cf2\n Setting subtract.equal = FALSE\n")
      subtract.equal <- FALSE
    }
    if(nrSamples == 1 && nrSamples2 == 1 && subtract.equal) {
      warning("samples numbers for both cf and cf2 equal to 1\n Setting subtract.equal = FALSE\n")
      subtract.equal <- FALSE
    }
    if(cf2$Time != Time) {
      stop("time extent in two loops does not agree... Aborting...!\n")
    }
    if(!real2 && smeared2) tcf2 <- cf2$sicf
    else if(!real2) tcf2 <- cf2$icf
    else if(smeared2) tcf2 <- cf2$scf
    else tcf2 <- cf2$cf

    ## compute vev2 now
    ## mean over all gauges and times
    vev2 <- 0.
    if(nrSamples2 == 1) vev2 <- mean(tcf2)
    else vev2 <- mean(tcf2[,sindex2,])
    if(verbose) message("vev2 = ", vev2, "\n")
    if(!subtract.vev2) vev2 <- 0.

    ## now we check using conf.index whether the data sets are matched
    ## we remove any non matched entries
    if(any(!(cf$conf.index %in% cf2$conf.index))) {
      missing.ii <- which(!(cf$conf.index %in% cf2$conf.index))
      tcf <- tcf[,,-missing.ii]
      cf$conf.index <- cf$conf.index[-missing.ii]
      warning(paste("removed config", missing.ii, "from data set cf, it could not be matched\n"))
    }
    if(any(!(cf2$conf.index %in% cf$conf.index))) {
      missing2.ii <- which(!(cf2$conf.index %in% cf$conf.index))
      tcf2 <- tcf2[,,-missing2.ii]
      cf2$conf.index <- cf2$conf.index[-missing.ii]
      warning(paste("removed config", missing.ii, "from data set cf2, it could not be matched\n"))
    }
    if(dim(tcf2)[3] != dim(tcf)[3]) {
      stop("number of gauges for the two loops does not agree... Aborting...!\n")
    }
    ## the unique matched configuration number index
    N <- dim(tcf2)[3]
    conf.index <- unique(cf2$conf.index, cf$conf.index)
    Cf <- array(NA, dim=c(N, Time/2+1))
    
    ## re-order data
    ## and average over samples, tcf and tcf2 have then dim(Time,N)
    if(nrSamples != 1) {
      mtcf <- tcf - vev
      tcf <- apply(mtcf[,sindex,], c(1,3), sum)
    }
    else {
      tcf <- tcf[,1,] - vev
    }
    if(nrSamples2 != 1) {
      mtcf2 <- tcf2 - vev2
      tcf2 <- apply(mtcf2[,sindex2,], c(1,3), sum)
    }
    else {
      tcf2 <- tcf2[,1,] - vev2
    }

    ## finally we correlate
    for(dt in c(0:(Time/2))) {
      ## here we do the time average (t and Time-1) in the same step
      Cf[,1+dt] <- apply(0.5*(tcf[i,]*tcf2[i2,] + sign*tcf2[i,]*tcf[i2,]), 2, mean)
      ## subtract product of equal samples
      if(subtract.equal) Cf[,1+dt] <- Cf[,1+dt] -
        apply(apply(0.5*(mtcf[i,sindex,]*mtcf2[i2,sindex2,] + sign*mtcf2[i,sindex2,]*mtcf[i2,sindex,]),
                    c(2,3), mean), 2, sum)
      ## shift the index array by 1 to the left
      i2 <- (i2) %% Time + 1
    }
    if(nrSamples2 == nrSamples) {
      if(subtract.equal) Cf <- Cf/nrSamples/(nrSamples-1)
      else Cf <- Cf/nrSamples/nrSamples
    }
    else {
      ## subtract.equal must be FALSE here
      Cf <- Cf/nrSamples/nrSamples2
    }
  }
  ret <- cf_orig(cf=Cf)
  ret <- cf_meta(.cf=ret, Time=Time, nrObs=1, nrStypes=1, symmetrised=TRUE)
  ret$nrSamples=nrSamples
  ret$nrSamples2=nrSamples2
  ret <- addConfIndex2cf(ret, conf.index)
  return(invisible(ret))
}

