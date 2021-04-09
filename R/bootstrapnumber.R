meanindexed <- function(data, indexvector) {
  return(invisible(mean(data[indexvector])))
}

sd.index <- function(data, indexvector) {
  return(invisible(sd(data[indexvector])))
}



#' Performs a Bootstrap with Blocking Analysis of a Timeseries
#' 
#' Performs a Bootstrap with Blocking Analysis of a Timeseries
#' 
#' the routine will compute the error, the error of the error and the
#' integrated autocorrelation time for different block size using a bootstrap
#' analysis. The blocksize is systematically increased starting from \code{1}
#' until \code{(length(data)-skip)/blocksize < 20}. Note that only data is kept
#' in exact multiples of the block length.
#' 
#' @param data a numerical vector containing the time series
#' @param skip integer value providing the warm up phase length.
#' @param boot.R number of bootstrap samples. See also \link[boot]{boot}, and
#' \link[boot]{tsboot}.
#' @param boot.l block length for blocked bootstrap.
#' @param tsboot.sim the \code{sim} parameter of \link[boot]{tsboot}.
#' @param pl logical, indicating whether or not to plot the result.
#' @return returns a data frame containing the mean value, the error
#' approximation, the estimate of the error of the error, the value of tau int
#' and the bias for all block sizes.
#' @author Carsten Urbach, \email{carsten.urbach@@liverpool.ac.uk}
#' @seealso for an alternative way to analyse such time series see
#' \code{\link{uwerr}} and \code{\link{computeacf}}
#' @keywords ts
#' @examples
#' 
#' data(plaq.sample)
#' plaq.boot <- bootstrap.analysis(plaq.sample, pl=TRUE)
#' 
#' @export bootstrap.analysis
bootstrap.analysis <- function(data, skip=0, boot.R=100,
                               tsboot.sim="geom", pl=FALSE, boot.l=2) {
  data <- data[skip:length(data)]
  data.mean = mean(data)
  error.naive = sd(data)/sqrt(length(data))

  message("mean value = ", data.mean, "\n")
  message("naive error = ", error.naive, "\n")
  
  data.boot <- boot::boot(data=data, statistic=meanindexed, R=boot.R, stype="i")
  data.boot.ci <- boot::boot.ci(data.boot, type = c("norm", "basic", "perc"))

  message("                  mean        -err           +err            stderr        bias\n")
  message("bootstrap      = ", data.boot$t0[1], "(", (data.boot.ci$normal[1,2]-data.boot$t0[1])/1.96
	, ",", -(data.boot$t0[1]-data.boot.ci$normal[1,3])/1.96, ")", sd(data.boot$t[,1]),
	mean(data.boot$t[,1])-data.boot$t0[1],"\n")
  Blocksize <-  numeric()
  Mean <- numeric()
  Error <- numeric()
  DError <- numeric()
  Tauint <- numeric()
  Bias <- numeric()
  Blocksize[1] <- 1
  Mean[1] <- data.boot$t0[1]
  Error[1] <- sd(data.boot$t[,1])
  DError[1] <- 0.
  Tauint[1] <- 0.
  Bias[1] <- 0.

  message("blocking analysis:\n")
  message("\t\t\t mean   \t stderr \t dstderr\t tau_int\t bias\n")
  j <- 1
  while((length(data))/boot.l > 20) {
    ndata <- block.ts(data, l=boot.l)
    j <- j+1
    data.tsboot <- boot::boot(ndata, statistic=meanindexed, R=boot.R)
    ## use the same seed ...
    set.seed(data.tsboot$seed)
    data.sdboot <- boot::boot(ndata, statistic=sd.index, R=boot.R)
    ##data.tsboot.ci <- boot.ci(data.tsboot, type = c("norm", "basic", "perc"))
    Blocksize[j] <- boot.l
    Mean[j] <- data.tsboot$t0[1]
    Error[j] <- sd(data.tsboot$t[,1])
    DError[j] <- sd(data.sdboot$t[,1])/sqrt(length(ndata))
    Tauint[j] <- sd(data.tsboot$t[,1])^2/error.naive^2/2
    Bias[j] <- data.tsboot$t0[1] - mean(data.tsboot$t[,1])

    message("blocklength =", boot.l, "\t", Mean[j], "\t", Error[j], "\t",
            DError[j], "\t", Error[j]^2/error.naive^2/2, "\t", Bias[j], "\n")
    if(boot.l < 32) {
      boot.l <- boot.l*2
    }
    else {
      boot.l <- boot.l+20
    }
  }
  df <- data.frame(Blocksize=Blocksize, Mean=Mean, Error=Error, DError=DError, Tauint=Tauint, Bias=Bias)
  if(pl) {
    plot(data.boot)
    new_window_if_appropriate()
    plotwitherror(df$Blocksize, df$Error, df$DError, xlab="l", ylab="Error")
  }
  return(invisible(df))
}
