computeacf <- function(tdata, W.max, Lambda=100) {
  N <- length(tdata)
  Gamma.tmp <- rep(0, times=2*W.max+W.max+1)
  dGamma <- rep(0., times=W.max+1)
  ## for t > W.max we set Gamma to 0 for simplicity
  Gamma.tmp[1:(W.max+1)] <- acf(tdata, lag.max = W.max, plot=FALSE)$acf

  ## now we determine the error using (E.11) from Luescher, hep-lat/0409106
  for(t in 0:W.max) {
    k <- c(max(1,(t-Lambda)):(t+Lambda))
    dGamma[t+1] <- sum((Gamma.tmp[(k+t+1)]+Gamma.tmp[(abs(k-t)+1)]-2*Gamma.tmp[t+1]*Gamma.tmp[(k+1)])^2);
    dGamma[t+1] <- sqrt(dGamma[t+1]/N)
  }
  Gamma <- Gamma.tmp[1:(W.max+1)]

  ## and we determine where to stop summing Gamma -> W
  ## for a more sophisticated approach see Wolff, hep-lat/0306017
  W <- 0
  for(t in 0:W.max) {
    W <- t
    if(Gamma[t+1] < dGamma[t+1]) break
  }

  ## compute integrated autocorrelation time
  tau <- 0.5 + sum(Gamma[2:(W+1)])
  ## Madras, Sokal approximation for the error, (E.14) in Luescher, hep-lat/0409106
  dtau <- sqrt((4*W + 2) * tau^2 / N)
  
  res <- list(lags = c(0:(W.max)), Gamma=Gamma, dGamma=dGamma,
              W.max=W.max, W=W, tdata=tdata, tau=tau, dtau=dtau)
  attr(res, "class") <- c("hadronacf", "list")
  return(invisible(res))
}

## generic function to plot an object of class "myGamma"
plot.hadronacf <- function (Gamma, col = "black", ...)
{
  ## this is to avoid a warning
  Gamma$dGamma[1] <- 0.001
  ## determine ylimits from data
  ylim <- c(min(Gamma$Gamma - 2 * Gamma$dGamma, na.rm = TRUE), 
            max(Gamma$Gamma + 2 * Gamma$dGamma, na.rm = TRUE))
  ## data points
  plot(Gamma$lags, Gamma$Gamma, ylim = ylim,
       col = col, xlab="t", ylab="rho[t]", ...)
  ## errors
  arrows(Gamma$lags, Gamma$Gamma - Gamma$dGamma, Gamma$lags, Gamma$Gamma + Gamma$dGamma,
         length = 0.01, angle = 90, code = 3, col = col)
  # cuttoff and baseline
  abline(h=0, col="red")
  abline(v=Gamma$W, col="red")
}

## generic summary function for an object of class "myGamma"
summary.hadronacf <- function(Gamma) {
  cat("Analysis based on Autocorrelation function\n")
  cat("cut-off parameter W:\t", Gamma$W, "\n")
  cat("tauint:\t\t\t", Gamma$tau, "\n")
  cat("dtauint:\t\t", Gamma$dtau, "\n")
  cat("data mean:\t\t", mean(Gamma$tdata), "\n")
  cat("data error (naive):\t", sqrt(var(Gamma$tdata)/length(Gamma$tdata)), "\n")
  cat("data error (corrected):\t", sqrt(2*Gamma$tau*var(Gamma$tdata)/length(Gamma$tdata)), "\n")
}

