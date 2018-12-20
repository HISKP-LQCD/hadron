meanindexed <- function(data, indexvector) {
  return(invisible(mean(data[indexvector])))
}

sd.index <- function(data, indexvector) {
  return(invisible(sd(data[indexvector])))
}

bootstrap.analysis <- function(data, skip=0, boot.R=100,
                               tsboot.sim="geom", pl=FALSE, boot.l=2) {
  data <- data[skip:length(data)]
  data.mean = mean(data)
  error.naive = sd(data)/sqrt(length(data))

  cat("mean value = ", data.mean, "\n")
  cat("naive error = ", error.naive, "\n")
  
  data.boot <- boot::boot(data=data, statistic=meanindexed, R=boot.R, stype="i")
  data.boot.ci <- boot::boot.ci(data.boot, type = c("norm", "basic", "perc"))

  cat("                  mean        -err           +err            stderr        bias\n")
  cat("bootstrap      = ", data.boot$t0[1], "(", (data.boot.ci$normal[1,2]-data.boot$t0[1])/1.96
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

  cat("blocking analysis:\n")
  cat("\t\t\t mean   \t stderr \t dstderr\t tau_int\t bias\n")
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

    cat("blocklength =", boot.l, "\t", Mean[j], "\t", Error[j], "\t",
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
    if(interactive() && (grepl(pattern="X11", x=names(dev.cur()), ignore.case=TRUE) || grepl(pattern="null", x=names(dev.cur()), ignore.case=TRUE))) {
      X11()
    }
    plotwitherror(df$Blocksize, df$Error, df$DError, xlab="l", ylab="Error")
  }
  return(invisible(df))
}
