mean.index <- function(data, indexvector) {
  return(invisible(mean(data[indexvector])))
}

bootstrap.analysis <- function(data, skip=0, boot.R=100,
                               tsboot.sim="geom", pl=F) {
  
  data.mean = mean(data[skip:length(data)])
  error.naive = sd(data[skip:length(data)])/sqrt(length(data)-skip)

  cat("mean value = ", data.mean, "\n")
  cat("naive error = ", error.naive, "\n")
  
  data.boot <- boot(data=data[skip:length(data)], statistic=mean.index, R=boot.R, stype="i")
  data.boot.ci <- boot.ci(data.boot, type = c("norm", "basic", "perc"))

  cat("                  mean        -err           +err            stderr        bias\n")
  cat("bootstrap      = ", data.boot$t0[1], "(", (data.boot.ci$normal[1,2]-data.boot$t0[1])/1.96
	, ",", -(data.boot$t0[1]-data.boot.ci$normal[1,3])/1.96, ")", sd(data.boot$t[,1]),
	mean(data.boot$t[,1])-data.boot$t0[1],"\n")
  Blocksize <-  numeric()
  Mean <- numeric()
  Error <- numeric()
  DError <- numeric()
  Tauint <- numeric()
  Blocksize[1] <- 1
  Mean[1] <- data.boot$t0[1]
  Error[1] <- sd(data.boot$t[,1])
  DError[1] <- 0.
  Tauint[1] <- 0.
  cat("blocking analysis:\n")
  boot.l <- 2
  cat("\t\t\t mean   \t stderr \t dstderr\t tau_int\n")
  j <- 1
  while((length(data)-skip)/boot.l > 20) {
    j <- j+1
    data.tsboot <- tsboot(data[skip:length(data)], statistic=mean, R=boot.R, l=boot.l,
	                  sim=tsboot.sim)
    data.tsboot.ci <- boot.ci(data.tsboot, type = c("norm", "basic", "perc"))
    cat("blocklength =", boot.l, "\t",
        data.tsboot$t0[1], "\t", sd(data.tsboot$t[,1]), "\t",
        sqrt(boot.l/(length(data)-skip))*sd(data.tsboot$t[,1]),
        "\t", sd(data.tsboot$t[,1])^2/error.naive^2/2, "\n")
    Blocksize[j] <- boot.l
    Mean[j] <- data.tsboot$t0[1]
    Error[j] <- sd(data.tsboot$t[,1])
    DError[j] <- sqrt(boot.l/(length(data)-skip))*sd(data.tsboot$t[,1])
    Tauint[j] <- sd(data.tsboot$t[,1])^2/error.naive^2/2
    
    if(boot.l < 32) {
      boot.l <- boot.l*2
    }
    else {
      boot.l <- boot.l+20
    }
  }
  df <- data.frame(Blocksize=Blocksize, Mean=Mean, Error=Error, DError=DError, Tauint=Tauint)
  if(pl) {
    plot(data.boot)
  }
  return(invisible(df))
}
