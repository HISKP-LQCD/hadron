bootstrap.cf <- function(cf, boot.R=400, boot.l=2, seed=1234) {
  if(!any(class(cf) == "cf")) {
    stop("bootstrap.cf requires an object of class cf as input! Aborting!\n")
  }
  cf$boot.samples <- TRUE
  cf$boot.R <- boot.R
  cf$boot.l <- boot.l
  cf$seed <- seed
  cf$cf0 <- apply(cf$cf, 2, mean)
  ## we set the seed for reproducability and correlation
  set.seed(seed)
  ## now we bootstrap the correlators
  cf$cf.tsboot <- tsboot(cf$cf, statistic = function(x){ return(apply(x,2,mean))},
                         R = boot.R, l=boot.l, sim="geom")
  return(invisible(cf))
}

plot.cf <- function(cf, boot.R=400, boot.l=2, ...) {
  if(!cf$boot.samples) {
    cf <- bootstrap.cf(cf, boot.R, boot.l)
  }
  Err <- apply(cf$cf.tsboot$t, 2, sd)
  plotwitherror(rep(c(0:(cf$Time/2)), times=cf$nrStypes*cf$nrObs), cf$cf0, Err, ...)
}
