## will first multiply with 
removeTemporal.cf <- function(cf, single.cf1, single.cf2, p1=c(0,0,0), p2=c(0,0,0), L, lat.disp=TRUE) {

  if(missing(cf)) {
    stop("cf must be provided to removeTemporal.cf at least\n")
  }
  if(missing(single.cf1)) {
    ## take the differences of C(t+1) and C(t)
    cf <- takeTimeDiff.cf(cf)
    return(invisible(cf))
  }

  if(!cf$boot.samples) {
    stop("please provide bootstrapped cf to removeTemporal.cf using the same configurations and seed for all!\n")
  }
  if(missing(L)) {
    L <- cf$Time/2
    warning("L was missing, set it to T/2\n")
  }
  mass1 <- list()
  mass2 <- list()
  
  if(missing(single.cf2)) {
    single.cf2 <- single.cf1
  }
  
  if(!cf$boot.samples || 
     cf$boot.R != single.cf1$boot.R || single.cf1$boot.R != single.cf2$boot.R ||
     cf$boot.l != single.cf1$boot.l || single.cf1$boot.l != single.cf2$boot.l) {
##     cf$seed != single.cf1$seed || single.cf1$seed != single.cf2$seed) {
    stop("please provide equally bootstrapped cfs to removeTemporal.cf using the same configurations and seed for all and the same boot.R and boot.l!\n")
  }
  
  if(inherits(single.cf1, "effectivemassfit")) {
    mass1$t0 <- single.cf1$opt.res$par[1]
    mass1$t <- single.cf1$massfit.tsboot[,1]
  }
  else if(inherits(single.cf1, "matrixfit")) {
    mass1$t0 <- single.cf1$opt.res$par[1]
    mass1$t <- single.cf1$opt.tsboot[1,]
  }
  else {
    stop("single.cf1 must be either of class effectivemassfit or matrixfit! Aborting...\n")
  }
  if(inherits(single.cf2, "effectivemassfit")) {
    mass2$t0 <- single.cf2$opt.res$par[1]
    mass2$t <- single.cf2$massfit.tsboot[,1]
  }
  else if(inherits(single.cf1, "matrixfit")) {
    mass2$t0 <- single.cf2$opt.res$par[1]
    mass2$t <- single.cf2$opt.tsboot[1,]
  }
  else {
    stop("single.cf2 must be either of class effectivemassfit or matrixfit! Aborting...\n")
  }
  ## use momenta p1 and p2 and lattice dispersion relation to shift energies
  if(any(p1 != 0)) {
    if(lat.disp) {
      pshift <- 2*sum(sin(pi*p1/L)^2)
      mass1$t0 <- acosh( cosh(mass1$t0) + pshift )
      mass1$t <- acosh( cosh(mass1$t) + pshift )
    }
    else {
      pshift <- sum((2*pi*p1/L)^2)
      mass1$t0 <- sqrt( mass1$t0^2 + pshift )
      mass1$t <- sqrt( mass1$t^2 + pshift )
    }
  }
  if(any(p2 != 0)) {
    if(lat.disp) {
      pshift <- 2*sum(sin(pi*p2/L)^2)
      mass2$t0 <- acosh( cosh(mass2$t0) + pshift )
      mass2$t <- acosh( cosh(mass2$t) + pshift )
    }
    else {
      pshift <- sum((2*pi*p2/L)^2)
      mass2$t0 <- sqrt( mass2$t0^2 + pshift )
      mass2$t <- sqrt( mass2$t^2 + pshift )
    }
  }
  ## multiply with the exponential correction factor
  Exptt <- exp((mass2$t0-mass1$t0)*c(0:(T/2)))
  if(!is.null(cf$cf)) {
    cf$cf <- cf$cf*t(array(Exptt, dim=dim(cf$cf)[c(2,1)]))
  }
  cf$cf.tsboot$t0 <- cf$cf.tsboot$t0*Exptt
  cf$cf0 <- cf$cf.tsboot$t0
  for(i in c(1:cf$boot.R)) {
    cf$cf.tsboot$t[i,] <- cf$cf.tsboot$t[i,]*exp((mass2$t[i]-mass1$t[i])*c(0:(T/2)))
  }
  ## take the differences of C(t+1) and C(t)
  cf <- takeTimeDiff.cf(cf)

  ## multiply with the exponetial inverse
  Exptt <- exp(-(mass2$t0-mass1$t0)*c(-1:(T/2-1)))
  if(!is.null(cf$cf)) {
    cf$cf <- cf$cf*t(array(Exptt, dim=dim(cf$cf)[c(2,1)]))
  }
  cf$cf.tsboot$t0 <- cf$cf.tsboot$t0*Exptt
  cf$cf0 <- cf$cf.tsboot$t0
  for(i in c(1:cf$boot.R)) {
    cf$cf.tsboot$t[i,] <- cf$cf.tsboot$t[i,]*exp(-(mass2$t[i]-mass1$t[i])*c(-1:(T/2-1)))
  }
  ## store masses in cf
  cf$mass1 <- mass1
  cf$mass2 <- mass2
  cf$weighted <- TRUE
  cf$weight.factor <- 1.
  return(invisible(cf))
}

takeTimeDiff.cf <- function(cf) {
  if(missing(cf)) {
    stop("takeTimeDiff: cf must be provided! Aborting...\n")
  }
  ## number of time slices (hopefully in units of T/2+1)
  T <- cf$Time
  Nt <- length(cf$cf0)
  nrObs <- floor(Nt/(T/2+1))
  ## the time indices to be subtracted
  tt0 <- c()
  for(i in c(1:nrObs)) {
    tt0 <- c(tt0, ((i-1)*(T/2+1)+1):((i-1)*(T/2+1)+T/2))
  }
  tt1 <- tt0+1

  ## take the differences, set the remaining points to NA
  cf$cf0[tt1] <- cf$cf0[tt0]-cf$cf0[tt1]
  cf$cf0[-tt1] <- NA
  if(!is.null(cf$cf)) {
    cf$cf[,tt1] <- cf$cf[,tt0]-cf$cf[,tt1]
    cf$cf[,-tt1] <- NA
  }
  ## now the bootstrap samples
  if(cf$boot.samples) {
    cf$cf.tsboot$t0[tt1] <- cf$cf.tsboot$t0[tt0]-cf$cf.tsboot$t0[tt1]
    cf$cf.tsboot$t0[-tt1] <- NA
    cf$cf.tsboot$t[,tt1] <- cf$cf.tsboot$t[,tt0]-cf$cf.tsboot$t[,tt1]
    cf$cf.tsboot$t[,-tt1] <- NA
  }
  ## save info
  cf$shifted <- TRUE
  ## return subtracted cf
  return(invisible(cf))
}
