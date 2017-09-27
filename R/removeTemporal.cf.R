## will first multiply with 
removeTemporal.cf <- function(cf, single.cf1, single.cf2,
                              p1=c(0,0,0), p2=c(0,0,0), L,
                              lat.disp=TRUE, weight.cosh=FALSE) {

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
  Time <- cf$Time
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
  ## use momenta p1 and p2 and lattice or continuum dispersion relation to shift energies
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
  c <- 0.
  if(weight.cosh) c <- 1.
  ## multiply with the exponential correction factor
  Exptt <- exp((mass2$t0-mass1$t0)*c(0:(Time/2))) + c*exp((mass2$t0-mass1$t0)*(Time-c(0:(Time/2))))
  if(!is.null(cf$cf)) {
    cf$cf <- cf$cf*t(array(Exptt, dim=dim(cf$cf)[c(2,1)]))
  }
  cf$cf.tsboot$t0 <- cf$cf.tsboot$t0*Exptt
  cf$cf0 <- cf$cf.tsboot$t0
  for(i in c(1:cf$boot.R)) {
    cf$cf.tsboot$t[i,] <- cf$cf.tsboot$t[i,]*
      (exp((mass2$t[i]-mass1$t[i])*c(0:(Time/2))) + c*exp((mass2$t[i]-mass1$t[i])*(Time-c(0:(Time/2)))))
  }
  ## take the differences of C(t+1) and C(t)
  cf <- takeTimeDiff.cf(cf)

  ## multiply with the exponetial inverse
  Exptt <- exp(-(mass2$t0-mass1$t0)*c(-1:(Time/2-1))) + c*exp(-(mass2$t0-mass1$t0)*(Time-c(-1:(Time/2-1))))
  if(!is.null(cf$cf)) {
    cf$cf <- cf$cf*t(array(Exptt, dim=dim(cf$cf)[c(2,1)]))
  }
  cf$cf.tsboot$t0 <- cf$cf.tsboot$t0*Exptt
  cf$cf0 <- cf$cf.tsboot$t0
  for(i in c(1:cf$boot.R)) {
    cf$cf.tsboot$t[i,] <- cf$cf.tsboot$t[i,]*
      (exp(-(mass2$t[i]-mass1$t[i])*c(-1:(Time/2-1))) +c*exp(-(mass2$t[i]-mass1$t[i])*(Time-c(-1:(Time/2-1)))) )
  }
  ## store masses in cf
  cf$mass1 <- mass1
  cf$mass2 <- mass2
  cf$weighted <- TRUE
  cf$weight.factor <- 1.
  cf$weight.cosh <- weight.cosh
  return(invisible(cf))
}

takeTimeDiff.cf <- function(cf, deltat = 1, forwardshift= FALSE) {
  if(missing(cf)) {
    stop("takeTimeDiff: cf must be provided! Aborting...\n")
  }
  ## number of time slices (hopefully in units of T/2+1)
  T <- cf$Time
  Nt <- dim(cf$cf)[2]
  
  nts <- cf$Time/2+1                                                                                                                          
  if( "symmetrised" %in% names(cf) ) {
    if(!cf$symmetrised){
      nts <- cf$Time
    }
  }

  nrObs <- floor(Nt/nts)
  ## the time indices to be subtracted
  tt0 <- c()
  for(i in c(1:nrObs)) {
    tt0 <- c(tt0, ((i-1)*(nts)+1):(i*(nts)-deltat))
  }
  tt1 <- tt0 + deltat

  ## the default is a type of backwards derivative: C'(t) = C(t-1) - C(t)
  ## alternatively, we can also do a forward derivative: C'(t) = C(t+1) - C(t)
  tlhs <- tt1
  trhs1 <- tt0
  trhs2 <- tt1
  if( forwardshift ){
    tlhs <- tt0
    trhs1 <- tt1
    trhs2 <- tt0
  }

  ## take the differences, set the remaining points to NA
  if(!is.null(cf$cf)) {
    cf$cf[,tlhs] <- cf$cf[,trhs1]-cf$cf[,trhs2]
    cf$cf[,-tlhs] <- NA
  }
  ## now the bootstrap samples
  if(cf$boot.samples) {
    cf$cf0[tlhs] <- cf$cf0[trhs1]-cf$cf0[trhs2]
    cf$cf0[-tlhs] <- NA

    cf$cf.tsboot$t0[tlhs] <- cf$cf.tsboot$t0[trhs1]-cf$cf.tsboot$t0[trhs2]
    cf$cf.tsboot$t0[-tlhs] <- NA
    cf$cf.tsboot$t[,tlhs] <- cf$cf.tsboot$t[,trhs1]-cf$cf.tsboot$t[,trhs2]
    cf$cf.tsboot$t[,-tlhs] <- NA
  }
  ## save info
  cf$shifted <- TRUE
  cf$deltat <- deltat
  cf$forwardshift <- forwardshift
  ## return subtracted cf
  return(invisible(cf))
}
