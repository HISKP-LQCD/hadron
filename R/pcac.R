
pcacsym <- function(data, t, T2, pa=FALSE) {
  # t and T2 in R counting (starting with 1)
  # so t must be larger than 1
  # use data from t and T2-t
  if(t < 2 || t > T2/2) {
    stop("t out of range in function pcacsym")
  }
  if(pa) {
    numerator1 <- 0.5*(0.5*(data[t+1]-data[2*T2+t+1])-0.5*(data[T2+2-(t+1)]-data[3*T2+2-(t+1)]))
    numerator2 <- 0.5*(0.5*(data[t-1]-data[2*T2+t-1])-0.5*(data[T2+2-(t-1)]-data[3*T2+2-(t-1)]))
  }
  else {
    numerator1 <- 0.5*(data[t+1]-data[T2+2-(t+1)])
    numerator2 <- 0.5*(data[t-1]-data[T2+2-(t-1)])
  }
  denumerator <- 0.5*(data[T2+t]+data[T2+T2+2-t])
  mass <- (numerator1-numerator2)/(4.*denumerator)
  rm(numerator1)
  rm(numerator2)
  rm(denumerator)
  return(invisible(mass))
}

pcacfit <- function(data, from, to, T2, pa=FALSE) {
  mass <- 0.
  for(t in from:to) {
    mass = mass + pcacsym(data, t, T2, pa)
  }
  mass = mass/(t-from+1)
  return(invisible(mass))
}



pcac <- function(psfilename, apfilename, pafilename, from=3, to=3, fit=F, skip=0, plotit=F, S=1.5, debug=F) {
  if(!missing(psfilename)) {
    psscar <- read.table(file=psfilename, col.names=c("t","ps"), header=F)
  }
  else {
    stop("Error! Option filename is mandatory!")
  }
  if(!missing(apfilename)) {
    axps <- read.table(file=apfilename, col.names=c("t","ap"), header=F)
    ap <- TRUE
  }
  else ap <- FALSE

  if(!missing(pafilename)) {
    if(ap) {
      psax <- read.table(file=pafilename, col.names=c("t","pa"), header=F)
      pa <- TRUE
    }
    else {
      axps <- read.table(file=pafilename, col.names=c("t","ap"), header=F)
      axps$ap <- -axps$ap
      pa <- FALSE
      ap = TRUE
    }
  }
  else pa <-  FALSE

  if(!pa && !ap) {
    stop("Error! neither axps nor psax file specified")
  }
  
  if(length(psscar$ps) != length(axps$ap)) {
    print("Files do not belong together, they have different sizes! Aborting...")
    return(0)
  }
  if(pa) {
    if(length(psscar$ps) != length(psax$pa)) {
      print("Files do not belong together, they have different sizes! Aborting...")
      return(0)
    }
  }

  T2 <- (max(psscar$t)-min(psscar$t)+1)
  L <- (to-from+1)
  Zps <- array(psscar$ps, dim=c(T2,length(psscar$ps)/T2))
  Zap <- array(axps$ap, dim=c(T2, length(axps$ap)/T2))
  if(pa) Zpa <- array(psax$pa, dim=c(T2, length(psax$pa)/T2))

  cat("Found", length(psscar$ps)/T2, "measurements, skipping", skip, " \n")

  rm(psscar)
  rm(axps)
  if(pa) rm(psax)
  # put ps and ap correlators behind eachother for each measurement.
  if(pa) {
    Z <- rbind(Zap, Zps, Zpa)
    rm(Zpa)
  }
  else {
    Z <- rbind(Zap, Zps)
  }
  rm(Zps)
  rm(Zap)
  # prepare result file
  result <- data.frame(t = array(0.,dim=c(L)), mass = array(0.,dim=c(L)), dmass = array(0.,dim=c(L)),
                       ddmass = array(0.,dim=c(L)), tauint = array(0.,dim=c(L)), dtauint = array(0.,dim=c(L)))
#                       chisqr = array(0.,dim=c(L)))
  i <- 1

  for(t in from:to) {
    if(fit) {
      mass <- uwerrderived(pcacfit, data=Z[,skip:(length(Z[1,]))], S=S, pl=debug, t, to, T2, pa)
    }
    else {
      mass <- uwerrderived(pcacsym, data=Z[,skip:(length(Z[1,]))], S=S, pl=debug, t, T2, pa)
    }
    result$t[i] <- t
    result$mass[i] <- mass$value[1]
    result$dmass[i] <- mass$dvalue[1]
    result$ddmass[i] <- mass$ddvalue[1]
    result$tauint[i] <- mass$tauint[1]
    result$dtauint[i] <- mass$dtauint[1]
    i=i+1
  }
  
  rm(mass)
  rm(Z)
  attr(result, "class") <- c("massfit", "data.frame")
  if(plotit == T) {
    plot(result)
  }
  return(invisible(result))
  
}
