readcmicor <- function(filename) {
  data <- read.table(filename, header=F,
                     colClasse=c("integer","integer","integer","numeric","numeric","integer"))
  attr(data, "class") <- c("cmicor", "data.frame")  
  return(invisible(data))
}

readcmidatafiles <- function(path="./", basename="onlinemeas", skip=1,
                             last.digits=4, verbose=FALSE) {
  ## read one file to determine parameters
  ofiles <- dir(path=path, pattern=paste(basename, "*", sep=""))

  ## sort input files using the last 4 characters of the filename
  e <- nchar(ofiles[1])
  s <- e-last.digits+1
  ## sort-index vector
  gaugeno <- as.integer(substr(ofiles,s,e))
  ii <- order(gaugeno)

  ## read single files
  j <- 1
  cmicor <- data.frame()
  for(i in ii) {
    if(verbose) {
      cat("Reading from file", ofiles[i], "\n")
    }
    cmicor <- rbind(cmicor, read.table(paste(path, ofiles[i], sep=""), skip=skip))
    j <- j+1
  }
  attr(cmicor, "class") <- c("cmicor", "data.frame")  
  return(invisible(cmicor))
}

extract.obs <- function(cmicor, vec.obs=c(1), ind.vec=c(1,2,3,4,5), sym.vec, verbose=FALSE) {
  if(missing(cmicor)) {
    stop("data missing in extract.obs\n")
  }
  
  nrObs <- max(cmicor[,ind.vec[1]])
  if(nrObs != length(unique(cmicor[,ind.vec[1]]))) {
    stop("extract.obs: data inconsistent, nrObs not equal to input data length\n")
  }
  nrStypes <- length(unique(cmicor[,ind.vec[2]]))
  Time <-  2*max(cmicor[,ind.vec[3]])
  Thalf <- max(cmicor[,ind.vec[3]])+1
  if(Thalf != length(unique(cmicor[,ind.vec[3]]))) {
    stop("extract.obs: data inconsistent, T not equal to input data length\n")
  }
  if(verbose) cat("nrObs=",nrObs, "nrStypes=",nrStypes, "T=", Time, "\n")

  nrObs <- length(vec.obs)
  data <- cmicor[cmicor[,ind.vec[1]] %in% vec.obs,]
  data[(data[,ind.vec[3]]!=0 & (data[,ind.vec[3]]!=(Thalf-1))),ind.vec[c(4,5)]] <-
      data[(data[,ind.vec[3]]!=0 & (data[,ind.vec[3]]!=(Thalf-1))),ind.vec[c(4,5)]]/2
  ## symmetrise or anti-symmetrise for given observable?
  isym.vec <- rep(+1, times= nrObs*Thalf*nrStypes)
  if(!missing(sym.vec)) {
    if(length(sym.vec) != nrObs) {
      stop("sym.vec was given, but had not the correct length")
    }
    for(i in c(1:nrObs)) {
      if(!sym.vec) {
        isym.vec[((i-1)*Thalf*nrStypes+1):((i)*Thalf*nrStypes)] <- -1
      }
    }
  }
  cf <- array((data[,ind.vec[4]] + isym.vec*data[,ind.vec[5]]),
              dim=c(nrObs*Thalf*nrStypes, length(data[,1])/(nrObs*Thalf*nrStypes)))

  return(invisible(list(cf=cf, Time=Time, nrStypes=nrStypes, nrObs=nrObs)))
}

readhlcor <- function(filename) {
  return(invisible(read.table(filename, header=F,
                              colClasse=c("integer", "integer","integer","integer","numeric","numeric","numeric","numeric","integer"))))
}

readoutputdata <- function(filename) {
  data <- read.table(filename, header=F)
  attr(data, "class") <- c("outputdata", "data.frame")  
  return(invisible(data))
}
