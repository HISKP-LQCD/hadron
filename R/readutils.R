readcmicor <- function(filename) {
  data <- read.table(filename, header=F,
                     colClasse=c("integer","integer","integer","numeric","numeric","integer"))
  attr(data, "class") <- c("cmicor", "data.frame")  
  return(invisible(data))
}

readcmidatafiles <- function(path=".", basename="onlinemeas", skip=1, ind.vec=c(1,2,3,4,5),
                             last.digits=4, verbose=FALSE, sym.vec) {
  ## read one file to determine parameters
  ofiles <- dir(path=path, pattern=paste(basename, "*", sep=""))
  tmp <- read.table(paste(path, ofiles[1], sep=""), skip=skip)
  nrObs <- max(tmp[,ind.vec[1]])
  if(nrObs != length(unique(tmp[,ind.vec[1]]))) {
    stop("data inconsistent, nrObs not equal to input data length\n")
  }
  nrStypes <- length(unique(tmp[,ind.vec[2]]))
  Time <-  2*max(tmp[,ind.vec[3]])
  Thalf <- max(tmp[,ind.vec[3]])+1
  if(Thalf != length(unique(tmp[,ind.vec[3]]))) {
    stop("data inconsistent, T not equal to input data length\n")
  }
  if(verbose) cat("nrObs=",nrObs, "nrStypes=",nrStypes, "T=", Time, "#files=", length(ofiles), "\n")

  ## memory for correlators
  cf <- array(0, dim=c(length(ofiles), nrObs*Thalf*nrStypes))

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

  ## sort input files using the last 4 characters of the filename
  e <- nchar(ofiles[1])
  s <- e-last.digits+1
  ## sort-index vector
  gaugeno <- as.integer(substr(ofiles,s,e))
  ii <- order(gaugeno)

  ## read single files
  j <- 1
  for(i in ii) {
    if(verbose) {
      cat("Reading from file", ofiles[i], "\n")
    }
    tmp <- read.table(paste(path, ofiles[i], sep=""), skip=skip)
    ## take care of the zeroth and last times, which don't need to be averaged
    tmp[(tmp[,ind.vec[3]]!=0 & (tmp[,ind.vec[3]]!=(Thalf-1))),] <-
      tmp[(tmp[,ind.vec[3]]!=0 & (tmp[,ind.vec[3]]!=(Thalf-1))),]/2
    cf[j,] <- (tmp[,ind.vec[4]] + isym.vec*tmp[,ind.vec[5]])
    j <- j+1
  }
  return(invisible(list(cf=cf, gaugeno=gaugeno[ii], Time=Time, nrStypes=nrStypes, nrObs=nrObs)))
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
