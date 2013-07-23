readcmicor <- function(filename, colClasses=c("integer","integer","integer","numeric","numeric","integer"),
                       skip=0) {
  data <- read.table(filename, header=F, skip=skip,
                     colClasses=colClasses)
  attr(data, "class") <- c("cmicor", "data.frame")  
  return(invisible(data))
}

getorderedfilelist <- function(path="./", basename="onlinemeas", last.digits=4) {

  ofiles <- dir(path=path, pattern=paste(basename, "*", sep=""))
  if(any(nchar(ofiles) != nchar(ofiles[1]))) {
    stop("we need all filenames to have the same length, aborting...\n")
  }

  ## sort input files using the last last.digits characters of each filename
  e <- nchar(ofiles[1])
  s <- e-last.digits+1
  ## sort-index vector
  gaugeno <- as.integer(substr(ofiles,s,e))
  ii <- order(gaugeno)
  return(invisible(ofiles[ii]))
}

readcmidatafiles <- function(files, excludelist=c(""), skip=1, verbose=FALSE,
                             colClasses=c("integer", "integer","integer","numeric","numeric")) {
  if(missing(files)) {
    stop("filelist missing, aborting...\n")
  }
  cmicor <- data.frame()
  ## read single files
  for(f in files) {
    if( !(f %in% excludelist) && file.exists(f)) {
      if(verbose) {
        cat("Reading from file", f, "\n")
      }
      cmicor <- rbind(cmicor, read.table(f, skip=skip, header=F,
                                         colClasses=colClasses))
    }
  }
  attr(cmicor, "class") <- c("cmicor", "data.frame")  
  return(invisible(cmicor))
}

readcmiloopfiles <- function(files, excludelist=c(""), skip=0, verbose=FALSE,
                             colClasses=c("integer", "integer","integer","integer",
                               "numeric","numeric","numeric","numeric")) {

  if(missing(files)) {
    stop("filelist missing, aborting...\n")
  }
  ldata <- data.frame()
  for(f in files) {
    if( !(f %in% excludelist) && file.exists(f)) {
      if(verbose) {
        cat("Reading from file", f, "\n")
      }
      ldata <- rbind(ldata,read.table(f, colClasses=colClasses, skip=skip))
    }
  }
  attr(ldata, "class") <- c("cmiloop", "data.frame")  
  return(invisible(ldata))
}

extract.loop <- function(cmiloop, obs=9, ind.vec=c(2,3,4,5,6,7,8), L) {
  ldata <- cmiloop[cmiloop[,ind.vec[1]] == obs,] 
  T <- max(ldata[,ind.vec[2]])
  nrSamples <- max(ldata[,ind.vec[3]])
  if(missing(L)) {
    L <- T/2
  }
  cf <- list(cf = array(ldata[,ind.vec[4]], dim=c(T, nrSamples, length(files)))/sqrt(L^3),
             icf = array(ldata[,ind.vec[5]], dim=c(T, nrSamples, length(files)))/sqrt(L^3),
             scf = array(ldata[,ind.vec[6]], dim=c(T, nrSamples, length(files)))/sqrt(L^3),
             sicf= array(ldata[,ind.vec[7]], dim=c(T, nrSamples, length(files)))/sqrt(L^3),
             Time=T, nrStypes=2, nrObs=1, nrSamples=nrSamples, obs=obs)
  return(invisible(cf))
}

extract.obs <- function(cmicor, vec.obs=c(1), ind.vec=c(1,2,3,4,5),
                        sym.vec, sign.vec, verbose=FALSE) {
  if(missing(cmicor)) {
    stop("data missing in extract.obs\n")
  }
  ## consistency check for data in file
  filenrObs <- max(cmicor[,ind.vec[1]])
  if(filenrObs != length(unique(cmicor[,ind.vec[1]]))) {
    stop("extract.obs: data inconsistent, nrObs in file not matching to input data length\n")
  }
  nrStypes <- length(unique(cmicor[,ind.vec[2]]))
  Time <-  2*max(cmicor[,ind.vec[3]])
  Thalf <- max(cmicor[,ind.vec[3]])+1
  if(Thalf != length(unique(cmicor[,ind.vec[3]]))) {
    stop("extract.obs: data inconsistent, T not equal to input data length\n")
  }
  if(verbose) cat("filenrObs=",nrObs, "nrStypes=",nrStypes, "T=", Time, "\n")

  nrObs <- length(vec.obs)
  data <- cmicor[cmicor[,ind.vec[1]] %in% vec.obs,]
  data[(data[,ind.vec[3]]!=0 & (data[,ind.vec[3]]!=(Thalf-1))),ind.vec[c(4,5)]] <-
      data[(data[,ind.vec[3]]!=0 & (data[,ind.vec[3]]!=(Thalf-1))),ind.vec[c(4,5)]]/2
  ## symmetrise or anti-symmetrise for given observable?
  isym.vec <- rep(+1, times= nrObs*Thalf*nrStypes)
  isign.vec <- rep(+1, times= nrObs*Thalf*nrStypes)
  if(!missing(sym.vec)) {
    if(length(sym.vec) != nrObs) {
      stop("sym.vec was given, but had not the correct length")
    }
    for(i in c(1:nrObs)) {
      if(!sym.vec[i]) {
        isym.vec[((i-1)*Thalf*nrStypes+1):((i)*Thalf*nrStypes)] <- -1
      }
    }
  }
  if(!missing(sign.vec)) {
    if(length(sign.vec) != nrObs) {
      stop("sign.vec was given, but does not have the correct length")
    }
    for(i in c(1:nrObs)) {
      if(sign.vec[i] < 0) {
        isign.vec[((i-1)*Thalf*nrStypes+1):((i)*Thalf*nrStypes)] <- -1
      }
    }
  }
  cf <- t(array(isign.vec*(data[,ind.vec[4]] + isym.vec*data[,ind.vec[5]]),
              dim=c(nrObs*Thalf*nrStypes, length(data[,1])/(nrObs*Thalf*nrStypes))))

  return(invisible(list(cf=cf, Time=Time, nrStypes=nrStypes, nrObs=nrObs)))
}

readhlcor <- function(filename) {
  return(invisible(read.table(filename, header=F,
                              colClasses=c("integer", "integer","integer","integer","numeric","numeric","numeric","numeric","integer"))))
}

readoutputdata <- function(filename) {
  data <- read.table(filename, header=F)
  attr(data, "class") <- c("outputdata", "data.frame")  
  return(invisible(data))
}

readbinarycf <- function(files, T=48, obs=5, endian="little",
                         excludelist=c(""), sym=TRUE, path="") {

  ## indices for averaging +-t
  i1 <- c(2:(T/2))+obs*T
  i2 <- c(T:(T/2+2))+obs*T
  sign <- +1
  if(!sym) sign <- -1

  Cf <- complex()
  for(f in files) {
    ifs <- paste(path, f, sep="")
    if( !(ifs %in% excludelist) && file.exists(ifs)) {
      to.read <- file(ifs, "rb")
      tmp <- readBin(to.read, complex(), n=(obs+1)*T, endian = endian)
      ## average +-t
      tmp[i1] <- 0.5*(tmp[i1] + sign * tmp[i2])
      
      Cf <- cbind(Cf, tmp[c(1:(T/2+1))+obs*T])
      close(to.read)
    }
  }
  cf <- list(cf=t(Re(Cf)), icf=t(Im(Cf)), Time=T, nrStypes=1, nrObs=1)
  return(invisible(cf))
}

readbinarydisc <- function(files, T=48, obs=5, endian="little",
                           excludelist=c(""), nrSamples=1, path="") {
  Cf <- complex()

  N <- length(files)/nrSamples
  if(nrSamples*N != length(files)) {
    stop("not the same number of samples per gauge")
  }
  for(f in files) {
    ifs <- paste(path,f,sep="")
    if( !(ifs %in% excludelist) && file.exists(ifs)) {
      to.read <- file(ifs, "rb")
      tmp <- readBin(to.read, complex(), n=(obs+1)*T, endian = endian)
      Cf <- cbind(Cf, tmp[c((obs*T+1):((obs+1)*T))])
      close(to.read)
    }
  }
  Cf <- array(Cf, dim=c(T, nrSamples, N))
  cf <- list(cf=Re(Cf), icf=Im(Cf), scf=NULL, sicf=NULL,
             Time=T, nrStypes=1, nrObs=1, nrSamples=nrSamples, obs=obs)
  return(invisible(cf))
}

readcmidisc <- function(files, obs=9, ind.vec=c(2,3,4,5,6,7,8),
                        excludelist=c(""), skip=0, L,
                        colClasses=c("integer", "integer","integer","integer",
                          "numeric","numeric","numeric","numeric")) {
  if(missing(files)) {
    stop("filelist missing, aborting...\n")
  }
  if(length(ind.vec) != 7) {
    stop("ind.vec must have length 7, aborting...\n")
  }
  ldata <- data.frame()
  for(f in files) {
    if( !(f %in% excludelist) && file.exists(f)) {
      tmp <- read.table(f, colClasses=colClasses, skip=skip)
      ldata <- rbind(ldata, tmp[tmp[,ind.vec[1]] %in% obs, ])
      T <- max(ldata[,ind.vec[2]])
      nrSamples <- max(ldata[,ind.vec[3]])
    }
  }
  if(missing(L)) L <- T/2
  cf <- list(cf = array(ldata[,ind.vec[4]], dim=c(T, nrSamples, length(files)))/sqrt(L^3),
             icf = array(ldata[,ind.vec[5]], dim=c(T, nrSamples, length(files)))/sqrt(L^3),
             scf = array(ldata[,ind.vec[6]], dim=c(T, nrSamples, length(files)))/sqrt(L^3),
             sicf= array(ldata[,ind.vec[7]], dim=c(T, nrSamples, length(files)))/sqrt(L^3),
             Time=T, nrStypes=2, nrObs=1, nrSamples=nrSamples, obs=obs)
  return(invisible(cf))
}
