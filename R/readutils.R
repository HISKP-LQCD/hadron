readcmicor <- function(filename, colClasses=c("integer","integer","integer","numeric","numeric","integer"),
                       skip=0) {
  data <- read.table(filename, header=F, skip=skip,
                     colClasses=colClasses)
  attr(data, "class") <- c("cmicor", class(data))  
  return(invisible(data))
}

getorderedfilelist <- function(path="./", basename="onlinemeas", last.digits=4) {
  ofiles <- Sys.glob( sprintf( "%s/%s*", path, basename ) ) 
  ii <- getorderedconfignumbers(path=path, basename=basename, last.digits=last.digits)
  return(invisible(ofiles[ii]))
}

getorderedconfignumbers <- function(path="./", basename="onlinemeas", last.digits=4) {
  ofiles <- Sys.glob( sprintf( "%s/%s*", path, basename ) ) 
  if(any(nchar(ofiles) != nchar(ofiles[1]))) {
    stop("getconfigurationnumbers: all filenames need to have the same length, aborting...\n")
  }

  ## sort input files using the last last.digits characters of each filename
  e <- nchar(ofiles[1])
  s <- e-last.digits+1
  ## sort-index vector
  gaugeno <- as.integer(substr(ofiles,s,e))
  return(invisible(order(gaugeno)))
}

readcmifiles <- function(files, excludelist=c(""), skip, verbose=FALSE,
                         colClasses, obs=NULL, obs.index) {
  if(missing(files)) {
    stop("readcmifiles: filelist missing, aborting...\n")
  }
  reduce <- FALSE
  if(!(is.null(obs) || missing(obs.index))) {
    reduce <- TRUE
  }
  tmpdata <- read.table(files[1], colClasses=colClasses, skip=skip)
  if(reduce) {
    tmpdata <- tmpdata[tmpdata[,obs.index] == obs,]
  }
  fLength <- length(tmpdata$V1)
  nFiles <- length(files)
  nCols <- length(tmpdata)
  ## we generate the full size data.frame first
  tmpdata[,] <- NA
  ldata <- tmpdata
  ldata[((nFiles-1)*fLength+1):(nFiles*fLength),] <- tmpdata
  # create a progress bar
  
  pb <- NULL
  if(!verbose) {
    pb <- txtProgressBar(min = 1, max = nFiles, style = 3)
  }
  for(i in c(1:nFiles)) {
    # update progress bar
    if(!verbose) {
      setTxtProgressBar(pb, i)
    }
    if( !(files[i] %in% excludelist) && file.exists(files[i])) {
      if(verbose) {
        cat("Reading from file", files[i], "\n")
      }
      ## read the data
      tmpdata <- read.table(files[i], colClasses=colClasses, skip=skip)
      if(reduce) {
        tmpdata <- tmpdata[tmpdata[,obs.index] == obs,]
      }
      ## sanity checks
      if(fLength < length(tmpdata$V1)) {
        warning(paste("readcmifiles: file ", files[i], " does not have the same length as the other files. We will cut and hope...\n"))
      }
      else if(fLength > length(tmpdata$V1)) {
        stop(paste("readcmifiles: file", files[i], "is too short. Aborting...\n"))
      }
      if(nCols != length(tmpdata)) {
        stop(paste("readcmifiles: file", files[i], "does not have the same number of columns as the other files. Aborting...\n"))
      }
      
      ldata[((i-1)*fLength+1):(i*fLength),] <- tmpdata
      rm(tmpdata)
    }
    else if(!file.exists(files[i])) {
      warning(paste("readcmifiles: dropped file", files[i], "because its missing\n"))
    }
  }
  if(!verbose) {
    close(pb)
  }
  ## remove NAs from missing files
  ldata <- na.omit(ldata)
  return(invisible(ldata))
}

readcmidatafiles <- function(files, excludelist=c(""), skip=1, verbose=FALSE,
                             colClasses=c("integer", "integer","integer","numeric","numeric"),
                             obs=NULL, obs.index=1) {

  data <- readcmifiles(files, excludelist=excludelist, skip=skip, verbose=verbose,
                       colClasses=colClasses, obs=obs, obs.index=obs.index)
  attr(data, "class") <- c("cmicor", class(data))  
  return(invisible(data))
}

readcmiloopfiles <- function(files, excludelist=c(""), skip=0, verbose=FALSE,
                             colClasses=c("integer", "integer","integer","integer",
                               "numeric","numeric","numeric","numeric"),
                             obs=NULL, obs.index=2) {
  data <- readcmifiles(files, excludelist=excludelist, skip=skip, verbose=verbose,
                       colClasses=colClasses, obs=obs, obs.index=obs.index)
  attr(data, "class") <- c("cmiloop", class(data))
  return(invisible(data))
}

extract.loop <- function(cmiloop, obs=9, ind.vec=c(2,3,4,5,6,7,8,1), L) {
  ldata <- cmiloop[cmiloop[,ind.vec[1]] == obs,] 
  T <- max(ldata[,ind.vec[2]])
  nrSamples <- max(ldata[,ind.vec[3]])
  if(missing(L)) {
    L <- T/2
  }
  cf <- list(cf = array(ldata[,ind.vec[4]], dim=c(T, nrSamples, length(ldata[,ind.vec[4]])/T/nrSamples))/sqrt(L^3),
             icf = array(ldata[,ind.vec[5]], dim=c(T, nrSamples, length(ldata[,ind.vec[5]])/T/nrSamples))/sqrt(L^3),
             scf = array(ldata[,ind.vec[6]], dim=c(T, nrSamples, length(ldata[,ind.vec[6]])/T/nrSamples))/sqrt(L^3),
             sicf= array(ldata[,ind.vec[7]], dim=c(T, nrSamples, length(ldata[,ind.vec[7]])/T/nrSamples))/sqrt(L^3),
             Time=T, nrStypes=2, nrObs=1, nrSamples=nrSamples, obs=obs, conf.index=unique(ldata[,ind.vec[8]]))
  return(invisible(cf))
}

extract.obs <- function(cmicor, vec.obs=c(1), ind.vec=c(1,2,3,4,5),
                        sym.vec, sign.vec, verbose=FALSE) {
  if(missing(cmicor)) {
    stop("extract.obs: data missing in extract.obs\n")
  }
  ## consistency check for data in the data
  filenrObs <- max(cmicor[,ind.vec[1]])
  if(filenrObs != length(unique(cmicor[,ind.vec[1]]))) {
    stop("extract.obs: data inconsistent, nrObs in the data does not match to input data length\n")
  }
  nrStypes <- length(unique(cmicor[,ind.vec[2]]))

  Time <-  2*max(cmicor[,ind.vec[3]])
  Thalf <- max(cmicor[,ind.vec[3]])+1
  if(Thalf != length(unique(cmicor[,ind.vec[3]]))) {
    stop("extract.obs: data inconsistent, T not equal to what was found in the input data\n")
  }
  if(verbose) cat("extract.obs: filenrObs=",nrObs, "nrStypes=",nrStypes, "T=", Time, "\n")

  nrObs <- length(vec.obs)
  data <- cmicor[cmicor[,ind.vec[1]] %in% vec.obs,]
  data[(data[,ind.vec[3]]!=0 & (data[,ind.vec[3]]!=(Thalf-1))),ind.vec[c(4,5)]] <-
      data[(data[,ind.vec[3]]!=0 & (data[,ind.vec[3]]!=(Thalf-1))),ind.vec[c(4,5)]]/2
  ## symmetrise or anti-symmetrise for given observable?
  isym.vec <- rep(+1, times= nrObs*Thalf*nrStypes)
  isign.vec <- rep(+1, times= nrObs*Thalf*nrStypes)
  if(!missing(sym.vec)) {
    if(length(sym.vec) != nrObs) {
      stop("extract.obs: sym.vec was given, but had not the correct length")
    }
    for(i in c(1:nrObs)) {
      if(!sym.vec[i]) {
        isym.vec[((i-1)*Thalf*nrStypes+1):((i)*Thalf*nrStypes)] <- -1
      }
    }
  }
  if(!missing(sign.vec)) {
    if(length(sign.vec) != nrObs) {
      stop("extract.obs: sign.vec was given, but does not have the correct length")
    }
    for(i in c(1:nrObs)) {
      if(sign.vec[i] < 0) {
        isign.vec[((i-1)*Thalf*nrStypes+1):((i)*Thalf*nrStypes)] <- -1
      }
    }
  }

  cf <- t(array(isign.vec*(data[,ind.vec[4]] + isym.vec*data[,ind.vec[5]]),
              dim=c(nrObs*Thalf*nrStypes, length(data[,1])/(nrObs*Thalf*nrStypes))))
  ret <- list(cf=cf, icf=NULL, Time=Time, nrStypes=nrStypes, nrObs=nrObs, boot.samples=FALSE)
  attr(ret, "class") <- c("cf", class(ret))
  return(invisible(ret))
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

readtextcf <- function(file, T=48, sym=TRUE, path="", skip=1, check.t=0, ind.vector=c(2,3)) {
  if(missing(file)) {
    stop("files must be given! Aborting...\n")
  }
  if(T < 1) {
    stop("T must be larger than 0 and integer, aborting...\n")
  }

  tmp <- read.table(paste(path, file, sep=""), skip=skip)
  if((length(ind.vector) < length(tmp)-1) || (max(ind.vector) > length(tmp))){
    stop("index vector too long or out of range\n")
  }
  
  if(check.t > 0 && max(tmp[[check.t]]) != T-1) {
    stop("T in function call does not match the one in the file, aborting...\n")
  }
  
  i1 <- c(2:(T/2))
  i2 <- c(T:(T/2+2))
  ii <- c(1:(T/2+1))
  sign <- +1
  if(!sym) sign <- -1

  tmp <- array(tmp[[ind.vector[1]]] + 1i*tmp[[ind.vector[2]]], dim=c(T, length(tmp[[ind.vector[1]]])/T))
  ## average +-t
  tmp[i1,] <- 0.5*(tmp[i1,] + sign * tmp[i2,])

  ret <- list(cf=t(Re(tmp[ii,])), icf=t(Im(tmp[ii,])), Time=T, nrStypes=1, nrObs=1, boot.samples=FALSE)
  attr(ret, "class") <- c("cf", class(ret))
  return(invisible(ret))
}

readbinarycf <- function(files, T=48, obs=5, Nop=1, endian="little",
                         op="aver", excludelist=c(""), sym=TRUE, path="") {

  if(missing(files)) {
    stop("files must be given! Aborting...\n")
  }
  if(Nop < 1) {
    stop("Nop must be larger than 0 and integer, aborting...\n")
  }
  if(obs < 0) {
    stop("obs must be a positive integer, aborting...\n")
  }
  if(T < 1) {
    stop("T must be larger than 0 and integer, aborting...\n")
  }
  if(endian != "little" && endian != "big") {
    stop("endian must be either little or big, abroting...\n")
  }
  
  ## indices for averaging +-t
  i1 <- c(2:(T/2))
  i2 <- c(T:(T/2+2))
  ii <- c(1:(Nop*T))+obs*T
  sign <- +1
  if(!sym) sign <- -1

  Cf <- complex()
  for(f in files) {
    ifs <- paste(path, f, sep="")
    if( !(ifs %in% excludelist) && file.exists(ifs)) {
      to.read <- file(ifs, "rb")
      tmp <- readBin(to.read, what=complex(), n=(obs+Nop)*T, endian = endian)[ii]

      ## we average the replica
      if(Nop > 1) {
        if(op == "aver") {
          tmp <- apply(array(tmp, dim=c(T, Nop)), 1, mean)
        }
        else {
          tmp <- apply(array(tmp, dim=c(T, Nop)), 1, sum)
        }
      }

      ## average +-t
      tmp[i1] <- 0.5*(tmp[i1] + sign * tmp[i2])

      Cf <- cbind(Cf, tmp[c(1:(T/2+1))])
      close(to.read)
    }
    else if(!file.exists(ifs)) {
      cat("file ", ifs, "does not exist...\n")
    }
  }
  ret <- list(cf=t(Re(Cf)), icf=t(Im(Cf)), Time=T, nrStypes=1, nrObs=1, boot.samples=FALSE)
  attr(ret, "class") <- c("cf", class(ret))
  return(invisible(ret))
}

readbinarydisc <- function(files, T=48, obs=5, endian="little",
                           excludelist=c(""), nrSamples=1, path="") {
  Cf <- complex()

  N <- length(files)/nrSamples
  if(nrSamples*N != length(files)) {
    stop("readbinarydisc: not the same number of samples per gauge")
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
                          "numeric","numeric","numeric","numeric"),
                        debug = FALSE) {
  if(missing(files)) {
    stop("readcmidisc: filelist missing, aborting...\n")
  }
  if(length(ind.vec) != 7) {
    stop("readcmidisc: ind.vec must have length 7, aborting...\n")
  }
  nFiles<-0
  for(f in files) {
    if( !(f %in% excludelist) && file.exists(f)) {
      nFiles<-nFiles+1
      lastFile<-f
    }
  }
  tmp <- read.table(lastFile, colClasses=colClasses, skip=skip)
  tmp <- tmp[tmp[,ind.vec[1]] %in% obs, ]
  nRows <- nrow(tmp)
  nCols <- ncol(tmp)
  ldata <- matrix(0,nFiles*nRows, nCols)

  n <- 0
  for(f in files) {
    if( !(f %in% excludelist) && file.exists(f)) {
      if(debug) print(f)
      tmp <- read.table(f, colClasses=colClasses, skip=skip)

      ldata[c((n*nRows+1):((n+1)*nRows)),] <- as.matrix(tmp[tmp[,ind.vec[1]] %in% obs, ])
      n <- n+1
    }
  }
  T <- max(ldata[,ind.vec[2]])
  nrSamples <- max(ldata[,ind.vec[3]])

  if(missing(L)) L <- T/2
  cf <- list(cf = array(ldata[,ind.vec[4]], dim=c(T, nrSamples, nFiles))/sqrt(L^3),
             icf = array(ldata[,ind.vec[5]], dim=c(T, nrSamples, nFiles))/sqrt(L^3),
             scf = array(ldata[,ind.vec[6]], dim=c(T, nrSamples, nFiles))/sqrt(L^3),
             sicf= array(ldata[,ind.vec[7]], dim=c(T, nrSamples, nFiles))/sqrt(L^3),
             Time=T, nrStypes=2, nrObs=1, nrSamples=nrSamples, obs=obs)
  return(invisible(cf))
}
