readcmicor <- function(filename, colClasses=c("integer","integer","integer","numeric","numeric","integer"),
                       skip=0) {
  data <- read.table(filename, header=F, skip=skip,
                     colClasses=colClasses)
  attr(data, "class") <- c("cmicor", class(data))  
  return(invisible(data))
}

getorderedfilelist <- function(path="./", basename="onlinemeas", last.digits=4, ending="") {
  ofiles <- Sys.glob( sprintf( "%s/%s*%s", path, basename,ending ) ) 
  ii <- getorderedconfigindices(path=path, basename=basename, last.digits=last.digits, ending=ending)
  return(invisible(ofiles[ii]))
}

getconfignumbers <- function(ofiles, basename="onlinemeas", last.digits=4, ending="") {
  if(any(nchar(ofiles) != nchar(ofiles[1]))) {
    stop("getconfigurationnumbers: all filenames need to have the same length, aborting...\n")
  }
  lending <- nchar(ending)
  
  ## sort input files using the last last.digits characters of each filename
  e <- nchar(ofiles[1]) - lending
  s <- e-last.digits+1
  ## sorted config numbers
  return(invisible(as.integer(substr(ofiles,s,e))))
}

getorderedconfigindices <- function(path="./", basename="onlinemeas", last.digits=4, ending="") {
  ofiles <- Sys.glob( sprintf( "%s/%s*%s", path, basename, ending ) ) 
  if(any(nchar(ofiles) != nchar(ofiles[1]))) {
    stop("getconfigurationnumbers: all filenames need to have the same length, aborting...\n")
  }
  lending <- nchar(ending)
  
  ## sort input files using the last last.digits characters of each filename
  e <- nchar(ofiles[1]) - lending
  s <- e-last.digits+1
  ## sort-index vector
  gaugeno <- as.integer(substr(ofiles,s,e))
  return(invisible(order(gaugeno)))
}

getorderedconfignumbers <- function(path="./", basename="onlinemeas", last.digits=4, ending="") {
  ofiles <- Sys.glob( sprintf( "%s/%s*%s", path, basename, ending ) ) 
  if(any(nchar(ofiles) != nchar(ofiles[1]))) {
    stop("getconfigurationnumbers: all filenames need to have the same length, aborting...\n")
  }
  lending <- nchar(ending)
  
  ## sort input files using the last last.digits characters of each filename
  e <- nchar(ofiles[1]) - lending
  s <- e-last.digits+1
  ## sort-index vector
  gaugeno <- as.integer(substr(ofiles,s,e))
  return(invisible(sort(gaugeno)))
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
    tmpdata <- tmpdata[tmpdata[,obs.index] %in% obs,]
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
        tmpdata <- tmpdata[tmpdata[,obs.index] %in% obs,]
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
                        sym.vec, sign.vec, verbose=FALSE, symmetrise=TRUE) {
  if(missing(cmicor)) {
    stop("extract.obs: data missing in extract.obs\n")
  }
  ## consistency check for data in the data
  if( !all( unique(vec.obs) %in% unique(cmicor[,ind.vec[1]]) )) {
    stop("extract.obs: vec.obs does not match or is not fully included in the observable list in the data\n")
  }
  nrObs <- length(vec.obs)
  nrStypes <- length(unique(cmicor[,ind.vec[2]]))
  Time <-  2*max(cmicor[,ind.vec[3]])
  Thalf <- max(cmicor[,ind.vec[3]])+1
  if(Thalf != length(unique(cmicor[,ind.vec[3]]))) {
    stop("extract.obs: data inconsistent, T not equal to what was found in the input data\n")
  }
  if(verbose) cat("extract.obs: nrObs=",nrObs, "nrStypes=",nrStypes, "T=", Time, "\n")

  data <- cmicor[cmicor[,ind.vec[1]] %in% vec.obs,]
  cf <- NULL
  
  if(symmetrise){
    ## we divide everything by 2 apart from t=0 and t=T/2
    data[(data[,ind.vec[3]]!=0 & (data[,ind.vec[3]]!=(Thalf-1))),ind.vec[c(4,5)]] <-
        data[(data[,ind.vec[3]]!=0 & (data[,ind.vec[3]]!=(Thalf-1))),ind.vec[c(4,5)]]/2
    ## symmetrise or anti-symmetrise for given observable?
    isym.vec <- rep(+1, times= nrObs*Thalf*nrStypes)
    isign.vec <- rep(+1, times= nrObs*Thalf*nrStypes)
    if(!missing(sym.vec)) {
      if(length(sym.vec) != nrObs) {
        stop("extract.obs: sym.vec was given, but does not have the correct length")
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
  }else{ # no symmetrisation
    cf <- t(array(0,dim=c( nrObs*Time*nrStypes, length(data[,1])/(nrObs*Thalf*nrStypes) ) ) )
    for(bw in c(0:1)){
      # select forward or backward correlator
      tmp <-  t(array(data=data[,ind.vec[4+bw]],dim=c(nrObs*Thalf*nrStypes, length(data[,1])/(nrObs*Thalf*nrStypes))))
      for(i in c(1:nrObs)){
        for(s in c(1:nrStypes)){
          # since we are not symmetrising, the individual correlators have T entries 
          lhs <- c((bw+1):(Thalf-bw)) + (i-1)*nrStypes*Time  + (s-1)*Time  + bw*(Thalf-1)
          # in the cmi format, the backward correlator is on time-slices 1 to T/2-1 (indices 2 to T/2)
          rhs <- c((bw+1):(Thalf-bw)) + (i-1)*nrStypes*Thalf + (s-1)*Thalf
          # but in "reverse" order (to ease averaging)
          if( bw == 1 ) rhs <- rev(rhs)
          cf[,lhs] <- tmp[,rhs]
        }
      }
    }
  }
  ret <- list(cf=cf, icf=NULL, Time=Time, nrStypes=nrStypes, nrObs=nrObs, boot.samples=FALSE, jackknife.samples=FALSE,
              symmetrised=symmetrise)
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

readtextcf <- function(file, T=48, sym=TRUE, path="", skip=1, check.t=0, ind.vector=c(2,3), symmetrise=TRUE,
                       sparsity=1, avg=1, Nmin=4, autotruncate=TRUE) {
  if(missing(file)) {
    stop("files must be given! Aborting...\n")
  }
  if(T < 1) {
    stop("T must be larger than 0 and integer, aborting...\n")
  }

  tmp <- read.table(paste(path, file, sep=""), skip=skip)
  if((length(ind.vector) < 2) || (max(ind.vector) > length(tmp)) || (min(ind.vector) < 1)){
    stop("index vector too short or out of range\n")
  }
  
  if(check.t > 0 && max(tmp[[check.t]]) != T-1) {
    stop("T in function call does not match the one in the file, aborting...\n")
  }

  if(length(tmp[[ind.vector[1]]]) %% T != 0) {
    stop("T does not devide the number of rows in file, aborting... check value of paramter skip to readtextcf!\n")
  }
  
  i1 <- c(2:(T/2))
  i2 <- c(T:(T/2+2))
  ii <- c(1:(T/2+1))
  sign <- +1
  if(!sym) sign <- -1

  tmp <- array(tmp[[ind.vector[1]]] + 1i*tmp[[ind.vector[2]]], dim=c(T, length(tmp[[ind.vector[1]]])/T))
  if( (sparsity > 1 | avg > 1) & (ncol(tmp) %% (sparsity*avg) != 0) ){
    if(autotruncate){
      cat(sprintf("sparsity=%d, avg=%d, ncol=%d\n",sparsity,avg,ncol(tmp)))
      cat("readtextcf: Sparsification and/or averaging requested, but their product does not divide the number of measurements!\n")
      cat("readtextcf: Reducing the number of total measurements to fit!\n")
      nmeas <- as.integer( (sparsity*avg)*floor( ncol(tmp)/(sparsity*avg) ))
      if( nmeas/(sparsity*avg) >= Nmin ){
        tmp <- tmp[,1:nmeas]
      } else {
        cat(sprintf("readtextcf: After sparsification and averaging, less than %d measurements remain, disabling sparsification and averaging!\n",Nmin))
        sparsity <- 1
        avg <- 1
      }
    } else {
      stop("readtextcf: Sparsification and/or averaging requested, but their product does not divide the number of measurements!\n")
    }
  }

  ## sparsify data
  if(sparsity > 1){
    sp.idx <- seq(from=1,to=ncol(tmp),by=sparsity)
    tmp <- tmp[,sp.idx]
  } 
  # average over 'avg' measurements sequentially
  if(avg > 1){
    tmp2 <- tmp
    tmp <- array(0, dim=c(T,ncol(tmp2)/avg))
    for( i in c(1:ncol(tmp)) ){
      tmp[,i] <- (1.0/avg)*apply(X=tmp2[,((i-1)*avg+1):(i*avg)],
                                 MARGIN=1,
                                 FUN=sum)
    }
  }
  
  ## average +-t
  if(symmetrise) {
    tmp[i1,] <- 0.5*(tmp[i1,] + sign * tmp[i2,])
  }else{
    ii <- c(1:T)
  }


  ret <- list(cf=t(Re(tmp[ii,])), icf=t(Im(tmp[ii,])), Time=T, nrStypes=1, nrObs=1, boot.samples=FALSE, jackknife.samples=FALSE,
              symmetrised=symmetrise)
  attr(ret, "class") <- c("cf", class(ret))
  return(invisible(ret))
}

readbinarycf <- function(files, T=48, obs=5, Nop=1, endian="little",
                         op="aver", excludelist=c(""), sym=TRUE, path="",
                         hdf5format=FALSE, hdf5name, hdf5index=c(1,2)) {

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
    stop("endian must be either little or big, aborting...\n")
  }
  if(hdf5format) {
    if(missing(hdf5name)) stop("hdf5name must be given, aborting...\n")
    h5avail <- require(rhdf5)
    if(!h5avail) stop("rhdf5 package not installed, aborting...\n")
    obs <- 1
    Nop <- 1
    if(length(hdf5index)<2) hdf5index <- c(hdf5index, hdf5index)
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
      tmp <- numeric()
      if(!hdf5format) {
        to.read <- file(ifs, "rb")
        tmp <- readBin(to.read, what=complex(), n=(obs+Nop)*T, endian = endian)[ii]
      }
      else {
        LS <- h5ls(ifs)
        if(as.integer(LS[LS$name == hdf5name,]$dim) != T) {
          stop(paste(hdf5name, "in file", ifs, "has not the correct length, aborting...\n"))
        }
        tmp <- h5read(file=ifs, name=hdf5name)
        tmp <- tmp[,hdf5index[1]] + 1i*tmp[,hdf5index[2]]
      }
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

      if(!hdf5format) {
        close(to.read)
      }
      else {
        if(exists("H5close")) H5close()
      }
    }
    else if(!file.exists(ifs)) {
      cat("file ", ifs, "does not exist...\n")
    }
  }
  ret <- list(cf=t(Re(Cf)), icf=t(Im(Cf)), Time=T, nrStypes=1, nrObs=1, boot.samples=FALSE, jackknife.samples=FALSE)
  attr(ret, "class") <- c("cf", class(ret))
  return(invisible(ret))
}


# read binary correlation functions sample by sample, return as a list of length "nosamples"
# where increasing indices refer to averaging over increasing numbers of samples
readbinarysamples <- function(files, T=48, nosamples=2, endian="little",
                              op="aver", excludelist=c(""), sym=TRUE, path="", ftype=double() ){

  if(missing(files)) {
    stop("files must be given! Aborting...\n")
  }
  if(T < 1) {
    stop("T must be larger than 0 and integer, aborting...\n")
  }
  if(endian != "little" && endian != "big") {
    stop("endian must be either little or big, abroting...\n")
  }
  
  Cf <- list()
  for( i in 1:nosamples ){
    Cf[[i]] <- ftype
  }
  ## indices for averaging +-t
  i1 <- c(2:(T/2))
  i2 <- c(T:(T/2+2))
  sign <- +1
  if(!sym) sign <- -1

  for(f in files){
    ifs <- paste(path, f, sep="")
    if( !(ifs %in% excludelist) && file.exists(ifs)){
      to.read <- file(ifs,"rb")
      tmp <- array(readBin(to.read, what=ftype, n=T*nosamples, endian=endian), dim=c(T, nosamples))
      close(to.read)
      # average over increasing numbers of samples
      for( i in 1:nosamples ){
        tmp2 <- ftype
        if( i == 1 ){
          tmp2 <- tmp[,1]
        } else {
          tmp2 <- apply(X=tmp[,1:i],MARGIN=1,FUN=mean)
        }
        # average over +- t
        tmp2[i1] <- 0.5 * ( tmp2[i1] + sign * tmp2[i2] )
        Cf[[i]] <- cbind(Cf[[i]],tmp2[1:(T/2+1)])
      }
    } else if(!file.exists(ifs)) {
      cat("file ", ifs, "does not exist...\n")
    }
  }

  ret <- list()
  for( i in 1:nosamples ){
    ret[[i]] <- list(cf=t(Re(Cf[[i]])), icf=t(Im(Cf[[i]])), Time=T, nrStypes=1, nrObs=1, boot.samples=FALSE, jackknife.samples=FALSE)
    attr(ret[[i]], "class") <- c("cf", class(ret[[i]]))
  }
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

readgradflow <- function(path, skip=0, basename="gradflow", col.names) {
  files <- getorderedfilelist(path=path, basename=basename, last.digits=6)
  # the trajectory numbers deduced from the filename
  gaugeno <- getconfignumbers(files, basename=basename, last.digits=6)
  files <- files[(skip+1):length(files)]
  if(length(files)==0) stop(sprintf("readgradflow: no tmlqcd gradient flow files found in path %s",path))

  tmpdata <- data.frame()
  if(missing(col.names)) {
    tmpdata <- read.table(file=files[1],colClasses="numeric",header=TRUE,stringsAsFactors=FALSE)
  }
  else {
    tmpdata <- read.table(file=files[1],colClasses="numeric",col.names=col.names,stringsAsFactors=FALSE)
  }
  ## add the trajectory number, if not present
  if(is.null(tmpdata$traj)) tmpdata$traj <- gaugeno[1]
  
  fLength <- length(tmpdata$t)
  nFiles <- length(files)
  nCols <- ncol(tmpdata)
  ## we generate the full size data.frame first
  tmpdata[,] <- NA
  ldata <- tmpdata
  ldata[((nFiles-1)*fLength+1):(nFiles*fLength),] <- tmpdata
  rm(tmpdata)
  
  pb <- NULL
  pb <- txtProgressBar(min = 1, max = length(files), style = 3)
  for( i in 1:length(files) ){
    setTxtProgressBar(pb, i)
    tmp <- data.frame()
    if(missing(col.names)) {
      tmp <- read.table(file=files[i], colClasses="numeric", header=TRUE, stringsAsFactors=FALSE)
    }
    else {
      tmp <- read.table(file=files[i], colClasses="numeric", col.names=col.names, stringsAsFactors=FALSE)
    }
    if(is.null(tmp$traj)) tmp$traj <- gaugeno[i]
    # the tmlqcd gradient flow routine has the tendency to crash, so we check if the files
    # are complete 
    if( dim( tmp )[1] != fLength ) {
      warning(sprintf("For file %s, number of rows is not correct: %d instead of %d\n",files[i],dim(tmp)[1],fLength) )
      ldata[((i-1)*fLength+1):(i*fLength),] <- NA
    } else {
      ldata[((i-1)*fLength+1):(i*fLength),] <- tmp
    }
  }
  close(pb)

  # keep only rows which contain all data
  ldata <- ldata[complete.cases(ldata),]

  # order by t as outermost index
  ldata <- ldata[order(ldata$t,ldata$traj),]
  rownames(ldata) <- NULL
  return(invisible(ldata))
}

