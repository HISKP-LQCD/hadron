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
                         colClasses, obs=NULL, obs.index, avg=1, stride=1) {
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
  # when stride is > 1, we only read a subset of files
  nFilesToLoad <- as.integer(nFiles/stride)
  nCols <- length(tmpdata)
  ## we generate the full size data.frame first
  tmpdata[,] <- NA
  ldata <- tmpdata
  ldata[((nFilesToLoad-1)*fLength+1):(nFilesToLoad*fLength),] <- tmpdata
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
    if( !(files[i] %in% excludelist) && file.exists(files[i]) && (i-1) %% stride == 0) {
      
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
      
      dat_idx_start <- ((i-1)/stride*fLength) + 1
      dat_idx_stop <- dat_idx_start+fLength-1
      ldata[dat_idx_start:dat_idx_stop,] <- tmpdata
      rm(tmpdata)
    }
    else if(!file.exists(files[i])) {
      warning(paste("readcmifiles: dropped file", files[i], "because it's missing\n"))
    }
  }
  if(!verbose) {
    close(pb)
  }

  # if we want to average over successive samples, we do this here
  if(avg > 1){
    for(i in seq(1,nFilesToLoad,by=avg)){
      out_idx_start <- (i-1)*fLength + 1
      out_idx_stop <- i*fLength
      for(j in seq(i+1,i+avg-1)){
        print(j)
        avg_idx_start <- (j-1)*fLength + 1
        avg_idx_stop <- j*fLength
        # add the next sample to the sample that we use as an output base
        ldata[out_idx_start:out_idx_stop,] <- ldata[out_idx_start:out_idx_stop,] +
                                              ldata[avg_idx_start:avg_idx_stop,]
        # invalidate the sample that we just added
        ldata[avg_idx_start:avg_idx_stop,] <- NA 
      }
      # take the average over the samples
      ldata[out_idx_start:out_idx_stop,] <- ldata[out_idx_start:out_idx_stop,]/avg
    }
  }
  ## remove NAs from missing files
  ldata <- na.omit(ldata)

  return(invisible(ldata))
}

readcmidatafiles <- function(files, excludelist=c(""), skip=1, verbose=FALSE,
                             colClasses=c("integer", "integer","integer","numeric","numeric"),
                             obs=NULL, obs.index=1, avg=1, stride=1) {

  data <- readcmifiles(files, excludelist=excludelist, skip=skip, verbose=verbose,
                       colClasses=colClasses, obs=obs, obs.index=obs.index, avg=avg, stride=stride)
  attr(data, "class") <- c("cmicor", class(data))  
  return(invisible(data))
}

readcmiloopfiles <- function(files, excludelist=c(""), skip=0, verbose=FALSE,
                             colClasses=c("integer", "integer","integer","integer",
                               "numeric","numeric","numeric","numeric"),
                             obs=NULL, obs.index=2) {
  data <- readcmifiles(files, excludelist=excludelist, skip=skip, verbose=verbose,
                       colClasses=colClasses, obs=obs, obs.index=obs.index, avg=1, stride=1)
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

  cf <- cf_meta(nrObs = 1, Time = T, nrStypes = 2)
  cf <- cf_orig(cf,
                cf = array(ldata[,ind.vec[4]], dim=c(T, nrSamples, length(ldata[,ind.vec[4]])/T/nrSamples))/sqrt(L^3),
                icf = array(ldata[,ind.vec[5]], dim=c(T, nrSamples, length(ldata[,ind.vec[5]])/T/nrSamples))/sqrt(L^3))
  cf <- cf_smeared(cf,
                   scf = array(ldata[,ind.vec[6]], dim=c(T, nrSamples, length(ldata[,ind.vec[6]])/T/nrSamples))/sqrt(L^3),
                   iscf =  array(ldata[,ind.vec[7]], dim=c(T, nrSamples, length(ldata[,ind.vec[7]])/T/nrSamples))/sqrt(L^3),
                   nrSamples = nrSamples,
                   obs = obs)

  # TODO: This should be set via a constructor.
  cf$conf.index <- unique(ldata[,ind.vec[8]])

  return (invisible(cf))
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
  
  ret <- cf_meta(nrObs = nrObs, Time = Time, nrStypes = nrStypes, symmetrised = symmetrise)
  ret <- cf_orig(ret, cf = cf, icf = NA)

  return (invisible(ret))
}

readhlcor <- function(filename) {
  return(invisible(read.table(filename, header=FALSE,
                              colClasses=c("integer", "integer","integer","integer","numeric","numeric","numeric","numeric","integer"))))
}

readoutputdata <- function(filename) {
  data <- read.table(filename, header=FALSE)
  attr(data, "class") <- c("outputdata", "data.frame")  
  return(invisible(data))
}

readtextcf <- function(file, T=48, sym=TRUE, path="", skip=1, check.t=0, ind.vector=c(2,3), symmetrise=TRUE,
                       stride=1, avg=1, Nmin=4, autotruncate=TRUE) {
  stopifnot(!missing(file))
  stopifnot(T >= 1)

  tmp <- read.table(paste(path, file, sep=""), skip=skip)
  stopifnot(!((length(ind.vector) < 2) || (max(ind.vector) > length(tmp)) || (min(ind.vector) < 1)))
  
  if(check.t > 0 && max(tmp[[check.t]]) != T-1) {
    stop("T in function call does not match the one in the file, aborting...\n")
  }

  if(length(tmp[[ind.vector[1]]]) %% T != 0) {
    stop("T does not devide the number of rows in file, aborting... check value of paramter skip to readtextcf!\n")
  }
  
  ii <- c(1:(T/2+1))

  tmp <- array(tmp[[ind.vector[1]]] + 1i*tmp[[ind.vector[2]]], dim=c(T, length(tmp[[ind.vector[1]]])/T))
  if( (stride > 1 | avg > 1) & (ncol(tmp) %% (stride*avg) != 0) ){
    if(autotruncate){
      cat(sprintf("stride=%d, avg=%d, ncol=%d\n",stride,avg,ncol(tmp)))
      cat("readtextcf: Sparsification and/or averaging requested, but their product does not divide the number of measurements!\n")
      cat("readtextcf: Reducing the number of total measurements to fit!\n")
      nmeas <- as.integer( (stride*avg)*floor( ncol(tmp)/(stride*avg) ))
      if( nmeas/(stride*avg) >= Nmin ){
        tmp <- tmp[,1:nmeas]
      } else {
        cat(sprintf("readtextcf: After sparsification and averaging, less than %d measurements remain, disabling sparsification and averaging!\n",Nmin))
        stride <- 1
        avg <- 1
      }
    } else {
      stop("readtextcf: Sparsification and/or averaging requested, but their product does not divide the number of measurements!\n")
    }
  }

  ## sparsify data
  if(stride > 1){
    sp.idx <- seq(from=1,to=ncol(tmp),by=stride)
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
  
  ret <- cf_meta(nrObs = 1, Time=T, nrStypes = 1)
  ret <- cf_orig(ret, cf = t(Re(tmp)), icf = t(Im(tmp)))

  if (symmetrise) {
    sign <- +1
    if (!sym) sign <- -1

    ret <- symmetrise.cf(ret, sign)
  }

  return (invisible(ret))
}

#' @title reader for Nissa text format correlation functions
#' @param file_baseames_to_read Character vector of file names without the
#'                              smearing combination suffixes (such as 'll', 'ls', 'sl', 'ss')
#'                              which will be added in the reading routine accordign to what was 
#'                              passed via `smear_combs_to_read`. An example would be
#'                              '0001/mes_contr_2pts', not the lack of the smearing suffix.
#' @param smear_combs_to_read Character vector containing the smearing cominations that are to be read.
#'                            These will be attached to the `file_basenames_to_read` in the reading routine.
#' @param Time Integer, time extent of the lattice.
#' @param combs_to_read Data frame containing the indices of the masses and r-paramter combinations to
#'                      be read as well as the name of the spin combination.
#'                      For a two-point function using the second and third mass (0-indexed), 
#'                      the (+^dag,+) r-combination and the pseudoscalar-pseudoscalar spin combination
#'                      would look as follows:
#'                      \tabular{rrrrr}{
#'                        m1_idx \tab m2_idx \tab r1_idx \tab r2_idx \tab spin_comb \cr
#'                        1      \tab 2      \tab 0      \tab 0      \tab "P5P5"
#'                      }
#' @param sym.vec Integer or numeric vector. Specifies whether the correlator at
#'                the given position is symmetric (+1.0) or anti-symmetric (-1.0 )
#'                under time reflection. This is passed to \code{symmetrise.cf}. This
#'                should be of sufficient length to cover all correlators that are
#'                going to be read (one number per row of \code{combs_to_read} and
#'                per entry of \code{smear_combs_to_read})
#' @param symmetrise Boolean, specifies whether averaging over backward and forward
#'                   correlators should be done after the correlator has been read in.
#' @param nts Integer, number of time slices to be read from the correlator files.
readnissatextcf <- function(file_basenames_to_read,
                            smear_combs_to_read,
                            Time,
                            combs_to_read,
                            nts = Time, 
                            sym.vec = c(1),
                            symmetrise = FALSE)
{
  tmp <- read_nissa_textcf_kernel(file_basenames_to_read,
                                  smear_combs_to_read,
                                  nts,
                                  combs_to_read)

  total_nts <- nts*length(smear_combs_to_read)*nrow(combs_to_read)

  realcols <- seq(1,2*total_nts,2)
  imagcols <- seq(2,2*total_nts,2)

  cf <- cf_meta(nrObs = nrow(combs_to_read), Time=Time, nrStypes = length(smear_combs_to_read), symmetrised = FALSE)
  cf <- cf_orig(cf, cf = tmp[,realcols,drop=FALSE], icf = tmp[,imagcols,drop=FALSE])

  if(symmetrise){
    # in some cases it makes sense to store only a subset of the time slices of a
    # correlation function. In this case, symmetrisation is not possible unless
    # the missing time slices are reconstructed or added manually somehow.
    if( nts != Time ){
      stop("The time extent and the number of time slices in the correlator do not agree, cannot symmetrise!")
    }
    cf <- symmetrise.cf(cf, sym.vec)
  }
  return (invisible(cf))
}


readbinarycf <- function(files, 
                         T, 
                         obs=5, 
                         Nop=1,
                         symmetrise=TRUE,
                         endian="little",
                         op="aver",
                         excludelist=c(""), 
                         sym=TRUE, 
                         path="",
                         hdf5format=FALSE, 
                         hdf5name, 
                         hdf5index=c(1,2)) {
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
    h5avail <- requireNamespace('rhdf5')
    if(!h5avail) stop("rhdf5 package not installed, aborting...\n")
    obs <- 1
    Nop <- 1
    if(length(hdf5index)<2) hdf5index <- c(hdf5index, hdf5index)
  }
  ii <- c(1:(Nop*T))+obs*T

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
        LS <- rhdf5::h5ls(ifs)
        if(as.integer(LS[LS$name == hdf5name,]$dim) != T) {
          stop(paste(hdf5name, "in file", ifs, "has not the correct length, aborting...\n"))
        }
        tmp <- rhdf5::h5read(file=ifs, name=hdf5name)
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

      Cf <- cbind(Cf, tmp[1:T])
      
      if(!hdf5format) {
        close(to.read)
      }
      else {
        if(exists("rhdf5::H5close")) rhdf5::H5close()
      }
    }
    else if(!file.exists(ifs)) {
      cat("file ", ifs, "does not exist...\n")
    }
  }

  ret <- cf_meta(nrObs = 1, Time=T, nrStypes = 1, symmetrised = FALSE)
  ret <- cf_orig(ret, cf = t(Re(Cf)), icf = t(Im(Cf)))

  if (symmetrise) {
    sign <- +1
    if (!sym) sign <- -1
    ret <- symmetrise.cf(ret, sign)
  }

  return (invisible(ret))
}


#' Read binary correlation function by sample
#'
#' Read binary correlation functions sample by sample, return as a list of
#' length `nosamples` where increasing indices refer to averaging over
#' increasing numbers of samples.
#'
#' @param files character vector. Paths to the file to read. As `path` is
#' prepended to each element, one can also just pass the filenames here.
#' @param T numeric. Time extent.
#' @param endian character, either `little` or `big`.
#' @param path character. Path that is prefixed to each of the paths given in
#' `files`.
#' @param excludelist character vector. Elements in `files` that are specified
#' in `excludelist` are skipped. The caller could also just pass
#' `setdiff(files, excludelist)`.
#' @param sym logical. Whether the read data shall be symmetrized in the end.
#' @param ftype numeric type. As the data is read in binary this type has to
#' match exactly the one in the file.
readbinarysamples <- function(files, T=48, nosamples=2, endian="little",
                              op="aver", excludelist=c(""), sym=TRUE, path="", ftype=double() ){

  if(missing(files)) {
    stop("files must be given! Aborting...\n")
  }
  stopifnot(length(files) > 0)

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
        Cf[[i]] <- cbind(Cf[[i]], tmp2)
      }
    } else if(!file.exists(ifs)) {
      cat("file ", ifs, "does not exist...\n")
    }
  }

  ret <- list()
  for (i in 1:nosamples) {
    ret[[i]] <- cf_meta(nrObs = 1, Time = T, nrStypes = 1, symmetrised = FALSE)
    ret[[i]] <- cf_orig(ret[[i]], cf = t(Re(Cf[[i]])), icf = t(Im(Cf[[i]])))

    sign <- +1
    if (!sym) sign <- -1
    ret[[i]] <- symmetrise.cf(ret[[i]], sign)
  }

  return (invisible(ret))
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

  cf <- cf_meta(Time = T)
  cf <- cf_orig(cf, cf = Re(Cf), icf = Im(Cf))
  cf <- cf_smeared(cf, scf = NA, iscf = NA, nrSamples = nrSamples, obs = obs)

  return (invisible(cf))
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

  cf <- cf_meta(nrObs = 1, Time = T, nrStypes = 2)
  cf <- cf_orig(cf,
                cf = array(ldata[, ind.vec[4]], dim=c(T, nrSamples, nFiles))/sqrt(L^3),
                icf = array(ldata[, ind.vec[5]], dim=c(T, nrSamples, nFiles))/sqrt(L^3))
  cf <- cf_smeared(cf,
                   scf = array(ldata[, ind.vec[6]], dim=c(T, nrSamples, nFiles))/sqrt(L^3),
                   iscf= array(ldata[, ind.vec[7]], dim=c(T, nrSamples, nFiles))/sqrt(L^3),
                   nrSamples = nrSamples,
                   obs = obs)

  return (invisible(cf))
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

