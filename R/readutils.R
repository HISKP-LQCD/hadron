#' @export
readcmicor <- function(filename, colClasses=c("integer","integer","integer","numeric","numeric","integer"),
                       skip=0) {
  data <- read.table(filename, header=F, skip=skip,
                     colClasses=colClasses)
  attr(data, "class") <- c("cmicor", class(data))  
  return(invisible(data))
}



#' Creates an ordered filelist from a basename and a path
#' 
#' These functions generate an ordered filelist and an order list of config
#' numbers by using a path and a basename and '*'.
#' 
#' All filenames are assumend to have equal length.
#' 
#' @param path the path to be searched
#' @param basename the basename of the files
#' @param last.digits the number of last characters in each filename to be used
#' for ordering the list.
#' @param ending the file extension after the digits.
#' @return returns the ordered list of strings.
#' @author Carsten Urbach, \email{curbach@@gmx.de}
#' @seealso \code{\link{readcmidatafiles}}, \code{\link{extract.obs}}
#' @keywords file
#' @export
#' @examples
#'
#' filelist <- getorderedfilelist(path=paste0(system.file(package="hadron"), "/extdata/"),
#'                                basename="testfile", last.digits=3, ending=".dat")
#' filelist
#' 
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

#' Creates an ordered vector of gauge config file numbers
#' 
#' These functions generate an ordered list of config
#' numbers by using a path and a basename and '*'.
#' 
#' All filenames are assumend to have equal length.
#' 
#' @param path the path to be searched
#' @param basename the basename of the files
#' @param last.digits the number of last characters in each filename to be used
#' for ordering the list.
#' @param ending the file extension after the digits.
#' @return returns the ordered list of gauge config numbers as a numeric vector.
#' @author Carsten Urbach, \email{curbach@@gmx.de}
#' @seealso \code{\link{readcmidatafiles}}, \code{\link{extract.obs}}
#' @keywords file
#'
#' @examples
#' confignumbers <- getorderedconfignumbers(path=paste0(system.file(package="hadron"), "/extdata/"),
#'                                basename="testfile", last.digits=3, ending=".dat")
#' confignumbers
#' @export
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



#' Read Single Data Files in Chris Michael Format
#' 
#' reads data from single files in Chris Michael format
#' 
#' These functions reads data from single data files. It is assumed that every
#' file has the same number of columns.
#' 
#' The cmi (Chris Michael) format for connected correlators comprises 6 colums
#' per file: 1) the observable type number (itype); 2) the operator type number
#' (iobs); 3) the time difference from source going from 0 to \eqn{Time/2} for
#' each operator type; 4) \eqn{c_1}{c1} correlator value at time value forward
#' in time; 5) \eqn{c_2}{c2} correlator value at time value backward in time;
#' 6) number of gauge configuration.
#' 
#' There are scripts shipped with the package converting the output written
#' into seperate files for each gauge configuration into the expected format.
#' They are called \code{puttogether.sh} and \code{puttogether_reverse.sh}
#' which will sort with increasing and with decreasing gauge configuration
#' number, respectively.
#' 
#' Note, that the normalisation of correlators needs multiplication by factor
#' of \eqn{0.5} (and possible \eqn{(2*\kappa)^2}{(2*k)^2} and \eqn{L^3} factors
#' dependent on your conventions).
#' 
#' The values of \code{itype} run from \code{1} to the total number of gamma
#' matrix combinations available. \code{iobs} equals \code{1} for local-local
#' correlators, \code{3} for local-smeared, \code{5} for smeared-local and
#' \code{7} for smeared-smeared
#' 
#' For charged mesons the order of gamma-matrix combinations is as follows:\cr
#' order PP PA AP AA 44 P4 4P A4 4A for pion like \eqn{P=\gamma_5}{P=g5}
#' \eqn{A=\gamma_4\gamma_5}{A=g4g5} \eqn{4=\gamma_4}{4=g4}\cr order 44 VV AA 4V
#' V4 4A A4 VA AV for rho-a1 like \eqn{4=\gamma_i\gamma_4}{4=gig4}
#' \eqn{V=\gamma_i}{V=gi} \eqn{A=\gamma_i\gamma_5}{A=gig5}\cr order BB SS -
#' total 20 \eqn{\gamma_i\gamma_4\gamma_5}{B=gig4g5} \eqn{S=I}\cr itype=21 is
#' conserved vector current at sink, \eqn{\gamma_5}{g5} at source
#' 
#' For neutral mesons the order of gamma-matrix combinations is as follows:\cr
#' order PP PA AP AA II PI IP AI IA for pion like \eqn{P=\gamma_5}{P=g5}
#' \eqn{A=\gamma_4\gamma_5}{A=g4g5} \eqn{I=1}{1=1}\cr order 44 VV BB 4V V4 4B
#' B4 VB BV for rho-b1 like \eqn{4=\gamma_i\gamma_4}{4=gig4}
#' \eqn{V=\gamma_i}{V=gi} \eqn{B=\gamma_i\gamma_4\gamma_5}{B=gig4g5}\cr order
#' XX AA - total 20 for a0-X like \eqn{A=\gamma_i\gamma_5}{A=gig5}
#' \eqn{X=\gamma_4}{X=g4}
#' 
#' For loops (disconnected contributions to neutral mesons) the convention is
#' as follows: files are assumed to have eight columns with gauge, gamma, t,
#' sample, ReTL, ImTL, ReTF, ImTF, where gamma is 1 to 16 as list of
#' (hermitian) gamma matrices: order g_5 g_1 g_2 g_3\cr -ig_4* g_5 g_1 g_2
#' g_3\cr -ig_5* i*g_5 g_1 g_2 g_3 ie 1,..\cr -ig_5g_4* -i*g_5 g_1 g_2 g_3 ie
#' g_4, g_5*row 2\cr (so P is 1; A4 is 5; S is 9; A_i is 10,11,12 etc)
#' 
#' t is t-value of trace (here spatial momentum is zero) sample is sample
#' number 1,...24 (or 96) ReTL is real part of trace at time t, with gamma
#' combination given and Local operator (F is Fuzzed == non-local) operator).
#' 
#' Normalisation is trace M^-1 with M=1+...
#' 
#' To make a disconnected correlator, one combines these traces for different t
#' (and different sample number) as a product. Note only Re Gamma=1 and Im
#' Gamma=gamma_5 have VEV's, see \code{\link{computeDisc}}
#' 
#' @aliases readcmicor readcmifiles readcmidatafiles readcmiloopfiles
#' @param files list of filenames to be read. Can be created using
#' \code{getorderedfilelist}.
#' @param skip Number of lines to be skipped at the beginning of each file
#' @param excludelist files to exclude from reading.
#' @param verbose Increases verbosity of the function.
#' @param colClasses The column data type classes, the \code{read.table}.
#' @param obs To reduce memory consumption it is possible to extract only one
#' of the observales. The column in which to match \code{obs} is to be given
#' with \code{obs.index}. This will only be affective if \code{obs} is not
#' \code{NULL}.
#' @param obs.index The column in which to match \code{obs} is to be given with
#' \code{obs.index}.
#' @param avg Integer. Average over successive number samples
#' @param stride Integer. Read only subset of files with corresponding stride.
#' @return \code{readcmicor} returns an object of class \code{cmicor}, read
#' from a single file.
#' 
#' \code{readcmidatafiles} returns an object of class \code{cmicor}, which is
#' an \code{rbind} of all \code{data.frame}s read from the single files in the
#' filelist.
#' 
#' \code{readcmiloopfiles} returns an object of class \code{cmiloop}, which is
#' an \code{rbind} of all \code{data.frame}s read from the single files in the
#' filelist.
#' @author Carsten Urbach, \email{curbach@@gmx.de}
#' @seealso \code{\link{getorderedfilelist}}, \code{\link{extract.obs}},
#' \code{\link{readcmidisc}}
#' @keywords file
#' @examples
#'
#' ## a running toy example
#' files <- paste0(system.file(package="hadron"), "/extdata/outprcvn.dddd.00.0000")
#' X <- readcmifiles(files, skip=0,
#'                   colClasses=c("integer", "integer","integer","numeric","numeric"))
#' X
#'
#' ## a more realistic example
#' \dontrun{filelist <- getorderedfilelist("ouptrc", last.digits=3, ending=".dat")}
#' \dontrun{cmicor <- readcmidatafiles(filelist, skip=1)}
#' 
#' @export readcmifiles
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
    pb <- txtProgressBar(min = 0, max = nFiles, style = 3)
  }
  for(i in c(1:nFiles)) {
    # update progress bar
    if(!verbose) {
      setTxtProgressBar(pb, i)
    }
    if( !(files[i] %in% excludelist) && file.exists(files[i]) && (i-1) %% stride == 0) {
      
      if(verbose) {
        message("Reading from file ", files[i], "\n")
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

#' @export
readcmidatafiles <- function(files, excludelist=c(""), skip=1, verbose=FALSE,
                             colClasses=c("integer", "integer","integer","numeric","numeric"),
                             obs=NULL, obs.index=1, avg=1, stride=1) {

  data <- readcmifiles(files, excludelist=excludelist, skip=skip, verbose=verbose,
                       colClasses=colClasses, obs=obs, obs.index=obs.index, avg=avg, stride=stride)
  attr(data, "class") <- c("cmicor", class(data))  
  return(invisible(data))
}

#' @export
readcmiloopfiles <- function(files, excludelist=c(""), skip=0, verbose=FALSE,
                             colClasses=c("integer", "integer","integer","integer",
                               "numeric","numeric","numeric","numeric"),
                             obs=NULL, obs.index=2) {
  data <- readcmifiles(files, excludelist=excludelist, skip=skip, verbose=verbose,
                       colClasses=colClasses, obs=obs, obs.index=obs.index, avg=1, stride=1)
  attr(data, "class") <- c("cmiloop", class(data))
  return(invisible(data))
}



#' Extract a single loop from an object of class \code{cmiloop}
#' 
#' Extracts all loop values from an object of class \code{cmiloop} for all
#' available times, samples and configurations.
#' 
#' 
#' @param cmiloop input object of class \code{cmiloop} generated for instance
#' with \code{readcmiloopfiles}.
#' @param obs the observable to extract
#' @param ind.vec index vector to be used during extraction with
#' \code{ind.vec[1]} the column with the observable number, \code{ind.vec[2]}
#' the time values, \code{ind.vec[3]} the sample numbers, \code{ind.vec[4]} the
#' real part of the local loop, \code{ind.vec[5]} the imaginary part of the
#' local loop, \code{ind.vec[6]} and \code{ind.vec[7]} the same for fuzzed (or
#' smeared) loops and \code{ind.vec[8]} for the configuraton number.
#' @param L The spatial lattice extent needed for normalisation. If not given
#' set to \code{Time/2}.
#' @return a list with elements as follows:
#' 
#' \code{cf}: real part of the local loop
#' 
#' \code{icf}: imaginary part of the local loop
#' 
#' \code{scf}: real part of the smeared loop
#' 
#' \code{iscf}: imaginary part of the smeared loop
#' 
#' \code{Time=Time}, \code{nrSamples}, \code{nrObs=1}, \code{nrStypes=2},
#' \code{obs=obs} and \code{conf.index}. The last is the list of configurations
#' corresponding to the loops.
#' @author Carsten Urbach, \email{curbach@@gmx.de}
#' @seealso \code{\link{readcmiloopfiles}}
#' @export extract.loop
extract.loop <- function(cmiloop, obs=9, ind.vec=c(2,3,4,5,6,7,8,1), L) {
  ldata <- cmiloop[cmiloop[,ind.vec[1]] == obs,] 
  Time <- max(ldata[,ind.vec[2]])
  nrSamples <- max(ldata[,ind.vec[3]])
  if(missing(L)) {
    L <- Time/2
  }

  cf <- cf_meta(nrObs = 1, Time = Time, nrStypes = 2)
  cf <- cf_orig(cf,
                cf = array(ldata[,ind.vec[4]], dim=c(Time, nrSamples, length(ldata[,ind.vec[4]])/Time/nrSamples))/sqrt(L^3),
                icf = array(ldata[,ind.vec[5]], dim=c(Time, nrSamples, length(ldata[,ind.vec[5]])/Time/nrSamples))/sqrt(L^3))
  cf <- cf_smeared(cf,
                   scf = array(ldata[,ind.vec[6]], dim=c(Time, nrSamples, length(ldata[,ind.vec[6]])/Time/nrSamples))/sqrt(L^3),
                   iscf =  array(ldata[,ind.vec[7]], dim=c(Time, nrSamples, length(ldata[,ind.vec[7]])/Time/nrSamples))/sqrt(L^3),
                   nrSamples = nrSamples,
                   obs = obs)

  # TODO: This should be set via a constructor.
  cf$conf.index <- unique(ldata[,ind.vec[8]])

  return (invisible(cf))
}



#' Extract One or More Gamma Combinations from am CMI Correlator
#' 
#' Extracts one or more gamma matrix combinations (observables) from a
#' correlator stored in cmi format
#' 
#' C(t) and C(-t) are averaged as indicated by \code{sym.vec}.
#' 
#' @param cmicor an correlator object in cmi format
#' @param vec.obs vector containing the numbers of observables to be extracted.
#' @param ind.vec Index vector indexing the column numbers in cmicor to be
#' used. The first must be the observable index, the second the smearing type
#' index, the third the time, the fourth C(+t) and the fifth C(-t).
#' 
#' Index vector indexing the column numbers in cmiloop to be used. The first
#' must be the observable index, the second the smearing type index, the third
#' the time, the fourth ReTL, the fifth ImTL, the sixth ReTF and the seventh
#' ImTF.
#' @param verbose Increases verbosity of the function.
#' @param sym.vec a vector of bools of length equal to the number of
#' observables indicating whether C(t) is symmetric in t, i.e. whether C(+t)
#' and C(-t) should be added or subtracted. If not given C(+t) and C(-t) will
#' be assumed to be symmetric.
#' @param sign.vec a sign vector of length equal to the number of observables
#' indicating whether the corresponding correlation function should be
#' multiplied by +-1.
#' @param symmetrise if set to \code{TRUE}, the correlation function will be
#' averaged for \code{t} and \code{Time-t}, with the sign depending on the value
#' of \code{sym}.  Note that currently the correlator with t-values larger than
#' \code{Time/2} will be discarded.
#' @return returns a list containing \item{cf}{ for \code{extract.obs}: array
#' containing the correlation function with dimension number of files times
#' (nrObs*nrStypes*(Time/2+1)). C(t) and C(-t) are averaged according to
#' \code{sym.vec}.
#' 
#' for \code{extract.loop}: ReTL } \item{icf}{ for \code{extract.loop} only:
#' ImTL } \item{scf}{ for \code{extract.loop} only: ReTF } \item{sicf}{ for
#' \code{extract.loop} only: ImTF } \item{Time}{ The time extent of the
#' correlation functions.  } \item{nrStypes}{ The number of smearing
#' combinations.  } \item{nrObs}{ The number of observables.  }
#' \item{nrSamples}{ for \code{extrac.loop} only: the number of samples found
#' in the files.  }
#' @author Carsten Urbach, \email{curbach@@gmx.de}
#' @seealso \code{\link{readcmicor}}, \code{\link{readcmidatafiles}},
#' @keywords ts
#' @examples
#' 
#' files <- paste0(system.file(package="hadron"), "/extdata/outprcvn.dddd.00.0000")
#' X <- readcmifiles(files, skip=0,
#'                   colClasses=c("integer", "integer","integer","numeric","numeric"))
#' Y <- extract.obs(X)
#' Y
#' 
#' @export extract.obs
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
    stop("extract.obs: data inconsistent, Time not equal to what was found in the input data\n")
  }
  if(verbose) message("extract.obs: nrObs=",nrObs, " nrStypes=",nrStypes, " Time=", Time, "\n")

  data <- cmicor[cmicor[,ind.vec[1]] %in% vec.obs,]
  cf <- NULL
  
  if(symmetrise){
    ## we divide everything by 2 apart from t=0 and t=Time/2
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
          # since we are not symmetrising, the individual correlators have Time entries 
          lhs <- c((bw+1):(Thalf-bw)) + (i-1)*nrStypes*Time  + (s-1)*Time  + bw*(Thalf-1)
          # in the cmi format, the backward correlator is on time-slices 1 to Time/2-1 (indices 2 to Time/2)
          rhs <- c((bw+1):(Thalf-bw)) + (i-1)*nrStypes*Thalf + (s-1)*Thalf
          # but in "reverse" order (to ease averaging)
          if( bw == 1 ) rhs <- rev(rhs)
          cf[,lhs] <- tmp[,rhs]
        }
      }
    }
  }
  
  ret <- cf_meta(nrObs = nrObs, Time = Time, nrStypes = nrStypes, symmetrised = symmetrise)
  ret <- cf_orig(ret, cf = cf, icf = NULL)

  return (invisible(ret))
}

#' readhlcor
#'
#' @param filename String. Filename of the heavy light correlator data file. The
#'    file is expected to have nine columns, the first four integer, the second four numeric
#'    and the last integer valued again.
#'
#' @return
#' Invisibly returns a \link{data.frame} object containing the file content.
#' @export
readhlcor <- function(filename) {
  return(invisible(read.table(filename, header=FALSE,
                              colClasses=c("integer", "integer","integer","integer","numeric","numeric","numeric","numeric","integer"))))
}



#' Read Data In output.data Format of tmLQCD
#' 
#' reads data from an output.data file written by tmLQCD
#' 
#' The data can be plotted directly using \dQuote{plot}.
#' 
#' @param filename filename of the data file
#' @return returns a data frame of class \dQuote{outputdata} containing the
#' data.
#' @author Carsten Urbach \email{curbach@@gmx.de}
#' @keywords file
#' @return
#' Returns an object of class `outputdata` derived from a data.frame
#' as generated by \link{read.table} applied to the input file.
#' 
#' @examples
#' 
#' plaq <- readoutputdata(paste0(system.file(package="hadron"), "/extdata/output.data"))
#' plot(plaq)
#' 
#' @export readoutputdata
readoutputdata <- function(filename) {
  data <- read.table(filename, header=FALSE)
  attr(data, "class") <- c("outputdata", "data.frame")  
  return(invisible(data))
}



#' Read correlator data from single file
#' 
#' Reads arbitrary number of samples for a complex correlation function from a
#' text file.
#' 
#' 
#' @param file filename of file to read from.
#' @param Time time extent of the correlation function
#' @param sym if \code{TRUE} average C(+t) and C(-t), otherwise C(+t) and
#' -C(-t). Averaging can be switched off using the \code{symmetrise} option.
#' @param skip number of lines to skip at beginning of file
#' @param check.t if set to an integer value larger than zero the function will
#' assume that in the corresponding column of the file the Euclidean time is
#' counted and it will check whether the maximum in this column is identical to
#' Time-1.
#' @param ind.vector index vector of length 2 with the indices of real and
#' imaginary values of correlator, respectivley.
#' @param symmetrise if set to \code{TRUE}, the correlation function will be
#' averaged for \code{t} and \code{Time-t}, with the sign depending on the value
#' of \code{sym}. Note that currently the correlator with t-values larger than
#' \code{Time/2} will be discarded.
#' @param path the path to the files.
#' @param autotruncate Boolean. Whether to autotruncate or not
#' @param avg Integer. Average over successive number samples
#' @param stride Integer. Read only subset of files with corresponding stride.
#' @param Nmin Integer. Minimal number of measurements that must remain after
#' sparsification and averaging. Default equals to 4.
#' @return returns a list with two arrays \code{cf} and \code{icf} with real
#' and imaginary parts of the correlator, and integers \code{Time},
#' \code{nrStypes=1} and \code{nrObs=1}. Both of the arrays have dimension
#' \code{c(N, (Time/2+1))}, where \code{N} is the number of measurements
#' (gauges).  \code{Time} is the time extent, \code{nrStypes} the number of
#' smearing levels and \code{nrObs} the number of operators, both of which are
#' currently fixed to 1.
#' @author Carsten Urbach, \email{curbach@@gmx.de}
#' @seealso \code{\link{readcmidatafiles}}, \code{\link{readbinarydisc}},
#' \code{\link{readcmidisc}}, \code{\link{readcmicor}},
#' \code{\link{readbinarycf}}
#' @keywords file
#' @export readtextcf
readtextcf <- function(file, Time=48, sym=TRUE, path="", skip=1, check.t=0, ind.vector=c(2,3), symmetrise=TRUE,
                       stride=1, avg=1, Nmin=4, autotruncate=TRUE) {
  stopifnot(!missing(file))
  stopifnot(Time >= 1)

  tmp <- read.table(paste(path, file, sep=""), skip=skip)
  stopifnot(!((length(ind.vector) < 2) || (max(ind.vector) > length(tmp)) || (min(ind.vector) < 1)))
  
  if(check.t > 0 && max(tmp[[check.t]]) != Time-1) {
    stop("Time in function call does not match the one in the file, aborting...\n")
  }

  if(length(tmp[[ind.vector[1]]]) %% Time != 0) {
    stop("Time does not devide the number of rows in file, aborting... check value of paramter skip to readtextcf!\n")
  }
  
  ii <- c(1:(Time/2+1))

  tmp <- array(tmp[[ind.vector[1]]] + 1i*tmp[[ind.vector[2]]], dim=c(Time, length(tmp[[ind.vector[1]]])/Time))
  if( (stride > 1 | avg > 1) & (ncol(tmp) %% (stride*avg) != 0) ){
    if(autotruncate){
      warning(sprintf("stride=%d, avg=%d, ncol=%d\n",stride,avg,ncol(tmp)))
      warning("readtextcf: Sparsification and/or averaging requested, but their product does not divide the number of measurements!\n")
      warning("readtextcf: Reducing the number of total measurements to fit!\n")
      nmeas <- as.integer( (stride*avg)*floor( ncol(tmp)/(stride*avg) ))
      if( nmeas/(stride*avg) >= Nmin ){
        tmp <- tmp[,1:nmeas]
      } else {
        warning(sprintf("readtextcf: After sparsification and averaging, less than %d measurements remain, disabling sparsification and averaging!\n",Nmin))
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
    tmp <- array(0, dim=c(Time,ncol(tmp2)/avg))
    for( i in c(1:ncol(tmp)) ){
      tmp[,i] <- (1.0/avg)*apply(X=tmp2[,((i-1)*avg+1):(i*avg)],
                                 MARGIN=1,
                                 FUN=sum)
    }
  }
  
  ret <- cf_meta(nrObs = 1, Time=Time, nrStypes = 1)
  ret <- cf_orig(ret, cf = t(Re(tmp)), icf = t(Im(tmp)))

  if (symmetrise) {
    sign <- +1
    if (!sym) sign <- -1

    ret <- symmetrise.cf(ret, sign)
  }

  return (invisible(ret))
}

#' @title reader for Nissa text format correlation functions
#' @param file_basenames_to_read Character vector of file names without the
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
#'
#' @return
#' Returns an object of class `cf`.
#' 
#' @export
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




#' read correlation function from binary files
#' 
#' Reads a correlation function from binary files, including hdf5 formatted
#' files.
#' 
#' It is assumend that each file contains at least \code{(obs+N)*Time} complex
#' doubles, where \code{Time} is the time extent, \code{obs} is the number of the
#' observable to read in and \code{Nop} the number of replicas for this
#' observable. It is assumed that complex is the fastest running index, next
#' time and then obs. The filelist is assumed to be ordered according to the
#' gauge configuration MC history.
#' 
#' @param files list of filenames to be read. Can be created using
#' \code{getorderedfilelist}. The filelist is assumed to be order according to
#' ascending gauge fields.
#' @param Time time extent of correlation functions.
#' @param obs each file may contain many correlation functions. With 'obs'
#' one choses which observable to read in. To be precise, in each file the
#' reading will start at point Time*obs*sizeof(complex\code{<double>}) and read
#' Nop*Time*sizeof(complex\code{<double>}).
#' @param symmetrise symmetrise the correlation function or not
#' @param Nop number of replicas for the correlator to read in.
#' @param endian the endianess of the binary file.
#' @param excludelist files to exclude from reading.
#' @param sym if \code{TRUE} average C(+t) and C(-t), otherwise C(+t) and
#' -C(-t).
#' @param op the N replicas can be either averaged (\code{op="aver"}) or summed
#' (\code{op="sum"}).
#' @param path path to be prepended to every filename.
#' @param hdf5format if \code{TRUE}, try to read from an hdf5 file.
#' @param hdf5name Name of the data set as a string.
#' @param hdf5index The data might be an array of size n x Time. \code{hdf5index}
#' is used to convert two columns of the data to a complex valued vector using
#' the first and second index for real and imaginary part, respectively. If
#' \code{hdf5index} has length smaller than 2 the first index is reused.
#' @return returns a list with two arrays \code{cf} and \code{icf} with real
#' and imaginary parts of the correlator, and integers \code{Time},
#' \code{nrStypes=1} and \code{nrObs=1}. Both of the arrays have dimension
#' \code{c(N, (Time/2+1))}, where \code{N} is the number of measurements
#' (gauges).  \code{Time} is the time extent, \code{nrStypes} the number of
#' smearing levels and \code{nrObs} the number of operators, both of which are
#' currently fixed to 1.
#' @author Carsten Urbach, \email{curbach@@gmx.de}
#' @seealso \code{\link{readcmidatafiles}}, \code{\link{readbinarydisc}},
#' \code{\link{readcmidisc}}, \code{\link{readcmicor}}
#' @keywords file
#' @examples
#'
#' X <- readbinarycf(path=paste0(system.file(package="hadron"), "/extdata/"),
#'                   files="C2_bin.dat", Time=64, obs=0)
#' X
#' X$cf
#' 
#' @export readbinarycf
readbinarycf <- function(files, 
                         Time, 
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
  if(Time < 1) {
    stop("Time must be larger than 0 and integer, aborting...\n")
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
  ii <- c(1:(Nop*Time))+obs*Time

  Cf <- complex()
  for(f in files) {
    ifs <- paste(path, f, sep="")
    if( !(ifs %in% excludelist) && file.exists(ifs)) {
      tmp <- numeric()
      if(!hdf5format) {
        to.read <- file(ifs, "rb")
        tmp <- readBin(to.read, what=complex(), n=(obs+Nop)*Time, endian = endian)[ii]
      }
      else {
        LS <- rhdf5::h5ls(ifs)
        if(as.integer(LS[LS$name == hdf5name,]$dim) != Time) {
          stop(paste(hdf5name, "in file", ifs, "has not the correct length, aborting...\n"))
        }
        tmp <- rhdf5::h5read(file=ifs, name=hdf5name)
        tmp <- tmp[,hdf5index[1]] + 1i*tmp[,hdf5index[2]]
      }
      ## we average the replica
      if(Nop > 1) {
        if(op == "aver") {
          tmp <- apply(array(tmp, dim=c(Time, Nop)), 1, mean)
        }
        else {
          tmp <- apply(array(tmp, dim=c(Time, Nop)), 1, sum)
        }
      }

      Cf <- cbind(Cf, tmp[1:Time])
      
      if(!hdf5format) {
        close(to.read)
      }
      else {
        if(exists("rhdf5::H5close")) rhdf5::H5close()
      }
    }
    else if(!file.exists(ifs)) {
      warning("file ", ifs, " does not exist...\n")
    }
  }

  ret <- cf_meta(nrObs = 1, Time=Time, nrStypes = 1, symmetrised = FALSE)
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
#' @param Time numeric. Time extent.
#' @param endian character, either `little` or `big`.
#' @param path character. Path that is prefixed to each of the paths given in
#' `files`.
#' @param excludelist character vector. Elements in `files` that are specified
#' in `excludelist` are skipped. The caller could also just pass
#' `setdiff(files, excludelist)`.
#' @param sym logical. Whether the read data shall be symmetrized in the end.
#' @param ftype numeric type. As the data is read in binary this type has to
#' match exactly the one in the file.
#' @param nosamples number of samples
#'
#' @return
#' Returns a \link{list} of `cf` objects.
#' 
#' @export
readbinarysamples <- function(files, Time=48, nosamples=2, endian="little",
                              excludelist=c(""), sym=TRUE, path="", ftype=double() ){

  if(missing(files)) {
    stop("files must be given! Aborting...\n")
  }
  stopifnot(length(files) > 0)

  if(Time < 1) {
    stop("Time must be larger than 0 and integer, aborting...\n")
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
      tmp <- array(readBin(to.read, what=ftype, n=Time*nosamples, endian=endian), dim=c(Time, nosamples))
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
      warning("file ", ifs, " does not exist...\n")
    }
  }

  ret <- list()
  for (i in 1:nosamples) {
    ret[[i]] <- cf_meta(nrObs = 1, Time = Time, nrStypes = 1, symmetrised = FALSE)
    ret[[i]] <- cf_orig(ret[[i]], cf = t(Re(Cf[[i]])), icf = t(Im(Cf[[i]])))

    sign <- +1
    if (!sym) sign <- -1
    ret[[i]] <- symmetrise.cf(ret[[i]], sign)
  }

  return (invisible(ret))
}




#' read disconnected loops from binary files
#' 
#' Reads disconnected loops from binary files.
#' 
#' It is assumend that each file contains O*Time complex doubles, where Time is the
#' time extent and O the number of observables in the file. It is assumed that
#' complex is the fastest running index, next time and then observables. The
#' different samples are assumend to be in different files. The file list is
#' assumed to be ordered with number of samples running fastest, and then
#' number of gauges.
#' 
#' @param files list of filenames to be read. Can be created for instance using
#' \code{getorderedfilelist}. The filelist is assumed to be ordered with number
#' of samples running fastest, and the next to fastest nubmer of gauges.
#' @param Time time extent of correlation functions.
#' @param obs each file may contain Time*obs correlation functions. With
#' \code{obs} one choses which observable to read in.
#' @param endian the endianess of the binary file.
#' @param excludelist files to exclude from reading.
#' @param nrSamples the number of samples
#' @param path path to be prepended to every filename.
#' @return returns a list with two arrays \code{cf} and \code{icf} with real
#' and imaginary parts of the loops, and integers \code{Time},
#' \code{nrStypes=1}, \code{nrSamples} and \code{nrObs=1}. Both of the arrays
#' have dimension \code{c(Time, N)}, where \code{N} is the number of measurements
#' (gauges) and \code{Time} the time extent, \code{nrStypes} the number of smearing
#' levels and \code{nrObs} the number of operators, both of which are currently fixed to 1.
#' @author Carsten Urbach, \email{curbach@@gmx.de}
#' @seealso \code{\link{readcmidatafiles}}, \code{\link{readbinarycf}},
#' \code{\link{readcmidisc}}, \code{\link{readcmicor}}
#' @keywords file
#' @examples
#'
#' ## running toy example
#' file <- paste0(system.file("extdata", package = "hadron"), "/C2_pi0.dat")
#' X <- readbinarydisc(files=file, Time=64, obs=0)
#' X$cf
#'
#' ## more realistic example
#' \dontrun{files <- character()}
#' \dontrun{for(i in seq(600,1744,8)) }
#' \dontrun{  files <- c(files, "C2_dis_u_conf", sprintf("%.04d", i), ".dat", sep="")}
#' \dontrun{cf <- readbinarydisc(files, obs=4, excludelist=c("C2_pi0_conf0632.dat"))}
#' 
#' @export readbinarydisc
readbinarydisc <- function(files, Time=48, obs=5, endian="little",
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
      tmp <- readBin(to.read, complex(), n=(obs+1)*Time, endian = endian)
      Cf <- cbind(Cf, tmp[c((obs*Time+1):((obs+1)*Time))])
      close(to.read)
    }
  }
  Cf <- array(Cf, dim=c(Time, nrSamples, N))

  cf <- cf_meta(Time = Time)
  cf <- cf_orig(cf, cf = Re(Cf), icf = Im(Cf))
  cf <- cf_smeared(cf, scf = NA, iscf = NA, nrSamples = nrSamples, obs = obs)

  return (invisible(cf))
}



#' reads disconnected loops in cmi format
#' 
#' reads disconnected loops in cmi (Chris Michael) format from a list of files.
#' 
#' 
#' @param files list of filenames to be read. Can be created using
#' \code{getorderedfilelist}.
#' @param obs index of operator to parse from files
#' @param ind.vec vector containing the index (column in file) of obs, t,
#' samples, Re(local), Im(local, Re(smeared), Im(smeared).
#' @param excludelist files to exclude from reading.
#' @param skip lines to skip at beginning of each file.
#' @param colClasses The column data type classes, the \code{read.table}.
#' @param L the spatial lattice extent, set to \code{Time/2} if missing.
#' @param debug setting debug to TRUE makes the routine more verbose by
#' spilling out separate filenames.
#' @return returns a list with four arrays \code{cf}, \code{icf} \code{scf} and
#' \code{sicf} containing real and imaginary parts of the local and smeared
#' loops, respectively, and integers \code{Time}, \code{nrStypes=2},
#' \code{nrSamples} and \code{nrObs=1}. The four arrays have dimension
#' \code{c(Time, S, N)}, where \code{S} is the nubmer of samples, \code{Time} is the
#' time extent and \code{N} is the number of measurements (gauges).
#' \code{Time} is the time extent, \code{nrStypes} the number of smearing
#' levels and \code{nrObs} the number of operators, which are currently fixed
#' to 1 and 2, respectively. \code{nrSamples} is the number of samples.
#' 
#' Note that the arrays are normalised by \code{1/sqrt(L^2)}.
#' 
#' The routine expects that all files have identical content. Otherwise the
#' routine will stop.
#' @author Carsten Urbach, \email{curbach@@gmx.de}
#' @seealso \code{\link{readcmidatafiles}}, \code{\link{readbinarycf}},
#' \code{\link{readbinarydisc}}, \code{\link{readcmicor}}
#' @keywords file
#' @examples
#'
#' # a running toy example
#' hpath <- system.file(package="hadron")
#' files <- paste0(hpath, "/extdata/newdisc.0.1373.0.006.k0v4.10")
#' X <- readcmidisc(files=files)
#' X
#'
#' ## a more realistic example
#' \dontrun{v4files <- character()}
#' \dontrun{for(i in seq(600,1744,8))}
#' \dontrun{  v4files <- }
#' \dontrun{   c(v4files, paste("disc.0.163265.0.006.k0v4.", sprintf("%.04d", i), sep=""))}
#' \dontrun{v4data <- readcmidisc(v4files)}
#' @export readcmidisc
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
  Time <- max(ldata[,ind.vec[2]])
  nrSamples <- max(ldata[,ind.vec[3]])

  if(missing(L)) L <- Time/2

  cf <- cf_meta(nrObs = 1, Time = Time, nrStypes = 2)
  cf <- cf_orig(cf,
                cf = array(ldata[, ind.vec[4]], dim=c(Time, nrSamples, nFiles))/sqrt(L^3),
                icf = array(ldata[, ind.vec[5]], dim=c(Time, nrSamples, nFiles))/sqrt(L^3))
  cf <- cf_smeared(cf,
                   scf = array(ldata[, ind.vec[6]], dim=c(Time, nrSamples, nFiles))/sqrt(L^3),
                   iscf= array(ldata[, ind.vec[7]], dim=c(Time, nrSamples, nFiles))/sqrt(L^3),
                   nrSamples = nrSamples,
                   obs = obs)

  return (invisible(cf))
}



#' Read Gradient Flow Output Files in tmLQCD format
#' 
#' given a pathname, reads all gradient flow output files in that directory
#' 
#' This function reads all tmLQCD gradient flow files in the given path and
#' returns a data frame which concatenates them all.
#' 
#' The single files are expected to be in the tmLQCD format which consists of a
#' header with the column names "traj t P Eplaq Esym tsqEplaq tsqEsym Wsym" and
#' the measurement for each flow time in rows. The columns can be ordered
#' arbitrarily as long as the header and the data are consistent.
#' 
#' @param path the path into which the function should descend
#' @param skip number of measurements to skip.
#' @param basename basename of the files to be read.
#' @param col.names column names of the columns in the files to be read. If not
#' given it will be infered from the files, if possible.
#' @return The function returns a data frame ordered first by the flow time and
#' then by the the trajectory number (so the trajectory number is the index
#' which runs fastest). The data frame has column names \itemize{ \item t -
#' flow time \item traj - trajectory number \item P - plaquette expectation
#' value (at flow time t) \item Eplaq - energy density from plaquette
#' definition (at flow time t) \item Esym - energy density from clover
#' definition (at flow time t) \item tsqEplaq - flow time squared multiplied by
#' plaquette energy density \item tsqEsym - flow time squared multiplied by
#' clover energy density \item Wsym - BMW 'w(t)' observable }.
#' @author Bartosz Kostrzewa, \email{bartosz.kostrzewa@@desy.de}
#' @keywords file
#' @examples
#' 
#' path <- system.file("extdata/", package="hadron")
#' raw.gf <- readgradflow(path)
#' 
#' @export readgradflow
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
  pb <- txtProgressBar(min = 0, max = length(files), style = 3)
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

