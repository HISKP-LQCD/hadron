\name{readbinarydisc}
\alias{readbinarydisc}
\title{read disconnected loops from binary files}
\description{
  Reads disconnected loops from binary files.
}
\usage{
readbinarydisc(files, T=48, obs=5, endian="little",
               excludelist=c(""), nrSamples, path="")
}
\arguments{
  \item{files}{
    list of filenames to be read. Can be created for instance using 
    \code{getorderedfilelist}. The filelist is assumed to be ordered
    with number of samples running fastest, and the next to fastest
    nubmer of gauges.
  }
  \item{T}{
    time extend of correlation functions.
  }
  \item{obs}{
    each file may contain T*obs correlation functions. With
    \code{obs} one choses which observable to read in.
  }
  \item{endian}{
    the endianess of the binary file.
  }
  \item{excludelist}{
    files to exclude from reading.
  }
  \item{nrSamples}{
    the number of samples
  }
  \item{path}{
    path to be prepended to every filename.
  }
}
\value{
  returns a list with two arrays \code{cf} and \code{icf} with real and
  imaginary parts of the loops, and integers \code{Time},
  \code{nrStypes=1}, \code{nrSamples} and \code{nrObs=1}. Both of the
  arrays have dimension  
  \code{c(T, N)}, where \code{N} is the number of measurements
  (gauges) and \code{T} the time extend.  \code{Time} is the time
  extend, \code{nrStypes} the number 
  of smearing levels and \code{nrObs} the number of operators, both of
  which are currently fixed to 1.
}
\details{
  It is assumend that each file contains O*T complex doubles,
  where T is the time extend and O the number of observables in the
  file. It is assumed that complex is the fastest running index, next
  time and then observables. The different samples are assumend to be in
  different files. The file list is assumed to be ordered with number of
  samples running fastest, and then number of gauges.
}
\seealso{
  \code{\link{readcmidatafiles}}, \code{\link{readbinarycf}},
  \code{\link{readcmidisc}}, \code{\link{readcmicor}}
}
\examples{
\dontrun{files <- character()}
\dontrun{for(i in seq(600,1744,8)) }
\dontrun{  files <- c(files, "C2_dis_u_conf", sprintf("\%.04d", i), ".dat", sep=""))}
\dontrun{cf <- readbinarydisc(files, obs=4, excludelist=c("C2_pi0_conf0632.dat"))}
}
\author{Carsten Urbach, \email{curbach@gmx.de}}
\keyword{file}
