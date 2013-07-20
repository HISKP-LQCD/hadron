\name{readbinarycf}
\alias{readbinarycf}
\title{read correlation function from binary files}
\description{
  Reads a correlation function from binary files.
}
\usage{
cf <- readbinarycf(files, T=48, obs=5, endian="little",
                         excludelist=c(""), sym=TRUE)
}
\arguments{
  \item{files}{
    list of filenames to be read. Can be created using
    \code{getorderedfilelist}. 
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
  \item{sym}{
    if \code{TRUE} average C(+t) and C(-t), otherwise C(+t) and -C(-t). 
  }
}
\value{
  returns a list with two arrays \code{cf} and \code{icf} with real and
  imaginary parts of the correlator, and integers \code{Time},
  \code{nrStypes=1} and \code{nrObs=1}. Both of the arrays have dimension 
  \code{c(N, (Time/2+1))}, where \code{N} is the number of measurements
  (gauges).  \code{Time} is the time extend, \code{nrStypes} the number
  of smearing levels and \code{nrObs} the number of operators, both of
  which are currently fixed to 1.
}
\details{
  It is assumend that each file contains O*T complex doubles,
  where T is the time extend and O the number of observables in the
  file. It is assumed that complex is the fastest running index, next
  time and then observables.
}
\seealso{
  \code{\link{readcmidatafiles}}, \code{\link{readbinarydisc}},
  \code{\link{readcmidisc}}, \code{\link{readcmicor}}
}
\examples{
\dontrun{cf <- readbinarycf(files, obs=4, excludelist=c("C2_pi0_conf0632.dat"))}
}
\author{Carsten Urbach, \email{curbach@gmx.de}}
\keyword{file}


