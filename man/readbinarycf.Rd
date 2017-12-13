\name{readbinarycf}
\alias{readbinarycf}
\title{read correlation function from binary files}
\description{
  Reads a correlation function from binary files, including hdf5
  formatted files.
}
\usage{
readbinarycf(files, T=48, obs=5, Nop=1, endian="little", op="aver", 
                   excludelist=c(""), sym=TRUE, path="",
                   hdf5format=FALSE, hdf5name, hdf5index=c(1,2))
}
\arguments{
  \item{files}{
    list of filenames to be read. Can be created using
    \code{getorderedfilelist}. The filelist is assumed to be order
    according to ascending gauge fields.
  }
  \item{T}{
    time extend of correlation functions.
  }
  \item{obs}{
    each file may contain many correlation functions. With
    \code{obs} one choses which observable to read in. To be precise, in
    each file the reading will start at point
    T*obs*sizeof(complex<double>) and read Nop*T*sizeof(complex<double>).
  }
  \item{Nop}{
    number of replicas for the correlator to read in.
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
  \item{op}{
    the N replicas can be either averaged (\code{op="aver"}) or summed
    (\code{op="sum"}).
  }
  \item{path}{
    path to be prepended to every filename.
  }
  \item{hdf5format}{
    if \code{TRUE}, try to read from an hdf5 file.
  }
  \item{hdf5name}{
    Name of the data set as a string.
  }
  \item{hdf5index}{
    The data might be an array of size n x T. \code{hdf5index} is used
    to convert two columns of the data to a complex valued vector using
    the first and second index for real and imaginary part,
    respectively. If \code{hdf5index} has length smaller than 2 the
    first index is reused. 
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
  It is assumend that each file contains at least \code{(obs+N)*T} complex
  doubles, where \code{T} is the time extend, \code{obs} is the number of the
  observable to read in and \code{Nop} the number of replicas for this
  observable. It is assumed that complex is the fastest running index,
  next time and then obs. The filelist is assumed to be ordered according to
  the gauge configuration MC history. 
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


