\name{readcmidisc}
\alias{readcmidisc}
\title{reads disconnected loops in cmi format}
\description{
  reads disconnected loops in cmi (Chris Michael) format from a list of files.
}
\usage{
readcmidisc(files, obs=9, ind.vec=c(2,3,4,5,6,7,8),
            excludelist=c(""), skip=0, L,
            colClasses=c("integer", "integer","integer","integer",
                          "numeric","numeric","numeric","numeric"),
            debug=FALSE)

}
\arguments{
  \item{files}{
    list of filenames to be read. Can be created using
    \code{getorderedfilelist}. 
  }
  \item{obs}{
    index of operator to parse from files
  }
  \item{ind.vec}{
    vector containing the index (column in file) of obs, t, samples,
    Re(local), Im(local, Re(smeared), Im(smeared).
  }
  \item{excludelist}{
    files to exclude from reading.
  }
  \item{skip}{
    lines to skip at beginning of each file.
  }
  \item{colClasses}{
    The column data type classes, the \code{read.table}.
  }
  \item{L}{
    the spatial lattice extend, set to \code{Time/2} if missing.
  }
  \item{debug}{
    setting debug to TRUE makes the routine more verbose by spilling out
    separate filenames.
  }
}
\value{
  returns a list with four arrays \code{cf}, \code{icf} \code{scf} and
  \code{sicf} containing real and imaginary parts of the local and
  smeared loops, respectively, and integers \code{Time}, 
  \code{nrStypes=2}, \code{nrSamples} and \code{nrObs=1}. The four
  arrays have dimension  
  \code{c(T, S, N)}, where \code{S} is the nubmer of samples, \code{T}
  is the time extend and \code{N} is the number of measurements
  (gauges).  \code{Time} is the time extend, \code{nrStypes} the number
  of smearing levels and \code{nrObs} the number of operators, which are
  currently fixed to 1 and 2, respectively. \code{nrSamples} is the number
  of samples.

  Note that the arrays are normalised by \code{1/sqrt(L^2)}.

  The routine expects that all files have identical content. Otherwise
  the routine will stop.
}
\seealso{
  \code{\link{readcmidatafiles}}, \code{\link{readbinarycf}},
  \code{\link{readbinarydisc}}, \code{\link{readcmicor}}
}
\examples{
\dontrun{v4files <- character()}
\dontrun{for(i in seq(600,1744,8)) v4files <- c(v4files, paste("disc.0.163265.0.006.k0v4.", sprintf("\%.04d", i), sep=""))}
\dontrun{v4data <- readcmidisc(v4files)}

}
\author{Carsten Urbach, \email{curbach@gmx.de}}
\keyword{file}
