\name{extract.loop}
\alias{extract.loop}
\title{Extract a single loop from an object of class \code{cmiloop}}
\description{
  Extracts all loop values from an object of class \code{cmiloop} for
  all available times, samples and configurations.
}
\usage{
extract.loop(cmiloop, obs=9, ind.vec=c(2,3,4,5,6,7,8,1), L)
}
\arguments{
  \item{cmiloop}{
    input object of class \code{cmiloop} generated for instance with
    \code{readcmiloopfiles}. 
  }
  \item{obs}{
    the observable to extract
  }
  \item{ind.vec}{
    index vector to be used during extraction with \code{ind.vec[1]} the
    column with the observable number, \code{ind.vec[2]} the time
    values, \code{ind.vec[3]} the sample numbers, \code{ind.vec[4]} the
    real part of the local loop, \code{ind.vec[5]} the
    imaginary part of the local loop, \code{ind.vec[6]} and
    \code{ind.vec[7]} the same for fuzzed (or smeared) loops and
    \code{ind.vec[8]} for the configuraton number.
  }
  \item{L}{
    The spatial lattice extend needed for normalisation. If not given
    set to \code{T/2}.
  }
}
\value{
  a list with elements as follows:

  \code{cf}: real part of the local loop

  \code{icf}: imaginary part of the local loop
  
  \code{scf}: real part of the smeared loop

  \code{iscf}: imaginary part of the smeared loop

  \code{Time=T}, \code{nrSamples}, \code{nrObs=1}, \code{nrStypes=2},
  \code{obs=obs} and \code{conf.index}. The last is the list of
  configurations corresponding to the loops.

}
\seealso{
  \code{\link{readcmiloopfiles}}
}
\author{Carsten Urbach, \email{curbach@gmx.de}}
