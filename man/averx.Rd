\name{averx}
\alias{averx}
\title{Compute <x> From Bare Data}
\description{
  Compute <x> from bare data by fitting to a constant
}
\usage{
averx(data3pt, data2pt, pionfit, boot.R=400, boot.l=2,
      piont1, piont2, useCov=FALSE, t1, t2, seed=123456,
      type="solve", force.optim=FALSE)
}
\arguments{
  \item{data3pt}{
    3pt correlator as read in with \code{\link{read.table}} from a file
    produced with the GWC code from routine momf2e
  }
  \item{data2pt}{
    2pt pseudo-scalar correlator as read in with \code{\link{read.table}} from a file
    produced with the GWC code from routine pssca.
  }
  \item{ind.vec}{
    index vector
  }
  \item{ind.vec2pt}{
    index vector
  }
  \item{skip}{
    number of measurements to be skipped from the beginning
  }
  \item{t1}{
    lower bound for the fitrange in time (t1,t2). Counting starts with 0.
  }
  \item{t2}{
    upper bound for the fitrange in time (t1,t2). Counting starts with 0.
  }
  \item{S}{
    passed to \code{uwerr}, see documentation of \code{\link{uwerr}}.
  }
  \item{mps}{
    value of the pseudo-scalar mass. If not given the value determined
    in averx is used
  }
  \item{par}{
    initial values for the fit of the pseudo-scalar mass. par[1] is the
    amplitude and par[2] is the mass for the cosh fit.
  }
  \item{method}{
    the type of error analysis to be used. Can be either \dQuote{uwerr},
    \dQuote{boot}, \dQuote{all} or \dQuote{no}. For \dQuote{no} (or any
    other string) no error analysis is performed. This might be helpful
    for a first impression and also to test different initial values for
    the fitting parameters. The latter is in particular needed for more
    than one state in the fit.
  }
  \item{nrep}{
    vector (N1, N2, ...) of replica length N1, N2. If missing it is
    assumed that there is only one ensemble. If there are two or more
    replica the parameter \code{skip} has no effect.
  }
}
\value{
  returns a list containing

  \item{averx}{
    value of <x>
  }
  \item{daverx}{
    error of <x>
  }
  \item{data}{
    average 3pt correlator and its error
  }
  \item{fit.uwerr}{
    return object of the uwerr call
  }
  \item{mps}{
    pseudo-scalar mass as given as input or as determined
  }
  \item{N}{
    number of measurements
  }
}
%\details{
%}
%\references{
%}
\seealso{
  \code{\link{uwerr}}, \code{\link{plot.averx}}
}
\examples{
library(hadron)
\dontrun{data3pt <- read.table("momf2e_op_d_0.dat")}
\dontrun{data2pt <- read.table("pssca_corr_0.dat")}
\dontrun{res3pt <- averx(data=data3pt, data2pt=data2pt, t1 = 6, t2 = 18, mps=0.13587)}
}
\author{Carsten Urbach, \email{carsten.urbach@liverpool.ac.uk}}
\keyword{optimize}
\keyword{ts}
