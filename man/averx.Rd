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
  \item{pionfit}{
    a fit to the 2pt function. If not provided, it will be performed internally.
  }
  \item{boot.R}{
    number of bootstrap samples.
  }
  \item{boot.l}{
    bootstrap block length.
  }
  \item{useCov}{
    whether or not to use a correlated fit.
  }
  \item{piont1}{
    t1 of fitrange to pion 2pt function
  }
  \item{piont2}{
    t2 of fitrange to pion 2pt function
  }
  \item{seed}{
    seed for PRNG.
  }
  \item{type}{
    type for the effective mass formula, see \link{bootstrap.effectivemass}.
  }
  \item{t1}{
    lower bound for the fitrange in time (t1,t2). Counting starts with 0.
  }
  \item{t2}{
    upper bound for the fitrange in time (t1,t2). Counting starts with 0.
  }
  \item{force.optim}{
    If set to `TRUE`, the usage of \link{optim} will be enforced.
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
