\name{onlinemeas}
\alias{onlinemeas}
\title{determines pion mass and pcac mass from online measured
  correlator of the HMC code}
\description{
  determines pion mass and pcac mass from online measured
  correlator of the HMC code
}
\usage{
onlinemeas(data, t1, t2, S=1.5, pl=FALSE, skip=0,
           iobs=1, ind.vec=c(1,3,4,5),
           boot.R=99, boot.l=10, tsboot.sim="geom",
           method="uwerr", fit.routine="optim", nrep,
           oldnorm = FALSE)

}
\arguments{
  \item{data}{
    data to be fitted to as e.g. the output of
    \code{\link{readcmicor}}. Currently only
    \code{\link{cmicor}} format is supported.
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
  \item{pl}{
    logical: if set to TRUE the function produces plots
  }
  \item{skip}{
    number of measurements to be discarded at the beginning of the time
    series. \code{skip} has no effect if two or more replica are used,
    see argument \code{nrep}.
  }
  \item{iobs}{
    if there are several operators available (local-local,
    local-smeared, etc.), then this labels these (for cmi format)
  }
  \item{ind.vec}{
    index vector indexing the column numbers in cmicor to be used
  }
  \item{boot.R}{
    number of bootstrap samples for bootstrap analysis
  }
  \item{boot.l}{
    average block size for blocking analysis with tsboot
  }
  \item{tsboot.sim}{
    The type of simulation required to generate the replicate
    time series. See \code{\link{tsboot}} for details.
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
  \item{fit.routine}{
    The fit routine to be used. Default is \dQuote{gsl}, which uses the
    gnu scientific library \dQuote{gsl_multifit_fdfsolver} solver to
    minimise the chisquare. All other values lead to the usage of R's
    \link{optim} function. The latter choice might be significantly
    slower.
  }
  \item{oldnorm}{
    If set to \dQuote{TRUE}, the old online measurement normalisation of
    \dQuote{tmLQCD} prior to version 5.2.0 is used in order to get correct
    values for the pion decay constant.
  }
}
\value{
  returns an object of \code{class} \code{ofit} with the following items

  \item{fitresult}{
    result from the fit as returned by \code{\link{optim}}
  }
  \item{fitresultpp}{
    Fit result of the PP correlator only
  }
  \item{t1}{
    lower bound for the fitrange in time (t1,t2). Counting starts with 0.
  }
  \item{t2}{
    upper bound for the fitrange in time (t1,t2). Counting starts with 0.
  }
  \item{N}{
    number of measurements found in the data
  }
  \item{Time}{
    Time extent found in the data
  }
  \item{fitdata}{
    \code{\link{data.frame}} containing the time values used in the fit,
    the averaged correlator and its error and the value of Chi for each
    time value
  }
  \item{uwerrresultmps}{
    the result of the time series analysis for the lowest mass as
    carried out by \code{\link{uwerr}} 
  }
  \item{uwerrresultmpcac}{
    the result of the time series analysis for the PCAC mass
    carried out by \code{\link{uwerr}}, see details
  }
  \item{effmass}{
    effective masses in the pion channel
  }
  \item{matrix.size}{
    size of the data matrix, copied from input
  }
  \item{boot}{
    object returned by the call to \code{\link{boot}} if \code{method}
    was set correspodingly. Otherwise \code{NULL}.
  }
  \item{tsboot}{
    object returned by the call to \code{\link{tsboot}} if \code{method}
    was set correspodingly. Otherwise \code{NULL}.
  }
  \item{method}{
    error analysis method as copied from input
  }
  \item{fit.routine}{
    \code{fit.routine} as copied from input
  }
  \item{nrep}{
    \code{nrep} as copied from input
  }
  \item{dpaopp}{
    \code{\link{data.frame}} containing the pcac masses computed not
    with a fit, but with the derivative method for all time values in
    between \code{t1} and \code{t2}
  }
}
\details{
  The online measurements in the HMC code compute the PP and PA
  correlation functions summed over spatial x for all t. We analyse
  these correlators in different ways:

  First, only the PP correlator is analysed and fitted by
  \eqn{p_1^2\cosh(-m(t-T/2))}{p1*p1*cosh(-m(t-T/2))} for \eqn{m} and
  \eqn{p_1}{p1}.
  
  Second, PP and PA correlators are fitted together with three parameters
  as \eqn{C_\mathrm{PP} = p_1^2\cosh(-m(t-T/2))}{C_PP =
    p1*p1*cosh(-m(t-T/2))} and \eqn{C_\mathrm{PA} =
    p_1p_2\cosh(-m(t-T/2))}{C_PA = p1*p2*cosh(-m(t-T/2))} in a
  simultaneous fit. \eqn{m} is then the pseudo scalar mass and the pcac
  mass is determined from
  \deqn{m_\mathrm{PCAC} = m_\mathrm{PS} \frac{p_2}{2p_1}}{%
    mpcac = mps p_2/(2 p_1)}

  Finally, the PCAC mass can also be determined computing
  \deqn{m_\mathrm{PCAC}(t) =%
    \frac{C_\mathrm{PA}(t+1)-C_\mathrm{PA}(t-1)}{4C_\mathrm{PP}(t)}}{%
    C_PA(t) = (C_PA(t+1)-C_PA(t-1))/(4 C_PP(t))}
  using the symmetric finite difference operator.
}
%\references{
%}
\seealso{
  \code{\link{readcmicor}}, \code{\link{uwerr}},
  \code{\link{variational}}
}
%\examples{
%}
\author{Carsten Urbach, \email{curbach@gmx.de}}
\keyword{optimize}
\keyword{ts}
