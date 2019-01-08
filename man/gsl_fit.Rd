\name{gsl_fit_correlator_matrix}
\alias{gsl_fit_correlator_matrix}
\title{Fit a correlator matrix}
\description{Fit a correlator matrix
}
\usage{gsl_fit_correlator_matrix(par, Thalf, x, y, err, tr, N,
no_masses=1, prec=c(1.e-10,1.e-4))
}
\arguments{
  \item{par}{initial parameter list}
  \item{Thalf}{half time extent of the lattice}
  \item{x}{explaining variable}
  \item{y}{data points}
  \item{err}{errors on data points}
  \item{tr}{length of \code{x}}
  \item{N}{matrix size. Can be 2, 4, or 6.}
  \item{no_masses}{the number of masses to fit}
  \item{prec}{vector of precision for the fit routine. First entry is
    the absolute precision, the second relative one}
}
\value{
  returns a list with entries:
  \item{par}{
    The set best fit parameters found.
  }
  \item{value}{
    The value of chi square corresponding to the best fit parameters
  }
  \item{convergence}{
    An integer code. '0' indicates successful convergence. Error codes
    are the once from Gnu Scientific Library.
  }
  \item{counts}{
    An integer giving the number of calls to the function to be
    minimised. 
  }
  \item{dof}{
    The degrees of freedom of the fit.
  }
}
\details{
  This function calls the multifit solver of the Gnu Scientific
  Library (GSL). The function to be minimised is specific for the functions
  \code{\link{pion}} and \code{rho}. For details see
  \code{\link{pion}} and \code{rho}. 
}
\references{
  See GSL webpage \code{http://www.gnu.org/software/gsl/}
}
\seealso{
  \code{\link{pion}}, \code{\link[stats]{optim}}
}
%\examples{
%
%}
\author{Carsten Urbach, \email{carsten.urbach@physik.hu-berlin.de}}
\keyword{optimize}
