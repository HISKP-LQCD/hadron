\name{variational}
\alias{variational}
\title{variational analysis of a correlator matrix}
\description{
  variational analysis of a correlator matrix
}
\usage{
variational(Cor, ta, tb, tmax, N, T1, matrix.size, no.masses) 
}
\arguments{
  \item{Cor}{
    correlator matrix
  }
  \item{ta}{
    first time value for eigenvector determination
  }
  \item{tb}{
    second time value for eigenvector determination
  }
  \item{tmax}{
    maximal time value to be considered in the analysis
  }
  \item{N}{
    matrix size. must be <= to the size of the correlator matrix
  }
  \item{T1}{
    maximal time value in correlator matrix
  }
  \item{matrix.size}{
    matrix size to be used for start value determination (e.g. for
    pionfit) 
  }
  \item{no.masses}{
    number of mass values for start value determination
  }
}
\value{
  returns a list with following entries
  \item{t}{
    the list of time values 
  }
  \item{res.values}{
    the mass values in a \code{tmax-ta} times \code{N} array for all
    analysed t-values
  }
  \item{par}{
    startvalues for a fit
  }
  \item{variational.masses}{
    the list of mass values as determined from the variational analysis 
  }
}
%\references{
%}
\seealso{
  pion
}
%\examples{

%}
\author{Carsten Urbach, \email{carsten.urbach@physik.hu-berlin.de}}
\keyword{variational}

