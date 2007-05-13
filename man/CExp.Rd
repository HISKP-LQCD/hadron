\name{CExp}
\alias{CExp}
\title{Cosh Or Sinh Build Out Of Two Exps}
\description{
  Evaluates
  \deqn{f(x) = \frac{1}{2}(\exp(-m(T-x))\pm\exp(-m x))}{%
  f(x) = 1/2 (exp(-m(T-x)) +/- exp(-m x))}
  for given mass \eqn{m}, vector \eqn{x} and time extent \eqn{T}. This
  form is better usable in \eqn{\chi^2}{chi^2} fitting than cosh or sinh.
}
\usage{
CExp(m, Time, x, sign=1.)
}
\arguments{
  \item{m}{
    mass value
  }
  \item{Time}{
    Time extent
  }
  \item{x}{
    vector of values on which to evaluate the function
  }
  \item{sign}{
    with sign=1 cosh is evaluated, with sign=-1 sinh
  }
}
\value{
  vector \eqn{f(x)}
}
\author{Carsten Urbach \email{carsten.urbach@liverpool.ac.uk}}
\keyword{math}
