\name{gevp}
\alias{gevp}
\title{solve GEVP for correlator matrix}
\description{
  solve GEVP for a real, symmetric correlator matrix
}
\usage{
gevp(cf, Time, t0=1, matrix.size=2, element.order=c(1,2,3,4),
for.tsboot=TRUE, sort.type="values")
}
\arguments{
  \item{cf}{
    correlation matrix preferably obtained with a call to
    \code{extrac.obs} (or at leas with the same structure) or an
    already averaged one.

    cf is supposed to be an array of \code{dim=c(N, n*(Time/2+1))}, where
    \code{N} is the number of observations and \code{n} is the number of
    single correlators in the matrix. E.g. for a 2x2 matrix \code{n}
    would be 4.
  }
  \item{Time}{
    time extend of the lattice.
  }
  \item{t0}{
    initial time value of the GEVP, must be in between 0 and
    \code{Time/2-2}. 
  }
  \item{matrix.size}{
    the size of the square matrix
  }
  \item{element.order}{
    specifies how to fit the \code{n} linearly ordered single
    correlators into the correlator
    matrix. \code{element.order=c(1,2,3,4)} leads to a matrix
    \code{matrix(cf[element.order], nrow=2)}.
  }
  \item{sort.type}{
    Sort the eigenvalues either in descending order, or by using the
    scalar product of the eigenvectors with the eigenvectors at
    \eqn{t=t_0+1}{t=t0+1}. Possible values are "values" or "vectors". 
  }
  \item{for.tsboot}{
    for internal use of \link{\code{bootstrap.gevp}}. Alters the returned
    values, see details. 
  }
}
\value{
  Returns a list with the sorted eigenvalues, sorted eigenvectors and
  sorted (reduced) amplitudes for all t > t0.

  In case \code{for.tsboot=TRUE} the same is returned as one long vector
  with first all eigenvalues concatenated, then all eigenvectors and
  then all (reduced) amplitudes concatenated. 
}
\details{
  The generalised eigenvalue problem\cr
  \eqn{
    C(t) v(t,t_0) = C(t_0)\lambda(t,t_0) v(t,t_0)
  }{
    C(t) v(t,t0) = C(t0)lambda(t,t0) v(t,t0)}\cr
  is solved by performing a Cholesky decomposition of \eqn{C(t_0)=L^t
    L}{C(t0)=t(L) L} and transforming the GEVP into a standard
  eigenvalue problem for all values of \eqn{t}. The matrices \eqn{C} are
  symmetrised for all \eqn{t}. So we solve for \eqn{\lambda}{lambda}\cr
  \eqn{(L^t)^{-1} C(t) L^{-1} w = \lambda w}{solve(t(L)) C(t) solve(L) w = lambda w}\cr
  with\cr
  \eqn{w = L v} or the wanted \eqn{v = L^{-1} w}.

  The amplitudes can be computed from\cr
  \eqn{
    A_i^{(n)}(t) = \sum_{j}C_{ij}(t) v_j^{(n)}(t,t_0)/(\sqrt{(v^{(n)},
      Cv^{(n)})(\exp(-mt)\pm \exp(-m(t-t)))})
  }
  and this is what the code returns up to the factor\cr
  \eqn{
    1/\sqrt{\exp(-mt)\pm \exp(-m(t-t))}
  }
}
\references{
  Michael, Christopher and Teasdale, I., Nucl.Phys.B215 (1983) 433, DOI:
  10.1016/0550-3213(83)90674-0\cr
  Blossier, B. et al., JHEP 0904 (2009) 094, DOI:
  10.1088/1126-6708/2009/04/094, arXiv:0902.1265
}
\seealso{
  \code{boostrap.gevp}, \code{extract.obs}
}
\examples{

}
\author{Carsten Urbach, \email{curbach@gmx.de}}
\keyword{GEVP}
