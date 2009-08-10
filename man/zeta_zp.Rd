\name{zetazp}
\alias{zetazp}
\title{Computes the running of Z_P from scale mu0 to scale
  mu2}
\description{
  Computes the running of the renomalisation constant \eqn{Z_P} from
  scale \eqn{\mu_0}{mu0} to scale \eqn{\mu_2}{mu2} in the renomalisation
  schema RI' for \eqn{N_f=2}{Nf=2} only. The running is done using
  perturbation theory up to \eqn{\alpha_s**3}{alpha_s^3} order. The
  corresponding values of \eqn{\alpha_s}{alpha_s} at the scales
  \eqn{\mu_0}{mu0} and \eqn{\mu_2}{mu2} are needed as input, see
  \code{\link{alphas}}.
  }
\usage{
zetazp(zp0, alpha0, alpha2, nl = 3)
}
\arguments{
  \item{zp0}{initial value of \eqn{Z_P}}
  \item{alpha0}{\eqn{\alpha_s}{alpha_s} at initial scale}
  \item{alpha2}{\eqn{\alpha_s}{alpha_s} at final scale}
  \item{nl}{order in PT, range 0 to 3}
}
\value{
  returns the value of Z_P at scale mu2 in the RI' scheme
}
%\details{
%  
%}
%\references{
%}
\seealso{
  \code{\link{alphas}}
}
\examples{
al2 <- alphas(mu = 3.0, nl = 3, lam0 = 0.250, Nc = 3, Nf = 2)
al0 <- alphas(mu = 2.0, nl = 3, lam0 = 0.250, Nc = 3, Nf = 2)
zetazp(zp0 = 0.6, alpha0 = al0, alpha2 = al2, nl = 3)
}
\author{Carsten Urbach, \email{curbach@gmx.de}}
