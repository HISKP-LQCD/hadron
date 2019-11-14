\name{alphas}
\alias{alphas}
\title{compute alpha strong at given scale}
\description{compute alpha strong (\eqn{\alpha_s}{alpha_s}) at given
  scale \eqn{\mu}{mu} up to N3LO in PT in the RI' renormalisation
  scheme.
}
\usage{
alphas(mu, nl=3, lam0=0.250, Nc=3., Nf=2., use.cimpl=TRUE)
}
\arguments{
  \item{mu}{the renormalisation scale \eqn{\mu}{mu} in GeV}
  \item{nl}{order in PT, range 0 to 3}
  \item{lam0}{\eqn{\Lambda_\mathrm{QCD}}{Lambda_QCD} in GeV}
  \item{Nc}{number of colours \eqn{N_c}{Nc}, defaults to 3}
  \item{Nf}{number of flavours \eqn{N_f}{Nf}, default is 2}
  \item{use.cimpl}{
    Use the C implementation instead of the R implementation, which
    might improve speed.
  }
}
\value{
  returns the value of alpha strong \eqn{\alpha_s}{alpha_s} at scale
  \eqn{\mu}{mu}
}
%\details{
%}
%\references{
%}
\seealso{
  \code{\link{zetazp}}
}
\examples{
alphas(mu=2.0, nl=3)
}
\author{Carsten Urbach, \email{curbach@gmx.de}}
