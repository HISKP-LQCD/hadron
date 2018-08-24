\name{gevp2cf}
\alias{gevp2cf}
\title{Extracts a principle correlator from a GEVEP}
\description{
  Extracts a principle correlator from a GEVP and converts it into an
  object of class \code{cf}
}
\usage{
gevp2cf(gevp, id=1)
}
\arguments{
  \item{gevp}{
    An object returned by \code{\link{bootstrap.gevp}}.
  }
  \item{id}{
    The index of the principal correlator to extract.
  }
}
\value{
  An object of class \code{cf}, which contains bootstrap samples
  already. So a call to \code{bootstrap.cf} is neither needed nor
  possible. It can be treated further by
  \code{\link{bootstrap.effectivemass}} or \code{\link{matrixfit}} to
  extract a mass value.
}
%\details{
%}
%\references{
%}
\seealso{
  \code{\link{gevp}}, \code{\link{matrixfit}}, \code{\link{bootstrap.effectivemass}}
}
\examples{
data(correlatormatrix)
## bootstrap the correlator matrix
correlatormatrix <- bootstrap.cf(correlatormatrix, boot.R=400, boot.l=1, seed=132435)
## solve the GEVP
t0 <- 4
correlatormatrix.gevp <- bootstrap.gevp(cf=correlatormatrix, t0=t0, element.order=c(1,2,3,4))
## extract the ground state and plot
pc1 <- gevp2cf(gevp=correlatormatrix.gevp, id=1)
plot(pc1, log="y")
}
\author{Carsten Urbach, \email{curbach@gmx.de}}
\keyword{GEVP}

