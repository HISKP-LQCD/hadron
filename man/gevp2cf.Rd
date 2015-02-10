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
\dontrun{## apply a GEVP analysis}
\dontrun{pion.cor.gevp <- bootstrap.gevp(pion.cor, t0=1)}
\dontrun{## extract the first principal correlator}
\dontrun{pion.pc1 <- gevp2cf(pion.cor.gevp, id=1)}
}
\author{Carsten Urbach, \email{curbach@gmx.de}}
\keyword{GEVP}

