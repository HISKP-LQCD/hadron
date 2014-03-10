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
    An object returned by \link{\code{bootstrap.gevp}}.
  }
  \item{id}{
    The index of the principal correlator to extract.
  }
}
\value{
  An object of class \code{cf}, which contains bootstrap samples
  already. So a call to \code{bootstrap.cf} is neither needed nor
  possible. It can be treated further by
  \link{\code{bootstrap.effectivemass}} or \link{\code{matrixfit}} to
  extract a mass value.
}
%\details{
%}
%\references{
%}
\seealso{
  \link{\code{gevp}}, \link{\code{matrixfit}}, \link{\code{bootstrap.effectivemass}}
}
\examples{
\dontrun{## apply a GEVP analysis}
\dontrun{pion.cor.gevp <- bootstrap.gevp(pion.cor, t0=1)}
\dontrun{## extract the first principal correlator}
\dontrun{pion.pc1 <- gevp2cf(pion.cor.gevp, id=1)}
}
\author{Carsten Urbach, \email{curbach@gmx.de}}
\keyword{GEVP}

