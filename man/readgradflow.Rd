\name{readgradflow}
\alias{readgradflow}
\title{Read Gradient Flow Output Files in tmLQCD format}
\description{
  given a pathname, reads all gradient flow output files in that directory
}
\usage{
raw.gf <- readgradflow(path)
}
\arguments{
  \item{path}{
    the path into which the function should descend
  }
}
\value{
  The function returns a data frame ordered first by the flow time
  and then by the the trajectory number (so the trajectory number is the
  index which runs fastest). The data frame has column names
  \itemize{
    \item t - flow time
    \item traj - trajectory number
    \item P - plaquette expectation value (at flow time t)
    \item Eplaq - energy density from plaquette definition (at flow time t)
    \item Esym - energy density from clover definition (at flow time t)
    \item tsqEplaq - flow time squared multiplied by plaquette energy density
    \item tsqEsym - flow time squared multiplied by clover energy density
    \item Wsym - BMW 'w(t)' observable
  }.
}
\details{
  This function reads all tmLQCD gradient flow files in the given path and returns
  a data frame which concatenates them all.

  The single files are expected to be in the tmLQCD format which consists of a header
  with the column names "traj t P Eplaq Esym tsqEplaq tsqEsym Wsym" and the
  measurement for each flow time in rows. The columns can be ordered arbitrarily
  as long as the header and the data are consistent.
}
\examples{
library(hadron)
\dontrun{raw.gf <- readgradflow(path)}
}
\author{Bartosz Kostrzewa, \email{bartosz.kostrzewa@desy.de}}  
\keyword{file}

