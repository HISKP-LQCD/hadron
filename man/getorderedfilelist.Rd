\name{getorderedfilelist}
\alias{getorderedfilelist}
\alias{getorderedconfignumbers}
\title{Creates an ordered filelist from a basename and a path}
\description{
  These functions generate an ordered filelist and an order list of
  config numbers by using a path and a
  basename and '*'.
}
\usage{
getorderedfilelist(path="./", basename="onlinemeas",
last.digits=4, ending="")

getorderedconfignumbers(path="./", basename="onlinemeas",
last.digits=4, ending="")
}
\arguments{
  \item{path}{
    the path to be searched
  }
  \item{basename}{
    the basename of the files
  }
  \item{last.digits}{
    the number of last characters in each filename to be used for
    ordering the list.
  }
  \item{ending}{
    the file extension after the digits.
  }
}
\value{
  returns the ordered list of strings.
}
\details{
  All filenames are assumend to have equal length.
}
\seealso{
  \code{\link{readcmidatafiles}}, \code{\link{extract.obs}}
}
\examples{
\dontrun{filelist <- getorderedfilelist("ouptrc")}
}
\author{Carsten Urbach, \email{curbach@gmx.de}}
\keyword{file}
