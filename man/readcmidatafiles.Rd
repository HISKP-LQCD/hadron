\name{readcmidatafiles}
\alias{readcmidatafiles}
\title{Read Single Data Files in Chris Michael Format}
\description{
  reads data from single files in Chris Michael format
}
\usage{
readcmidatafiles(files, skip=1, verbose=FALSE) 
}
\arguments{
  \item{files}{
    list of filenames to be read. Can be created using
    \code{getorderedfilelist}. 
  }
  \item{skip}{
    Number of lines to be skipped at the beginning of each file
  }
  \item{ind.vec}{
    Index vector indexing the column numbers in cmicor to be used. The
    first must be the observable index, the second the smearing type
    index, the third the time, the fourth C(+t) and the fifth C(-t).
  }
  \item{verbose}{
    Increases verbosity of the function.
  }
  \item{sym.vec}{
    a vector of bools of length equal to the number of observables
    indicating whether C(+t) and C(-t) should be added or subtracted. If
    not given C(+t) and C(-t) will be averaged for all observables.
  }
}
\value{
  returns an object of class \code{cmicor}, which is an \code{rbind} of
  all \code{data.frame}s read from the single files.
}
\details{
  This function reads data from single data files. The data in the files is assumed
  to have 6 columns, specified by \code{ind.vector}: 1) the observable
  type number (itype); 2) the operator 
  type number (iobs); 3) the time difference from source going from 0 to \eqn{T/2}
  for each operator type; 4) \eqn{c_1}{c1} correlator value at time value forward in
  time; 5) \eqn{c_2}{c2} correlator value at time value backward in time;
}
\seealso{
  \code{\link{readcmicor}}, \code{\link{getorderedfilelist}},
    \code{\link{extract.obs}}
}
\examples{
library(hadron)
\dontrun{filelist <- getorderedfilelist("ouptrc")}
\dontrun{cf <- readcmidatafiles(filelist, skip=1)}
}
\author{Carsten Urbach, \email{curbach@gmx.de}}
\keyword{file}
