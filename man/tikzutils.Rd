\name{tikzutils}
\alias{tikzutils}
\alias{tikz.init}
\alias{tikz.finalize}
\title{Convenience Functions for \code{tikzDevice}}
\description{
  initialize and finalize a \code{tikzDevice} and carry out optional post-processing
}
\usage{
tikzfiles <- tikz.init(basename="basename",...)
plot(...)
tikz.finalize(tikzfiles=tikzfiles,crop=TRUE,margins=0,clean=TRUE, standAlone=TRUE)
}
\arguments{
  \item{basename}{
    the base of the files which will be used by \code{tikzDevice}, e.g. "basename" -> "basename.pdf", etc.
  }
  \item{...}{
    optional arguments which are passed to \code{tikz}, see \code{\link{tikz}}
  }
  \item{tikzfiles}{
    a list with members $pdf, $tex, $aux and $log, returned by \code{tikz.init}
    which must be passed to \code{tikz.finalize}
  }
  \item{standAlone}{
    A logical value indicating whether the output file should
    be suitable for direct processing by LaTeX. A value of
    \code{FALSE} indicates that the file is intended for inclusion
    in a larger document.
  }
  \item{crop}{
    boolean indicating whether \code{pdfcrop} should be called on the resulting pdf
    ( existence of \code{pdfcrop} is checked before the command is called ),
    default TRUE
  }
  \item{margins}{
    margins argument for pdfcrop command, should be passed as a string consisting
    of one or multiple numbers (e.g. "10" or "10.5 7.5 6.2 10"), default 0
  }
  \item{clean}{
    boolean indicating whether temporary files, e.g. "basename.tex", 
    "basename.aux" and "basename.log" should be deleted after the pdf
    has been generated, default TRUE
  }
}
\value{
  \code{tikz.init} returns a list with character vector members, $pdf, $tex, $aux $log containing
  the corresponding filenames
}
\examples{
library(hadron)
\dontrun{tikzfiles <- tikz.init("plotname",width=3,height=4)
plot(...)}
tikz.finalize(tikzfiles=tikzfiles,clean=FALSE)}
}
\author{Bartosz Kostrzewa, \email{bartosz.kostrzewa@desy.de}}  
\keyword{file}

