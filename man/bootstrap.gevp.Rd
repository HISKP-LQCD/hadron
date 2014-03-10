\name{bootstrap.gevp}
\alias{bootstrap.gevp}
\title{perform a bootstrap analysis of a GEVP}
\description{
  perform a bootstrap analysis of a GEVP for a real, symmetric correlator matrix
}
\usage{
bootstrap.gevp(cf, t0, boot.R=400, boot.l=2, matrix.size=2, element.order=c(1,2,3,4), seed=1234)
}
\arguments{
  \item{cf}{
    correlation matrix obtained with a call to \code{extrac.obs}.
  }
  \item{t0}{
    initial time value of the GEVP
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
  \item{boot.R}{
    number of bootstrap samples. Used only if cf was not yet
    bootstrapped using \code{bootstrap.cf}.
  }
  \item{boot.l}{
    block size for autocorrelation analysis
  }
  \item{seed}{
    seed for the random number generation used for boostrapping.
  }

}
\value{
  Returns an object of class \code{gevp} with member objects:

  \code{cf}:\cr
  The input data, if needed bootstrapped with
  \link{\code{bootstrap.cf}}.

  \code{res.gevp}:\cr
  The object returned from the call to \link{\code{gevp}}. For the
  format see \link{\code{gevp}}. 

  \code{gevp.tsboot}:\cr
  The bootstrap samples of the GEVP. For the format see
  \link{\code{gevp}}. 
}
\references{
  Michael, Christopher and Teasdale, I., Nucl.Phys.B215 (1983) 433, DOI:
  10.1016/0550-3213(83)90674-0\cr
  Blossier, B. et al., JHEP 0904 (2009) 094, DOI:
  10.1088/1126-6708/2009/04/094, arXiv:0902.1265
}
\seealso{
  \code{gevp}, \code{extract.obs}, \code{bootstrap.cf}
}
\examples{
\dontrun{pion.cor.gevp <- bootstrap.gevp(pion.cor, t0=1)}
\dontrun{## extract the first principal correlator}
\dontrun{pion.pc1 <- gevp2cf(pion.cor.gevp, id=1)}
\dontrun{## which can now be treated like a bootstrapped correlation function}
\dontrun{pion.pc1.effectivemass <- bootstrap.effectivemass(cf=pion.pc1, type="acosh")}
\dontrun{pion.pc1.effectivemass <- fit.effectivemass(pion.pc1.effectivemass, t1=12, t2=23, useCov=FALSE)}
\dontrun{summary(pion.pc1.effectivemass)}
\dontrun{plot(pion.pc1.effectivemass, xlab=c("t/a"), ylab=c("aM"))}
}
\author{Carsten Urbach, \email{curbach@gmx.de}}
\keyword{GEVP}

