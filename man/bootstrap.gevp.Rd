\name{bootstrap.gevp}
\alias{bootstrap.gevp}
\title{perform a bootstrap analysis of a GEVP}
\description{
  perform a bootstrap analysis of a GEVP for a real, symmetric correlator matrix
}
\usage{
bootstrap.gevp(cf, t0 = 1, element.order = 1:cf$nrObs, sort.type = "vectors", sort.t0 = TRUE)
}
\arguments{
  \item{cf}{
    correlation matrix obtained with a call to \code{extrac.obs}.
  }
  \item{t0}{
    initial time value of the GEVP, must be in between 0 and
    \code{Time/2-2}. Default is 1.
  }
  \item{matrix.size}{
    the size of the square matrix
  }
  \item{element.order}{
    specifies how to fit the \code{n} linearly ordered single
    correlators into the correlator
    matrix. \code{element.order=c(1,2,3,4)} leads to a matrix
    \code{matrix(cf[element.order], nrow=2)}.
    Double indexing is allowed.
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
  \item{sort.type}{
    Sort the eigenvalues either in descending order, or by using the
    scalar product of the eigenvectors with the eigenvectors at
    \eqn{t=t_0+1}{t=t0+1}. Possible values are "values", "vectors" and
    "det". The last one represents a time consuming, but in principle
    better version of sorting by vectors. 
  }
  \item{sort.t0}{
    for \code{sort.type} "vectors" use \eqn{t_0}{t0} as reference or \eqn{t-1}{t-1}.
  }
}
\details{
  Say something on "det" sorting method.
}
\value{
  Returns an object of class \code{gevp} with member objects:

  \code{cf}:\cr
  The input data, if needed bootstrapped with
  \code{\link{bootstrap.cf}}.

  \code{res.gevp}:\cr
  The object returned from the call to \code{\link{gevp}}. For the
  format see \code{\link{gevp}}. 

  \code{gevp.tsboot}:\cr
  The bootstrap samples of the GEVP. For the format see
  \code{\link{gevp}}. 
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
data(correlatormatrix)
## bootstrap the correlator matrix
correlatormatrix <- bootstrap.cf(correlatormatrix, boot.R=400, boot.l=1, seed=132435)
## solve the GEVP
t0 <- 4
correlatormatrix.gevp <- bootstrap.gevp(cf=correlatormatrix, t0=t0, element.order=c(1,2,3,4))
## extract the ground state and plot
pc1 <- gevp2cf(gevp=correlatormatrix.gevp, id=1)
plot(pc1, log="y")
## determine the corresponding effective masses
pc1.effectivemass <- bootstrap.effectivemass(cf=pc1)
pc1.effectivemass <- fit.effectivemass(cf=pc1.effectivemass, t1=5, t2=20)
## summary and plot
summary(pc1.effectivemass)
plot(pc1.effectivemass)

## we can also use matrixfit with a special model for a principal
## correlators
pc1.matrixfit <- matrixfit(pc1, t1=2, t2=24, fit.method="lm", model="pc", useCov=FALSE)
summary(pc1.matrixfit)
plot(pc1.matrixfit)

## the same can be achieved using bootstrap.nlsfit
model <- function(par, x, t0) return(exp(-par[1]*(x-t0))*(par[3]+(1-par[3])*exp(-par[2]*(x-t0))))
ii <- c(2:4, 6:25)
fitres <- parametric.nlsfit(fn=model, par.guess=c(0.5, 1, .9), y=pc1$cf0[ii], dy=pc1$tsboot.se[ii], x=ii-1, boot.R=pc1$boot.R, t0=t0)
summary(fitres)
plot(fitres, log="y")
}
\author{Carsten Urbach, \email{curbach@gmx.de}}
\keyword{GEVP}

