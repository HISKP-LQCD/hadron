\name{bootstrap.nlsfit}
\alias{bootstrap.nlsfit}
\alias{summary.bootstrapfit}
\alias{print.bootstrapfit}
\alias{plot.bootstrapfit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  bootstrap an NLS fit
}
\description{
  Performs and bootstraps a non-linear least-squares fit to data with y and x errors.
}
\usage{
bootstrap.nlsfit(fn, gr=NULL, dfn=NULL, par.guess, errormodel = "yerrors", sim =
"parametric", boot.R, y, dy, x, dx = NULL, bsamples, useCov=FALSE, CovMatrix, use.minpack.lm=TRUE, parallel=FALSE, ...)

summary(object, digits=2, ...)

print(x, digits=2, ...)

plot(x, ..., xlim, ylim, rep=FALSE, col.line="black",
col.band="gray", lty=c(1), lwd=c(1), xlab="x", ylab="y", supports=1000)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fn}{
    the (non-linear) function to be fitted to the data. Its first
    argument must be the fit parameters named \code{par}. The second
    must be \code{x}, the explaining variable.
  }
  \item{gr}{
    \code{gr=d(fn) / d(par)} is a function to return the gradient of \code{fn}. It must return an
    array with \code{length(x)} rows and \code{length(par)} columns.
  }
  \item{dfn}{
    \code{dfn=d(fn) / dx} is the canonical derivative of \code{fn} by \code{x} and only relevant
    if x-errors are provided.
  }
  \item{par.guess}{
    initial guess values for the fit parameters.
  }
  \item{errormodel}{
    can be either "yerrors" for errors only on the y-axis, or "xyerrors"
    for errors on both axes.
  }
  \item{sim}{
    can be "parametric" for a parametric bootstrap approach. Using the
    values \code{dy} (and \code{dx} or \code{CovMatrix}, if applicable), bootstrap samples
    are generated assuming a Gaussian error distribution.
  }
  \item{boot.R}{
    the number of bootstrap samples to use.
  }
  \item{y}{
    the data as a one-dimensional numerical vector to be described by
    the fit function. 
  }
  \item{dy}{
    errors of \code{y} as a one-dimensional numerical vector.
  }
  \item{x}{
    values of the explaining variable in form of a one-dimensional
    numerical vector.
  }
  \item{dx}{
    errors of \code{x} as a one-dimensional numerical vector.
  }
  \item{bsamples}{
    bootstrap samples of \code{y} (and \code{x}, if applicable). Must be
    provided as array of dimensions \code{c(boot.R, n)} with
    \code{n} equals to \code{length(y)} in case of 'yerrors'
    and For 'xyerrors' to \code{length(y) + length(x)}.
  }
  \item{useCov}{
    perform a fully correlated fit by deterining the variance-covariance
    matrix from the provided bootstrap samples or using the covariance matrix \code{CovMatrix} directly.
  }
  \item{CovMatrix}{
    complete variance-covariance matrix of dimensions \code{c(length(y), length(y))}
    or \code{c(length(y)+length(x), length(y)+length(x))} depending on the errormodel.
  }
  \item{use.minpack.lm}{
    use the \code{minpack.lm} library if available. This is usually faster
    than the default \code{optim} but somtimes also less stable.
  }
  \item{parallel}{
    parallelise over bootstrap samples.
  }
  \item{\dots}{
    additional parameters to be passed on to the fit function \code{fn}.
  }
  \item{object}{
    object returned by \code{bootstrap.nlsfit}
  }
  \item{digits}{
    number of significant digits to print in summary or print.
  }
  \item{xlim}{
    x limits of the plot.
  }
  \item{ylim}{
    y limits of the plot.
  }
  \item{rep}{
    If set to \code{TRUE}, operate like "replot" in gnuplot. Allows
    adding points with error bars to the current plot. Switches the underlying
    plotting routine from \code{plot} to \code{points}.
  }
  \item{col.line}{
    line colour.
  }
  \item{col.band}{
    error band colour.
  }
  \item{opacity.band}{
    error band opacity.
  }
  \item{lwd}{
    line width for fitted curve.
  }
  \item{lty}{
    line type of fitted curve.
  }
  \item{supports}{
    number of supporting points for plotting the function.
  }
  \item{...}{
    arguments to be passed on to \code{fn} and \code{gr} in case of
    \code{bootstrap.nlsfit} or additional general parameters passed on
    to \code{summary}, \code{print} or \code{plot}.
  }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  returns a list of class 'bootstrapfit'. It returns all input
  parameters and adds in addition the following:
  \item{t0}{the one dimensional numerical vector of length
    \code{npar+1}. \code{npar} is the number of fit parameters. In case
    of 'yerrors' this equals \code{length(par.guess)}. For 'xyerrors'
    this equals \code{length(par.guess) + length(x)}. \code{t0} contains
    the best fit parameters
    obtained on the original data. The last element in \code{t0} is the
    chisquare value.}
  \item{t}{an array of dimensions \code{(npar+1, boot.R)} with
    \code{npar} as in \code{t0}. The rows contain the individual
    bootstrap observations.}
  \item{bsamples}{the bootstrap samples used as an array of dimensions
    \code{(length(y), boot.R)} or \code{(length(y)+length(x), boot.R)}
    depending on the error model with \code{npar} as in \code{t0}. }
  \item{Qval}{the p-value of the fit on the original data}
  \item{chisqr}{the residual chisqr value.}
  \item{dof}{the residual degrees of freedom of the fit.}
  \item{nx}{the number of x-values.}
  \item{tofn}{
    the original \code{...} list of parameters to be passed on to the
    fit function
  }
}
%%\references{
%% ~put references to the literature/web site here ~
%%}
\author{
Carsten Urbach, \email{curbach@gmx.de}; Johann Ostmeyer, \email{ostmeyer@hiskp.uni-bonn.de}
}
%%\note{
%%  ~~further notes~~
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{optim}}
}
\examples{
value <- c(0.1, 0.2, 0.3)
dvalue <- c(0.01, 0.01, 0.015)
x <- c(1,2,3)
dx <- c(0.1, 0.1, 0.1)
boot.R <- 1500
## with xy-errors and no correlated fit
fitres <- bootstrap.nlsfit(fn=function(par, x) par[1] + par[2]*x, par.guess=c(1,1), errormodel="xyerrors", boot.R=boot.R, y=value, dy=dvalue, x=x, dx=dx, useCov=FALSE)
summary(fitres)
## with covariance matrix and y-errors
bootstrapsamples <- array(NA, dim=c(boot.R+1, length(value)))
bootstrapsamples[1,] <- value
## some indices
crr <- c(1:(boot.R+1))
rr <- c(2:(boot.R+1))
bootstrapsamples[rr, 1] <- rnorm(n=boot.R, mean = value[1], sd = dvalue[1])
bootstrapsamples[rr, 2] <- rnorm(n=boot.R, mean = value[2], sd = dvalue[2])
bootstrapsamples[rr, 3] <- rnorm(n=boot.R, mean = value[3], sd = dvalue[3])

## the fit without correlation and y-errors only
fitres <- bootstrap.nlsfit(fn=function(par, x) par[1] + par[2]*x, par.guess=c(1,1), errormodel="yerrors", boot.R=boot.R, y=value, dy=dvalue, x=x, useCov=TRUE, bsamples=bootstrapsamples[rr,])
summary(fitres)
}
\keyword{ models }
