\name{bootstrap.nlsfit}
\alias{bootstrap.nlsfit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  bootstrap an NLS fit
}
\description{
  performs and bootstraps an NLS fit to data with y and x errors.
}
\usage{
bootstrap.nlsfit(fn, par.guess, errormodel = "yerrors", sim = "parametric", boot.R, y, dy, x, dx = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fn}{
    the (non-linear) function to be fitted to the data. Its first
    argument must be the fit parameters named \code{par}. The second
    must be \code{x}, the explaining variable.
}
  \item{par.guess}{
    Initial guess values for the fit parameters.
}
  \item{errormodel}{
    can be either "yerrors" for errors only on the y-axis, or "xyerrors"
    for errors on both axis.
}
  \item{sim}{
    can be "parametric" for a parametric bootstrap approach. Using the
    values \code{dy} (and \code{dx}, if applicable), bootstrap samples
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
    Values of the explaining variable in form of a one-dimensional
    numerical vector.
}
  \item{dx}{
    errors of \code{x} as a one-dimensional numerical vector.
}
  \item{\dots}{
    additional parameters to be passed on to the fit function \code{fn}.
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
  \item{t}{an array of dimensions \code{(npar+1, boot.R+1)} with
    \code{npar} as in \code{t0}. The first
    row is identical to \code{t0}. The other rows contain the individual
    bootstrap observations.}
  \item{bsamples}{the bootstrap samples used as an array of dimensions
    \code{(npar+1, boot.R+1)} with \code{npar} as in \code{t0}. }
  \item{Qval}{the p-value of the fit on the original data}
  \item{chisqr}{the residual chisqr value.}
  \item{dof}{the residual degrees of freedom of the fit.}
  \item{nx}{the number of x-values.}
}
%%\references{
%% ~put references to the literature/web site here ~
%%}
\author{
Carsten Urbach, \email{curbach@gmx.de}
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
fitres <- bootstrap.nlsfit(fn=function(par, x) par[1] + par[2]*x, par.guess=c(1,1), errormodel="xyerrors", boot.R=boot.R, y=value, dy=dvalue, x=x, dx=dx)
summary(fitres)
}
\keyword{ models }
