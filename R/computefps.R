#' Computes the pseudoscalar decay constant for the twisted mass case from the
#' pseudoscalar amplitude and mass
#' 
#' From a mass and amplitude determination (using \code{\link{matrixfit}} or
#' \code{\link{fit.effectivemass}}, \code{\link{bootstrap.gevp}} and
#' \code{\link{gevp2amplitude}} the pseudoscalar decay constant is determined
#' for the case of Wilson twisted mass fermions from the pseudoscalar amplitude
#' and mass
#' 
#' The pseudoscalar decay constant is computed from\cr \deqn{f_\mathrm{PS} =
#' 2\kappa(\mu_1+\mu_2)\frac{PP}{\sqrt{2}\sqrt{m_\mathrm{PS}}^3}}{% fps = 2
#' kappa (mu1+mu2) PP/sqrt(2)/sqrt(mps)^3} for \code{normalisation="cmi"} or
#' \deqn{f_\mathrm{PS} =
#' (\mu_1+\mu_2)\frac{PP}{\sqrt{2}\sqrt{m_\mathrm{PS}}^3}}{% fps = (mu1+mu2)
#' PP/sqrt(2)/sqrt(mps)^3} expecting physical normalisation of the
#' amplitudes.\cr When \code{disprel="lattice"},\cr
#' \deqn{\sqrt{m_{\mathrm{PS}}^3}}{% sqrt(mps^3)} is replaced with
#' \deqn{\sqrt{m_{\mathrm{PS}}} \sinh{m_{\mathrm{PS}}}}{% sqrt(mps)*sinh(mps))}
#' which can reduce lattice artefacts for heavy meson masses.
#' 
#' @param mfit An object of type \code{matrixfit} or \code{gevp.amplitude}
#' generated with \code{\link{matrixfit}} or \code{\link{gevp2amplitude}},
#' respectively.
#' @param PP If \code{mfit} is missing this must contain the value for the
#' pseudoscalar amplitude.
#' @param mass If \code{mfit} is missing this must contain the value for the
#' pseudoscalar mass.
#' @param mu1,mu2 The values for the twisted quark masses involved in the
#' pseudoscalar meson. If \code{mu2} is missing it will be assumed to be equal
#' to \code{mu1}.
#' @param Kappa The \eqn{\kappa}{kappa}-value of the run, needed only if
#' \code{normalisation="cmi"}.
#' @param normalisation normalisation of the correlators. If set to "cmi" the
#' \eqn{\kappa}{kappa} value must be specified.
#' @param disprel One of "continuum" or "lattice". Indicates whether the
#' formula for the decay constant should take into account the lattice
#' dispersion relation for the meson. Theoretically this can reduce lattice
#' artefacts for heavy mesons.
#' @param boot.fit If set to \code{FALSE}, the computation is not bootstrapped,
#' even if the \code{matrixfit} or \code{gevp.amplitude} contain bootstrap
#' samples.  This is a useful time-saver if error information is not strictly
#' necessary.  Of course, this affects the return values related to the
#' bootstrap, which are set to \code{NA}.
#' @return If \code{mfit} ist missing the value of fps will printed to stdout
#' and returned as a simple numerical value.
#' 
#' If \code{mfit} is available, this object will be returned but with
#' additional objects added: \code{fps}, \code{fps.tsboot}, \code{mu1,mu2},
#' \code{normalistaion} and \code{Kappa} if applicable.
#' @author Carsten Urbach, \email{curbach@@gmx.de}
#' @seealso \code{\link{matrixfit}}, \code{\link{gevp2amplitude}},
#' @keywords GEVP optimise ts
#'
#' @examples
#'
#' cfnew <- extractSingleCor.cf(correlatormatrix, id=1)
#' cfnew <- bootstrap.cf(cfnew, boot.R=99, boot.l=1)
#' cfnew.fit <- matrixfit(cf=cfnew, t1=12, t2=20, parlist=array(c(1,1),
#'                        dim=c(2,1)), sym.vec=c("cosh"), neg.vec=c(1))
#' cfnew.fps <- computefps(mfit=cfnew.fit, mu1=0.004, normalisation="new")
#' summary(cfnew.fps)
#' 
#' @export computefps
computefps <- function(mfit, PP, mass, mu1, mu2, Kappa, normalisation="cmi", disprel="continuum", boot.fit=TRUE) {
  if(missing(mu1)) {
    stop("computefps: mu1  must be specified! Aborting...\n")
  }
  if(missing(Kappa)) {
    Kappa <- NA
  }
  if(missing(mu2)) mu2 <- mu1
  if(missing(mfit)) {
    denominator <- sqrt(mass^3)
    if(disprel == "lattice") denominator <- sinh(mass)*sqrt(mass)
    
    if(normalisation != "cmi") return((mu1+mu2)*abs(PP)/denominator)
    else return(sqrt(2)*Kappa*(mu1+mu2)*abs(PP)/denominator)
  }
  else if(any(class(mfit) == "matrixfit")) {
    k <- Kappa
    if(normalisation != "cmi") k <- sqrt(0.5)
    else  mfit$kappa <- Kappa
    
    denominator <- sqrt(mfit$opt.res$par[1]^3)
    if(disprel == "lattice") denominator <- sinh(mfit$opt.res$par[1])*sqrt(mfit$opt.res$par[1])
    mfit$fps <- 2*k*(mu1+mu2)/sqrt(2)*abs(mfit$opt.res$par[2])/denominator

    if(boot.fit) {
      denominator <- sqrt(mfit$opt.tsboot[1,]^3)
      if(disprel == "lattice") denominator <- sinh(mfit$opt.tsboot[1,])*sqrt(mfit$opt.tsboot[1,])
      mfit$fps.tsboot <- sqrt(2.)*k*(mu1+mu2)*abs(mfit$opt.tsboot[2,])/denominator
    } else {
      mfit$fps.tsboot <- NA
    }

    mfit$mu1 <- mu1
    mfit$mu2 <- mu2
    mfit$normalisation <- normalisation
    mfit$disprel <- disprel
    return(invisible(mfit))
  }
  else if(inherits(mfit, "gevp.amplitude")) {
    k <- Kappa
    if(normalisation != "cmi") k <- sqrt(0.5)
    else  mfit$kappa <- Kappa

    denominator <- sqrt(mfit$m0^3)
    if(disprel == "lattice") denominator <- sinh(mfit$m0)*sqrt(mfit$m0)
    mfit$fps <- sqrt(2)*k*(mu1+mu2)*abs(mfit$meanAmplitude)/denominator

    if(boot.fit) {
      denominator <- sqrt(mfit$m0.tsboot^3)
      if(disprel == "lattice") denominator <- sinh(mfit$m0.tsboot)*sqrt(mfit$m0.tsboot)
      mfit$fps.tsboot <- sqrt(2)*k*(mu1+mu2)*abs(mfit$meanAmplitude.tsboot[,1])/denominator
    } else {
      mfit$fps.tsboot <- NA
    }
    
    mfit$mu1 <- mu1
    mfit$mu2 <- mu2
    mfit$normalisation <- normalisation
    mfit$disprel <- disprel
    return(invisible(mfit))
  }
  else{
    stop("computefps: expecting object of class matrixfit, gevp.amplitude or at least an amplitude\nAborting...!\n")
  }
}



#' Computes the pseudoscalar decay constant for the Osterwalder Seiler case
#' from the pseudoscalar amplitude and mass
#' 
#' From a mass and amplitude determination (using \code{\link{matrixfit}}) the
#' pseudoscalar decay constant is determined for the case of Osterwalder Seiler
#' (OS) fermions from the AS and SS amplitude (in the twisted basis), ZA and
#' the OS pion mass.
#' 
#' The pseudoscalar decay constant is computed from\cr
#' \deqn{f_\mathrm{PS}^\mathrm{OS} = Z_A \sqrt{2}\kappa\frac{\langle 0|
#' A|\pi\rangle}{m_\mathrm{PS}}}{% fpsOS = sqrt(2) kappa ZA <0|A|pi>/mps } for
#' \code{normalisation="cmi"} or \deqn{f_\mathrm{PS}^\mathrm{OS} = Z_A
#' \frac{\langle 0| A|\pi\rangle}{m_\mathrm{PS}}}{% fpsOS = ZA <0|A|pi>/mps}
#' expecting physical normalisation of the amplitudes.\cr
#' 
#' @param mfit An object of type \code{matrixfit} generated with
#' \code{\link{matrixfit}}. The correlation matrix (SS, SA, AS, AA) must have
#' been analysed, where the correlators are in the twisted basis.
#' @param ZA The value of the renormalisation constant \eqn{Z_A}{ZA}.
#' @param dZA The value of the (normally distributed) error of the
#' renormalisation constant \eqn{Z_A}{ZA}.
#' @param ZAboot Bootstrap samples for \eqn{Z_A}{ZA}. If they are provided,
#' they are used for computing fps, if not, bootstrap samples are generated
#' from \code{dZA}. If both are missing, the error of \eqn{Z_A}{ZA} is not
#' taken into account.
#' @param Kappa The \eqn{\kappa}{kappa}-value of the run, needed only if
#' \code{normalisation="cmi"}.
#' @param normalisation normalisation of the correlators. If set to "cmi" the
#' \eqn{\kappa}{kappa} value must be specified.
#' @param boot.fit If set to \code{FALSE}, the computation is not bootstrapped,
#' even if the \code{matrixfit} or \code{gevp.amplitude} contain bootstrap
#' samples.  This is a useful time-saver if error information is not strictly
#' necessary.  Of course, this affects the return values related to the
#' bootstrap, which are set to \code{NA}.
#' @return If \code{mfit} is available, this object will be returned but with
#' additional objects added: \code{fpsOS}, \code{fpsOS.tsboot},
#' \code{normalistaion}, \code{ZA}, \code{ZAboot} and \code{kappa} if
#' applicable.
#' @author Carsten Urbach, \email{curbach@@gmx.de}
#' @seealso \code{\link{matrixfit}}
#' @keywords optimise ts
#' @export computefpsOS
computefpsOS <- function(mfit, Kappa=sqrt(0.5), normalisation="cmi", boot.fit=TRUE, ZA=1, ZAboot, dZA) {
  if(any(class(mfit) == "matrixfit")) {
    k <- Kappa
    if(normalisation != "cmi") k <- sqrt(0.5)
    else  mfit$kappa <- Kappa

    mfit$ZA <- ZA
    zab <- rep(1, times=mfit$boot.R)
    if(!missing(ZAboot) || !missing(dZA)) {
      if(missing(ZAboot)){
        zab <- rnorm(n=mfit$boot.R, mean=ZA, sd=dZA)
      }
      else {
        zab <- ZAboot
        if(length(ZAboot) != mfit$boot.R) {
          dZA <- sd(ZAboot)
          zab <- rnorm(n=mfit$boot.R, mean=ZA, sd=dZA)
        }
      }
    }
    mfit$ZAboot <- zab
    
    mfit$fpsOS <- ZA*sqrt(2)*k*mfit$opt.res$par[3]*sqrt(mfit$opt.res$par[1])/mfit$opt.res$par[1]    
    if(boot.fit) {
      mfit$fpsOS.tsboot <- zab*sqrt(2)*k*mfit$opt.tsboot[3,]*sqrt(mfit$opt.tsboot[1,])/mfit$opt.tsboot[1,]
    }
    else {
      mfit$fps.tsboot <- NA
    }
    mfit$normalisationOS <- normalisation
    return(invisible(mfit))
  }
  else {
    stop("computefps: expecting object of class matrixfit\nAborting...!\n")
  }
}
