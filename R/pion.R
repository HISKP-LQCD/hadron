#' Fit Pseudo Scalar Sector
#' 
#' fits one or several cosh and sinh to data for extracting the pseudo-scalar
#' meson mass, excited states in the pseudo-scalar channel and determines
#' decays constant and \eqn{m_\mathrm{PCAC}}{m_paca}
#' 
#' This is so far only for twisted mass lattice QCD in the twisted basis.
#' 
#' A call to \code{pion} returns an object, for which several generic functions
#' are available, such as \code{summary} and \code{plot}. The function should
#' always be used by also considering the output of \code{plot}. In particular,
#' in case of the gamma matrix method it should always be checked whether the
#' autocorrelation function could be integrated in a sensible way. For this
#' reason the autocorrelation function is plotted as well as the integrated
#' autocorrelation time as a function of the cut-off value for the integral.
#' Obviously, the autocorrelation function should be zero at the cut-off and
#' for larger values only fluctuate around zero. In other words, the
#' autocorrelation time should reach a plateau. Changing the parameter \code{S}
#' (see above) allows to influence the automatic determination of the cut-off
#' value (see \code{\link{uwerr}} for details).
#' 
#' Another useful check is also to devide the full ensemble into two (or more)
#' sub-ensembles (using the parameter \code{nrep}). The summary function will
#' then also compute a Q-value, which should be larger than \eqn{0.1} at least
#' for the error to be trustworthy. Otherwise the autocorrlation times are
#' presumably too large to determine the error reliably.
#' 
#' Finally it is also usefull to check the value of \eqn{\chi}{chi} per data
#' point. Such a plot is produced when \code{plot} is used.
#' 
#' The data must be ordered as in the output of \code{\link{readcmicor}}, see
#' \code{help(readcmicor)} for details.
#' 
#' The expected order of gamma matrices and operators (local-local,
#' local-fuzzed, fuzzed-local and fuzzed-fuzzed) (fuzzed = non-local) is as
#' follows for all charged mesons:
#' 
#' 1) the 4 operators for each type must be sorted like local-local,
#' local-fuzzed, fuzzed-local and fuzzed-fuzzed. (fuzzed=non-local)
#' 
#' 2) The 20 available types must be in the following order:\cr order PP PA AP
#' AA 44 P4 4P A4 4A for pion like \eqn{P=\gamma_5}{P=g5}
#' \eqn{A=\gamma_4\gamma_5}{A=g4g5} \eqn{4=\gamma_4}{4=g4}\cr order 44 VV AA 4V
#' V4 4A A4 VA AV for rho-a1 like \eqn{4=\gamma_i\gamma_4}{4=gig4}
#' \eqn{V=\gamma_i}{V=gi} \eqn{A=\gamma_i\gamma_5}{A=gig5}\cr order BB SS -
#' total 20 \eqn{\gamma_i\gamma_4\gamma_5}{B=gig4g5}\cr \eqn{S=I}
#' 
#' In this routine only PP PA AP AA 44 P4 4P A4 4A are used. See also
#' \code{cmicor}!
#' 
#' (cases with space index "i" are summed over i=1,2,3) (best choice is weaker
#' coupling at sink - ie PA rather than AP\cr order of magnitudes \eqn{P > 4 >
#' A} (4 mixes A)\cr order of magnitudes \eqn{4\sim A > V}{4 ~ A > V} (A mixes
#' V))\cr
#' 
#' itype=21 is conserved vector current at sink, \eqn{\gamma_5}{g5} at source (
#' iobs is LV 1, FV 5 )
#' 
#' pion will perform a fit of of the following matrix \tabular{lcccccc}{ \tab
#' PL \tab PF \tab AL \tab AF \tab 4L \tab 4F \cr PL \tab p1 p1 cosh \tab p1 p2
#' cosh \tab p1 p3 sinh \tab p1 p4 sinh \tab p1 p5 sinh \tab p1 p6 sinh \cr PF
#' \tab p2 p1 cosh \tab p2 p2 cosh \tab p2 p3 sinh \tab p2 p4 sinh \tab p2 p5
#' sinh \tab p2 p6 sinh \cr AL \tab p3 p1 sinh \tab p3 p2 sinh \tab p3 p3 cosh
#' \tab p3 p4 cosh \tab p3 p5 cosh \tab p3 p6 cosh \cr AF \tab p4 p1 sinh \tab
#' p4 p2 sinh \tab p4 p3 cosh \tab p4 p4 cosh \tab p4 p5 cosh \tab p4 p6 cosh
#' \cr 4L \tab p5 p1 sinh \tab p5 p2 sinh \tab p5 p3 cosh \tab p5 p4 cosh \tab
#' p5 p5 cosh \tab p5 p6 cosh \cr 4F \tab p6 p1 sinh \tab p6 p2 sinh \tab p6 p3
#' cosh \tab p6 p4 cosh \tab p6 p5 cosh \tab p6 p6 cosh \cr } for coupling
#' parameter \eqn{p_1} to \eqn{p_6} and a mass \eqn{m_\mathrm{PS}}{mps}
#' entering \eqn{\cosh}{cosh} (and \eqn{\sinh}{sinh} in the same way) as
#' \eqn{\cosh(-m_\mathrm{PS}(T/2-t))}{cosh(-mps(T/2-t))}. The values of \eqn{t}
#' are running from \code{t1} to \code{t2} as specified by the user.
#' 
#' Corresponding to the value of \code{matrix.size} only a submatrix of the
#' matrix given above is fitted. Moreover, if \code{no.masses} larger than one
#' additional masses and coupling parameters are introduced.
#' 
#' The pion decay constant is computed from the relation \deqn{f_\mathrm{PS} =
#' 4\kappa\mu\frac{p_1}{\sqrt{m_\mathrm{PS}}^3}}{% fps = 4 kappa mu
#' p_1/sqrt(mps)^3} and the PCAC mass (for \code{matrix.size} > 2) from
#' \deqn{m_\mathrm{PCAC} = m_\mathrm{PS} \frac{p_3}{2p_1}}{% mpcac = mps p_3/(2
#' p_1)} when the \code{\link{summary.pionfit}} function is called.  For
#' \code{matrix.size} > 4 also \eqn{Z_V} is computed as
#' \deqn{Z_\mathrm{V}=\frac{2\mu}{m_\mathrm{PS}}\frac{p_1}{p_5}}{% Z_V = 2 mu
#' p_1/(mps p_5)}
#' 
#' @param cmicor data to be fitted to as e.g. the output of
#' \code{\link{readcmicor}}
#' @param mu the value of the bare quark twisted mass
#' @param kappa the value of the hopping parameter
#' @param t1 lower bound for the fitrange in time (t1,t2). Counting starts with
#' 0.
#' @param t2 upper bound for the fitrange in time (t1,t2). Counting starts with
#' 0.
#' @param S passed to \code{uwerr}, see documentation of \code{\link{uwerr}}.
#' @param pl logical: if set to TRUE the function produces plots
#' @param skip number of measurements to be discarded at the beginning of the
#' time series. \code{skip} has no effect if two or more replica are used, see
#' argument \code{nrep}.
#' @param variational list of parameters used for the variational analysis
#' @param ind.vec index vector indexing the column numbers in cmicor to be used
#' @param no.masses number of masses to be extracted. This argument will set
#' the number of exponentials to be used in the fit.
#' 
#' Note that this is not yet stable for \code{no.masses > 1}.
#' @param matrix.size matrix size to be used in the fit. Can be currently set
#' to 2,4 and 6.
#' @param boot.R number of bootstrap samples for bootstrap analysis
#' @param boot.l average block size for blocking analysis with tsboot
#' @param tsboot.sim The type of simulation required to generate the replicate
#' time series. See \code{\link{tsboot}} for details.
#' @param method the type of error analysis to be used. Can be either
#' \dQuote{uwerr}, \dQuote{boot}, \dQuote{all} or \dQuote{no}. For \dQuote{no}
#' (or any other string) no error analysis is performed. This might be helpful
#' for a first impression and also to test different initial values for the
#' fitting parameters. The latter is in particular needed for more than one
#' state in the fit.
#' @param mass.guess numerical vector of mass-values to be used as initial
#' values in the fit. If given, it must have at least \code{no.masses} entries.
#' @param par.guess numerical vector of couling parameter values to be used as
#' initial values in the fit. The order is \eqn{P_L}, \eqn{P_F}, \eqn{A_L},
#' \eqn{A_F}, \eqn{4_L}, \eqn{4_F}, for notation see below. If given, it must
#' have at least \code{no.masses} times \code{matrix.size} entries.
#' @param nrep vector (N1, N2, ...) of replica length N1, N2. If missing it is
#' assumed that there is only one ensemble. If there are two or more replica
#' the parameter \code{skip} has no effect.
#' @param fit.routine The fit routine to be used. Default is \dQuote{gsl},
#' which uses the gnu scientific library \dQuote{gsl_multifit_fdfsolver} solver
#' to minimise the chisquare. All other values lead to the usage of R's
#' \link{optim} function. The latter choice might be significantly slower.
#' @return returns an object of \code{class} \code{pionfit} with the following
#' items
#' 
#' \item{fitresult}{ result from the fit as returned by \code{\link{optim}} }
#' \item{t1}{ lower bound for the fitrange in time (t1,t2). Counting starts
#' with 0.  } \item{t2}{ upper bound for the fitrange in time (t1,t2). Counting
#' starts with 0.  } \item{N}{ number of measurements found in the data }
#' \item{Time}{ Time extent found in the data } \item{fitdata}{
#' \code{\link{data.frame}} containing the time values used in the fit, the
#' averaged correlator and its error and the value of Chi for each time value }
#' \item{uwerrresultmps}{ the result of the time series analysis for the lowest
#' mass as carried out by \code{\link{uwerr}} } \item{uwerrresultmps2}{ the
#' result of the time series analysis for the second lowest mass as carried out
#' by \code{\link{uwerr}} if no.masses larger than 1.  }
#' \item{uwerrresultmps3}{ the result of the time series analysis for the
#' second lowest mass as carried out by \code{\link{uwerr}}, if no.masses
#' larger than 2.  } \item{uwerrresultfps}{ the result of the time series
#' analysis for the pseudo-scalar decay constant as carried out by
#' \code{\link{uwerr}}.  } \item{uwerrresultmpcac}{ the result of the time
#' series analysis for the PCAC mass as carried out by \code{\link{uwerr}} if
#' no.masses larger than 1.  } \item{uwerrresultzv}{ the gamma method analysis
#' for \eqn{Z_V}{Z_V}, if matrix size equals 6 } \item{effmass}{ effective
#' masses } \item{mu}{ the value of the bare quark twisted mass } \item{kappa}{
#' the value of the hopping parameter } \item{variational.masses}{ mass values
#' as determined by the variational analysis } \item{no.masses}{ no.masses
#' determined, copied from input } \item{matrix.size}{ size of the data matrix,
#' copied from input } \item{boot}{ object returned by the call to
#' \code{\link{boot}} if \code{method} was set correspodingly. Otherwise
#' \code{NULL}.  } \item{tsboot}{ object returned by the call to
#' \code{\link{tsboot}} if \code{method} was set correspodingly. Otherwise
#' \code{NULL}.  } \item{var.res}{ the full result of the variational analysis,
#' as returned by a call to \code{\link{variational}} } \item{method}{ error
#' analysis method as copied from input } \item{fit.routine}{
#' \code{fit.routine} as copied from input } \item{nrep}{ \code{nrep} as copied
#' from input }
#' @author Carsten Urbach, \email{curbach@@gmx.de}
#' @seealso \code{\link{readcmicor}}, \code{\link{uwerr}},
#' \code{\link{variational}}
#' @keywords optimize ts
#' @examples
#' 
#' library(hadron)
#' \dontrun{cmicor <- readcmicor("pion.dat")}
#' \dontrun{pionfit <- pion(cmicor, kappa=0.160856, mu=0.004, t1=10, t2=23,}
#' \dontrun{         no.masses=1, matrix.size=4)}
#' \dontrun{summary(pionfit)}
#' 
pion <- function(cmicor, mu=0.1, kappa=0.156, t1, t2, S=1.5, pl=FALSE, skip=0,
                variational=list(ta=3, tb=4, N=6), ind.vec=c(1,3,4,5),
                no.masses=1, matrix.size=2, boot.R=99, boot.l=10, tsboot.sim="geom",
                method="uwerr", fit.routine="optim", mass.guess, par.guess, nrep) {
  
  if(missing(cmicor)) {
    stop("Error! Data is missing!")
  }
  if(missing(t1) || missing(t2)) {
    stop("Error! t1 and t2 must be specified!")
  }
  if(missing(mass.guess)) {
    mass.guess <- c(0.2, 1., 3.)
  }
  else {
    if(length(mass.guess) < no.masses) {
      stop("mass.guess has not the correct length!")
    }
  }
  if(missing(par.guess)) {
    par.guess <- c(1.,0.8,0.1,0.1,0.1,0.1, 0.1,0.1,0.1,0.1,0.1,0.1, 0.1,0.1,0.1,0.1,0.1,0.1) 
  }
  else{
    if(length(par.guess) < no.masses*matrix.size) {
      stop("par.guess has not the correct length!")
    }
  }
  par <- numeric()
  length(par) <- no.masses*(matrix.size+1)
  for(i in 1:no.masses) {
    par[i*(matrix.size+1)] <- mass.guess[i]
    par[(c(1:matrix.size)+(i-1)*(matrix.size+1))] = par.guess[(c(1:matrix.size)+(i-1)*(matrix.size))]
  }
  
  Time <-  2*max(cmicor[,ind.vec[2]])
  Thalf <- max(cmicor[,ind.vec[2]])
  T1 <- Thalf+1
  t1p1 <- (t1+1)
  t2p1 <- (t2+1)
  nrObs <- max(cmicor[,ind.vec[1]])
  Skip <- (skip*(T1)*nrObs*4+1)
  Length <- length(cmicor[,ind.vec[3]])
  if(missing(nrep)) {
    nrep <- c(length(cmicor[((Skip):Length),ind.vec[3]])/(nrObs*(T1)*4))
  }
  else {
    skip <- 0
    if(sum(nrep) != length(cmicor[((Skip):Length),ind.vec[3]])/(nrObs*(T1)*4)) {
      stop("sum of replica differs from total no of measurements!")
    }
  }

  Z <- array(cmicor[((Skip):Length),ind.vec[3]], 
             dim=c(nrObs*(T1)*4,(length(cmicor[((Skip):Length),ind.vec[3]])/(nrObs*(T1)*4))))
  # negative times
  W <- array(cmicor[((Skip):Length),ind.vec[4]], 
             dim=c(nrObs*(T1)*4,(length(cmicor[((Skip):Length),ind.vec[4]])/(nrObs*(T1)*4))))

  rm(cmicor)
  W <- arrangeCor.pion(T1=T1, W=W, Z=Z)
  rm(Z)
  
#  options(show.error.messages = FALSE)
  pion.eff.ll <- effectivemass(from=(t1+1), to=(t2+1), Time, W[1:T1,] , pl=FALSE, S=1.5, nrep=nrep)
  pion.eff.lf <- effectivemass(from=(t1+1), to=(t2+1), Time, W[(T1+1):(2*T1),] , pl=FALSE, S=1.5, nrep=nrep)
  pion.eff.ff <- effectivemass(from=(t1+1), to=(t2+1), Time, W[(3*T1+1):(4*T1),] , pl=FALSE, S=1.5, nrep=nrep)
  options(show.error.messages = TRUE)
  

  pion.eff <- data.frame(t=pion.eff.ll$t, mll=pion.eff.ll$mass, dmll=pion.eff.ll$dmass,
                         mlf=pion.eff.lf$mass, dmlf=pion.eff.lf$dmass,
                         mff=pion.eff.ff$mass, dmff=pion.eff.ff$dmass)

#  pion.eff <- NULL
  
  Cor <- rep(0., times=9*4*T1)
  E <- rep(0., times=9*4*T1)
  
  for(i in 1:(9*4*T1)) {
    Cor[i] <- mean(W[(i),])
    tmpe <- try(uwerrprimary(W[(i),], pl=F, nrep=nrep)$dvalue, TRUE)
    if(!inherits(tmpe, "try-error")) E[i] = tmpe
    else {
      warning("error of correlator replaced by naive estimate!\n", call.=F)
      E[i] = sd(W[(i),])/sqrt(length(W[(i),]))
    }
  }

  res.var <- variational(Cor=Cor, N=variational$N, ta=variational$ta, tb=variational$tb, tmax = T1,
                         T1=T1, matrix.size=matrix.size, no.masses = no.masses)
  N <- max(matrix.size,variational$N)
#  par <- res.var$par
#  cat(res.var$par, "\n")
  variational.masses <-  res.var$variational.masses

  # Index vector of data to be used in the analysis
  ii <- c((t1p1):(t2p1), (t1p1+T1):(t2p1+T1), (t1p1+3*T1):(t2p1+3*T1))
  if(matrix.size > 2) {
    ii <- c(ii, (t1p1+4*T1):(t2p1+4*T1), (t1p1+5*T1):(t2p1+5*T1),
            (t1p1+6*T1):(t2p1+6*T1), (t1p1+7*T1):(t2p1+7*T1),
            (t1p1+12*T1):(t2p1+12*T1), (t1p1+13*T1):(t2p1+13*T1),
            (t1p1+15*T1):(t2p1+15*T1))
  }
  if(matrix.size > 4) {
    ii <- c(ii, (t1p1+16*T1):(t2p1+16*T1), (t1p1+17*T1):(t2p1+17*T1),
            (t1p1+19*T1):(t2p1+19*T1),
            (t1p1+20*T1):(t2p1+20*T1), (t1p1+21*T1):(t2p1+21*T1),
            (t1p1+22*T1):(t2p1+22*T1), (t1p1+23*T1):(t2p1+23*T1),
            (t1p1+28*T1):(t2p1+28*T1), (t1p1+29*T1):(t2p1+29*T1),
            (t1p1+30*T1):(t2p1+30*T1), (t1p1+31*T1):(t2p1+31*T1))
  }

  #BFGS
  if(no.masses == 1) {
    pionfit <- optim(par, ChiSqr.1mass, method="BFGS", control=list(trace=0),Thalf=Thalf,
                     x=c((t1):(t2)), y=Cor[ii], err=E[ii], tr = (t2-t1+1), N=matrix.size)
    fit.mass <- abs(pionfit$par[matrix.size+1])
    fit.fpi <- 2*kappa*2*mu/sqrt(2)*abs(pionfit$par[1])/sqrt(fit.mass^3)
    if(matrix.size > 2) {
      fit.pcac <- abs(pionfit$par[5])*pionfit$par[3]/pionfit$par[1]/2.
    }
    if(matrix.size > 4) {
      fit.zv <- 2.*mu/abs(pionfit$par[5])*pionfit$par[1]/pionfit$par[5]
    }
  }
  else if(no.masses == 2) {
    pionfit <- optim(par, ChiSqr.2mass, method="BFGS", control=list(trace=0),Thalf=Thalf,
                     x=c((t1):(t2)), y=Cor[ii], err=E[ii], tr = (t2-t1+1), N=matrix.size)
    fit.mass <- sort(abs(pionfit$par[c((matrix.size+1),(2*matrix.size+2))]))
  }
  else if(no.masses > 2) {
    pionfit <- optim(par, ChiSqr.3mass, method="BFGS", control=list(trace=0),Thalf=Thalf,
                     x=c((t1):(t2)), y=Cor[ii], err=E[ii], tr = (t2-t1+1), N=matrix.size)
    fit.mass <- sort(abs(pionfit$par[c((matrix.size+1),(2*matrix.size+2),(3*matrix.size+3))]))
  }
  if(pionfit$convergence!=0) {
    warning("optim did not converge for pionfit! ", pionfit$convergence)
  }

  fit.dof <- (t2-t1+1)*3-length(pionfit$par)
  fit.chisqr <- pionfit$value

  if(pl) {
    plot.effmass(m=fit.mass, ll=pion.eff.ll, lf=pion.eff.lf, ff=pion.eff.ff)
  }

  fit.uwerrm <- NULL
  fit.uwerrf <- NULL
  fit.uwerrpcac <- NULL
  fit.uwerrzv <- NULL
  fit.uwerrm2 <- NULL
  fit.uwerrm3 <- NULL
  fit.boot <- NULL
  fit.tsboot <- NULL
  if(method == "uwerr" || method == "all") {
    fit.uwerrm <- uwerr(f=fitmasses.pion, data=t(W[ii,]), S=S, pl=pl, nrep=nrep,
                        Time=Time, t1=t1, t2=t2, Err=E[ii], par=par, N=matrix.size,
                        no.masses=no.masses, fit.routine=fit.routine)

    fit.uwerrf <- uwerr(f=fitf.pion, data=t(W[ii,]), S=S, pl=pl, nrep=nrep,
                        Time=Time, t1=t1, t2=t2, Err=E[ii], par=par, N=matrix.size, no.masses=no.masses,
                        kappa=kappa, mu=mu, fit.routine=fit.routine)
    if(matrix.size > 2) {
      fit.uwerrpcac <- uwerr(f=fitmpcac.pion, data=t(W[ii,]), S=S, pl=pl, nrep=nrep,
                             Time=Time, t1=t1, t2=t2, Err=E[ii], par=par, N=matrix.size, no.masses=no.masses,
                             kappa=kappa, mu=mu, fit.routine=fit.routine)
    }
    if(matrix.size > 4) {
      fit.uwerrzv <- uwerr(f=fitzv.pion, data=t(W[ii,]), S=S, pl=pl, nrep=nrep, Time=Time, t1=t1, t2=t2,
                           Err=E[ii], par=par, N=matrix.size, no.masses=no.masses, kappa=kappa,
                           mu=mu, fit.routine=fit.routine)
    }
    
    if(no.masses == 2) {
      fit.uwerrm2 <- uwerr(f=fitmasses.pion, data=t(W[ii,]), S=S, pl=pl, nrep=nrep,
                           Time=Time, t1=t1, t2=t2, Err=E[ii], par=par, N=matrix.size,
                           no.masses=no.masses, no=2, fit.routine=fit.routine)
    }
    if(no.masses > 2) {
      fit.uwerrm3 <- uwerr(f=fitmasses.pion, data=t(W[ii,]), S=S, pl=pl, nrep=nrep,
                           Time=Time, t1=t1, t2=t2, Err=E[ii], par=par, N=matrix.size,
                           no.masses=no.masses, no=3, fit.routine=fit.routine)
    }
  }
  if(method == "boot" || method == "all") {
    fit.boot <- boot::boot(data=t(W[ii,]), statistic=fit.pion.boot, R=boot.R, stype="i",
                           Time=Time, t1=t1, t2=t2, Err=E[ii], par=par, N=matrix.size, no.masses=no.masses,
                           kappa=kappa, mu=mu, fit.routine=fit.routine)

    fit.tsboot <- boot::tsboot(tseries=t(W[ii,]), statistic=fit.pion.boot, R=boot.R, l=boot.l, sim=tsboot.sim,
                               Time=Time, t1=t1, t2=t2, Err=E[ii], par=par, N=matrix.size, no.masses=no.masses,
                               kappa=kappa, mu=mu, fit.routine=fit.routine)
  }

  
  Chi <- rep(0., times=9*4*T1)
  Fit <- rep(0., times=9*4*T1)

  jj <-  c(t1p1:t2p1)
  Fit[jj] <- pionfit$par[1]^2*CExp(m=fit.mass[1], Time=2*Thalf, x=jj-1)
  Fit[jj+T1] <- pionfit$par[1]*pionfit$par[2]*CExp(m=fit.mass[1], Time=2*Thalf, x=jj-1)
  Fit[jj+2*T1] <- pionfit$par[1]*pionfit$par[2]*CExp(m=fit.mass[1], Time=2*Thalf, x=jj-1)
  Fit[jj+3*T1] <- pionfit$par[2]*pionfit$par[2]*CExp(m=fit.mass[1], Time=2*Thalf, x=jj-1)
  if(matrix.size > 2) {
    Fit[jj+4*T1] <- pionfit$par[1]*pionfit$par[3]*CExp(m=fit.mass[1], Time=2*Thalf, x=jj-1, sign=-1.)
    Fit[jj+5*T1] <- pionfit$par[1]*pionfit$par[4]*CExp(m=fit.mass[1], Time=2*Thalf, x=jj-1, sign=-1.)
    Fit[jj+6*T1] <- pionfit$par[2]*pionfit$par[3]*CExp(m=fit.mass[1], Time=2*Thalf, x=jj-1, sign=-1.)
    Fit[jj+7*T1] <- pionfit$par[2]*pionfit$par[4]*CExp(m=fit.mass[1], Time=2*Thalf, x=jj-1, sign=-1.)

    Fit[jj+12*T1] <- pionfit$par[3]*pionfit$par[3]*CExp(m=fit.mass[1], Time=2*Thalf, x=jj-1)
    Fit[jj+13*T1] <- pionfit$par[3]*pionfit$par[4]*CExp(m=fit.mass[1], Time=2*Thalf, x=jj-1)
    Fit[jj+14*T1] <- pionfit$par[3]*pionfit$par[4]*CExp(m=fit.mass[1], Time=2*Thalf, x=jj-1)
    Fit[jj+15*T1] <- pionfit$par[4]*pionfit$par[4]*CExp(m=fit.mass[1], Time=2*Thalf, x=jj-1)
  }
  if(matrix.size > 4) {
    Fit[jj+16*T1] <- pionfit$par[5]*pionfit$par[5]*CExp(m=fit.mass[1], Time=2*Thalf, x=jj-1)
    Fit[jj+17*T1] <- pionfit$par[5]*pionfit$par[6]*CExp(m=fit.mass[1], Time=2*Thalf, x=jj-1)
    Fit[jj+18*T1] <- pionfit$par[5]*pionfit$par[6]*CExp(m=fit.mass[1], Time=2*Thalf, x=jj-1)
    Fit[jj+19*T1] <- pionfit$par[6]*pionfit$par[6]*CExp(m=fit.mass[1], Time=2*Thalf, x=jj-1)

    Fit[jj+20*T1] <- pionfit$par[1]*pionfit$par[5]*CExp(m=fit.mass[1], Time=2*Thalf, x=jj-1, sign=-1.)
    Fit[jj+21*T1] <- pionfit$par[1]*pionfit$par[6]*CExp(m=fit.mass[1], Time=2*Thalf, x=jj-1, sign=-1.)
    Fit[jj+22*T1] <- pionfit$par[2]*pionfit$par[5]*CExp(m=fit.mass[1], Time=2*Thalf, x=jj-1, sign=-1.)
    Fit[jj+23*T1] <- pionfit$par[2]*pionfit$par[6]*CExp(m=fit.mass[1], Time=2*Thalf, x=jj-1, sign=-1.)
    
    Fit[jj+28*T1] <- pionfit$par[3]*pionfit$par[5]*CExp(m=fit.mass[1], Time=2*Thalf, x=jj-1)
    Fit[jj+29*T1] <- pionfit$par[3]*pionfit$par[6]*CExp(m=fit.mass[1], Time=2*Thalf, x=jj-1)
    Fit[jj+30*T1] <- pionfit$par[4]*pionfit$par[5]*CExp(m=fit.mass[1], Time=2*Thalf, x=jj-1)
    Fit[jj+31*T1] <- pionfit$par[4]*pionfit$par[6]*CExp(m=fit.mass[1], Time=2*Thalf, x=jj-1)
  }
  
  Chi[ii] <- (Fit[ii]-Cor[ii])/E[ii]
  
  res <- list(fitresult=pionfit, t1=t1, t2=t2, N=length(W[1,]), Time=Time,
              fitdata=data.frame(t=(jj-1), Fit=Fit[ii], Cor=Cor[ii], Err=E[ii], Chi=Chi[ii]),
              uwerrresultmps=fit.uwerrm, uwerrresultmps2=fit.uwerrm2, uwerrresultmps3=fit.uwerrm3,
              uwerrresultfps=fit.uwerrf, uwerrresultmpcac=fit.uwerrpcac, uwerrresultzv=fit.uwerrzv,
              boot=fit.boot, tsboot=fit.tsboot, method=method,
              effmass=pion.eff, kappa=kappa, mu=mu, fit.routine=fit.routine,
              variational.masses=variational.masses, no.masses=no.masses,
              matrix.size = matrix.size, nrep=nrep, res.var=res.var)
  attr(res, "class") <- c("pionfit", "list")  
  return(invisible(res))
}

arrangeCor.pion <- function(T1, W, Z) {

  for(i in 1:(T1)) {
    two <- 2.
    if(i==1 || i==(T1)) {
      # Take care of zeros in the correlators when summing t and T-t+1
      two <- 1.
    }

    # PP for LL, (LF + FL)/2, FF -> symmetric cosh
    W[i,] <- (W[i,] + Z[i,])/two
    W[(i+T1),] <- (W[(i+T1),] + W[(i+2*T1),] + Z[(i+T1),] + Z[(i+2*T1),])/two/2.
    W[(i+2*T1),] <- W[(i+T1),]
    W[(i+3*T1),] <- (W[(i+3*T1),] + Z[(i+3*T1),])/two
    # PA for LL, LF, FL, FF -> antisymmetric -sinh
    W[(i+4*T1),] <- (W[(i+4*T1),] - Z[(i+4*T1),])/two
    W[(i+5*T1),] <- (W[(i+5*T1),] - Z[(i+5*T1),])/two
    W[(i+6*T1),] <- (W[(i+6*T1),] - Z[(i+6*T1),])/two
    W[(i+7*T1),] <- (W[(i+7*T1),] - Z[(i+7*T1),])/two
    # use again PA -> -sinh
    W[(i+8*T1),] <- W[(i+4*T1),]
    W[(i+9*T1),] <- W[(i+5*T1),]
    W[(i+10*T1),] <- W[(i+6*T1),]
    W[(i+11*T1),] <- W[(i+7*T1),]
    # AA for LL, (LF + FL)/2, FF -> symmetric cosh
    # if needed, bug fix minus here
    W[(i+12*T1),] <- (W[(i+12*T1),] + Z[(i+12*T1),])/two
    W[(i+13*T1),] <- (W[(i+13*T1),] + W[(i+14*T1),] + Z[(i+13*T1),] + Z[(i+14*T1),])/two/2.
    W[(i+14*T1),] <-  W[(i+13*T1),]
    W[(i+15*T1),] <- (W[(i+15*T1),] + Z[(i+15*T1),])/two
    # 44 for LL, LF + FL/2, FF -> -cosh
    W[(i+16*T1),] <- -(W[(i+16*T1),] + Z[(i+16*T1),])/two
    W[(i+17*T1),] <- -(W[(i+17*T1),] + W[(i+17*T1),] + Z[(i+18*T1),] + Z[(i+18*T1),])/two/2.
    W[(i+18*T1),] <-   W[(i+17*T1),]
    W[(i+19*T1),] <- -(W[(i+19*T1),] + Z[(i+19*T1),])/two
    # P4 for LL, LF, FL, FF -> antisymmetric -sinh
    W[(i+20*T1),] <- (W[(i+20*T1),] - Z[(i+20*T1),])/two
    W[(i+21*T1),] <- (W[(i+21*T1),] - Z[(i+21*T1),])/two
    W[(i+22*T1),] <- (W[(i+22*T1),] - Z[(i+22*T1),])/two
    W[(i+23*T1),] <- (W[(i+23*T1),] - Z[(i+23*T1),])/two
    # 4P use again P4 -> -sinh
    W[(i+24*T1),] <- W[(i+20*T1),]
    W[(i+25*T1),] <- W[(i+21*T1),]
    W[(i+26*T1),] <- W[(i+22*T1),]
    W[(i+27*T1),] <- W[(i+23*T1),]
    # A4 use 4A -> cosh
    # if needed, bug fix minus here
    W[(i+28*T1),] <- (W[(i+32*T1),] + Z[(i+32*T1),])/two
    W[(i+29*T1),] <- (W[(i+34*T1),] + Z[(i+34*T1),])/two
    W[(i+30*T1),] <- (W[(i+33*T1),] + Z[(i+33*T1),])/two
    W[(i+31*T1),] <- (W[(i+35*T1),] + Z[(i+35*T1),])/two
    # 4A for LL, LF, FL, FF -> cosh
    W[(i+32*T1),] <- W[(i+28*T1),]
    W[(i+33*T1),] <- W[(i+29*T1),]
    W[(i+34*T1),] <- W[(i+30*T1),]
    W[(i+35*T1),] <- W[(i+31*T1),]
  }
  return(invisible(W))
}

fitmasses.pion <- function(Cor, Err, t1, t2, Time, par=c(1.,0.1,0.12),
                         N=2, no.masses=1, no=1, kludge=FALSE, fit.routine="optim") {
  Thalf <- Time/2
  T1 <- Thalf+1
  t1p1 <- (t1+1)
  t2p1 <- (t2+1)
  tr <- (t2-t1+1)

  if(no.masses == 1) {
    fit <- optim(par, ChiSqr.1mass, method="BFGS", Thalf=Thalf,
                 x=c((t1):(t2)), y=Cor, err=Err, tr=tr, N=N)

    return(abs(fit$par[N+1]))
  }
  else if (no.masses == 2) {
    fit <- optim(par, ChiSqr.2mass, method="BFGS", Thalf=Thalf,
                 x=c((t1):(t2)), y=Cor, err=Err, tr=tr, N=N, kludge=kludge)

    return(sort(abs(fit$par[c((N+1),(2*N+2))]))[no])
  }
  else if (no.masses == 3) {
    fit <- optim(par, ChiSqr.3mass, method="BFGS", Thalf=Thalf,
                 x=c((t1):(t2)), y=Cor, err=Err, tr=tr, N=N, kludge=kludge)
    return(sort(abs(fit$par[c((N+1),(2*N+2),(3*N+3))]))[no])
  }
}

fitf.pion <- function(Cor, Err, t1, t2, Time, par=c(1.,0.1,0.12),
                         N=2, no.masses=1, no=1, kappa, mu, kludge=FALSE, fit.routine="optim") {
  Thalf <- Time/2
  T1 <- Thalf+1
  t1p1 <- (t1+1)
  t2p1 <- (t2+1)
  tr <- (t2-t1+1)

  if(no.masses == 1) {
    fit <- optim(par, ChiSqr.1mass, method="BFGS", Thalf=Thalf,
                 x=c((t1):(t2)), y=Cor, err=Err, tr=tr, N=N)
    sort.ind <- c(1)
  }
  else if (no.masses == 2) {
    fit <- optim(par, ChiSqr.2mass, method="BFGS", Thalf=Thalf,
                 x=c((t1):(t2)), y=Cor, err=Err, tr=tr, N=N, kludge=kludge)
    sort.ind <- order(fit$par[c((N+1),(2*N+2))])
  }
  else if (no.masses == 3) {
    fit <- optim(par, ChiSqr.3mass, method="BFGS", Thalf=Thalf,
                 x=c((t1):(t2)), y=Cor, err=Err, tr=tr, N=N, kludge=kludge)
    sort.ind <- order(fit$par[c((N+1),(2*N+2),(3*N+3))])
  }
#  fit.fpi <- 2*kappa*2*mu/sqrt(2)*abs(fit$par[(sort.ind[1]-1)*(N+1)+1])/sqrt(abs(fit$par[sort.ind[1]*(N+1)])^3)
  fit.fpi <- abs(fit$par[(sort.ind[1]-1)*(N+1)+1])/sqrt(abs(fit$par[sort.ind[1]*(N+1)])^3)
  return(fit.fpi)
}

fitmpcac.pion <- function(Cor, Err, t1, t2, Time, par=c(1.,0.1,0.12),
                          N=2, no.masses=1, no=1, kappa, mu, kludge=FALSE,
                          fit.routine="optim") {
  Thalf <- Time/2
  T1 <- Thalf+1
  t1p1 <- (t1+1)
  t2p1 <- (t2+1)
  tr <- (t2-t1+1)

  if(no.masses == 1) {
    fit <- optim(par, ChiSqr.1mass, method="BFGS", Thalf=Thalf,
                 x=c((t1):(t2)), y=Cor, err=Err, tr=tr, N=N)

    sort.ind <- c(1)
  }
  else if (no.masses == 2) {
    fit <- optim(par, ChiSqr.2mass, method="BFGS", Thalf=Thalf,
                 x=c((t1):(t2)), y=Cor, err=Err, tr=tr, N=N, kludge=kludge)
    sort.ind <- order(fit$par[c((N+1),(2*N+2))])
  }
  else if (no.masses == 3) {
    fit <- optim(par, ChiSqr.3mass, method="BFGS", Thalf=Thalf,
                 x=c((t1):(t2)), y=Cor, err=Err, tr=tr, N=N, kludge=kludge)
    sort.ind <- order(fit$par[c((N+1),(2*N+2),(3*N+3))])
  }
  fit.pcac <- abs(fit$par[sort.ind[1]*(N+1)])*fit$par[(sort.ind[1]-1)*(N+1)+3]/fit$par[(sort.ind[1]-1)*(N+1)+1]/2.
  return(fit.pcac)
}

fitzv.pion <- function(Cor, Err, t1, t2, Time, par=c(1.,0.1,0.12),
                       N=2, no.masses=1, no=1, kappa, mu, kludge=FALSE,
                       fit.routine="optim") {
  Thalf <- Time/2
  T1 <- Thalf+1
  t1p1 <- (t1+1)
  t2p1 <- (t2+1)
  tr <- (t2-t1+1)

  if(no.masses == 1) {
    fit <- optim(par, ChiSqr.1mass, method="BFGS", Thalf=Thalf,
                 x=c((t1):(t2)), y=Cor, err=Err, tr=tr, N=N)
    sort.ind <- c(1)
  }
  else if (no.masses == 2) {
    fit <- optim(par, ChiSqr.2mass, method="BFGS", Thalf=Thalf,
                 x=c((t1):(t2)), y=Cor, err=Err, tr=tr, N=N, kludge=kludge)
    sort.ind <- order(fit$par[c((N+1),(2*N+2))])
  }
  else if (no.masses == 3) {
    fit <- optim(par, ChiSqr.3mass, method="BFGS", Thalf=Thalf,
                 x=c((t1):(t2)), y=Cor, err=Err, tr=tr, N=N, kludge=kludge)
    sort.ind <- order(fit$par[c((N+1),(2*N+2),(3*N+3))])
  }
  fit.zv <- 2*mu/abs(fit$par[sort.ind[1]*(N+1)])*fit$par[(sort.ind[1]-1)*(N+1)+1]/fit$par[(sort.ind[1]-1)*(N+1)+5]
  return(fit.zv)
}


fit.pion.boot <- function(Z, d, Err, t1, t2, Time, par=c(1.,0.1,0.12),
                          N=2, no.masses=1, kludge=FALSE, kappa, mu,
                          fit.routine="optim") {
  Thalf <- Time/2
  T1 <- Thalf+1
  t1p1 <- (t1+1)
  t2p1 <- (t2+1)
  tr <- (t2-t1+1)
  Cor <- rep(0., times=length(Z[1,]))
  if(!missing(d)) {
    for(i in 1:length(Z[1,])) {
      Cor[i] = mean(Z[d,(i)])
    }
  }
  else {
    for(i in 1:length(Z[1,])) {
      Cor[i] = mean(Z[,(i)])
    }
  }

  if(no.masses == 1) {
    fit <- optim(par, ChiSqr.1mass, method="BFGS", Thalf=Thalf,
                 x=c((t1):(t2)), y=Cor, err=Err, tr=tr, N=N)
    sort.ind <- c(1)
    fit.fpi <- 2*kappa*2*mu/sqrt(2)*abs(fit$par[(sort.ind[1]-1)*(N+1)+1])/sqrt(abs(fit$par[sort.ind[1]*(N+1)])^3)
    fit.pcac <- abs(fit$par[sort.ind[1]*(N+1)])*fit$par[(sort.ind[1]-1)*(N+1)+3]/fit$par[(sort.ind[1]-1)*(N+1)+1]/2.
    fit.zv <- 2.*mu/abs(fit$par[sort.ind[1]*(N+1)])*fit$par[(sort.ind[1]-1)*(N+1)+1]/fit$par[(sort.ind[1]-1)*(N+1)+5]
    if(N > 4) {
      return(c(abs(fit$par[N+1]), fit.fpi, fit$par[c(1:N)],
               fit.pcac, fit.zv,
               fit$value))
    }
    else if(N > 2) {
      return(c(abs(fit$par[N+1]), fit.fpi, fit$par[c(1:N)],
               fit.pcac,
               fit$value))
    }
    else{
      return(c(abs(fit$par[N+1]), fit.fpi, fit$par[c(1:N)],
               fit$value))
    }
  }
  else if (no.masses == 2) {
    fit <- optim(par, ChiSqr.2mass, method="BFGS", Thalf=Thalf,
                 x=c((t1):(t2)), y=Cor, err=Err, tr=tr, N=N, kludge=kludge)
    sort.ind <- order(fit$par[c((N+1),(2*N+2))])
    fit.fpi <- 2*kappa*2*mu/sqrt(2)*abs(fit$par[(sort.ind[1]-1)*(N+1)+1])/sqrt(abs(fit$par[sort.ind[1]*(N+1)])^3)
    fit.pcac <- abs(fit$par[sort.ind[1]*(N+1)])*fit$par[(sort.ind[1]-1)*(N+1)+3]/fit$par[(sort.ind[1]-1)*(N+1)+1]/2.
    fit.zv <- 2.*mu/abs(fit$par[sort.ind[1]*(N+1)])*fit$par[(sort.ind[1]-1)*(N+1)+1]/fit$par[(sort.ind[1]-1)*(N+1)+5]
    if(N > 4) {
      return(c(abs(fit$par[sort.ind[1]*(N+1)]), fit.fpi,
               fit$par[c(((sort.ind[1]-1)*(N+1)+1):((sort.ind[1])*(N+1)-1))],
               fit.pcac, fit.zv,
               abs(fit$par[sort.ind[2]*(N+1)]),
               fit$par[c(((sort.ind[2]-1)*(N+1)+1):((sort.ind[2])*(N+1)-1))],
               fit$value))
    }
    if(N > 2) {
      return(c(abs(fit$par[sort.ind[1]*(N+1)]), fit.fpi,
               fit$par[c(((sort.ind[1]-1)*(N+1)+1):((sort.ind[1])*(N+1)-1))],
               fit.pcac,
               abs(fit$par[sort.ind[2]*(N+1)]),
               fit$par[c(((sort.ind[2]-1)*(N+1)+1):((sort.ind[2])*(N+1)-1))],
               fit$value))
    }
    else {
      return(c(abs(fit$par[sort.ind[1]*(N+1)]), fit.fpi,
               fit$par[c(((sort.ind[1]-1)*(N+1)+1):((sort.ind[1])*(N+1)-1))],
               abs(fit$par[sort.ind[2]*(N+1)]),
               fit$par[c(((sort.ind[2]-1)*(N+1)+1):((sort.ind[2])*(N+1)-1))],
               fit$value))
    }
  }
  else if (no.masses == 3) {
    fit <- optim(par, ChiSqr.3mass, method="BFGS", Thalf=Thalf,
                 x=c((t1):(t2)), y=Cor, err=Err, tr=tr, N=N, kludge=kludge)
    sort.ind <- order(fit$par[c((N+1),(2*N+2),(3*N+3))])
    fit.fpi <- 2*kappa*2*mu/sqrt(2)*abs(fit$par[(sort.ind[1]-1)*(N+1)+1])/sqrt(abs(fit$par[sort.ind[1]*(N+1)])^3)
    fit.pcac <- abs(fit$par[sort.ind[1]*(N+1)])*fit$par[(sort.ind[1]-1)*(N+1)+3]/fit$par[(sort.ind[1]-1)*(N+1)+1]/2.
    fit.zv <- 2.*mu/abs(fit$par[sort.ind[1]*(N+1)])*fit$par[(sort.ind[1]-1)*(N+1)+1]/fit$par[(sort.ind[1]-1)*(N+1)+5]
    if(N > 4) {
      return(c(abs(fit$par[sort.ind[1]*(N+1)]), fit.fpi,
               fit$par[c(((sort.ind[1]-1)*(N+1)+1):((sort.ind[1])*(N+1)-1))],
               fit.pcac, fit.zv,
               abs(fit$par[sort.ind[2]*(N+1)]),
               fit$par[c(((sort.ind[2]-1)*(N+1)+1):((sort.ind[2])*(N+1)-1))],
               abs(fit$par[sort.ind[3]*(N+1)]),
               fit$par[c(((sort.ind[3]-1)*(N+1)+1):((sort.ind[3])*(N+1)-1))],
               fit$value))
    }
    if(N > 2) {
      return(c(abs(fit$par[sort.ind[1]*(N+1)]), fit.fpi,
               fit$par[c(((sort.ind[1]-1)*(N+1)+1):((sort.ind[1])*(N+1)-1))],
               fit.pcac,
               abs(fit$par[sort.ind[2]*(N+1)]),
               fit$par[c(((sort.ind[2]-1)*(N+1)+1):((sort.ind[2])*(N+1)-1))],
               abs(fit$par[sort.ind[3]*(N+1)]),
               fit$par[c(((sort.ind[3]-1)*(N+1)+1):((sort.ind[3])*(N+1)-1))],
               fit$value))
    }
    else {
      return(c(abs(fit$par[sort.ind[1]*(N+1)]), fit.fpi,
               fit$par[c(((sort.ind[1]-1)*(N+1)+1):((sort.ind[1])*(N+1)-1))],
               abs(fit$par[sort.ind[2]*(N+1)]),
               fit$par[c(((sort.ind[2]-1)*(N+1)+1):((sort.ind[2])*(N+1)-1))],
               abs(fit$par[sort.ind[3]*(N+1)]),
               fit$par[c(((sort.ind[3]-1)*(N+1)+1):((sort.ind[3])*(N+1)-1))],
               fit$value))
    }
  }
}


