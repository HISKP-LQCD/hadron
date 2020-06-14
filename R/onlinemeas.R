#' determines pion mass and pcac mass from online measured correlator of the
#' HMC code
#' 
#' determines pion mass and pcac mass from online measured correlator of the
#' HMC code
#' 
#' The online measurements in the HMC code compute the PP and PA correlation
#' functions summed over spatial x for all t. We analyse these correlators in
#' different ways:
#' 
#' First, only the PP correlator is analysed and fitted by
#' \eqn{p_1^2\cosh(-m(t-T/2))}{p1*p1*cosh(-m(t-T/2))} for \eqn{m} and
#' \eqn{p_1}{p1}.
#' 
#' Second, PP and PA correlators are fitted together with three parameters as
#' \eqn{C_\mathrm{PP} = p_1^2\cosh(-m(t-T/2))}{C_PP = p1*p1*cosh(-m(t-T/2))}
#' and \eqn{C_\mathrm{PA} = }{C_PA = p1*p2*cosh(-m(t-T/2))}\eqn{
#' p_1p_2\cosh(-m(t-T/2))}{C_PA = p1*p2*cosh(-m(t-T/2))} in a simultaneous fit.
#' \eqn{m} is then the pseudo scalar mass and the pcac mass is determined from
#' \deqn{m_\mathrm{PCAC} = m_\mathrm{PS} \frac{p_2}{2p_1}}{% mpcac = mps p_2/(2
#' p_1)}
#' 
#' Finally, the PCAC mass can also be determined computing
#' \deqn{m_\mathrm{PCAC}(t) =% }{% C_PA(t) = (C_PA(t+1)-C_PA(t-1))/(4
#' C_PP(t))}\deqn{
#' \frac{C_\mathrm{PA}(t+1)-C_\mathrm{PA}(t-1)}{4C_\mathrm{PP}(t)}}{% C_PA(t) =
#' (C_PA(t+1)-C_PA(t-1))/(4 C_PP(t))} using the symmetric finite difference
#' operator.
#' 
#' @param data data to be fitted to as e.g. the output of
#' \code{\link{readcmicor}}. Currently only \code{cmicor} format is supported.
#' @param t1 lower bound for the fitrange in time (t1,t2). Counting starts with
#' 0.
#' @param t2 upper bound for the fitrange in time (t1,t2). Counting starts with
#' 0.
#' @param mu twisted mass parameter.
#' @param kappa hopping parameter.
#' @param S passed to \code{uwerr}, see documentation of \code{\link{uwerr}}.
#' @param pl logical: if set to TRUE the function produces plots
#' @param stat_range range of data to be included in the analysis.
#' @param skip number of measurements to be discarded at the beginning of the
#' time series. \code{skip} has no effect if two or more replica are used, see
#' argument \code{nrep}.
#' @param iobs if there are several operators available (local-local,
#' local-smeared, etc.), then this labels these (for cmi format)
#' @param ind.vec index vector indexing the column numbers in cmicor to be used
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
#' @param nrep vector (N1, N2, ...) of replica length N1, N2. If missing it is
#' assumed that there is only one ensemble. If there are two or more replica
#' the parameter \code{skip} has no effect.
#' @param fit.routine The fit routine to be used. Default is \dQuote{gsl},
#' which uses the gnu scientific library \dQuote{gsl_multifit_fdfsolver} solver
#' to minimise the chisquare. All other values lead to the usage of R's
#' \link{optim} function. The latter choice might be significantly slower.
#' @param oldnorm If set to \dQuote{TRUE}, the old online measurement
#' normalisation of \dQuote{tmLQCD} prior to version 5.2.0 is used in order to
#' get correct values for the pion decay constant.
#' @return returns an object of \code{class} \code{ofit} with the following
#' items
#' 
#' \item{fitresult}{ result from the fit as returned by \code{\link{optim}} }
#' \item{fitresultpp}{ Fit result of the PP correlator only } \item{t1}{ lower
#' bound for the fitrange in time (t1,t2). Counting starts with 0.  }
#' \item{t2}{ upper bound for the fitrange in time (t1,t2). Counting starts
#' with 0.  } \item{N}{ number of measurements found in the data } \item{Time}{
#' Time extent found in the data } \item{fitdata}{ \code{\link{data.frame}}
#' containing the time values used in the fit, the averaged correlator and its
#' error and the value of Chi for each time value } \item{uwerrresultmps}{ the
#' result of the time series analysis for the lowest mass as carried out by
#' \code{\link{uwerr}} } \item{uwerrresultmpcac}{ the result of the time series
#' analysis for the PCAC mass carried out by \code{\link{uwerr}}, see details }
#' \item{effmass}{ effective masses in the pion channel } \item{matrix.size}{
#' size of the data matrix, copied from input } \item{boot}{ object returned by
#' the call to \code{\link{boot}} if \code{method} was set correspodingly.
#' Otherwise \code{NULL}.  } \item{tsboot}{ object returned by the call to
#' \code{\link{tsboot}} if \code{method} was set correspodingly. Otherwise
#' \code{NULL}.  } \item{method}{ error analysis method as copied from input }
#' \item{fit.routine}{ \code{fit.routine} as copied from input } \item{nrep}{
#' \code{nrep} as copied from input } \item{dpaopp}{ \code{\link{data.frame}}
#' containing the pcac masses computed not with a fit, but with the derivative
#' method for all time values in between \code{t1} and \code{t2} }
#' @author Carsten Urbach, \email{curbach@@gmx.de}
#' @seealso \code{\link{readcmicor}}, \code{\link{uwerr}},
#' @keywords optimize ts
#' @export onlinemeas
onlinemeas <- function(data, t1, t2, 
                       stat_range, 
                       S=1.5, pl=FALSE, skip=0,
                       iobs=1, ind.vec=c(1,3,4,5), mu=0.1, kappa=0.125,
                       boot.R=99, boot.l=10, tsboot.sim="geom",
                       method="uwerr", fit.routine="optim", nrep,
                       oldnorm = FALSE) {
  
  if(missing(data)) {
    stop("Error! Data is missing!")
  }
  if(missing(t1) || missing(t2)) {
    stop("Error! t1 and t2 must be specified!")
  }
  if( missing(stat_range) ){
    stat_range <- c(skip,nrow(data))
  }
  par <- numeric(2)
  sign <- +1.
  data <- data[data$V1<3,]
  data <- data[data$V2==iobs,]

  Time <-  2*max(data[,ind.vec[2]])
  Thalf <- max(data[,ind.vec[2]])
  T1 <- Thalf+1
  t1p1 <- (t1+1)
  t2p1 <- (t2+1)
  nrObs <- 2
  nrType <- 1
  Skip <- (skip*(T1)*nrType*nrObs+1)
  Length <- length(data[,ind.vec[3]])
  message("time =", Time, " Thalf =", Thalf, "\n")
#  Thalf <- Thalf + 1
  if(missing(nrep)) {
    nrep <- c(length(data[((Skip):Length),ind.vec[3]])/(nrObs*(T1)*nrType))
  }
  else {
    skip <- 0
    if(sum(nrep) != length(data[((Skip):Length),ind.vec[3]])/(nrObs*(T1)*nrType)) {
      stop("sum of replica differs from total no of measurements!")
    }
  }

  Z <- array(data[((Skip):Length),ind.vec[3]], 
             dim=c(nrObs*(T1)*nrType,(length(data[((Skip):Length),ind.vec[3]])/(nrObs*(T1)*nrType))))
  # negative times
  W <- array(data[((Skip):Length),ind.vec[4]], 
             dim=c(nrObs*(T1)*nrType,(length(data[((Skip):Length),ind.vec[4]])/(nrObs*(T1)*nrType))))

  rm(data)
  W <- getCor(T1=T1, W=W, Z=Z, type=c("cosh","sinh"))
  if(oldnorm) {
    W <- W*(2*kappa)^2
  }
  rm(Z)
#  print(t(W))

  # some effective pion masses
  eff <- effectivemass(from=(t1+1), to=(t2+1), Time, W[1:T1,] , pl=FALSE, S=1.5, nrep=nrep)

  mass.eff <- data.frame(t=eff$t, m=eff$mass, dm=eff$dmass)

  Cor <- rep(0., times=nrObs*T1)
  E <- rep(0., times=nrObs*T1)
  
  for(i in 1:(T1*nrObs)) {
    Cor[i] <- mean(W[(i),])
    tmpe <- try(uwerrprimary(W[(i),], pl=F, nrep=nrep)$dvalue, TRUE)
    if(!inherits(tmpe, "try-error")) E[i] = tmpe
    else {
      warning("error of correlator replaced by naive estimate!\n", call.=F)
      E[i] = sd(W[(i),])/sqrt(length(W[(i),]))
    }
  }

  # now the PP correlator first _only_ 
  par[2] <- eff$mass[1]
  par[1] <- sqrt(abs(Cor[(t1+1)]*exp(par[2]*(t1+1))))
  
  # Index vector of data to be used in the analysis
  ii <- c((t1p1):(t2p1))

  #BFGS
  massfit <- optim(par, ChiSqr.singleCor, method="BFGS", control=list(trace=0),Thalf=Thalf,
                   x=c((t1):(t2)), y=Cor[ii], err=E[ii], tr = (t2-t1+1), sign=sign)
  sfit.mass <- abs(massfit$par[2])
  sfit.fpi <- 2*kappa*2*mu/sqrt(2)*abs(massfit$par[1])/sqrt(sfit.mass^3)
  if(massfit$convergence!=0) {
    warning("optim did not converge for massfit! ", massfit$convergence)
  }
  sfit.fpi <- 2*kappa*2*mu/sqrt(2)*abs(massfit$par[1])/sqrt(sfit.mass^3)
#    sfit.fpi <- 2*mu/sqrt(2)*abs(massfit$par[1])/sqrt(sfit.mass^3)
  message("mpi =", sfit.mass, " fpi =",sfit.fpi, "\n")
  
  sfit.dof <- (t2-t1+1)-length(massfit$par)
  sfit.chisqr <- massfit$value

  sfit.uwerrm <- NULL
  sfit.boot <- NULL
  sfit.tsboot <- NULL
  if(method == "uwerr" || method == "all") {
    sfit.uwerrm <- uwerr(f=fitmass, data=t(W[ii,]), S=S, pl=pl, nrep=nrep,
                        Time=Time, t1=t1, t2=t2, Err=E[ii], par=par, sign=sign,
                        fit.routine=fit.routine)
  }
  if(method == "boot" || method == "all") {
    sfit.boot <- boot::boot(data=t(W[ii,]), statistic=getfit.boot, R=boot.R, stype="i",
                     Time=Time, t1=t1, t2=t2, Err=E[ii], par=par, sign=sign,
                     fit.routine=fit.routine)

    sfit.tsboot <- boot::tsboot(tseries=t(W[ii,]), statistic=getfit.boot, R=boot.R, l=boot.l, sim=tsboot.sim,
                         Time=Time, t1=t1, t2=t2, Err=E[ii], par=par, sign=sign,
                         fit.routine=fit.routine)
  }

  
  # now do pcac mass, i.e. PP and PA
  # Index vector of data to be used in the analysis
  ii <- c((t1p1):(t2p1), (t1p1+T1):(t2p1+T1))

  par2 <- c(1, 0.1, 0.1)
  par2[3] <- eff$mass[1]
  par2[1] <- sqrt(abs(Cor[(t1+1)]*exp(par[2]*(t1+1))))
  par2[2] <- 0.
                                        #BFGS
  pcacfit <- optim(par2, ChiSqr.pcac, method="BFGS", control=list(trace=0),Thalf=Thalf,
                   x=c((t1):(t2)), y=Cor[ii], err=E[ii], tr = (t2-t1+1))
  fit.pcac <- 0.5*abs(pcacfit$par[3])*pcacfit$par[2]/pcacfit$par[1]

  if(pcacfit$convergence!=0) {
    warning("optim did not converge for pcacfit! ", massfit$convergence)
  }
  fit.dof <- (t2-t1+1)-length(massfit$par)
  fit.chisqr <- massfit$value

  fit.uwerrm <- NULL
  fit.uwerrpcac <- NULL
  fit.uwerrfpi <- NULL
  fit.boot <- NULL
  fit.tsboot <- NULL
  if(method == "uwerr" || method == "all") {
    fit.uwerrm <- uwerr(f=fitmass.online, data=t(W[ii,]), S=S, pl=pl, nrep=nrep,
                      Time=Time, t1=t1, t2=t2, Err=E[ii], par=par2,
                      fit.routine=fit.routine)
    fit.uwerrpcac <- uwerr(f=fitmpcac.online, data=t(W[ii,]), S=S, pl=pl, nrep=nrep,
                           Time=Time, t1=t1, t2=t2, Err=E[ii], par=par2,
                           fit.routine=fit.routine)
    fit.uwerrfpi <- uwerr(f=fitf.online, data=t(W[ii,]), S=S, pl=pl, nrep=nrep,
                      Time=Time, t1=t1, t2=t2, Err=E[ii], par=par2,
                      fit.routine=fit.routine)
  }
  if(method == "boot" || method == "all") {
    fit.boot <- boot::boot(data=t(W[ii,]), statistic=fit.online.boot, R=boot.R, stype="i",
                           Time=Time, t1=t1, t2=t2, Err=E[ii], par=par2,
                           fit.routine=fit.routine)

    fit.tsboot <- boot::tsboot(tseries=t(W[ii,]), statistic=fit.online.boot, R=boot.R, l=boot.l, sim=tsboot.sim,
                               Time=Time, t1=t1, t2=t2, Err=E[ii], par=par2, 
                               fit.routine=fit.routine)
  }

  Chi <- rep(0., times=2*T1)
  Fit <- rep(0., times=2*T1)

  jj <-  c(t1p1:t2p1)
  Fit[jj] <- pcacfit$par[1]^2*CExp(m=pcacfit$par[3], Time=2*Thalf, x=jj-1, sign=+1.)
  Fit[jj+T1] <- pcacfit$par[1]*pcacfit$par[2]*CExp(m=pcacfit$par[3], Time=2*Thalf, x=jj-1, sign=-1.)
  Chi[ii] <- (Fit[ii]-Cor[ii])/E[ii]

  dpaopp <- data.frame(t = array(0.,dim=c((t2-t1+1))), mass = array(0.,dim=c((t2-t1+1))),
                        dmass = array(0.,dim=c((t2-t1+1))), ddmass = array(0.,dim=c((t2-t1+1))),
                        tauint = array(0.,dim=c((t2-t1+1))), dtauint = array(0.,dim=c((t2-t1+1))))
  i <- 1
  for(t in t1p1:t2p1) {
    mass <- uwerrderived(pcacsym.online, data=t(W), S=S, pl=FALSE, t=t, T1=T1)
    dpaopp$t[i] <- t-1
    dpaopp$mass[i] <- mass$value[1]
    dpaopp$dmass[i] <- mass$dvalue[1]
    dpaopp$ddmass[i] <- mass$ddvalue[1]
    dpaopp$tauint[i] <- mass$tauint[1]
    dpaopp$dtauint[i] <- mass$dtauint[1]
    i=i+1
  }

  MChist.dpaopp <- rep(0., times=length(W[1,]))
  for(i in 1:length(W[1,])) {
    MChist.dpaopp[i] <- 0.
    for(t in t1p1:t2p1) {
      MChist.dpaopp[i] <- MChist.dpaopp[i] + pcacsym.online(data=W[,i], t=t, T1=T1)
    }
    # (t2p1-t1p1+1) values for this fitrange!
    MChist.dpaopp[i] <- MChist.dpaopp[i]/(t2p1-t1p1+1)
  }
  
  res <- list(fitresult=pcacfit, fitresultpp=massfit, t1=t1, t2=t2, N=length(W[1,]), Time=Time,
              fitdata=data.frame(t=(jj-1), Fit=Fit[ii], Cor=Cor[ii], Err=E[ii], Chi=Chi[ii]),
              uwerrresultmps=fit.uwerrm, 
              uwerrresultmpcac=fit.uwerrpcac, 
              uwerrresultfps=fit.uwerrfpi, 
              boot=fit.boot, tsboot=fit.tsboot, method=method, skip=skip,
              effmass=mass.eff, fit.routine=fit.routine, dpaopp=dpaopp, MChist.dpaopp=MChist.dpaopp,
              iobs=iobs, mu=mu, kappa=kappa,
              nrep=nrep, matrix.size=2)
  attr(res, "class") <- c("ofit", "list")  
  return(invisible(res))
}

fitmpcac.online <- function(Cor, Err, t1, t2, Time, par=c(1.,0.1,0.12),
                            N=2, no.masses=1, no=1, kappa, mu,
                            fit.routine="optim") {
  Thalf <- Time/2
  T1 <- Thalf+1
  t1p1 <- (t1+1)
  t2p1 <- (t2+1)
  tr <- (t2-t1+1)

  if(no.masses == 1) {
    fit <- optim(par, ChiSqr.pcac, method="BFGS", Thalf=Thalf,
                 x=c((t1):(t2)), y=Cor, err=Err, tr=tr)

    sort.ind <- c(1)
  }
  return(0.5*abs(fit$par[3])*fit$par[2]/fit$par[1])
}

fitmass.online <- function(Cor, Err, t1, t2, Time, par=c(1.,0.1,0.12),
                            N=2, no.masses=1, no=1, kappa, mu,
                            fit.routine="optim") {
  Thalf <- Time/2
  T1 <- Thalf+1
  t1p1 <- (t1+1)
  t2p1 <- (t2+1)
  tr <- (t2-t1+1)

  if(no.masses == 1) {
    fit <- optim(par, ChiSqr.pcac, method="BFGS", Thalf=Thalf,
                 x=c((t1):(t2)), y=Cor, err=Err, tr=tr)

    sort.ind <- c(1)
  }
  return(abs(fit$par[3]))
}

fitf.online <- function(Cor, Err, t1, t2, Time, par=c(1.,0.1,0.12),
                        N=2, no.masses=1, no=1, kappa, mu,
                        fit.routine="optim") {
  Thalf <- Time/2
  T1 <- Thalf+1
  t1p1 <- (t1+1)
  t2p1 <- (t2+1)
  tr <- (t2-t1+1)

  if(no.masses == 1) {
    fit <- optim(par, ChiSqr.pcac, method="BFGS", Thalf=Thalf,
                 x=c((t1):(t2)), y=Cor, err=Err, tr=tr)

    sort.ind <- c(1)
  }
  return(abs(fit$par[1])/sqrt(abs(fit$par[3]^3)))
}


pcacsym.online <- function(data, t, T1) {
  # PA[t+1]-PA[t-1]/4PP[t]
  return((data[t+1+T1]-data[t-1+T1])/(4.*data[t]))
}

fit.online.boot <- function(Z, d, Err, t1, t2, Time, par=c(1.,0.1,0.12),
                          kappa, mu, fit.routine="optim") {
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

  fit <- optim(par, ChiSqr.pcac, method="BFGS", Thalf=Thalf,
               x=c((t1):(t2)), y=Cor, err=Err, tr=tr)
  #fit.fpi <- 2*kappa*2*mu/sqrt(2)*abs(fit$par[(sort.ind[1]-1)*(N+1)+1])/sqrt(abs(fit$par[sort.ind[1]*(N+1)])^3)
  fit.pcac <- abs(fit$par[3])*fit$par[2]/fit$par[1]/2.
  #fit.zv <- 2.*mu/abs(fit$par[sort.ind[1]*(N+1)])*fit$par[(sort.ind[1]-1)*(N+1)+1]/fit$par[(sort.ind[1]-1)*(N+1)+5]
  return(c(abs(fit$par[3]), fit.pcac, fit$par[1:2]))
}
