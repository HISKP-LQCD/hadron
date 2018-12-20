## sum_i (a_i cosh(m_i*t))
## a_i are the amplitudes, m_i the masses and t is a vector of times
sum.cosh <- function(masses, amplitudes, t){
  colSums(amplitudes*cosh(masses%o%t))
}

## extract the effective mass from the sum of coshs
cosh.to.effmass <- function(masses, amplitudes, t, Thalf, type){
  if(type == "solve"){
    sapply(t, FUN=function(x) {invcosh(sum.cosh(masses, amplitudes, x-Thalf)/sum.cosh(masses, amplitudes, x+1-Thalf), timeextent=2*Thalf, t=x+1)})
  }else if(type == "acosh"){
    acosh((sum.cosh(masses, amplitudes, t-Thalf-1)+sum.cosh(masses, amplitudes, t-Thalf+1))/(2*sum.cosh(masses, amplitudes, t-Thalf)))
  }
  ## the other types might be included as well...
}

## this fit-function works completely analogous to fit.effecitvemass
## but it fits the correlator directly to a sum of cosh-terms
fit.cosh <- function(effMass, cf, t1, t2, useCov=FALSE, m.init, par, n.cosh=2, adjust.n.cosh=FALSE, every, ...){
  if(missing(effMass) && missing(cf)){
    stop("effMass or cf has to be provided!\n")
  }
  if(!missing(effMass)){
    stopifnot(inherits(effMass, 'effectivemass'))
    cf = effMass$cf
    Thalf = effMass$Time/2
    use.effmass <- TRUE
  }else{
    stopifnot(inherits(cf, 'cf_meta'))
    stopifnot(inherits(cf, 'cf_boot'))
    Thalf = cf$Time/2
    use.effmass <- FALSE
  }

  stopifnot(t1+2*n.cosh < t2)
  stopifnot(0 <= t1)

  cf.save <- cf$cf.tsboot$t
  cf0.save <- cf$cf0
  
  if(missing(par)){
    ## some initial values, probably not optimal...
    if(missing(m.init)){
      if(use.effmass){
        m.init = effMass$effMass[t2]
      }else{
        m.init = log(cf0.save[t2]/cf0.save[t2+1])
        if(is.na(m.init)){
          stop("Please provide m.init or par. The fit does not converge stably otherwise.\n")
        }
      }
    }
    if(n.cosh == 2){
      C1 = cf0.save[t1+1]
      C2 = cf0.save[t2+1]
      CMiddle = cf0.save[t1+2]
      a0 = C2/cosh(m.init*(t2-Thalf))
      C1 = C1-a0*cosh(m.init*(t1-Thalf))
      CMiddle = CMiddle-a0*cosh(m.init*(t1+1-Thalf))
      masses = c(m.init, acosh(C1/CMiddle))
      amplitudes = c(a0, C1/cosh(masses[[2]]*(t1-Thalf)))
      if(any(is.na(masses)) || any(is.na(amplitudes))){
        cat("The higher cosh-term cannot be resolved properly.\nA plateau fit is recommended.\n")
        if(adjust.n.cosh){
          cat("Changing number of cosh-terms. Now n.cosh=1.\n")
          n.cosh = 1
          masses = c(m.init)
          amplitudes = c(a0)
        }else{
          masses = c(0:(n.cosh-1)) + m.init
          amplitudes = cf0.save[t1+1]/n.cosh/cosh(masses*(t1-Thalf))
        }
      }
    }else{
      masses = c(0:(n.cosh-1)) + m.init
      amplitudes = cf0.save[t1+1]/n.cosh/cosh(masses*(t1-Thalf))
    }
  }else{
    masses = par[1:n.cosh]
    amplitudes = par[(n.cosh+1):(2*n.cosh)]
  }

  ## extract the relevant time slices
  ii <- c((t1+1):(t2+1))
  if(!missing(every)){
    ii <- ii[which(ii%%every == 0)]
  }
  if(any(is.na(cf0.save[ii]))){
    ii.na <- which(is.na(cf0.save[ii]))
    ii <- ii[-ii.na]
  }
  tt = ii-1-Thalf

  ## function to fit
  sum.cosh.fit <- function(par, x, ...) {
    sum.cosh(par[1:n.cosh], abs(par[(n.cosh+1):(2*n.cosh)]), x)
  }
  ## corresponding Jacobian
  sum.cosh.jac <- function(par, x, ...) {
    df.dm <- x * t(abs(par[(n.cosh+1):(2*n.cosh)])* sinh(par[1:n.cosh]%o%x)) 
    df.da <- t(sign(par[(n.cosh+1):(2*n.cosh)]) * cosh(par[1:n.cosh]%o%x))
    return(cbind(df.dm, df.da))
  }

  dy <- cf$tsboot.se[ii]
  boot.R <-cf$boot.R
  if (useCov) {
    fit.res <- bootstrap.nlsfit(fn = sum.cosh.fit,
                                par.guess = c(masses, amplitudes),
                                y = cf0.save[ii],
                                x = tt,
                                bsamples = cf.save[,ii],
                                gr = sum.cosh.jac,
                                CovMatrix = cov(cf.save[,ii]),
                                ...)
  } else {
    fit.res <- bootstrap.nlsfit(fn = sum.cosh.fit,
                                par.guess = c(masses, amplitudes),
                                y = cf0.save[ii],
                                x = tt,
                                bsamples = cf.save[,ii],
                                gr = sum.cosh.jac,
                                ...)
  }

  opt.res <- fit.res$t0[1:(2*n.cosh)]
  opt.res[1:n.cosh] <- abs(opt.res[1:n.cosh])
  opt.res[(n.cosh+1):(2*n.cosh)] <- abs(opt.res[(n.cosh+1):(2*n.cosh)][order(opt.res[1:n.cosh])])
  opt.res[1:n.cosh] <- opt.res[1:n.cosh][order(opt.res[1:n.cosh])]

    ## now we bootstrap the fit
  massfit.tsboot <- abs(fit.res$t)
  if(n.cosh >= 2){
    massfit.tsboot[,(n.cosh+1):(2*n.cosh)] <- t(apply(massfit.tsboot, FUN=function(res) {res[(n.cosh+1):(2*n.cosh)][order(res[1:n.cosh])]}, MARGIN=1))
    massfit.tsboot[,1:n.cosh] <- t(apply(massfit.tsboot[,1:n.cosh], FUN=function(res) {res[order(res)]}, MARGIN=1))
  }

  if(!use.effmass){
    effMass <- list()
    effMass$cf <- cf
  }
  effMass$coshfit <- list()
  effMass$coshfit$use.effmass <- use.effmass
  effMass$coshfit$t0 <- opt.res
  effMass$coshfit$t <- massfit.tsboot
  effMass$coshfit$se <- fit.res$se
  effMass$coshfit$n.cosh <- n.cosh
  effMass$t1 <- t1
  effMass$coshfit$t1 <- t1
  effMass$t2 <- t2
  effMass$coshfit$t2 <- t2
  effMass$coshfit$ii <- ii
  effMass$coshfit$useCov <- useCov
  effMass$coshfit$invCovMatrix <- fit.res$invCovMatrix
  effMass$dof <- length(ii)-2*n.cosh
  effMass$chisqr <- fit.res$chisqr
  effMass$coshfit$dof <- length(ii)-2*n.cosh
  effMass$coshfit$chisqr <- fit.res$chisqr
  effMass$coshfit$Qval <- fit.res$Qval
  attr(effMass, "class") <- c("coshfit", class(effMass))
  return(invisible(effMass))
}

plot.coshfit <- function(effMass, col.fitline="black", plot.mass=TRUE, plot.corr=FALSE, ...) {
  stopifnot(inherits(effMass, 'coshfit'))
  if(!inherits(effMass, 'effectivemass')){
    plot.mass = FALSE
    if(missing(plot.corr)){
      plot.corr = TRUE
    }
  }

  t <- c(effMass$t1:effMass$t2)
  Thalf <- effMass$cf$Time/2
  n.cosh <- effMass$coshfit$n.cosh

  pcol <- col2rgb(col.fitline,alpha=TRUE)/255
  pcol[4] <- 0.65
  pcol <- rgb(red=pcol[1],green=pcol[2],blue=pcol[3],alpha=pcol[4])

  if(plot.mass){
    t.all <- effMass$t.idx
    op <- options()
    options(warn=-1)
    plotwitherror(x=t.all-1, y=effMass$effMass[t.all], dy=effMass$deffMass[t.all], ...)
    options(op)

    if(!is.null(effMass$coshfit)){
      Y <- cosh.to.effmass(effMass$coshfit$t0[1:n.cosh], effMass$coshfit$t0[(n.cosh+1):(2*n.cosh)], t, Thalf, type=effMass$type)
      Y.boot <- apply(effMass$coshfit$t, FUN=function(x) {cosh.to.effmass(x[1:n.cosh], x[(n.cosh+1):(2*n.cosh)], t, Thalf, type=effMass$type)}, MARGIN=1)
      se <- apply(Y.boot, MARGIN=1, FUN=sd, na.rm=TRUE)

      ## plot it
      polyval <- c(Y+se, rev(Y-se))
      polygon(x=c(t, rev(t)), y=polyval, col=pcol, lty=0, lwd=0.001, border=pcol)

      ## plot the fitted curve on top
      lines(x=t, y=Y, col=col.fitline, lwd=1.3)
    }
  }

  if(plot.corr){
    plot(effMass$cf, ...)

    if(!is.null(effMass$coshfit)){
      Y <- sum.cosh(effMass$coshfit$t0[1:n.cosh], effMass$coshfit$t0[(n.cosh+1):(2*n.cosh)], t-Thalf)
      Y.boot <- apply(effMass$coshfit$t, FUN=function(x) {sum.cosh(x[1:n.cosh], x[(n.cosh+1):(2*n.cosh)], t-Thalf)}, MARGIN=1)
      se <- apply(Y.boot, MARGIN=1, FUN=sd, na.rm=TRUE)

      ## plot it
      polyval <- c(Y+se, rev(Y-se))
      polygon(x=c(t, rev(t)), y=polyval, col=pcol, lty=0, lwd=0.001, border=pcol)

      ## plot the fitted curve on top
      lines(x=t, y=Y, col=col.fitline, lwd=1.3)
    }
  }
}

summary.coshfit <- function(effMass, verbose=FALSE) {
  cat("\n ** Result of", effMass$coshfit$n.cosh, "-fold cosh-fit **\n\n")
  cat("no. measurements\t=\t", length(effMass$cf$cf[,1]), "\n")
  cat("boot.R\t=\t", effMass$cf$boot.R, "\n")
  cat("boot.l\t=\t", effMass$cf$boot.l, "\n")
  cat("Time extend\t=\t", effMass$cf$Time, "\n")
  cat("NA count in fitted bootstrap samples:\t", length(which(is.na(effMass$cf$cf.tsboot$t[,effMass$coshfit$ii]))),
      "(",100*length(which(is.na(effMass$cf$cf.tsboot$t[,effMass$coshfit$ii])))/ length(effMass$cf$cf.tsboot$t[,effMass$coshfit$ii]), "%)\n")
  cat("time range from", effMass$coshfit$t1, " to ", effMass$coshfit$t2, "\n")
  if(verbose) {
    cat("\nmasses:\n")
    print(data.frame(m = effMass$coshfit$t0[1:effMass$coshfit$n.cosh], dm = effMass$coshfit$se[1:effMass$coshfit$n.cosh]))
    cat("\namplitudes:\n")
    print(data.frame(a = effMass$coshfit$t0[(effMass$coshfit$n.cosh+1):(2*effMass$coshfit$n.cosh)], da = effMass$coshfit$se[(effMass$coshfit$n.cosh+1):(2*effMass$coshfit$n.cosh)]))
    cat("\ncorrelations:\n")
    fit.cor <- cor(effMass$coshfit$t[,1:(2*effMass$coshfit$n.cosh)], effMass$coshfit$t[,1:(2*effMass$coshfit$n.cosh)], use="na.or.complete")
    print(fit.cor)
    cat("\n")
  }else{
    cat("m\t=\t", effMass$coshfit$t0[1], "\n")
    cat("dm\t=\t", effMass$coshfit$se[1], "\n")
  }
  cat("chisqr\t=\t", effMass$coshfit$chisqr, "\n")
  cat("dof\t=\t", effMass$coshfit$dof, "\n")
  cat("chisqr/dof=\t",
      effMass$coshfit$chisqr/effMass$coshfit$dof, "\n")
  cat("Quality of the fit (p-value):",   effMass$coshfit$Qval, "\n")

}
