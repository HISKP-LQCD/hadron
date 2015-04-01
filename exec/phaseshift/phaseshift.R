## omega_lm
## this needs qtilde as input!
omegalm <- function(l=0, m=0, q, gamma=1, dvec=c(0,0,0)) {
  return( LuescherZeta(qsq = q^2, l=l, m=m, gamma=gamma, dvec=dvec)/(pi^(3/2)*sqrt(2*l+1)*q^(l+1)*gamma) )
}

## center of mass energy from q^2 and Mpi
Ecmofqsq <- function(q, Mpi) {
  return(2*acosh(cosh(Mpi) +  2*sin(q/2)^2 ))
}

## energy level in MF from Ecm
EofEcm <- function(Ecm, dvec, L) {
  return(acosh( cosh(Ecm) + 2*sum(sin(pi*dvec/L)^2) ))
}

## computes the boost factor as a function of q^2
gammaofqsq <- function(q, dvec, Mpi, L) {
  Ecm <- Ecmofqsq(q, Mpi)
  ## E/Ecm
  return(EofEcm(Ecm, dvec, L)/Ecm)
}

## computes the energy level given q^2
Eofqsq <- function(q, dvec, Mpi, L) {
  Ecm <- Ecmofqsq(q, Mpi)
  return(EofEcm(Ecm, dvec, L))
}

prepdetEqCMscan <- function(q, Mpi, L) {
  W <- list()
  W$dvec <- c(0,0,0)
  W$gamma <- gammaofqsq(q, dvec=W$dvec, Mpi=Mpi, L=L)
  W$qt <- L*q/2/pi
  W$L <- L
  W$Mpi <- Mpi
  W$q <- q

  W$w00 <- Re(omegalm(l=0, m=0, q=W$qt, gamma=W$gamma))
  return(W)
}

detEqCMscan <- function(par, W) {
  return(par[1] + 0.5*par[2]*W$q^2 - W$q*W$w00)
}

detEqCM <- function(q, par, Mpi, L) {
  gamma <- rep(1., times=(length(q)))
  qt <- L*q/2/pi
  return(par[1] + 0.5*par[2]*q^2 - q*Re(omegalm(l=0, m=0, q=qt, gamma=gamma)))
}

prepdetEqMF1scan <- function(q, Mpi, L) {
  W <- list()
  W$dvec <- c(0,0,1)
  W$gamma <- gammaofqsq(q, dvec=W$dvec, Mpi=Mpi, L=L)
  W$qt <- L*q/2/pi
  W$L <- L
  W$Mpi <- Mpi
  W$q <- q

  W$w00 <- Re(omegalm(l=0, m=0, q=W$qt, gamma=W$gamma, dvec=W$dvec))
  W$w20 <- Re(omegalm(l=2, m=0, q=W$qt, gamma=W$gamma, dvec=W$dvec))
  W$w40 <- Re(omegalm(l=4, m=0, q=W$qt, gamma=W$gamma, dvec=W$dvec))

  return(W)
}

detEqMF1scan <- function(par, W) {
  cd0 <- par[1] + 0.5*par[2]*W$q^2
  cd2 <- par[3]
  return( (cd0 - W$q*W$w00)*
         (cd2 - W$q^5*(W$w00 + 10/7*W$w20 + 18/7*W$w40))
         -W$q^6*5*W$w20^2
         )
}

detEqMF1 <- function(q, par, Mpi, L) {
  dvec <- c(0,0,1)

  gamma <- gammaofqsq(q, dvec=dvec, Mpi=Mpi, L=L)

  qt <- L*q/2/pi
  w00 <- Re(omegalm(l=0, m=0, q=qt, gamma=gamma, dvec=dvec))
  w20 <- Re(omegalm(l=2, m=0, q=qt, gamma=gamma, dvec=dvec))
  w40 <- Re(omegalm(l=4, m=0, q=qt, gamma=gamma, dvec=dvec))

  cd0 <- par[1] + 0.5*par[2]*q^2
  cd2 <- par[3]
  return( (cd0 - q*w00)*
         (cd2 - q^5*(w00 + 10/7*w20 + 18/7*w40))
         -q^6*5*w20^2
         )
}

prepdetEqMF2scan <- function(q, Mpi, L) {
  W <- list()
  W$dvec <- c(1,1,0)
  W$gamma <- gammaofqsq(q, dvec=W$dvec, Mpi=Mpi, L=L)
  W$qt <- L*q/2/pi
  W$L <- L
  W$Mpi <- Mpi
  W$q <- q

  W$w00 <- Re(omegalm(l=0, m=0, q=W$qt, gamma=W$gamma, dvec=W$dvec))
  W$w20 <- Re(omegalm(l=2, m=0, q=W$qt, gamma=W$gamma, dvec=W$dvec))
  W$w40 <- Re(omegalm(l=4, m=0, q=W$qt, gamma=W$gamma, dvec=W$dvec))
  W$w22 <- omegalm(l=2, m=2, q=W$qt, gamma=W$gamma, dvec=W$dvec)
  W$w44 <- Re(omegalm(l=4, m=4, q=W$qt, gamma=W$gamma, dvec=W$dvec))
  W$w42 <- omegalm(l=4, m=2, q=W$qt, gamma=W$gamma, dvec=W$dvec)

  return(W)
}

detEqMF2scan <- function(par, W) {
  cd0 <- (par[1] + 0.5*par[2]*W$q^2)/W$q
  cd2 <- par[3]/W$q^5

  return(Re(q^11*(10/49*(10*(-2*cd0 + 2*W$w00 + 7*W$w20)*W$w22^2 +
                         3*sqrt(15)*(4*cd0 - 4*W$w00 - 7*W$w20)*W$w22*W$w42 + 27*(W$w00 - cd0)*W$w42^2)
                  +(-cd2 + W$w00 + 2/7*(5*W$w20 + 9*W$w40))*(10*W$w22^2 + 1/7*(cd0 - W$w00)*(7*cd2 - 7*W$w00 + 10*W$w20 - 3*W$w40 + 3*sqrt(70)*W$w44)) +
                  5/7*W$w20*(20*W$w22^2 - 6*sqrt(15)*W$w22*W$w42 + W$w20*(7*cd2 - 7*W$w00 + 10*W$w20 - 3*W$w40 + 3*sqrt(70)*W$w44))
                  )
            ))
}

detEqMF2 <- function(q, par, Mpi, L) {
  dvec=c(1,1,0)
  gamma <- gammaofqsq(q, dvec=dvec, Mpi=Mpi, L=L)

  qt <- L*q/2/pi
  w00 <- Re(omegalm(l=0, m=0, q=qt, gamma=gamma, dvec=dvec))
  w20 <- Re(omegalm(l=2, m=0, q=qt, gamma=gamma, dvec=dvec))
  w22 <- omegalm(l=2, m=2, q=qt, gamma=gamma, dvec=dvec)
  w40 <- Re(omegalm(l=4, m=0, q=qt, gamma=gamma, dvec=dvec))
  w44 <- Re(omegalm(l=4, m=4, q=qt, gamma=gamma, dvec=dvec))
  w42 <- omegalm(l=4, m=2, q=qt, gamma=gamma, dvec=dvec)

  cd0 <- (par[1] + 0.5*par[2]*q^2)/q
  cd2 <- par[3]/q^5
  
  return(Re(q^11*(10/49*(10*(-2*cd0 + 2*w00 + 7*w20)*w22^2 +
                      3*sqrt(15)*(4*cd0 - 4*w00 - 7*w20)*w22*w42 + 27*(w00 - cd0)*w42^2)
               +(-cd2 + w00 + 2/7*(5*w20 + 9*w40))*(10*w22^2 + 1/7*(cd0 - w00)*(7*cd2 - 7*w00 + 10*w20 - 3*w40 + 3*sqrt(70)*w44)) +
               5/7*w20*(20*w22^2 - 6*sqrt(15)*w22*w42 + w20*(7*cd2 - 7*w00 + 10*w20 - 3*w40 + 3*sqrt(70)*w44))
               )
         ))
}

prepdetEqMF3scan <- function(q, Mpi, L) {
  W <- list()
  W$dvec <- c(1,1,1)
  W$gamma <- gammaofqsq(q, dvec=W$dvec, Mpi=Mpi, L=L)
  W$qt <- L*q/2/pi
  W$L <- L
  W$Mpi <- Mpi
  W$q <- q

  W$w00 <- Re(omegalm(l=0, m=0, q=W$qt, gamma=W$gamma, dvec=W$dvec))
  W$w40 <- Re(omegalm(l=4, m=0, q=W$qt, gamma=W$gamma, dvec=W$dvec))
  W$w22 <- omegalm(l=2, m=2, q=W$qt, gamma=W$gamma, dvec=W$dvec)
  W$w42 <- omegalm(l=4, m=2, q=W$qt, gamma=W$gamma, dvec=W$dvec)

  return(W)
}

detEqMF3scan <- function(par, W) {
  cd0 <- (par[1] + 0.5*par[2]*W$q^2)
  cd2 <- par[3]

  return(Re(( cd0 - W$q*W$w00 )*
            ( cd2 - W$q^5*(W$w00 - 12/7*W$w40 - 12*sqrt(10)/7*1i*W$w42 - 10*sqrt(6)/7*1i*W$w22)) +
            W$q^6*30*W$w22^2
            )
         )
}

detEqMF3 <- function(q, par, Mpi, L) {
  dvec=c(1,1,1)
  gamma <- gammaofqsq(q, dvec=dvec, Mpi=Mpi, L=L)

  qt <- L*q/2/pi
  w00 <- Re(omegalm(l=0, m=0, q=qt, gamma=gamma, dvec=dvec))
  w40 <- Re(omegalm(l=4, m=0, q=qt, gamma=gamma, dvec=dvec))
  w22 <- omegalm(l=2, m=2, q=qt, gamma=gamma, dvec=dvec)
  w42 <- omegalm(l=4, m=2, q=qt, gamma=gamma, dvec=dvec)

  cd0 <- (par[1] + 0.5*par[2]*q^2)
  cd2 <- par[3]
  return(Re(( cd0 - q*w00 )*
            ( cd2 - q^5*(w00 - 12/7*w40 - 12*sqrt(10)/7*1i*w42 - 10*sqrt(6)/7*1i*w22)) +
            q^6*30*w22^2
            )
         )
}

findSignChanges <- function(fn, makeplot=FALSE, par, W, no=3, threshold=100) {
  res <- fn(par=par, W=W)
  if(makeplot) {
    if(interactive()) X11()
    plot(W$q^2, res, type="l", ylim=c(-threshold, threshold))
    abline(a=0, b=0, col="red")
  }
  
  ii <- integer(0)
  j <- 1
  for(i in c(1:(length(W$q)-1))) {
    ## find all sign-changes 
    if(res[i]/abs(res[i]) != res[i+1]/abs(res[i+1])) {
      ## remove the poles
      if((res[i] < 0 && res[i-1] < res[i]) ||
         (res[i] > 0 && res[i-1] > res[i])) {
        ii[j] <- i
        j <- j+1
        if(j == no+1) break
      }
    }
  }
  return(ii)
}

findZeros <- function(fn, q, ii, makeplot=FALSE, tol=1.e-12, ...) {

  ## now determine the root more precisely
  zeros <- numeric(0.)
  for(j in c(1:length(ii))) {
    z <- uniroot(fn, interval=c(q[ii[j]], q[ii[j]+1]), tol = tol, ...)
    if(z$f.root < 1) {
      zeros[j] <- z$root
##      cat("root", j, ":", z$root, z$f.root, z$estim.prec, "\n")
    }
  }
  if(makeplot) {
    abline(v=zeros^2, col="blue")
  }

  return(zeros)
}

