source("analyse_pipi.R")

compute.weights <- function(err, pvalues) {
  return(pvalues^2 * min(err, na.rm=TRUE)^2/err^2)
}

compute.boots <- function(res, index=1, piononly=FALSE) {
  pvalues <- (1-2*abs(res[1,,5]-0.5)) * (1-2*abs(res[1,,6]-0.5))
  if(piononly) pvalues <- (1-2*abs(res[1,,6]-0.5))
  err <- apply(res[,,index], 2, sd, na.rm=TRUE)
  w <- compute.weights(err, pvalues)
  return(apply(res[,,index], 1, weighted.quantile, prob=c(0.5), w=w, na.rm=TRUE))
}


estimate.error <- function(res, index=1, prob=c(0.1573, 0.8427), main) {
  pvalues <- as.vector((1-2*abs(res[1,,5]-0.5)) * (1-2*abs(res[1,,6]-0.5)))

  err <- apply(res[,,index], 2, sd, na.rm=TRUE)
  w <- compute.weights(err, pvalues)
  x <- weighted.quantile(as.vector(res[1,,index]), prob=c(0.5), w=w, na.rm=TRUE)

  weighted.hist(x=as.vector(res[1,,index]), w=w, main=main, na.rm=TRUE)
  ## statistical
  x[2] <- sd(apply(res[,,index], 1, weighted.quantile, prob=c(0.5), w=w, na.rm=TRUE), na.rm=TRUE)
  ## systematic error
  ## lower and upper
  x[c(3:4)] <- weighted.quantile(as.vector(res[1,,index]), w=w, prob=prob, na.rm=TRUE)-x[1]
  return(x)
}

compile.ratio.sldata <- function(ens, path="./") {

  filelist <- Sys.glob(paste(path, "data.", ens, "*.Rdata", sep=""))
  N <- length(filelist)
  cat("processing", N, "data files for ensemble", ens, "\n")
  load(filelist[1])
  R <- pipi.data$boot.R
  rr <- c(2:(R+1))
  L <- pipi.data$L
  res <- array(0, dim=c(R+1, N, 9))
  for(i in c(1:N)) {
    load(filelist[i])

    ## deltaE
    res[1,i,1] <- pipi.data$opt.res$par[2]
    res[rr,i,1] <- pipi.data$opt.tsboot[,2]
    ## a0
    res[1,i,2] <- pipi.data$a0
    res[rr,i,2] <- pipi.data$opt.tsboot[,5]
    ## mpia0
    res[1,i,3] <- pipi.data$mpia0
    res[rr,i,3] <- pipi.data$opt.tsboot[,4]
    ## mpi
    res[1,i,4] <- pipi.data$pion.effmass$opt.res$par[1]
    res[rr,i,4] <- pipi.data$pion.effmass$massfit.tsboot[,1]
    ## Qval
    res[1,i,5] <- pipi.data$Qval
    res[rr,i,5] <- 1-pchisq(pipi.data$opt.tsboot[,3], pipi.data$dof)
    ## pion Qval
    res[1,i,6] <- pipi.data$pionQval
    res[rr,i,6] <- 1-pchisq(pipi.data$pion.effmass$massfit.tsboot[,2], pipi.data$pion.effmass$dof)

    qtilde <- compute.qtildesq(E=res[1,i,1]+2*res[1,i,4], L=L, dvec=c(0,0,0), mpi=res[1,i,4])
    ## qcotdelta
    res[1,i,7] <- 2.*Re( LuescherZeta(qtilde$qtsq, l=0, m=0, gamma=qtilde$gammaboost, dvec=c(0,0,0)) )/(qtilde$gammaboost*L*sqrt(pi))
    ## Ecm
    res[1,i,8] <- qtilde$Ecm
    ## q^2
    res[1,i,9] <- qtilde$q^2
    qtildeboot <- compute.qtildesq(E=res[rr,i,1]+2*res[rr,i,4], L=L, dvec=c(0,0,0), mpi=res[rr,i,4])
    res[rr,i,7] <- 2.*Re( LuescherZeta(qtildeboot$qtsq, l=0, m=0, gamma=qtildeboot$gammaboost, dvec=c(0,0,0)) )/(qtildeboot$gammaboost*L*sqrt(pi))
    res[rr,i,8] <- qtildeboot$Ecm
    res[rr,i,9] <- qtildeboot$q^2
  }

  save(res, L, file=paste(path, "res-ratio.", ens, ".Rdata", sep=""))
  return(invisible(res))
}

## extract scattering length from Luescher formula using
## energy shift values determined from effective mass fits
## to 2pt and 4pt functions

compile.efm.sldata <- function(path="./", data, L, ens) {
  R <- dim(data)[1]
  Nd <- dim(data)[c(2:3)]

  res <- array(0, dim=c(R, Nd[1]*Nd[2], 6))
  
  for(i in c(1:Nd[1])) {
    for(j in c(1:Nd[2])) {
      ## deltaE
      res[,(i-1)*Nd[2]+j,1] <- data[,i,j,6] - 2*data[,i,j,5]
      ## mpi
      res[,(i-1)*Nd[2]+j,4] <- data[,i,j,5]
      ## Qval
      res[,(i-1)*Nd[2]+j,5] <- data[,i,j,8]
      ## pion Qval
      res[,(i-1)*Nd[2]+j,6] <- data[,i,j,7]
    }
  }
  for(i in c(1:R)) {
    for(j in c(1:(Nd[1]*Nd[2]))) {
      if(res[i,j,1] < 0) res[i,j,2] <- NA
      else res[i,j,2] <- uniroot(deltaEvL, c(0.,-6.), deltaE=res[i,j,1], L=L, m=res[i,j,4])$root
      res[i,j,3] <- res[i,j,2]*res[i,j,4]
    }
  }
  save(res, file=paste(path, "res-efm.", ens, ".Rdata", sep=""))
  return(invisible(res))
}


compute.error.piL <- function(res, prob=c(0.1573, 0.8427)) { 
  deltaE <- estimate.error(res, index=1, prob=prob, main="deltaE")
  a0 <- estimate.error(res, index=2, prob=prob, main="a0")
  mpia0 <- estimate.error(res, index=3, prob=prob, main="mpia0")
  qcotdelta <- estimate.error(res, index=7, prob=prob, main="qcotdelta")
  Ecm <- estimate.error(res, index=8, prob=prob, main="Ecm")
  qsq <- estimate.error(res, index=9, prob=prob, main="q^2")
  return(list(deltaE=deltaE, a0=a0, mpia0=mpia0, qcotdelta=qcotdelta, Ecm=Ecm, qsq=qsq))
}
