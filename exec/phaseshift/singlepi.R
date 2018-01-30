singlepi <- function(t1, t2, T, redofit=TRUE, t.threshold=5, p="p0", path="./", srcpath="./", useCov=TRUE, verbose=TRUE, ending=".dat", boot.R=400, boot.l=1, ens, ind.vector=c(2,3)) {

  ## read pion data from file
  if(!file.exists(paste(path,"pidata.", p, ".Rdata", sep=""))) {
    pidata <- bootstrap.cf(readtextcf(paste("pi_corr_", p, ending, sep=""), T=T, check.t=1, path=srcpath, ind.vector=ind.vector), boot.R=boot.R, boot.l=boot.l, seed=seed)
    save(pidata, file=paste(path,"pidata.", p, ".Rdata", sep=""))
  }
  load(file=paste(path, "pidata.", p, ".Rdata", sep=""))

  if(redofit) {
    cat("...determining pion mass from effective mass and matrix fit\n")
    t1 <- t1
    t2 <- t2
    for(ta in c(t1:t2)) {
      if(t2-ta < t.threshold) break
      pion.effectivemass <- fit.effectivemass(bootstrap.effectivemass(cf=pidata, type="solve"), t1=ta, t2=t2-1, useCov=useCov)
      if(verbose) {
        summary(pion.effectivemass)
        plot(pion.effectivemass, xlab=c("t/a"), ylab=c("Meff"), main=paste("pion.", p, ".effectivemass.", ta, ".", t2, "p=", pion.effectivemass$Qval, sep=""))
      }
      ## save all of them
      meta.data <- list(ens=ens, T=T, srcpath=srcpath, t1=ta, t2=t2)
      save(meta.data, pion.effectivemass, file=paste(path, "pion.", p, ".effectivemass.", ta, ".", t2, ".Rdata", sep=""))
      pion.matrixfit <- matrixfit(pidata, t1=ta, t2=t2, useCov=useCov)
      save(pion.matrixfit, file=paste(path, "pion.", p, ".matrixfit.", ta, ".", t2, ".Rdata", sep=""))
    }
  }
}

pion.sys <- function(p = "p0", rep.outliers=FALSE) {

  compute.weights <- function(pvalues, err) {
    return(pvalues^2*min(err)^2/err^2)
  }

  filelist <- Sys.glob(paste("pion.", p, ".effectivemass.*.Rdata", sep=""))
  load(filelist[1])
  boot.R <- pion.effectivemass$boot.R

  res <- array(0., dim=c(boot.R+1, length(filelist), 2))
  for(i in c(1:length(filelist))) {
    load(filelist[i])
    if(pion.effectivemass$boot.R != boot.R) {
      stop("inconsistent boot.R!\n")
    }
    res[1,i,1] <- pion.effectivemass$opt.res$par[1]
    res[c(2:(boot.R+1)), i, 1] <- pion.effectivemass$massfit.tsboot[,1]
    res[1,i,2] <- pion.effectivemass$Qval
    res[c(2:(boot.R+1)), i, 2] <- 1-pchisq(pion.effectivemass$massfit.tsboot[,2], pion.effectivemass$dof)
  }

  ## this is the naive statistical uncertainty
  err <- apply(res[,,1], 2, sd, na.rm=TRUE)
  ## weights
  pvalues <- 1-2*abs(res[1,,2]-0.5)
  ## weights
  w <- compute.weights(pvalues, err)
  ## value
  mpi <- weighted.quantile(res[1,,1], w=w, prob=c(0.5), na.rm=TRUE)
  ## statistical error
  dmpi <- sd(apply(res[,,1], 1, weighted.quantile, w=w, prob=c(0.5), na.rm=TRUE), na.rm=TRUE)
  ## systematic error
  ## lower
  smpim <- weighted.quantile(res[1,,1], w=w, prob=c(0.1573), na.rm=TRUE)-mpi
  ## upper
  smpip <- weighted.quantile(res[1,,1], w=w, prob=c(0.8427), na.rm=TRUE)-mpi

  return(invisible(list(res=res, mpi=data.frame(mpi=mpi, dmpi=dmpi, smpip=smpip, smpim=smpim))))
}
