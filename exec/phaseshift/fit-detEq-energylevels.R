dirlist <- c("A40.32", "A40.24", "A40.20")

boot.R <- 1500
## cut in chisq
cut <- 0.01

compute.weights <- function(err, pvalues) {
##  return(1./pvalues^2)
  return(1./(pvalues^2 * err^2))
}

Epipi <- data.frame()
Epi <- data.frame()
reslist <- list()
maxcomb <- 0

boot.data <- array(0., dim=c(boot.R+1,2))
bdata <- array(0., dim=c(boot.R+1,2))
counter <- 1
rr <- c(2:(boot.R+1))

if(!file.exists("Energies.Rdata")) {
  for(k in c(1:length(dirlist))) {
    L <- 32
    frame <- 1
    level <- 1
    if(dirlist[k] == "A40.24") L <- 24
    if(dirlist[k] == "A40.20") L <- 20
    filelist <- Sys.glob(paste(dirlist[k], "/res*.Rdata", sep=""))
    for(i in c(1:length(filelist))) {
      ##cat("read from ", filelist[i], "\n")
      Str <- strsplit(filelist[i], ".", fixed=TRUE)[[1]]
      pc <- Str[length(Str)-2]
      level <- as.integer(substr(pc, nchar(pc), nchar(pc)))
      pc <- Str[length(Str)-1]
      frame <- as.integer(substr(pc, nchar(pc), nchar(pc))) + 1
      load(filelist[i])
      if(length(res[,1,1,1]) != boot.R+1) {
        stop(paste("inconsistent value of boot.R for", filelist[i], "\n"))
      }
      maxcomb <- max(maxcomb, dim(res)[2]*dim(res)[3])
      reslist <- append(reslist, list(res))
      
      ii <- which(dim(res) == 1)-1
      jj <- c(2,3)
      if(length(ii) > 0) {
        jj <- jj[-ii]
      }
      if(length(jj) == 0) {
        stop("only one fit range combination?\n")
      }
      
      ## this is the statistical uncertainty
      err <- as.vector(apply(res[,,,1], jj, sd, na.rm=TRUE))
      ## p-values
      pvalues <- as.vector(abs(res[1,,,7]-0.5) * abs(res[1,,,8]-0.5))
      ## q^2/mpi^2
      w <- compute.weights(err, pvalues*max(err))
      qsq <- weighted.median(as.vector(res[1,,,1]), w=w, na.rm=TRUE)
      
      ## this is the naive statistical uncertainty
      err <- as.vector(apply(res[,,,6], jj, sd, na.rm=TRUE))
      ## weights
      w <- compute.weights(err, pvalues*max(err))
      E <- weighted.quantile(as.vector(res[1,,,6]), prob=c(0.5), w=w, na.rm=TRUE)
      bdata[1,1] <- E
      bdata[rr,1] <- apply(res[-1,,,6], 1, weighted.quantile, prob=c(0.5), w=w, na.rm=TRUE)
      
      ## statistical error
      E[2] <- sd(bdata[,1], na.rm=TRUE)
      ## systematic error
      ## lower / upper
      E[c(3:4)] <- weighted.quantile(as.vector(res[1,,,6]), w=w, prob=c(0.1573, 0.8427), na.rm=TRUE)-E[1]
      
      Epipi <- rbind(Epipi, data.frame(E=E[1], dE=E[2], dmE=E[3], dpE=E[4], frame=frame, level=level, counter=counter, L=L, qsq=qsq))
      
      ## this is the statistical uncertainty
      err <- as.vector(apply(res[,,,5], jj, sd, na.rm=TRUE))
      ## q^2/mpi^2
      w <- compute.weights(err, pvalues*max(err))
      E <- weighted.quantile(as.vector(res[1,,,5]), prob=c(0.5), w=w, na.rm=TRUE)
      bdata[1,2] <- E
      bdata[rr,2] <- apply(res[-1,,,5], 1, weighted.quantile, prob=c(0.5), w=w, na.rm=TRUE)
      if(counter == 1) boot.data <- bdata
      else boot.data <- cbind(boot.data, bdata)
      ##qsqovmpisq.statd
      E[2] <- sd(bdata[,2], na.rm=TRUE)
      E[c(3,4)] <- weighted.quantile(as.vector(res[1,,,5]), w=w, prob=c(0.1573, 0.8427), na.rm=TRUE)-E[1]
      
      Epi <- rbind(Epi, data.frame(E=E[1], dE=E[2], dmE=E[3], dpE=E[4], frame=frame, level=level, counter=counter, L=L, qsq=qsq))
      counter <- counter + 1
      rm(w, res)
    }
  }
  ## inverse variance-covariance matrix
  M <- invertCovMatrix(boot.data[,seq(1,(length(boot.data[1,])-1),by=2)], boot.samples=TRUE)
  save(Epipi, Epi, boot.data, M, file="Energies.Rdata")
}
load("Energies.Rdata")
print(Epipi)
print(Epi)

#source("~/daten/workdir/hadron/exec/phaseshift/phaseshift.R")
source("phaseshift.R")

getEvalues <- function(par, Wlist, Lvalues, nolist, framelist, scanfnlist, fnlist) {
  evalues <- c()
  j <- 1
  for(f in framelist) {
    for(L in Lvalues) {
      if(nolist[[j]] > 0) {
        ii <- findSignChanges(fn=scanfnlist[[f]], par=par, makeplot=FALSE, W=Wlist[[j]], no=nolist[[j]])
        ##cat(f, L, ii, "\n")

        if(length(ii) != nolist[[j]]) return(NA)
        zeros <- findZeros(fn=fnlist[[f]], q=Wlist[[j]]$q, ii=ii, L=L, Mpi=Wlist[[j]]$Mpi, par=par, makeplot=FALSE)
        if(length(zeros) != nolist[[j]]) return(NA)
        evalues <- c(evalues, Eofqsq(zeros, dvec=Wlist[[j]]$dvec, Mpi=Wlist[[j]]$Mpi, L=L))
      }
      j <- j+1
    }
  }
  return(evalues)
}

chisqr <- function(par, M, energydata, Wlist, Lvalues, nolist, framelist, scanfnlist, fnlist) {
  xv <- getEvalues(par, Wlist, Lvalues, nolist, framelist, scanfnlist, fnlist) - energydata
  return(xv %*% M %*% xv)
}

chi <- function(par, LM, energydata, Wlist, Lvalues, nolist, framelist, scanfnlist, fnlist) {
  #xv <- LM %*% (getEvalues(par, Wlist, Lvalues, nolist, framelist, scanfnlist, fnlist) - energydata)
  #cat(xv, "\n")
  #return(xv)
  return(LM %*% (getEvalues(par, Wlist, Lvalues, nolist, framelist, scanfnlist, fnlist) - energydata))
}


Lvalues <- unique(Epipi$L)
q <- seq(0.00001, .6, 0.0005)
lm.avail <- require(minpack.lm)
framelist <- c(1:4)
scanfnlist <- list(detEqCMA1scan, detEqMF1A1scan, detEqMF2A1scan, detEqMF3A1scan)
fnlist <- list(detEqCMA1, detEqMF1A1, detEqMF2A1, detEqMF3A1)

par <- c(-1, 5, -100)

boot.par <- array(0., dim=c(boot.R+1,length(par)+1))

for(r in c(1, rr)) {
#for(r in c(1: 10)) {

  energydata <- c()
  iiEdata <- c()
  j <- 1

  ## generate the omega_lm for the moving frames
  ## we store them in a list
  Wlist <- list()
  nolist <- list()
  for(f in framelist) {
    for(L in Lvalues) {
      if(any(Epipi[Epipi$L==L,]$frame == f)) {
        jj <- which(Epi$frame==f & Epi$L==L)
        if(f == 1) Wlist[[j]] <- prepdetEqCMscan(q=q, L=L, Mpi=boot.data[r,2*jj[1]])
        if(f == 2) Wlist[[j]] <- prepdetEqMF1scan(q=q, L=L, Mpi=boot.data[r,2*jj[1]])
        if(f == 3) Wlist[[j]] <- prepdetEqMF2scan(q=q, L=L, Mpi=boot.data[r,2*jj[1]])
        if(f == 4) Wlist[[j]] <- prepdetEqMF3scan(q=q, L=L, Mpi=boot.data[r,2*jj[1]])
        nolist[[j]] <- max(Epipi[Epipi$frame==f & Epipi$L==L,]$level)
        iiEdata <- c(iiEdata, jj)
      }
      else {
        Wlist[[j]] <- NULL
        nolist[[j]] <- 0
      }
      j <- j+1
    }
  }
  energydata <- boot.data[r,(2*iiEdata-1)]
  cat(energydata, "\n")

  M <- invertCovMatrix(boot.data[,(2*iiEdata-1)], boot.samples=TRUE)
                                        #M <- diag(1./Epipi$dE[iiEdata]^2)
  LM <- chol(M)
  

  if(!lm.avail) {
    opt.res <- nls.lm(par, fn=chi, control = nls.lm.control(ftol=1.e-12, ptol=1.e-12, nprint=10), LM=LM, energydata=boot.data[r,(2*iiEdata-1)], Wlist=Wlist, Lvalues=Lvalues, nolist=nolist, framelist=framelist, scanfnlist=scanfnlist, fnlist=fnlist)
    boot.par[r, length(par)+1] <- opt.res$rsstrace[length(opt.res$rsstrace)]
  }
  else {
    opt.res <- optim(par, chisqr, method = c("BFGS"), control = list(parscale=1./par, trace=5), M=M, energydata=boot.data[r,(2*iiEdata-1)], Wlist=Wlist, Lvalues=Lvalues, nolist=nolist, framelist=framelist, scanfnlist=scanfnlist, fnlist=fnlist)
    ##opt.res <- optim(opt.res$par, chisqr, method = c("BFGS"), control = list(parscale=1./opt.res$par), M=M, energydata=energydata, Wlist=Wlist, Lvalues=Lvalues, nolist=nolist, framelist=framelist, scanfnlist=scanfnlist, fnlist=fnlist)
    boot.par[r, length(par)+1] <- opt.res$value
  }
  boot.par[r, c(1:length(par))] <- opt.res$par
  cat(r, "\n")
}
evalues <- getEvalues(boot.par[1,c(1:length(par))], Wlist, Lvalues, nolist, framelist, scanfnlist, fnlist)

cat(evalues, "\n")
cat(energydata, "\n")
cat((evalues-energydata)/Epipi$dE[iiEdata], "\n")

exit()

ii <- which(Qsq$qsq < cut)
rQcot <- Qcot[Qsq$qsq < cut, ]
rQsq <- Qsq[Qsq$qsq < cut, ]



## fit on the original data
res <- optim(par, chisqr, method = c("BFGS"), x=rQsq$qsq, y=rQcot$qcotdelta, M=M)
res <- optim(res$par, chisqr, method = c("BFGS"), x=rQsq$qsq, y=rQcot$qcotdelta, M=M, control = list(parscale=1./res$par))

## compute the statistical error
cat("cut in qsq", cut, "\n")
cat("Computing statistical uncertainty\n")
boot.par <- array(0., dim=c(1500,length(par)+1))
for(r in c(1:boot.R)) {
  boot.res <- optim(res$par, chisqr, method = c("BFGS"), x=boot.data[r,2*ii], y=boot.data[r,2*ii-1], M=M, control = list(parscale=1./res$par))
  boot.par[r,c(1:length(par))] <- boot.res$par
  boot.par[r,length(par)+1] <- boot.res$value
}
cat("a = ", 1./res$par[length(par)-1], "+-", sd(1./boot.par[,length(par)-1], na.rm=TRUE), " (in lattice units)\n")
cat("r = ", res$par[length(par)]*2, "+-", sd(boot.par[,length(par)]*2, na.rm=TRUE), " (in lattice units)\n")

## compute the systematic error
cat("Computing systematic uncertainty\n")
syst.par <- array(0., dim=c(4*maxcomb, length(par)+1))
sreslist <- reslist[ii]

for(r in c(1:(4*maxcomb))) {
  rQsq <- data.frame()
  rQcot <- data.frame()
  for(i in c(1:length(ii))) {
    a <- sample.int(dim(sreslist[[i]])[2], size=1)
    b <- sample.int(dim(sreslist[[i]])[3], size=1)
    qcotdelta <- sreslist[[i]][1,a,b,3]
    
    rQcot <- rbind(rQcot, data.frame(qcotdelta=qcotdelta, i=i, betaval=k, counter=counter))
    
    qsq <- sreslist[[i]][1,a,b,1]
    
    rQsq <- rbind(rQsq, data.frame(qsq=qsq, i=i, betaval=k, counter=counter))
  }
  if(any(is.na(c(rQsq$qsq, rQcot$qcotdelta)))) {
    syst.par[r,] <- NA
  }
  else {
    syst.res <- optim(res$par, chisqr, method = c("BFGS"), x=rQsq$qsq, y=rQcot$qcotdelta, M=M, control = list(parscale=1./res$par))
    syst.par[r,c(1:length(par))] <- syst.res$par
    syst.par[r,length(par)+1] <- syst.res$value
  }
}

cat("a = ", 1./res$par[length(par)-1],
    " - ", -quantile(1./syst.par[,length(par)-1], prob=c(0.1573), na.rm=TRUE)+1./res$par[length(par)-1],
    " + ", quantile(1./syst.par[,length(par)-1], prob=c(0.8427), na.rm=TRUE)-1./res$par[length(par)-1], "\n")

cat("r = ", res$par[length(par)]*2,
    " - ", 2*(-quantile(syst.par[,length(par)], prob=c(0.1573), na.rm=TRUE)+res$par[length(par)]),
    " + ", 2*quantile(syst.par[,length(par)], prob=c(0.8427), na.rm=TRUE)-res$par[length(par)], "\n")

## generate statistical error band
x <- seq(0,cut,0.0001)

fn <- function(x, par) {
  return(par[1] + par[2]*x)
}

fnx <- fn(x=x, par=res$par[c(length(par)-1, length(par))])
bfnx <- array(apply(boot.par[,5:6], 1, fn, x=x), dim=c(length(x),dim(boot.par)[1]))
sdfnx <- apply(bfnx, 1, sd, na.rm=TRUE)

## generate systematic error band
bfnx <- array(apply(syst.par[,5:6], 1, fn, x=x), dim=c(length(x),dim(syst.par)[1]))
sysfnxm <- apply(bfnx, 1, quantile, prob=c(0.1573),  na.rm=TRUE)
sysfnxp <- apply(bfnx, 1, quantile, prob=c(0.8427),  na.rm=TRUE)

## write to file
write.table(data.frame(x, fnx, sdfnx, sysfnxm, sysfnxp), file="errorband.dat", row.names=FALSE, col.names=FALSE, quote=FALSE)
