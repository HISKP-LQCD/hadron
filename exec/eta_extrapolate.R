read.etadata <- function(path) {

  meta <- read.table(paste(path,"eta_trick/L_2x2.blocking.masses.cosh.state.0.all", sep=""))
  metap <- read.table(paste(path,"eta_trick/L_2x2.blocking.masses.cosh.state.1.all", sep=""))
  
  mpi <- read.table(paste(path,"pion+/LF_4x4.blocking.masses.cosh.state.0.all", sep=""))
  if(file.exists(paste(path,"kaon/LF_8x8.masses.cosh.state.0.all", sep=""))) {
    mK <- read.table(paste(path,"kaon/LF_8x8.masses.cosh.state.0.all", sep=""))
  }
  else if(file.exists(paste(path,"kaon/resampling.dat", sep=""))) {
    mK <- read.table(paste(path,"kaon/resampling.dat", sep=""))
  }
  r0 <- read.table(paste(path,"r0/r0_chiral.resampled", sep=""))
  xb <- array(0, dim=c(1000, 5))
  xb[,1] <- r0$V2
  xb[,2] <- meta$V2
  xb[,3] <- metap$V2
  xb[,4] <- mpi$V2
  xb[,5] <- mK$V2
  return(invisible(xb))
}

R0 = 0.45
##R0 = 0.5
dR0 = 0.02
MPI0 = 0.1349766
MK = 0.497648
fm2GeV = 0.197327
xfit <- seq(0, 1.4, 0.01)


ensembles <- dir("/dsk/lattice02-0/ottnad/data/",
                 pattern = "^[ABD]",full.names=FALSE, ignore.case = FALSE)
Nens <- length(ensembles)
r0 <- c(rep(5.231, times=10),
        rep(5.710, times=6),
        rep(7.538, times=4))
dr0 <- c(rep(0.038, times=10),
        rep(0.041, times=6),
        rep(0.058, times=4))


if(eta.reread) {
  bres <- array(0, dim=c(Nens-5, 1000, 5))
  dres <- array(0, dim=c(Nens-5, 5))
  ens <- rep("", times=Nens-5)
  j <- 1

  for(i in 1:Nens) {
    if(!(i %in% c(4,7,13,18,21))) {
      path <- paste("/dsk/lattice02-0/ottnad/data/", ensembles[i], "/analyse/", sep="")
      cat(ensembles[i], i, "\n")
      bres[j,,c(1:5)] <- read.etadata(path)
      dres[j,] <- apply(bres[j,,], 2, sd)
      ens[j] <- ensembles[i]
      j <- j+1
    }
  }
}

reR0 <- c(R0, rnorm(boot.R-1, R0, dR0))
pp <- (reR0*MPI0/fm2GeV)^2
pK <- (reR0*MK/fm2GeV)^2

indexA <- c(1:8)
indexB <- c(9:13)
indexD <- c(14:16)

for(i in indexA) bres[,i,1] <- bres[,indexA[1],1]
for(i in indexB) bres[,i,1] <- bres[,indexB[1],1]
for(i in indexD) bres[,i,1] <- bres[,indexD[1],1]

result.frame <- data.frame(ens = c(ens, "phys"),
                           r0mpssq = rep(0, times=Nens-5+1), dr0mpssq = rep(0, times=Nens-5+1),
                           meovmk=rep(0, times=Nens-5+1), dmeovmk=rep(0, times=Nens-5+1),
                           gmo=rep(0, times=Nens-5+1), dgmo=rep(0, times=Nens-5+1),
                           r0metabar=rep(0, times=Nens-5+1), dr0metabar=rep(0, times=Nens-5+1),
                           dr0metabar2=rep(0, times=Nens-5+1),
                           r0meta=rep(0, times=Nens-5+1), dr0meta=rep(0, times=Nens-5+1),
                           r0metap=rep(0, times=Nens-5+1), dr0metap=rep(0, times=Nens-5+1)
                           )

boot.R <- length(bres[1,,1])
## extrapolate (meta/mK)^2 [(r0 mpi)^2]
eta1res <- array(0, dim=c(boot.R,2))
dx <- apply((bres[,,2]/bres[,,5])^2, 1, sd)
result.frame$r0mpssq <- c((bres[,1,1]*bres[,1,4])^2, pp[1])
result.frame$dr0mpssq <- c(apply((bres[,,1]*bres[,,4])^2, 1, sd), 0)
result.frame$meovmk <- c((bres[,1,2]/bres[,1,5]), 0)
result.frame$dmeovmk <- c(apply((bres[,,2]/bres[,,5]), 1, sd), 0)
  
for(i in 1:boot.R) {
  x <- (bres[,i,2]/bres[,i,5])^2
  t <- (bres[,i,1]*bres[,i,4])^2
  res <- lm(x ~ t, weights=1./dx^2)
  eta1res[i,1] <- res$coef[1]
  eta1res[i,2] <- res$coef[2]
}
result.frame$meovmk[Nens-4] <- sqrt(eta1res[1,1] + pp[1]*eta1res[1,2])
result.frame$dmeovmk[Nens-4] <- sd(sqrt(eta1res[2:boot.R,1] + pp[2:boot.R]*eta1res[2:boot.R,2]))
yfit <- sqrt(eta1res[1,1] + eta1res[1,2]*xfit)
dyfit <- array(0, dim=c(boot.R, length(xfit)))
for(i in 1:boot.R) {
  dyfit[i,] <- sqrt(eta1res[i,1] + eta1res[i,2]*xfit)
}
write.table(data.frame(x=xfit, y=yfit, dy=apply(dyfit, 2, sd)),
            file="meovmkfit.dat", quote=FALSE, col.name=FALSE, row.names=FALSE)


## extrapolate GMO as a function of [(r0 mpi)^2]
eta2res <- array(0, dim=c(boot.R,2))
dx <- apply(3*bres[,,2]^2/(4*bres[,,5]^2-bres[,,4]^2), 1, sd)
result.frame$gmo <- c(3*bres[,1,2]^2/(4*bres[,1,5]^2-bres[,1,4]^2), 0)
result.frame$dgmo <- c(dx, 0)
  
for(i in 1:boot.R) {
  x <- 3*bres[,i,2]^2/(4*bres[,i,5]^2-bres[,i,4]^2)
  t <- (bres[,i,1]*bres[,i,4])^2
  res <- lm(x ~ t, weights=1./dx^2)
  eta2res[i,1] <- res$coef[1]
  eta2res[i,2] <- res$coef[2]
}
result.frame$gmo[Nens-4] <- eta2res[1,1] + pp[1]*eta2res[1,2]
result.frame$dgmo[Nens-4] <- sd(eta2res[2:boot.R,1] + pp[2:boot.R]*eta2res[2:boot.R,2])
yfit <- eta2res[1,1] + eta2res[1,2]*xfit
for(i in 1:boot.R) {
  dyfit[i,] <- eta2res[i,1] + eta2res[i,2]*xfit
}
write.table(data.frame(x=xfit, y=yfit, dy=apply(dyfit, 2, sd)),
            file="gmofit.dat", quote=FALSE, col.name=FALSE, row.names=FALSE)



## extrapolate metap^2 [(r0 mpi)^2]
etapres <- array(0, dim=c(boot.R,2))
dx <- apply((bres[,,1]*bres[,,3])^2, 1, sd)
result.frame$r0metap <- c((bres[,1,1]*bres[,1,3]), 0)
## here we don't include the r0 error in the error
## for presentation reasons...
result.frame$dr0metap <- c(apply((bres[,,1]*bres[,,3]), 1, sd), 0)

result.frame$r0meta <- c((bres[,1,1]*bres[,1,2]), 0)
## here we don't include the r0 error in the error
## for presentation reasons...
result.frame$dr0meta <- c(apply((bres[,,1]*bres[,,2]), 1, sd), 0)

for(i in 1:boot.R) {
  x <- (bres[,i,1]*bres[,i,3])^2
  t <- (bres[,i,1]*bres[,i,4])^2
  res <- lm(x ~ t, weights=1./dx^2)
  etapres[i,1] <- res$coef[1]
  etapres[i,2] <- res$coef[2]
}
result.frame$r0metap[Nens-4] <- sqrt(etapres[1,1] + pp[1]*etapres[1,2])
result.frame$dr0metap[Nens-4] <- sd(sqrt(etapres[2:boot.R,1] + pp[2:boot.R]*etapres[2:boot.R,2]))
yfit <- sqrt(etapres[1,1] + etapres[1,2]*xfit)
for(i in 1:boot.R) {
  dyfit[i,] <- sqrt(etapres[i,1] + etapres[i,2]*xfit)
}
write.table(data.frame(x=xfit, y=yfit, dy=apply(dyfit, 2, sd)),
            file="etapfit.dat", quote=FALSE, col.name=FALSE, row.names=FALSE)

## compute dMeta^2/dMK^2
Detaj <- array(0, dim=c(1000,2))
for(i in 1:boot.R) {
  ## A100, A100s
  Detaj[i,1] <- (bres[1,i,2]^2 -bres[2,i,2]^2)/(bres[1,i,5]^2 -bres[2,i,5]^2)
  ## A80, A80s
  Detaj[i,2] <- (bres[7,i,2]^2 -bres[8,i,2]^2)/(bres[7,i,5]^2 -bres[8,i,5]^2)
}
Deta <- apply(Detaj[,], 1, mean)

## fit (aMK)^2 = aMKchi^2 + s*aMPs^2
## B-ensembles
ii <- c(9, 10, 11, 12, 13)
## A-ensembles
#ii <- c(1,3,5,6,7)
err <- apply(bres[ii,,5]^2, 1, sd)
kslope <- rep(0, times=boot.R)
for(i in 1:boot.R) {
  x <- bres[ii,i,5]^2
  t <- bres[ii,i,4]^2
  res <- lm (x~t, weights=1./err^2)
  kslope[i] <- res$coef[2]
}

## adjust intercept (r0*MKchi)^2to go through physical pion mass
kchi <- (pK - kslope*pp)
## compute differences in squared kaon masses
deltaKsq <- array(0, dim=c(Nens-5, 1000))
for(i in 1:(Nens-5)) {
  deltaKsq[i,] <- (bres[i,,1]^2*(bres[i,,5]^2 - kslope*bres[i,,4]^2) - kchi)
}
Metabar <- array(0, dim=c(Nens-5, 1000))
Metabarplot <- array(0, dim=c(Nens-5, 1000))
for(i in 1:(Nens-5)) {
  Metabar[i,] <- bres[i,,1]^2*bres[i,,2]^2 - Deta*deltaKsq[i,]
  Metabarplot[i,] <- bres[i,,1]^2*bres[i,,2]^2 - Deta[1]*deltaKsq[i,1]
}
dMetabar <- apply(Metabar, 1, sd)

## extrapolate metap^2 [(r0 mpi)^2]
eta3res <- array(0, dim=c(boot.R,2))
dx <- dMetabar
result.frame$r0metabar <- c(sqrt(Metabar[,1]), 0)
## here we don't include the r0 error in the error
result.frame$dr0metabar <- c(apply(sqrt(Metabar)/bres[,,1], 1, sd)*bres[,1,1], 0)
result.frame$dr0metabar2 <- c(apply(sqrt(Metabarplot)/bres[,,1], 1, sd)*bres[,1,1], 0)
  
for(i in 1:boot.R) {
  x <- Metabar[,i]
  t <- (bres[,i,1]*bres[,i,4])^2
  res <- lm(x ~ t, weights=1./dx^2)
  eta3res[i,1] <- res$coef[1]
  eta3res[i,2] <- res$coef[2]
}
result.frame$r0metabar[Nens-4] <- sqrt(eta3res[1,1] + pp[1]*eta3res[1,2])
result.frame$dr0metabar[Nens-4] <- sd(sqrt(eta3res[2:boot.R,1] + pp[2:boot.R]*eta3res[2:boot.R,2]))
result.frame$dr0metabar2[Nens-4] <- sd(sqrt(eta3res[2:boot.R,1] + pp[2:boot.R]*eta3res[2:boot.R,2]))
yfit <- sqrt(eta3res[1,1] + eta3res[1,2]*xfit)
for(i in 1:boot.R) {
  dyfit[i,] <- sqrt(eta3res[i,1] + eta3res[i,2]*xfit)
}
write.table(data.frame(x=xfit, y=yfit, dy=apply(dyfit, 2, sd)),
            file="etafit.dat", quote=FALSE, col.name=FALSE, row.names=FALSE)

result <- data.frame(metaovmK = sqrt(eta1res[1,1] + pp[1]*eta1res[1,2]),
                     dmetaovmK = sd(sqrt(eta1res[2:boot.R,1] + pp[2:boot.R]*eta1res[2:boot.R,2])),
                     meta1 = sqrt(eta1res[1,1] + pp[1]*eta1res[1,2])*MK,
                     dmeta1 =  sd(sqrt(eta1res[2:boot.R,1] + pp[2:boot.R]*eta1res[2:boot.R,2])*MK),
                     GMO = eta2res[1,1] + pp[1]*eta2res[1,2],
                     dGMO = sd(eta2res[2:boot.R,1] + pp[2:boot.R]*eta2res[2:boot.R,2]),
                     metaGMO = sqrt((eta2res[1,1] + pp[1]*eta2res[1,2])/3*(4*MK^2-MPI0^2)),
                     dmetaGMO = sd(sqrt((eta2res[2:boot.R,1] + pp[2:boot.R]*eta2res[2:boot.R,2])/3*(4*MK^2-MPI0^2))),
                     meta3 = sqrt(eta3res[1,1] + pp[1]*eta3res[1,2])/R0*fm2GeV,
                     dmeta3 = sd(sqrt(eta3res[2:boot.R,1] + pp[2:boot.R]*eta3res[2:boot.R,2])/reR0[2:boot.R]*fm2GeV),
                     metap = sqrt(etapres[1,1] + pp[1]*etapres[1,2])/R0*fm2GeV,
                     dmetap = sd(sqrt(etapres[2:boot.R,1] + pp[2:boot.R]*etapres[2:boot.R,2])/reR0[2:boot.R]*fm2GeV)
                     )

write.table(result.frame, file="result.dat", quote=FALSE, col.names=FALSE, row.names=FALSE)

ploteta <- function(result.frame) {
  plotwitherror(result.frame$r0mpssq[indexA], result.frame$r0metabar[indexA], result.frame$dr0metabar2[indexA], xlim=c(0,1.4), ylim=c(0,3.0), col="red")
  plotwitherror(result.frame$r0mpssq[indexB], result.frame$r0metabar[indexB], result.frame$dr0metabar2[indexB], col="blue", rep=TRUE)
  plotwitherror(result.frame$r0mpssq[indexD], result.frame$r0metabar[indexD], result.frame$dr0metabar2[indexD], col="green", rep=TRUE)
  plotwitherror(result.frame$r0mpssq[indexA], result.frame$r0metap[indexA], result.frame$dr0metap[indexA], col="red", rep=TRUE)
  plotwitherror(result.frame$r0mpssq[indexB], result.frame$r0metap[indexB], result.frame$dr0metap[indexB], col="blue", rep=TRUE)
  plotwitherror(result.frame$r0mpssq[indexD], result.frame$r0metap[indexD], result.frame$dr0metap[indexD], col="green", rep=TRUE)
  plotwitherror(result.frame$r0mpssq[17], result.frame$r0metabar[17], result.frame$dr0metabar2[17], col="black", rep=TRUE)
  plotwitherror(result.frame$r0mpssq[17], result.frame$r0metap[17], result.frame$dr0metap[17], col="black", rep=TRUE)
}
