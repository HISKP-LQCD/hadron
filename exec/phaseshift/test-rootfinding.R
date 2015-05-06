source("phaseshift.R")

Mpi <- 0.1414706
L <- 32
no <- 3

par <- c(1./-1.2921363062, 14.2906242026, 1./-49.8148003974)

q <- sqrt(c(seq(0.000001, .1, 0.00008), seq(0.100001, .25, 0.0002)))

if(!file.exists("WCM.Rdata")) {
  WCM <- prepdetEqCMscan(q=q, L=L, Mpi=Mpi)
  save(WCM, file="WCM.Rdata")
}
load("WCM.Rdata")
if(!file.exists("WMF1.Rdata")) {
  WMF1 <- prepdetEqMF1scan(q=q, L=L, Mpi=Mpi)
  save(WMF1, file="WMF1.Rdata")
}
load("WMF1.Rdata")
if(!file.exists("WMF2.Rdata")) {
  WMF2 <- prepdetEqMF2scan(q=q, L=L, Mpi=Mpi)
  save(WMF2, file="WMF2.Rdata")
}
load("WMF2.Rdata")
if(!file.exists("WMF3.Rdata")) {
  WMF3 <- prepdetEqMF3scan(q=q, L=L, Mpi=Mpi)
  save(WMF3, file="WMF3.Rdata")
}
load("WMF3.Rdata")
rm(q)

cat(par, "\n")
cat("L = ", L, "\n")
cat("Mpi = ", Mpi, "\n")
ii <- findSignChanges(fn=detEqCMA1scan, par=par, makeplot=TRUE, W=WCM, no=no, threshold=10)
zerosCM <- findZeros(fn=detEqCMA1, q=WCM$q, ii=ii, L=L, Mpi=Mpi, par=par, makeplot=TRUE)
ECM <- Eofqsq(zerosCM, dvec=c(0,0,0), Mpi=Mpi,L=L)
cat("CM A1", zerosCM, "\n")

ii <- findSignChanges(fn=detEqCMEscan, par=par, makeplot=TRUE, W=WCM, no=no, threshold=1)
zerosCME <- findZeros(fn=detEqCME, q=WCM$q, ii=ii, L=L, Mpi=Mpi, par=par, makeplot=TRUE)
ECM <- Eofqsq(zerosCME, dvec=c(0,0,0), Mpi=Mpi,L=L)
cat("CM E", zerosCME, "\n")
cat("energies", ECM, "\n")

ii <- findSignChanges(fn=detEqCMT2scan, par=par, makeplot=TRUE, W=WCM, no=no, threshold=1)
zerosCMT2 <- findZeros(fn=detEqCMT2, q=WCM$q, ii=ii, L=L, Mpi=Mpi, par=par, makeplot=TRUE)
ECM <- Eofqsq(zerosCME, dvec=c(0,0,0), Mpi=Mpi,L=L)
cat("CM T2", zerosCMT2, "\n")


ii <- findSignChanges(fn=detEqMF1A1scan, par=par, makeplot=TRUE, W=WMF1, no=no, threshold=1)
zerosMF1 <- findZeros(fn=detEqMF1A1, q=WMF1$q, ii=ii, L=L, Mpi=Mpi, par=par, makeplot=TRUE)
EMF1 <- Eofqsq(zerosMF1, dvec=c(0,0,1), Mpi=Mpi,L=L)
cat("MF1", zerosMF1, "\n")

ii <- findSignChanges(fn=detEqMF2A1scan, par=par, makeplot=TRUE, W=WMF2, no=no, threshold=.01)
zerosMF2 <- findZeros(fn=detEqMF2A1, q=WMF2$q, ii=ii, L=L, Mpi=Mpi, par=par, makeplot=TRUE)
EMF2 <- Eofqsq(zerosMF2, dvec=c(1,1,0), Mpi=Mpi, L=L)
cat("MF2", zerosMF2, "\n")

ii <- findSignChanges(fn=detEqMF3A1scan, par=par, makeplot=TRUE, W=WMF3, no=no, threshold=1)
zerosMF3 <- findZeros(fn=detEqMF3A1, q=WMF3$q, ii=ii, L=L, Mpi=Mpi, par=par, makeplot=TRUE)
EMF3 <- Eofqsq(zerosMF3, dvec=c(1,1,1), Mpi=Mpi, L=L)
cat("MF3", zerosMF3, "\n")
