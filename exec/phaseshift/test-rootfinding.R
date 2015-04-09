source("phaseshift.R")

Mpi <- 0.142
L <- 24
no <- 3

par <- c(-1,20,100)

q <- seq(0.00001, .6, 0.0005)

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

ii <- findSignChanges(fn=detEqCMA1scan, par=par, makeplot=TRUE, W=WCM, no=no, threshold=10)
zerosCM <- findZeros(fn=detEqCMA1, q=q, ii=ii, L=L, Mpi=Mpi, par=par, makeplot=TRUE)
ECM <- Eofqsq(zerosCM, dvec=c(0,0,0), Mpi=Mpi,L=L)

ii <- findSignChanges(fn=detEqMF1A1scan, par=par, makeplot=TRUE, W=WMF1, no=no, threshold=1000)
zerosMF1 <- findZeros(fn=detEqMF1A1, q=q, ii=ii, L=L, Mpi=Mpi, par=par, makeplot=TRUE)
EMF1 <- Eofqsq(zerosMF1, dvec=c(0,0,1), Mpi=Mpi,L=L)

ii <- findSignChanges(fn=detEqMF2A1scan, par=par, makeplot=TRUE, W=WMF2, no=no, threshold=10000)
zerosMF2 <- findZeros(fn=detEqMF2A1, q=q, ii=ii, L=L, Mpi=Mpi, par=par, makeplot=TRUE)
EMF2 <- Eofqsq(zerosMF2, dvec=c(1,1,0), Mpi=Mpi, L=L)

ii <- findSignChanges(fn=detEqMF3A1scan, par=par, makeplot=TRUE, W=WMF3, no=no, threshold=500)
zerosMF3 <- findZeros(fn=detEqMF3A1, q=q, ii=ii, L=L, Mpi=Mpi, par=par, makeplot=TRUE)
EMF3 <- Eofqsq(zerosMF3, dvec=c(1,1,1), Mpi=Mpi, L=L)

