

boot.R <- 1500
rr <- c(2:(boot.R+1))

Lvalues <- c(32, 24, 20)
## too coarse in the first part! q <- sqrt(c(seq(0.000001, .1, 0.00015), seq(0.100001, .15, 0.00025)))                                                                                                                                        
q <- sqrt(c(seq(0.000001, .1, 0.00009), seq(0.100001, .15, 0.00025)))
lm.avail <- require(minpack.lm)
framelist <- c(1:4)
irreplist <- c(1:3)

doSamples <- 10

## cut in chisq                                                                                                                                                                                                                               
cut <- 3.
lat.disp <- TRUE

if(file.exists("parameters.R")) {
  source("parameters.R")
}

load(file=paste("fit-detEq-result-cut", cut, ".Rdata", sep=""))

load("Energies.Rdata")

qsq <- numeric()
## determine the rough q^2 values for applying a cut                                                                                                                                                                                          
for(i in c(1:length(Epipi$E))) {
  dvec=c(0,0,0)
  if(Epipi$frame[i] == 2) dvec=c(0,0,1)
  if(Epipi$frame[i] == 3) dvec=c(0,1,1)
  if(Epipi$frame[i] == 4) dvec=c(1,1,1)
  if(lat.disp) qtsq <- compute.qtildesq(Epipi$E[i], dvec=dvec, L=Epipi$L[i], mpi=Epi$E[i])
  else qtsq <- compute.qtildesq.contdisp(Epipi$E[i], dvec=dvec, L=Epipi$L[i], mpi=Epi$E[i])
  qsq[i] <- qtsq$q^2
}
Epipi$qsq <- qsq
Epi$qsq <- qsq

cat("applying cut at q^2/Mpi^2 =", cut, "\n")
k <- which(Epi$qsq/Epi$E^2 <= cut)
Epipi <- Epipi[k,]
Epi <- Epi[k,]
kk <- rep(0, times=2*length(k))
l <- c(1:length(k))
kk[2*l-1] <- 2*k-1
kk[2*l] <- 2*k
boot.data <- boot.data[,kk]
print(cbind(Epipi, Epi$E))



## remove some obvious outliers
ind1 <- which(Epipi$irrep[iiEdata] == 1)
indn1 <- which(Epipi$irrep[iiEdata] != 1)

##kj <- which(bootqcotdelta[c(1:doSamples),ind1[17]] > -0.3)
##bootq[kj,ind1[17]] <- NA
##bootqcotdelta[kj,ind1[17]] <- NA

##jk <- which(bootqcotdelta[c(1:doSamples),ind1[16]] < -0.5)
##bootq[jk, ind1[16]] <- NA
##bootqcotdelta[jk, ind1[16]] <- NA


## now prepare the plot
x <- seq(0.00001,0.08,0.0001)
y <- atan(sqrt(x)/(boot.par[1,1] + boot.par[1,2]*x/2))
y2 <- atan(x^2*sqrt(x)/(boot.par[1,3]))

y <- array(0, dim=c(doSamples, length(x)))
y2 <- array(0, dim=c(doSamples, length(x)))
for(i in c(1:doSamples)) {
  y[i,] <- atan(sqrt(x)/(boot.par[i,1] + boot.par[i,2]*x/2))
  y2[i,] <- atan(x^2*sqrt(x)/(boot.par[i,3]))
}
## error band                                                                                                                                                                                                                                 
dy <- apply(y, 2, sd)
dy2 <- apply(y2, 2, sd)


tikzfiles <- tikz.init(paste("delta-vs-qsq.cut", cut, sep=""),width=6,height=5)

plot(NA, ylim=c(-0.8, 0.05), xlim=c(0, 0.08), xlab=c("$a^2q^2$"), ylab=c("$\\delta_\\ell$"))
polygon(x=c(x, rev(x)), y=c(y[1,]+dy, rev(y[1,]-dy)), col="gray", lty=0, lwd=0.001, border="gray")
lines(x,y[1,])
polygon(x=c(x, rev(x)), y=c(y2[1,]+dy2, rev(y2[1,]-dy2)), col="gray", lty=0, lwd=0.001, border="gray")
lines(x,y2[1,], lty=2)
legend(x=0.065, y=-0.4, legend=c("$\\ell = 0$ fit", "$\\ell = 2$ fit"), bty="n", lty=c(1,2))

#legend(x=0.045, y=-0.45, legend=c("$\\ell = 0$"), bty="n")                                                                                                                                                                                   
#legend(x=0.05, y=0.05, legend=c("$\\ell = 2$"), bty="n")

pchlist <- c(21,22,23,24)
collist <- c("red", "blue", "darkgreen")

for(l in c(1:length(Lvalues))) {
  for(f in c(1:length(framelist))) {
    ind1 <- which(Epipi$irrep[iiEdata] == 1 & Epipi$L[iiEdata] == Lvalues[l] & Epipi$frame[iiEdata] == f)
    indn1 <- which(Epipi$irrep[iiEdata] != 1 & Epipi$L[iiEdata] == Lvalues[l] & Epipi$frame[iiEdata] == f)

    for(p in c(1,2)) {
      if(p==1) {
        ii <- ind1
        bg <- collist
        pow <- 1
      }
      else {
        ii <- indn1
        bg <- rep("white", times=length(collist))
        pow <- 5
      }
      if(length(ii) > 1) {
        dx <- apply(bootq[c(1:doSamples),ii]^2, 2, sd, na.rm=TRUE)
        dy <- apply(atan(bootq[c(1:doSamples),ii]^pow/bootqcotdelta[c(1:doSamples),ii]), 2, sd, na.rm=TRUE)
      }
      else {
        dx <- sd(bootq[c(1:doSamples),ii]^2, na.rm=TRUE)
        dy <- sd(atan(bootq[c(1:doSamples),ii]^pow/bootqcotdelta[c(1:doSamples),ii]), na.rm=TRUE)
      }
      if(!any(is.na(bootq[1,ii])) && !any(is.na(bootqcotdelta[1,ii]))) {
        plotwitherror(bootq[1,ii]^2, atan(bootq[1,ii]^pow/bootqcotdelta[1,ii]), dx=dx,
                      dy=dy, rep=TRUE, col=collist[l], bg=bg[l], pch=pchlist[f])
      }
      else {
        plotwitherror(mean(bootq[c(1:doSamples),ii], na.rm=TRUE)^2, atan(mean(bootq[c(1:doSamples),ii], na.rm=TRUE)^pow/mean(bootqcotdelta[c(1:doSamples),ii], na.rm=TRUE)), dx=dx,
                      dy=dy, rep=TRUE, col=collist[l], bg=bg[l], pch=pchlist[f])
      }
    }
  }
}
legend("topright", legend=c("filled symbols: $\\delta_0$","open symbols: $\\delta_2$"), bty="n")
legend(x=0.015, y=-0.65, legend=c("$L=32$", "$L=24$", "$L=20$"), bty="n", text.col=collist)
legend("bottomleft", legend=c("$\\vec d = [0,0,0]$","$\\vec d = [0,0,1]$","$\\vec d = [1,1,0]$","$\\vec d = [1,1,1]$"), bty="n", pch=pchlist)

tikz.finalize(tikzfiles=tikzfiles,clean=TRUE)

