Lvalues <- c(32,24,20)
ulim <- c(0.8, 1.0, 1.2)
load("Energies.Rdata")

for(L in Lvalues) {
  i <- which(Lvalues == L)
  tikzfiles <- tikz.init(paste("energylevels", L, sep=""),width=5,height=5)
  
  ii <- which(Epipi$L == L)
  mirrep <- max(Epipi[ii,]$irrep)
  mlevel <- max(Epipi[ii,]$level)
  mframe <- max(Epipi[ii,]$frame)
  N <- length(Epipi[Epipi$L==L & Epipi$level==1,]$E)
  
  xv <- (Epipi[ii,]$irrep-1)*mframe + Epipi[ii,]$frame
  xv[xv==9] <- 6
  E <- Epipi[Epipi$L==L,]$E
  dE <- Epipi[Epipi$L==L,]$dE
  Mpi <- Epi[Epipi$L==L,]$E[1]
  p1 <- array(c(0,0,0, 0,0,1, 0,1,1, 1,1,1, 0,0,1, 0,1,1), dim=c(3,6))
  p2 <- array(c(0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,1, 0,1,1), dim=c(3,6)) 
  
  p12 <- array(c(0,0,1, 1,0,0,  0,1,0, 1,0,0, 0,1,1), dim=c(3,5))
  p22 <- array(c(0,0,1, -1,0,1, 0,0,1, 0,1,1, 0,1,1), dim=c(3,5))
  
  plot(NA, type="n", ylab="$aE$", xlab="", ylim=c(0,ulim[i]), xlim=c(0.5,N+0.5), xaxt="n")
  rect(xleft=xv-0.2, xright=xv, ytop=E+dE, ybottom=E-dE, col="gray")
  legend("topright", legend=c("$E_{\\pi\\pi}$"), bty="n", fill=c("gray"))
  
  labels <- c("A1, $\\vec d=[0,0,0]$", "A1, $\\vec d=[0,0,1]$", "A1, $\\vec d=[0,1,1]$", "A1, $\\vec d=[1,1,1]$", "E, $\\vec d=[0,0,0]$", "T2, $\\vec d=[0,0,0]$")
  axis(1, at=c(1:N), labels=rep("", times=N))
  text(x=c(1:N), y=-0.09, labels=labels, srt = 90, adj = 1, xpd = TRUE,cex=0.8,srt=60)
  legend("topleft", legend=c(paste("A40 L =", L)), bty="n", cex=1.2)
  
  disp <- c()
  disp2 <- c()
  for(i in c(1:6)) {
    disp[i] <- acosh( cosh(Mpi) + 2*sum(sin(pi*p1[,i]/L)^2)) + acosh( cosh(Mpi) + 2*sum(sin(pi*p2[,i]/L)^2))
  }
  for(i in c(1:5)) {
    disp2[i] <- acosh( cosh(Mpi) + 2*sum(sin(pi*p12[,i]/L)^2)) + acosh( cosh(Mpi) + 2*sum(sin(pi*p22[,i]/L)^2))
  }
  arrows(y0=disp, y1=disp, x0=c(1:6)-0.2, x1=c(1:6)+0.2, length=0, col="blue")
  arrows(y0=disp2, y1=disp2, x0=c(1:5)-0.2, x1=c(1:5)+0.2, length=0, col="blue")
  legend("bottomright", legend=c("$E_{\\pi}(\\vec p_1)+E_{\\pi}(\\vec p_2)$"), bty="n", lty=c(1), col=c("blue"))
  

  tikz.finalize(tikzfiles=tikzfiles,clean=TRUE)
}

cut <- 3.
load(paste("fit-detEq-result-cut", cut, ".Rdata", sep=""))
qsq <- numeric()
## determine the rough q^2 values for applying a cut
for(i in c(1:length(Epipi$E))) {
  dvec=c(0,0,0)
  if(Epipi$frame[i] == 2) dvec=c(0,0,1)
  if(Epipi$frame[i] == 3) dvec=c(0,1,1)
  if(Epipi$frame[i] == 4) dvec=c(1,1,1)
  qtsq <- compute.qtildesq(Epipi$E[i], dvec=dvec, L=Epipi$L[i], mpi=Epi$E[i])
  qsq[i] <- qtsq$q^2
}
Epipi$qsq <- qsq
Epi$qsq <- qsq

k <- which(Epi$qsq/Epi$E^2 <= cut)
Epipicut <- Epipi[k,]
Epicut <- Epi[k,]

Epipicut <- Epipicut[iiEdata,]
Epicut <- Epicut[iiEdata,]

for(L in Lvalues) {
  i <- which(Lvalues == L)
  tikzfiles <- tikz.init(paste("energylevels-pfit", L, "-cut", cut, sep=""),width=5,height=6)
  op <- options(cex=1.3)

  ii <- which(Epipi$L == L)
  mirrep <- max(Epipi[ii,]$irrep)
  mlevel <- max(Epipi[ii,]$level)
  mframe <- max(Epipi[ii,]$frame)
  N <- length(Epipi[Epipi$L==L & Epipi$level==1,]$E)
  
  xv <- (Epipi[ii,]$irrep-1)*mframe + Epipi[ii,]$frame
  xv[xv==9] <- 6
  E <- Epipi[Epipi$L==L,]$E
  dE <- Epipi[Epipi$L==L,]$dE
  Mpi <- Epi[Epipi$L==L,]$E[1]

  a <- 0.2
  plot(NA, type="n", ylab="$aE$", xlab="", ylim=c(a,ulim[i]), xlim=c(0.5,N+0.5), xaxt="n")
  rect(xleft=xv-0.2, xright=xv, ytop=E+dE, ybottom=E-dE, col="gray")

  legend("topright", legend=c("data $E_{\\pi\\pi}$", "fitted $E_{\\pi\\pi}$"), bty="n", fill=c("gray", "red"), border=c("black", "red"))
  legend("bottomleft", legend=c(paste("$q^2/M_\\pi^2\\leq$", cut, sep="")), bty="n", cex=1.3)
  
  labels <- c("A1, $\\vec d=[0,0,0]$", "A1, $\\vec d=[0,0,1]$", "A1, $\\vec d=[0,1,1]$", "A1, $\\vec d=[1,1,1]$", "E, $\\vec d=[0,0,0]$", "T2, $\\vec d=[0,0,0]$")
  axis(1, at=c(1:N), labels=rep("", times=N))

  if(L==32) a <- a-0.03/(ulim[i]-a)
  if(L==24) a <- a-0.04/(ulim[i]-a)
  if(L==20) a <- a-0.07/(ulim[i]-a)
  
  text(x=c(1:N), y=a, labels=labels, srt = 90, adj = 1, xpd = TRUE,cex=0.8,srt=60)
  legend("topleft", legend=c(paste("A40 ensemble, L =", L)), bty="n", cex=1.2)

  p1 <- array(c(0,0,0, 0,0,1, 0,1,1, 1,1,1, 0,0,1, 0,1,1), dim=c(3,6))
  p2 <- array(c(0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,1, 0,1,1), dim=c(3,6)) 
  
  p12 <- array(c(0,0,1, 1,0,0,  0,1,0, 1,0,0, 0,1,1), dim=c(3,5))
  p22 <- array(c(0,0,1, -1,0,1, 0,0,1, 0,1,1, 0,1,1), dim=c(3,5))

  disp <- c()
  disp2 <- c()
  for(i in c(1:6)) {
    disp[i] <- acosh( cosh(Mpi) + 2*sum(sin(pi*p1[,i]/L)^2)) + acosh( cosh(Mpi) + 2*sum(sin(pi*p2[,i]/L)^2))
  }
  for(i in c(1:5)) {
    disp2[i] <- acosh( cosh(Mpi) + 2*sum(sin(pi*p12[,i]/L)^2)) + acosh( cosh(Mpi) + 2*sum(sin(pi*p22[,i]/L)^2))
  }
  arrows(y0=disp, y1=disp, x0=c(1:6)-0.2, x1=c(1:6)+0.2, length=0, col="blue")
  arrows(y0=disp2, y1=disp2, x0=c(1:5)-0.2, x1=c(1:5)+0.2, length=0, col="blue")
  legend("bottomright", legend=c("$E_{\\pi}(\\vec p_1)+E_{\\pi}(\\vec p_2)$"), bty="n", lty=c(1), col=c("blue"))

  ## add the fitted values to the plot
  ii <- which(Epipicut$L == L)
  mirrep <- max(Epipicut[ii,]$irrep)
  mlevel <- max(Epipicut[ii,]$level)
  mframe <- max(Epipicut[ii,]$frame)
  N <- length(Epipicut[Epipicut$L==L & Epipicut$level==1,]$E)
  
  xv <- (Epipicut[ii,]$irrep-1)*mframe + Epipicut[ii,]$frame
  xv[xv==9] <- 6
  E <- bootevalues[1,ii]
  dE <- apply(bootevalues[c(1:10),ii], 2, sd)

  rect(xleft=xv, xright=xv+0.2, ytop=E+dE, ybottom=E-dE, col="red", border="red")


  tikz.finalize(tikzfiles=tikzfiles,clean=TRUE, crop=FALSE)
  options(op)
}
