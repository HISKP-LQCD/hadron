
averx <- function(data3pt, data2pt, ind.vec=c(1,2), ind.vec2pt=c(1,2), skip=0, t1, t2, mps, par=c(0.6, 0.15), S=1.5, method="uwerr") {

  Time <- max(data3pt[,ind.vec[1]])+1
  Thalf <- Time/2
  T1 <- Thalf+1
  Length <- min(length(data3pt[,ind.vec[1]]), length(data2pt[,ind.vec2pt[1]]))
  nrObs <- 1
  Skip <- (skip*(Time)*nrObs+1)
  
  W <- -array(data3pt[((Skip):Length),ind.vec[2]], 
             dim=c(nrObs*(Time),(length(data3pt[((Skip):Length),ind.vec[2]])/(nrObs*(Time)))))
  W2pt <- array(data2pt[((Skip):Length),ind.vec2pt[2]], 
                dim=c(nrObs*(Time),(length(data2pt[((Skip):Length),ind.vec2pt[2]])/(nrObs*(Time)))))

  Z <- array(0., dim=c(nrObs*(T1),(length(data3pt[((Skip):Length),ind.vec[2]])/(nrObs*(Time)))))
  Z2pt <- array(0., dim=c(nrObs*(T1),(length(data2pt[((Skip):Length),ind.vec2pt[2]])/(nrObs*(Time)))))

  Cor <- rep(0, times=nrObs*T1)
  Err <- rep(0, times=nrObs*T1)
  Cor2pt <- rep(0, times=nrObs*T1)
  Err2pt <- rep(0, times=nrObs*T1)
  
  Z[1,] <- W[1,]
  Z2pt[1,] <- W2pt[1,]
  for(i in 2:T1) {
    Z2pt[i,] <- (W2pt[i,] + W2pt[(Time-i+2),])*0.5
  }
  for(i in 2:T1) {
    Z[i,] <- (W[i,] + W[(Time-i+2),])*0.5
  }
  ii <- c((t1+1):(t2+1))
  Zall <- rbind(Z[ii,], Z2pt[ii,], Z2pt[Thalf,])

  for(i in 1:T1) {
    Cor2pt[i] <- mean(Z2pt[i,])
    Err2pt[i] <- uwerrprimary(Z2pt[i,])$dvalue
  }
  for(i in 1:T1) {
    Cor[i] <- mean(Z[i,])/Cor2pt[T1]
    Err[i] <- uwerrderived(f=get.ratio, data=rbind(Z[i,],Z2pt[T1,]), S=S)$dvalue
  }
  rm(Z, Z2pt, W, W2pt)
  
  averx.fit <- optimise(f=chisqr.averx, c(0.,1.), Cor=Cor[ii], Err=Err[ii])
  blub <-  get.averx(Cor=c(Cor[ii]*Cor2pt[T1],Cor2pt[ii],Cor2pt[T1]),
                     Err=c(Err[ii], Err2pt[ii]), Time=Time, t1=t1, t2=t2, par=par)
  if(missing(mps)) {
    mps = averx.fit$minimum/blub
  }
  averx.uwerr <- NULL
  if(method == "uwerr") {
    averx.uwerr <- uwerrderived(f=get.averx, data=Zall, S=S, Err=c(Err[ii], Err2pt[ii]), Time=Time, t1=t1, t2=t2, par=par)
  }
  res <- list(averx=averx.fit$minimum/mps, daverx=averx.uwerr$dvalue,
              data=data.frame(Cor=Cor, Err=Err), fit.uwerr=averx.uwerr,
              mps=mps, t1=t1, t2=t2, N=averx.uwerr$N)
  attr(res, "class") <- c("averx", "list")  
  return(invisible(res))
}

chisqr.averx <- function(x, Cor, Err) {
  return( sum(((x-Cor)/Err)^2) )
}

ChiSqr.mpi <- function(par, Time, x, y, err) {
  Sumall = sum(((y
    - par[1]*par[1]*(CExp(m=abs(par[2]), Time=Time, x=x)))/err)^2)
  return(Sumall)
}

get.ratio <- function(data) {
  return(data[1]/data[2])
}

get.averx <- function(Cor, Err, Time, t1, t2, par) {
  Half <- (length(Cor)-1)/2
  ii1 <- c(1:Half)
  ii2 <- c((Half+1):(2*Half))
  # 2pt function at T/2
  ii3 <- c(length(Cor))
  #fit 3pt function to a constant
  fitit <- optimise(f=chisqr.averx, c(0.,1.), Cor=Cor[ii1]/Cor[ii3], Err=Err[ii1])
  #get mpi
  fit.mpi <- optim(par=par, ChiSqr.mpi, method="BFGS", control=list(trace=0),Time=Time,
                    x=c((t1):(t2)), y=Cor[ii2], err=Err[ii2])
  # cat(fitit$minimum, abs(fit.mpi$par[2]), "\n")
  return(fitit$minimum/abs(fit.mpi$par[2]))
}

summary.averx <- function(averx) {
  cat("<x>      =", averx$averx, "\n")
  cat("error    =", averx$daverx, "\n")
  if(!is.null(averx$fit.uwerr)) {
    cat("tauint   =", averx$fit.uwerr$tauint, "\n")
    cat("dtauint  =", averx$fit.uwerr$dtauint, "\n")
    cat("Wopt     =", averx$fit.uwerr$Wopt, "\n")
  }
  cat("no.meas  =", averx$fit.uwerr$N, "\n")
}

print.averx <- function(averx) {
  summary.averx(averx)
}


