cont_extr.mpsfixed <- function(ref.points, betas) {

 data3.8 <- read.table("b3.8/chiralplots.dat",
                       col.name=list("mu", "mps", "dmps", "fps", "dfps", "L", "mpcac", "dmpcac", "mN", "dmN", "mDp", "dmDp", "mDpp", "dmDpp"))
 data3.9 <- read.table("b3.9/chiralplots.dat",
                       col.name=list("mu", "mps", "dmps", "fps", "dfps", "L", "mpcac", "dmpcac"))
 data4.05 <- read.table("b4.05/chiralplots.dat",
                        col.name=list("mu", "mps", "dmps", "fps", "dfps", "L", "mpcac", "dmpcac"))
 datar0 <- read.table("r0.dat",
                      col.names=list("beta", "r0", "dr0"))

 # fixed physical situation in r0 mps
 #ref.points <- c(0.70, 0.80, 0.90, 1.00, 1.10)
 N <- length(ref.points)
 load("idx.Rdata")
 intpoints.poly <- array(0, dim=c(3,N))
 rpsq <- ref.points^2

 fit <- optim(par=c(1,1,-1), fn=chisqr.poly, x=data3.8$mps^2, y=data3.8$fps,
              err=data3.8$dfps)
 # now we get the interpolated points from the fitted function
 # fixing r0*mps to the reference values.
 for(i in 1:N) {
   intpoints.poly[1,i] <-  eval.poly(fit$par, rpsq[i]/datar0$r0[1]^2)
 }
 fit <- optim(par=c(1,1,-1), fn=chisqr.poly, x=data3.9$mps^2, y=data3.9$fps,
              err=data3.9$dfps)
 for(i in 1:N) {
   intpoints.poly[2,i] <-  eval.poly(fit$par, rpsq[i]/datar0$r0[2]^2)
 }
 fit <- optim(par=c(1,1,-1), fn=chisqr.poly, x=data4.05$mps^2, y=data4.05$fps,
              err=data4.05$dfps)
 for(i in 1:N) {
   intpoints.poly[3,i] <-  eval.poly(fit$par, rpsq[i]/datar0$r0[3]^2)
 }
   
 

 # here we do a simple linear interpolation to the reference points
 # keeping r0*mps fixed to the reference values
 intpoints <- array(0, dim=c(3,N))
 intpoints[1,] <- get.intpoints(r0=datar0$r0[1], mps=data3.8$mps, fps=data3.8$fps,
                               ref.points=ref.points, sort.idx=intpoints.idx[1,,])
 intpoints[2,] <- get.intpoints(r0=datar0$r0[2], mps=data3.9$mps, fps=data3.9$fps,
                               ref.points=ref.points, sort.idx=intpoints.idx[2,,])
 intpoints[3,] <- get.intpoints(r0=datar0$r0[3], mps=data4.05$mps, fps=data4.05$fps,
                                ref.points=ref.points, sort.idx=intpoints.idx[3,,])

 load(file="../chiralfits/b3.8/bootsamples3.8.RData")
 load(file="../chiralfits/b3.9/bootsamples3.9.RData")
 load(file="../chiralfits/b4.05/bootsamples4.05.RData")


 intpoints.boot <- array(0, dim=c(1000, 3, N))
 for(i in 1:1000) {
   intpoints.boot[i,1,] <- get.intpoints(r0=datar0$r0[1],
                                        mps=bootsamples3.8[i,1,], fps=bootsamples3.8[i,2,],
                                        ref.points=ref.points, sort.idx=intpoints.idx[1,,])
   intpoints.boot[i,2,] <- get.intpoints(r0=datar0$r0[2],
                                        mps=bootsamples3.9[i,1,], fps=bootsamples3.9[i,2,],
                                        ref.points=ref.points, sort.idx=intpoints.idx[2,,])
   intpoints.boot[i,3,] <- get.intpoints(r0=datar0$r0[3],
                                        mps=bootsamples4.05[i,1,], fps=bootsamples4.05[i,2,],
                                        ref.points=ref.points, sort.idx=intpoints.idx[3,,])
 }

 #intpoints.boot[553,2,2] <- intpoints.boot[55,2,2]
 #intpoints.boot[553,2,1] <- intpoints.boot[55,2,1]
 
 intpoints.error <-  array(0, dim=c(3,N))
 for(i in 1:N) {
   for(j in 1:3) {
     intpoints.error[j,i] <- sd(intpoints.boot[,j,i])
   }
 }

 result <- matrix(nrow=3, ncol=(2+2*N))
 for(i in 1:3) {
   result[i,1] <- datar0$r0[i]
   result[i,2] <- datar0$dr0[i]
 }
 for(i in 1:3) {
   for(j in 1:N) {
     result[i,(2*j+1)] <- result[i,1]*intpoints[i, j]
     result[i,(2*(j+1))] <-  sqrt((result[i,1]*intpoints.error[i, j])^2 + (intpoints[i, j]*result[i,2])^2)
   }
 }
 cont <- array(0, dim=c(2,N))
 cont.boot <- array(0, dim=c(1000,2,N))
 ii <- betas
 lb <- length(betas)
 if(lb > 2) {
   for(i in 1:N) {
     cont[,i] <- optim(par=c(0,1), fn=chisqr.cont, x=1./datar0$r0[ii]^2, y=datar0$r0[ii]*intpoints[ii,i],
                       err=sqrt((datar0$r0[ii]*intpoints.error[ii,i])^2+(datar0$dr0[ii]*intpoints[ii,i])^2))$par
   }
   for(j in 1:1000) {
     for(i in 1:N) {
       cont.boot[j,,i] <- optim(par=c(0,1), fn=chisqr.cont, x=1./datar0$r0[ii]^2,
                                y=datar0$r0[ii]*intpoints.boot[j,ii,i],
                                err=sqrt((datar0$r0[ii]*intpoints.error[ii,i])^2+(datar0$dr0[ii]*intpoints[ii,i])^2))$par
     } 
   }
 }
 else {
   for(i in 1:N) {
     cont[2,i] <- sum(datar0$r0[ii]*intpoints[ii,i])/lb
     cont[1,i] <- 0
     #cont[,i] <- get.par(x=1./datar0$r0[ii]^2,
     #                    y=datar0$r0[ii]*intpoints[ii,i])
   }
   for(j in 1:1000) {
     for(i in 1:N) {
       cont.boot[j,2,i] <- sum( datar0$r0[ii]*intpoints.boot[j,ii,i])/lb
       cont.boot[j,1,i] <- 0
       #cont.boot[j,,i] <- get.par(x=1./datar0$r0[ii]^2,
       #                  y=datar0$r0[ii]*intpoints.boot[j,ii,i])
     }
   }
 }
 cont.error <- array(0, dim=c(2,N))
 for(i in 1:N) {
   for(j in 1:2) {
     cont.error[j,i] <- sd(cont.boot[,j,i])
   }
 }
 cont.res <- data.frame(no=c(1:N), r0fps=cont[2,], dr0fps=cont.error[2,])
 return(invisible(
                  list(data3.8=data3.8, data3.9=data3.9, data4.05=data4.05,
                       datar0=datar0, refpoints=ref.points,
                       intpoints=intpoints, intpoints.boot=intpoints.boot,
                       intpoints.error=intpoints.error, result=result,
                       intpoints.idx=intpoints.idx, intpoints.poly=intpoints.poly,
                       cont=cont, cont.boot=cont.boot, cont.error=cont.error,
                       cont.res=cont.res
                  )))
}

cont_extr.fpsfixed <- function(ref.points, betas) {

 data3.8 <- read.table("b3.8/chiralplots.dat",
                       col.name=list("mu", "mps", "dmps", "fps", "dfps", "L", "mpcac", "dmpcac", "mN", "dmN", "mDp", "dmDp", "mDpp", "dmDpp"))
 data3.9 <- read.table("b3.9/chiralplots.dat",
                       col.name=list("mu", "mps", "dmps", "fps", "dfps", "L", "mpcac", "dmpcac"))
 data4.05 <- read.table("b4.05/chiralplots.dat",
                        col.name=list("mu", "mps", "dmps", "fps", "dfps", "L", "mpcac", "dmpcac"))
 datar0 <- read.table("r0.dat",
                      col.names=list("beta", "r0", "dr0"))

 N <- length(ref.points)
 load("idx.Rdata")
 rpsq <- ref.points^2

 # here we do a simple linear interpolation to the reference points
 # keeping r0*mps fixed to the reference values
 intpoints <- array(0, dim=c(3,N))
 intpoints[1,] <- get.intpoints.f(r0=datar0$r0[1], mps=data3.8$mps, fps=data3.8$fps,
                                  ref.points=ref.points, sort.idx=intpoints.idx[1,,])
 intpoints[2,] <- get.intpoints.f(r0=datar0$r0[2], mps=data3.9$mps, fps=data3.9$fps,
                                  ref.points=ref.points, sort.idx=intpoints.idx[2,,])
 intpoints[3,] <- get.intpoints.f(r0=datar0$r0[3], mps=data4.05$mps, fps=data4.05$fps,
                                  ref.points=ref.points, sort.idx=intpoints.idx[3,,])

 load(file="../chiralfits/b3.8/bootsamples3.8.RData")
 load(file="../chiralfits/b3.9/bootsamples3.9.RData")
 load(file="../chiralfits/b4.05/bootsamples4.05.RData")

 
 intpoints.boot <- array(0, dim=c(1000, 3, N))
 for(i in 1:1000) {
   intpoints.boot[i,1,] <- get.intpoints.f(r0=datar0$r0[1],
                                           mps=bootsamples3.8[i,1,], fps=bootsamples3.8[i,2,],
                                           ref.points=ref.points, sort.idx=intpoints.idx[1,,])
   intpoints.boot[i,2,] <- get.intpoints.f(r0=datar0$r0[2],
                                           mps=bootsamples3.9[i,1,], fps=bootsamples3.9[i,2,],
                                           ref.points=ref.points, sort.idx=intpoints.idx[2,,])
   intpoints.boot[i,3,] <- get.intpoints.f(r0=datar0$r0[3],
                                           mps=bootsamples4.05[i,1,], fps=bootsamples4.05[i,2,],
                                           ref.points=ref.points, sort.idx=intpoints.idx[3,,])
 }

 intpoints.error <-  array(0, dim=c(3,N))
 for(i in 1:N) {
   for(j in 1:3) {
     intpoints.error[j,i] <- sd(intpoints.boot[,j,i])
   }
 }

 result <- matrix(nrow=3, ncol=(2+2*N))
 for(i in 1:3) {
   result[i,1] <- datar0$r0[i]
   result[i,2] <- datar0$dr0[i]
 }
 for(i in 1:3) {
   for(j in 1:N) {
     result[i,(2*j+1)] <- result[i,1]*intpoints[i, j]
     result[i,(2*(j+1))] <-  sqrt((result[i,1]*intpoints.error[i, j])^2 + (intpoints[i, j]*result[i,2])^2)
   }
 }
 cont <- array(0, dim=c(2,N))
 cont.boot <- array(0, dim=c(1000,2,N))
 ii <- betas
 lb <- length(betas)
 if(lb > 2) {
   for(i in 1:N) {
     cont[,i] <- optim(par=c(0,1), fn=chisqr.cont, x=1./datar0$r0[ii]^2, y=datar0$r0[ii]*intpoints[ii,i],
                       err=sqrt((datar0$r0[ii]*intpoints.error[ii,i])^2+(datar0$dr0[ii]*intpoints[ii,i])^2))$par
   }
   for(j in 1:1000) {
     for(i in 1:N) {
       cont.boot[j,,i] <- optim(par=c(0,1), fn=chisqr.cont, x=1./datar0$r0[ii]^2,
                                y=datar0$r0[ii]*intpoints.boot[j,ii,i],
                                err=sqrt((datar0$r0[ii]*intpoints.error[ii,i])^2+(datar0$dr0[ii]*intpoints[ii,i])^2))$par
     } 
   }
 }
 else {
   for(i in 1:N) {
     cont[2,i] <- sum(datar0$r0[ii]*intpoints[ii,i])/lb
     cont[1,i] <- 0
     #cont[,i] <- get.par(x=1./datar0$r0[ii]^2,
     #                    y=datar0$r0[ii]*intpoints[ii,i])
   }
   for(j in 1:1000) {
     for(i in 1:N) {
       cont.boot[j,2,i] <- sum( datar0$r0[ii]*intpoints.boot[j,ii,i])/lb
       cont.boot[j,1,i] <- 0
       #cont.boot[j,,i] <- get.par(x=1./datar0$r0[ii]^2,
       #                  y=datar0$r0[ii]*intpoints.boot[j,ii,i])
     }
   }
 }
 cont.error <- array(0, dim=c(2,N))
 for(i in 1:N) {
   for(j in 1:2) {
     cont.error[j,i] <- sd(cont.boot[,j,i])
   }
 }
 cont.res <- data.frame(no=c(1:N), r0fps=cont[2,], dr0fps=cont.error[2,])
 return(invisible(
                  list(data3.8=data3.8, data3.9=data3.9, data4.05=data4.05,
                       datar0=datar0, refpoints=ref.points,
                       intpoints=intpoints, intpoints.boot=intpoints.boot,
                       intpoints.error=intpoints.error, result=result,
                       intpoints.idx=intpoints.idx,
                       cont=cont, cont.boot=cont.boot, cont.error=cont.error,
                       cont.res=cont.res
                  )))
}

get.idx <- function(r0, mps, fps, ref.points) {

  rpsq <- ref.points^2
  mpssq <- (r0*mps)^2
  sort.idx <- array(0, dim=c(length(ref.points),2))
  for(i in 1:length(ref.points)) {
    sort.idx[i,] <- order(abs(mpssq-rpsq[i]))[1:2]
    cat(sort.idx[i,], "\n")
    cat(rpsq[i], " from ", mpssq, "\n")
  }
  cat("\n")  
  return(sort.idx)
}

get.intpoints <- function(r0, mps, fps, ref.points, sort.idx) {

  rpsq <- ref.points^2
  mpssq <- (r0*mps)^2
  intpoints <- rep(0, times=length(ref.points))
  
  for(i in 1:length(ref.points)) {
    sort.idx2 <- order(abs(mpssq-rpsq[i]))
    a <- (fps[sort.idx2[1]]-fps[sort.idx2[2]])/(mpssq[sort.idx2[1]]-mpssq[sort.idx2[2]])
    b <- fps[sort.idx2[1]] - a*mpssq[sort.idx2[1]]
    #a <- (fps[sort.idx[i,1]]-fps[sort.idx[i,2]])/(mpssq[sort.idx[i,1]]-mpssq[sort.idx[i,2]])
    #b <- fps[sort.idx[i,1]] - a*mpssq[sort.idx[i,1]]
    intpoints[i] <- (a*rpsq[i] + b)
  }


#  cat(r0, intpoints, "\n")
  return(intpoints)
}

get.intpoints.f <- function(r0, mps, fps, ref.points, sort.idx) {

  mpssq <- (mps)^2
  r0fps <- r0*fps
  intpoints <- rep(0, times=length(ref.points))
  
  for(i in 1:length(ref.points)) {
    sort.idx2 <- order(abs(r0fps-ref.points[i]))
    a <- (r0fps[sort.idx2[1]]-r0fps[sort.idx2[2]])/(mpssq[sort.idx2[1]]-mpssq[sort.idx2[2]])
    b <- r0fps[sort.idx2[1]] - a*mpssq[sort.idx2[1]]
    intpoints[i] <- sqrt((ref.points[i] - b)/a)
  }
#  cat(r0, intpoints, "\n")
  return(intpoints)
}

get.par <-  function(x, y) {
  par <- c(1,1)
  par[1] <- (y[1]-y[2])/(x[1]-x[2])
  par[2] <- y[1] - par[1]*x[1]
  return(par)
}

chisqr.cont <- function(par, x, y, err) {
  return(sum((y - par[1]*x - par[2])^2/err^2))
}

chisqr.poly <- function(par, x, y, err) {
  deg <- length(par)-1
  vec.pt <- c(0:deg)
  res <- 0.
  for(i in 1:length(x)) {
    vec.x <- (rep(x[i], times=(deg+1)))^vec.pt
    res = res + ((y[i]-sum(par*vec.x))^2/err[i]^2)
  }
  return(res)
}

eval.poly <- function(par, x) {
  deg <- length(par)-1
  vec.x <- (rep(x, times=(deg+1)))^c(0:deg)
  return(sum(par*vec.x))
}
