avercycle <- function(cmicor, cycle.l, ind.vec=c(1,3,4,5,6)) {
  Time <-  2*max(cmicor[,ind.vec[2]])
  Thalf <- max(cmicor[,ind.vec[2]])
  T1 <- Thalf+1
  nrObs <- max(cmicor[,ind.vec[1]])
  Skip <- 0
  Length <- length(cmicor[,ind.vec[3]])
  nrep <- c(length(cmicor[((Skip):Length),ind.vec[3]])/(nrObs*(T1)*4))
  cat("total no of measurements", nrep, "\n")
  start.no <- min(cmicor[,6])
  end.no <- max(cmicor[,6])
  cat("first measurement no", start.no, "\n")
  cat("last  measurement no", end.no, "\n")
  gauge.no <- cmicor[seq(from=1, to=(nrep*(nrObs*(T1)*4)), by=(nrObs*(T1)*4)), 6]
  cycle.no <- floor((end.no -start.no)/cycle.l)
  if((end.no -start.no)/cycle.l >  cycle.no) cycle.no <- cycle.no+1
  cat("no of cycles is", cycle.no, "\n")
  cat("average cycle length", nrep/cycle.no, "\n")
  cycle.ind <- array(0, dim=c(3,cycle.no))
  s <- start.no
  no.thiscycle <- 0
  c.no <- 1
  cycle.ind[2, 1] <- 1
  for(i in 1:nrep) {
    if(gauge.no[i] <= (s + cycle.l)) {
      no.thiscycle <- no.thiscycle+1
    }
    else {
      s <- s+cycle.l
      cycle.ind[1, c.no] <- no.thiscycle
      cycle.ind[3, c.no] <- i-1
      no.thiscycle <- 1
#      cat(c.no, cycle.ind[,c.no], "\n")
      c.no <- c.no + 1
      cycle.ind[2, c.no] <- i
    }
    if(i == nrep) {
      cycle.ind[1, c.no] <- no.thiscycle
      cycle.ind[3, c.no] <- i
#      cat(c.no, cycle.ind[,c.no], "\n")
    }
  }
  cat("sum of measurements used", sum(cycle.ind[1,]), "(check!) \n")
  newcor <- array(0, dim=c(cycle.no*(nrObs*(T1)*4), 6))
  ii <- c(1:3,6)
  
  for(i in 1:cycle.no) {
    for(j in 1:4) {
      newcor[(i-1)*(nrObs*(T1)*4)+c(1:(nrObs*(T1)*4)), ii[j]] <-
        cmicor[(cycle.ind[2, i]-1)*(nrObs*(T1)*4) + c(1:(nrObs*(T1)*4)), ii[j]]
    }
    for(j in cycle.ind[2, i]:cycle.ind[3, i]) {
#    jj <- c(cycle.ind[2,i]:cycle.ind[3, i])-1
      for(k in 4:5) {
#      for(p in 1:(nrObs*(T1)*4)) {
#        newcor[(i-1)*(nrObs*(T1)*4)+p,k] <- sum(cmicor[p+jj*(nrObs*(T1)*4), k])
        newcor[(i-1)*(nrObs*(T1)*4)+c(1:(nrObs*(T1)*4)), k] <-
          sum(newcor[(i-1)*(nrObs*(T1)*4)+c(1:(nrObs*(T1)*4)), k],
              cmicor[(j-1)*(nrObs*(T1)*4) + c(1:(nrObs*(T1)*4)), k]/nrep*cycle.no)
        
      }
    }
    cat(".")
  }
  gc(reset=TRUE)
  return(invisible(newcor))

}
