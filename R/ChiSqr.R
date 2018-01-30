ChiSqr.singleCor <- function(par, Thalf, x, y, err, tr, sign) {
  ii <- c(1:tr)
  return( sum(((y[ii] - par[1]*par[1]*( CExp(m=par[2], Time=2*Thalf, x=x, sign=sign) ))/err[ii])^2) )
}

ChiSqr.pcac <- function(par, Thalf, x, y, err, tr) {
  ii <- c(1:tr)
  return( sum(((y[ii] - par[1]*par[1]*( CExp(m=par[3], Time=2*Thalf, x=x, sign=+1.) ))/err[ii])^2)
         + sum(((y[ii+tr] - par[1]*par[2]*( CExp(m=par[3], Time=2*Thalf, x=x, sign=-1.) ))/err[ii+tr])^2))
}

ChiSqr.cst <- function(par, y, M) {
  return( sum( (y - par[1]) %*% M %*% (y - par[1]) ) )
}

ChiSqr.smeared <- function(par, Thalf, x, y, err, tr) {
  ii <- c(1:tr)
  cv1 <- CExp(m = abs(par[3]), Time = 2*Thalf, x=x)
  Sumall = (sum(((y[ii] - par[1]*par[2]*(cv1))/err[ii])^2)
            + sum(((y[ii+tr] - par[2]*par[2]*(cv1))/err[ii+tr])^2))
  return(Sumall)
}

dChiSqrdpar.smeared <- function(par, Thalf, x, y, err, tr) {
  ii <- c(1:tr)
  Sumall <- c(0.,0.,0.)
  m1 <- abs(par[3])
  cv1 <- CExp(m=m1, Time=2*Thalf, x=x)
  dcv1 <- dCExpdm(m=m1, Time=2*Thalf, x=x)
  Sumall[1] <- sum( 2*(y[ii] - par[1]*par[2]*cv1)*(-par[2]*cv1) / err[ii]^2 )
  Sumall[2] <- sum( 2*(y[ii] - par[1]*par[2]*cv1)*(-par[1]*cv1) / err[ii]^2 ) +
    sum( 2*(y[ii+tr] - par[2]^2*cv1)*(-2*par[2]*cv1) / err[ii+tr]^2 )
  Sumall[3] <- sum( 2*(y[ii] - par[1]*par[2]*cv1)*(-par[1]*par[2]*dcv1) /  err[ii]^2) +
    sum( 2*(y[ii+tr] - par[2]^2*cv1)*(-par[2]^2*dcv1) /  err[ii+tr]^2)
  return(Sumall)
}

ChiSqr.1mass <- function(par, Thalf, x, y, err, tr, N=2) {
  # index of mass
  # l <- length(par)/no.masses
  ii <- c(1:tr)
  Sumall <- 0.
  m1 <- abs(par[N+1])
  cv1 <- CExp(m=m1, Time=2*Thalf, x=x)
  sv1 <- CExp(m=m1, Time=2*Thalf, x=x, sign=-1.)
  if(N > 1) {
    # PP or 44
    Sumall = Sumall + (sum(((y[ii]
      - par[1]*par[1]*(cv1))/err[ii])^2)
    + sum(((y[ii+tr]
            - par[1]*par[2]*(cv1))/err[ii+tr])^2)
    + sum(((y[ii+2*tr]
            - par[2]*par[2]*(cv1))/err[ii+2*tr])^2))
    
  }
  if(N > 2) {
    # PA or 4A (sinh!)
    Sumall = Sumall + (sum(((y[ii+3*tr]
      - par[1]*par[3]*(sv1))/err[ii+3*tr])^2)
    + sum(((y[ii+4*tr]
            - par[1]*par[4]*(sv1))/err[ii+4*tr])^2)
    + sum(((y[ii+5*tr]
            - par[2]*par[3]*(sv1))/err[ii+5*tr])^2)
    + sum(((y[ii+6*tr]
            - par[2]*par[4]*(sv1))/err[ii+6*tr])^2))
    
    # AA
    Sumall = Sumall + (sum(((y[ii+7*tr]
      - par[3]*par[3]*(cv1))/err[ii+7*tr])^2)
    + sum(((y[ii+8*tr]
            - par[3]*par[4]*(cv1))/err[ii+8*tr])^2)
    + sum(((y[ii+9*tr]
            - par[4]*par[4]*(cv1))/err[ii+9*tr])^2))

  }
  if(N > 4) {
    # 44 or VV
    Sumall = Sumall + (sum(((y[ii+10*tr]
      - par[5]*par[5]*(cv1))/err[ii+10*tr])^2)
    + sum(((y[ii+11*tr]
            - par[5]*par[6]*(cv1))/err[ii+11*tr])^2)
    + sum(((y[ii+12*tr]
            - par[6]*par[6]*(cv1))/err[ii+12*tr])^2))

    # P4 or 4V (sinh!)
    Sumall = Sumall + (sum(((y[ii+13*tr]
      - par[1]*par[5]*(sv1))/err[ii+13*tr])^2)
    + sum(((y[ii+14*tr]
            - par[1]*par[6]*(sv1))/err[ii+14*tr])^2)
    + sum(((y[ii+15*tr]
            - par[2]*par[5]*(sv1))/err[ii+15*tr])^2)
    + sum(((y[ii+16*tr]
            - par[2]*par[6]*(sv1))/err[ii+16*tr])^2))

    # 4A or AV cosh
    Sumall = Sumall + (sum(((y[ii+17*tr]
      - par[3]*par[5]*(cv1))/err[ii+17*tr])^2)
    + sum(((y[ii+18*tr]
            - par[3]*par[6]*(cv1))/err[ii+18*tr])^2)
    + sum(((y[ii+19*tr]
            - par[4]*par[5]*(cv1))/err[ii+19*tr])^2)
    + sum(((y[ii+20*tr]
            - par[4]*par[6]*(cv1))/err[ii+20*tr])^2))
  }
  return(Sumall)
}


ChiSqr.2mass <- function(par, Thalf, x, y, err, tr, N=2, kludge=TRUE) {
  # index of mass
  l <- length(par)/2
  ii <- c(1:tr)
  m1 <- abs(par[N+1])
  m2 <- abs(par[2*N+2])
  Sumall <- 0.
  cv1 <- CExp(m=m1, Time=2*Thalf, x=x)
  cv2 <- CExp(m=m2, Time=2*Thalf, x=x)
  sv1 <- CExp(m=m1, Time=2*Thalf, x=x, sign=-1.)
  sv2 <- CExp(m=m2, Time=2*Thalf, x=x, sign=-1.)
  if(N > 1) {
    # PP or 44
    # there is never any a1 in gamma_i gamma_4
    Sumall = Sumall + (sum(((y[ii]
      - par[1]^2 *(cv1)
      - par[N+2]^2 *(cv2))/err[ii])^2)
    + sum(((y[ii+tr]
            - par[1]*par[2]*(cv1)
            - par[N+2]*par[N+3]*(cv2))/err[ii+tr])^2)
    + sum(((y[ii+2*tr]
            - par[2]^2 *(cv1)
            - par[N+3]^2 *(cv2))/err[ii+2*tr])^2))
    
  }
  if(N > 2) {
    # PA or 4A (sinh!)
    Sumall = Sumall + (sum(((y[ii+3*tr]
      - par[1]*par[3]*(sv1)
      - par[N+2]*par[N+4]*(sv2))/err[ii+3*tr])^2)
      + sum(((y[ii+4*tr]
              - par[1]*par[4]*(sv1)
              - par[N+2]*par[N+5]*(sv2))/err[ii+4*tr])^2)
      + sum(((y[ii+5*tr]
              - par[2]*par[3]*(sv1)
              - par[N+3]*par[N+4]*(sv2))/err[ii+5*tr])^2)
      + sum(((y[ii+6*tr]
              - par[2]*par[4]*(sv1)
              - par[N+3]*par[N+5]*(sv2))/err[ii+6*tr])^2))
      
    # AA (no a1 in A at maximal twist)
    Sumall = Sumall + (sum(((y[ii+7*tr]
      - par[3]*par[3]*(cv1)
      - par[N+4]*par[N+4]*(cv2))/err[ii+7*tr])^2)
      + sum(((y[ii+8*tr]
              - par[3]*par[4]*(cv1)
              - par[N+4]*par[N+5]*(cv2))/err[ii+8*tr])^2)
      + sum(((y[ii+9*tr]
              - par[4]*par[4]*(cv1)
              - par[N+5]*par[N+5]*(cv2))/err[ii+9*tr])^2))

  }
  if(N > 4) {
    # 44 or VV
    Sumall = Sumall + (sum(((y[ii+10*tr]
      - par[5]*par[5]*(cv1)
      - par[N+6]*par[N+6]*(cv2))/err[ii+10*tr])^2)
      + sum(((y[ii+11*tr]
              - par[5]*par[6]*(cv1)
              - par[N+6]*par[N+7]*(cv2))/err[ii+11*tr])^2)
      + sum(((y[ii+12*tr]
              - par[6]*par[6]*(cv1)
              - par[N+7]*par[N+7]*(cv2))/err[ii+12*tr])^2))

    # P4 or 4V (sinh!)
    Sumall = Sumall + (sum(((y[ii+13*tr]
      - par[1]*par[5]*(sv1)
      - par[N+2]*par[N+6]*(sv2))/err[ii+13*tr])^2)
      + sum(((y[ii+14*tr]
              - par[1]*par[6]*(sv1)
              - par[N+2]*par[N+7]*(sv2))/err[ii+14*tr])^2)
      + sum(((y[ii+15*tr]
              - par[2]*par[5]*(sv1)
              - par[N+3]*par[N+6]*(sv2))/err[ii+15*tr])^2)
      + sum(((y[ii+16*tr]
              - par[2]*par[6]*(sv1)
              - par[N+3]*par[N+7]*(sv2))/err[ii+16*tr])^2))

    # 4A or AV cosh
    Sumall = Sumall + (sum(((y[ii+17*tr]
      - par[3]*par[5]*(cv1)
      - par[N+4]*par[N+6]*(cv2))/err[ii+17*tr])^2)
      + sum(((y[ii+18*tr]
              - par[3]*par[6]*(cv1)
              - par[N+4]*par[N+7]*(cv2))/err[ii+18*tr])^2)
      + sum(((y[ii+19*tr]
              - par[4]*par[5]*(cv1)
              - par[N+5]*par[N+6]*(cv2))/err[ii+19*tr])^2)
      + sum(((y[ii+20*tr]
              - par[4]*par[6]*(cv1)
              - par[N+5]*par[N+7]*(cv2))/err[ii+20*tr])^2))
  }
  if(kludge) {
    Sumall = Sumall + 100000*(par[N+2]^2 + par[N+3]^2)
  }
  return(Sumall)
}


ChiSqr.3mass <- function(par, Thalf, x, y, err, tr, N=3, kludge=TRUE) {
  # index of mass
  l <- length(par)/2
  ii <- c(1:tr)
  m1 <- abs(par[N+1])
  m2 <- abs(par[2*N+2])
  m3 <- abs(par[3*N+3])
  cv1 <- CExp(m=m1, Time=2*Thalf, x=x)
  cv2 <- CExp(m=m2, Time=2*Thalf, x=x)
  cv3 <- CExp(m=m3, Time=2*Thalf, x=x)
  sv1 <- CExp(m=m1, Time=2*Thalf, x=x, sign=-1.)
  sv2 <- CExp(m=m2, Time=2*Thalf, x=x, sign=-1.)
  sv3 <- CExp(m=m3, Time=2*Thalf, x=x, sign=-1.)
  Sumall <- 0.
  if(N > 1) {
    # PP or 44
    Sumall = Sumall + sum(((y[ii]
      - par[1]^2 *(cv1)
      - par[N+2]^2 *(cv2)
      - par[2*N+3]^2 *(cv3))/err[ii])^2)
    
    Sumall = Sumall + sum(((y[ii+tr]
      - par[1]*par[2] *(cv1)
      - par[N+2]*par[N+3] *(cv2)
      - par[2*N+3]*par[2*N+4] *(cv3))/err[ii+tr])^2)

    Sumall = Sumall + sum(((y[ii+2*tr]
      - par[2]^2 *(cv1)
      - par[N+3]^2 *(cv2)
      - par[2*N+4]^2 *(cv3))/err[ii+2*tr])^2)
    
  }
  if(N > 2) {
    # PA or 4A (sinh!)
    Sumall = Sumall + sum(((y[ii+3*tr]
      - par[1]*par[3]*(sv1)
      - par[N+2]*par[N+4]*(sv2)
      - par[2*N+3]*par[2*N+5]*(sv3))/err[ii+3*tr])^2)
    Sumall = Sumall + sum(((y[ii+4*tr]
      - par[1]*par[4]*(sv1)
      - par[N+2]*par[N+5]*(sv2)
      - par[2*N+3]*par[2*N+6]*(sv3))/err[ii+4*tr])^2)
    Sumall = Sumall + sum(((y[ii+5*tr]
      - par[2]*par[3]*(sv1)
      - par[N+3]*par[N+4]*(sv2)
      - par[2*N+4]*par[2*N+5]*(sv3))/err[ii+5*tr])^2)
    Sumall = Sumall + sum(((y[ii+6*tr]
      - par[2]*par[4]*(sv1)
      - par[N+3]*par[N+5]*(sv2)
      - par[2*N+4]*par[2*N+6]*(sv3))/err[ii+6*tr])^2)
    
    # AA (no a1 in A)
    Sumall = Sumall + sum(((y[ii+7*tr]
      - par[3]^2*(cv1)
      - par[N+4]^2*(cv2)
      - par[2*N+5]^2*(cv3))/err[ii+7*tr])^2)
    Sumall = Sumall + sum(((y[ii+8*tr]
      - par[3]*par[4]*(cv1)
      - par[N+4]*par[N+5]*(cv2)
      - par[2*N+5]*par[2*N+6]*(cv3))/err[ii+8*tr])^2)
    Sumall = Sumall + sum(((y[ii+9*tr]
      - par[4]^2 *(cv1)
      - par[N+5]^2 *(cv2)
      - par[2*N+6]^2 *(cv3))/err[ii+9*tr])^2)

  }
  if(N > 4) {
    # 44 or VV
    Sumall = Sumall + sum(((y[ii+10*tr]
      - par[5]^2 *(cv1)
      - par[N+6]^2 *(cv2)
      - par[2*N+7]^2 *(cv3))/err[ii+10*tr])^2)
    Sumall = Sumall + sum(((y[ii+11*tr]
      - par[5]*par[6]*(cv1)
      - par[N+6]*par[N+7]*(cv2)
      - par[2*N+7]*par[2*N+8]*(cv3))/err[ii+11*tr])^2)
    Sumall = Sumall + sum(((y[ii+12*tr]
      - par[6]^2 *(cv1)
      - par[N+7]^2 *(cv2)
      - par[2*N+8]^2 *(cv3))/err[ii+12*tr])^2)    

    # P4 or 4V (sinh!)
    Sumall = Sumall + sum(((y[ii+13*tr]
      - par[1]*par[5]*(sv1)
      - par[N+2]*par[N+6]*(sv2)
      - par[2*N+3]*par[2*N+7]*(sv3))/err[ii+13*tr])^2)
    Sumall = Sumall + sum(((y[ii+14*tr]
      - par[1]*par[6]*(sv1)
      - par[N+2]*par[N+7]*(sv2)
      - par[2*N+3]*par[2*N+8]*(sv3))/err[ii+14*tr])^2)
    Sumall = Sumall + sum(((y[ii+15*tr]
      - par[2]*par[5]*(sv1)
      - par[N+3]*par[N+6]*(sv2)
      - par[2*N+4]*par[2*N+7]*(sv3))/err[ii+15*tr])^2)
    Sumall = Sumall + sum(((y[ii+16*tr]
      - par[2]*par[6]*(sv1)
      - par[N+3]*par[N+7]*(sv2)
      - par[2*N+4]*par[2*N+8]*(sv3))/err[ii+16*tr])^2)

    # 4A or AV cosh
    Sumall = Sumall + sum(((y[ii+17*tr]
      - par[3]*par[5]*(cv1)
      - par[N+4]*par[N+6]*(cv2)
      - par[2*N+5]*par[2*N+7]*(cv3))/err[ii+17*tr])^2)
    Sumall = Sumall + sum(((y[ii+18*tr]
      - par[3]*par[6]*(cv1)
      - par[N+4]*par[N+7]*(cv2)
      - par[2*N+5]*par[2*N+8]*(cv3))/err[ii+18*tr])^2)
    Sumall = Sumall + sum(((y[ii+19*tr]
      - par[4]*par[5]*(cv1)
      - par[N+5]*par[N+6]*(cv2)
      - par[2*N+6]*par[2*N+7]*(cv3))/err[ii+19*tr])^2)
    Sumall = Sumall + sum(((y[ii+20*tr]
      - par[4]*par[6]*(cv1)
      - par[N+5]*par[N+7]*(cv2)
      - par[2*N+6]*par[2*N+8]*(cv3))/err[ii+20*tr])^2)
  }
  if(kludge) {
    Sumall = Sumall + (par[N+2]*100000)^2 + (par[N+3]*100000)^2
  }
  return(Sumall)
}
