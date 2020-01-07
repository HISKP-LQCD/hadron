One <- diag(1+0.*1i, 4, 4)

Gamma4 <- matrix(c(0+0*1i,  0+0*1i, -1+0*1i,  0+0*1i,
                   0+0*1i,  0+0*1i,  0+0*1i, -1+0*1i,
                   -1+0*1i, 0+0*1i,  0+0*1i,  0+0*1i,
                   0+0*1i, -1+0*1i,  0+0*1i,  0+0*1i),
                 nrow=4, ncol=4)


Gamma5 <- matrix(c(1+0*1i,  0+0*1i,  0+0*1i,  0+0*1i,
                   0+0*1i,  1+0*1i,  0+0*1i,  0+0*1i,
                   0+0*1i,  0+0*1i, -1+0*1i,  0+0*1i,
                   0+0*1i,  0+0*1i,  0+0*1i, -1+0*1i),
                 nrow=4, ncol=4)



#' proton
#' 
#' proton
#' 
#' 
#' @param data Proton correlator matrix
#' @param bc String. Boundary conditions, default 'antiperiodic'
#' @param twistangle Numeric. Angle in twisted boundary conditions.
proton <- function(data, bc="antiperiodic", twistangle=0) {

  Time <- max(data[,1])+1
  Thalf <- Time/2

  datasum <- cbind(data$V2+1i*data$V3, data$V4+1i*data$V5, data$V6+1i*data$V7, data$V8+1i*data$V9)

  # forward in time we need to project with (One+Gamma4)
  # backward in time with (One-Gamma4)
  #
  # in case of twisted mass:
  # u^T Cg_5 d  u_i needs twist rotation on index i
  # so rotate with exp(\pmi \omega \gamma5/5) as
  # appropriate
  #
  # in case of antiperiodic BC in time:
  # in the dynamical update code the APBC are implemented as phase factors
  # at each t, so multiply with
  # exp(3*i*t*pi/T)
  # where the 3 comes from the three involved quarks.

  if(bc == "antiperiodic") {
    for(t in 1:Time) {
      datasum[(4*t-3):(4*t),] = exp(1i*3.*t*pi/Time)*datasum[(4*t-3):(4*t),]
    }
  }

  s <- sin(twistangle/2)
  c <- cos(twistangle/2)
  Cor <- rep(0., times=Thalf)
  Cor2 <- rep(0., times=Thalf)
  Cor[1] <- Re(sum(diag((One + Gamma4) %*% (c*One+1i*s*Gamma5)
                        %*% datasum[(1:4),] %*% (c*One+1i*s*Gamma5)
                        %*% (One + Gamma4))))/4
  Cor2[1] <- Re(sum(diag((One - Gamma4) %*% (c*One+1i*s*Gamma5)
                        %*% datasum[(1:4),] %*% (c*One+1i*s*Gamma5)
                        %*% (One - Gamma4))))/4
  for(t in 2:(Thalf+1)) {
    Cor[t] <- Re(sum(diag(
                          0.5* (((One + Gamma4) %*% (c*One+1i*s*Gamma5)
                                %*% datasum[((4*t-3):(4*t)),] %*% (c*One+1i*s*Gamma5)
                                %*% (One + Gamma4))
                                + ((One - Gamma4) %*% (c*One+1i*s*Gamma5)
                                %*% datasum[((4*(Time-t+2)-3):(4*(Time-t+2))),]
                                %*% (c*One+1i*s*Gamma5) %*% (One - Gamma4))
                                )
                          )))/4

    Cor2[t] <- Re(sum(diag(
                          0.5* ((One - Gamma4) %*% (c*One+1i*s*Gamma5)
                                %*% datasum[((4*t-3):(4*t)),] %*% (c*One+1i*s*Gamma5)
                                %*% (One - Gamma4)
                                + (One + Gamma4) %*% (c*One+1i*s*Gamma5)
                                %*% datasum[((4*(Time-t+2)-3):(4*(Time-t+2))),]
                                %*% (c*One+1i*s*Gamma5) %*% (One + Gamma4)
                                )
                          )))/4
    
  }
  return(data.frame(Cor, Cor2))
}
