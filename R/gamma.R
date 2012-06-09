gamma <- array(complex(16), dim=c(4,4,4))

## gamma0
gamma[1,3,1] <- -1.
gamma[1,4,2] <- -1.
gamma[1,1,3] <- -1.
gamma[1,2,4] <- -1.

## gamma1
gamma[2,4,1] <- -1i
gamma[2,3,2] <- -1i
gamma[2,2,3] <- +1i
gamma[2,1,4] <- +1i

## gamma2
gamma[3,4,1] <- -1.
gamma[3,3,2] <- +1.
gamma[3,2,3] <- +1.
gamma[3,1,4] <- -1.

## gamma3
gamma[4,3,1] <- -1i
gamma[4,4,2] <- +1i
gamma[4,1,3] <- +1i
gamma[4,2,4] <- -1i

sigmamunu <- array(complex(16), dim=c(4,4,4,4))

tmatrix <- array(0., dim=c(4,4))
tmatrix[1,1] <- 1
tmatrix[1,2] <- 10
tmatrix[2,1] <- 100
tmatrix[2,2] <- 1000

tmatrix[3,3] <- 10000
tmatrix[3,4] <- 100000
tmatrix[4,3] <- 1000000
tmatrix[4,4] <- 10000000

for(mu in 1:4) {
  for(nu in 1:4) {
    sigmamunu[mu, nu,,] <- gamma[mu,,] %*% gamma[nu,,]
  }
}

Tr <- function(x) {
  return(sum(diag(x)))
}
