files <- Sys.glob( "C2_pi+-*" )
pion.cor <- readbinarycf(files, obs=0, path="./")
files <- Sys.glob( "C4_1*" )
pipi1 <- readbinarycf(files, obs=0, path="./")
files <- Sys.glob( "C4_2*" )
pipi2 <- readbinarycf(files, obs=0, path="./")

files <- Sys.glob( "C4_3*" )
pipi3 <- readbinarycf(files, obs=0, path="./")
pipi3 <- mul.cf(pipi3, a=2.)
rm(files)
boot.R <- 400
boot.l <- 1
Thalf <- pion.cor$Time/2

pipi.cor <- pipi1 + pipi2 - pipi3
pipi.cor <- bootstrap.cf(pipi.cor, boot.R=boot.R, boot.l=boot.l)
pion.cor <- bootstrap.cf(pion.cor, boot.R=boot.R, boot.l=boot.l)

pion.effmass <- bootstrap.effectivemass(pion.cor, boot.R=boot.R, boot.l=boot.l, type="acosh")
pion.effmass <- fit.effectivemass(pion.effmass, t1=10, t2=23)

compRpipi <- function(c4, c2, Thalf) {
  tt <- c(1:Thalf)
  return((c4[tt] - c4[tt+1])/(c2[tt]^2 - c2[tt+1]^2))
}
## this is R(t+1/2)
Rpipi <- compRpipi(c4=pipi.cor$cf0, c2=pion.cor$cf0, Thalf=Thalf)

Rpipi.tsboot <- array(NA, dim=c(boot.R, Thalf))
for(i in c(1:boot.R)) {
  Rpipi.tsboot[i,] <- compRpipi(pipi.cor$cf.tsboot$t[i,], pion.cor$cf.tsboot$t[i,], Thalf=Thalf)
}

dRpipi <- apply(Rpipi.tsboot, 2, sd)


## fit formel: R(t+1/2) = A*(cosh(dE*t') +sinh(dE*t')*coth(2*Mpi*t'))
## t' = t+1/2-T/2
Rfn <- function(par, tp, m) {
  return(par[1]*(cosh(par[2]*tp) + sinh(par[2]*tp)/tanh(2*m*tp)))
}

fitfn <- function(par, y, t, m, Thalf, M) {
  tp <- t - Thalf 
  z <- Rfn(par, tp, m)
  return(sum((z-y) %*% M %*% (z-y)))
}
tt <- c(11:23)
tphys <- tt-0.5
M <- diag(1./dRpipi[tt]^2)
M <- invertCovMatrix(Rpipi.tsboot[,tt], boot.samples=TRUE)
par <- c(1.6,0.01)
opt.res <- optim(par, fn = fitfn, method="BFGS", M=M, Thalf=Thalf, t=tt, y=Rpipi[tt], m=pion.effmass$opt.res$par[1])
opt.res <- optim(opt.res$par, fn = fitfn, method="BFGS", M=M, Thalf=Thalf, t=tphys, y=Rpipi[tt],
                 m=pion.effmass$opt.res$par[1], control=list(parscale=1./opt.res$par))

par <- opt.res$par
opt.tsboot <- array(NA, dim=c(boot.R,3))
for(i in 1:boot.R) {
  opt <- optim(par, fn = fitfn,
               control=list(parscale=1/par), t=tphys, Thalf=Thalf,
               method="BFGS", M=M, y = Rpipi.tsboot[i,tt], m=pion.effmass$massfit.tsboot[i,1])
  opt.tsboot[i, c(1,2)] <- opt$par
  opt.tsboot[i, 3] <- opt$value
}

plotwitherror(c(1:Thalf)-0.5, Rpipi, dRpipi, ylim=c(1.6,1.9), xlim=c(5,22))

x <- seq(tphys[1], tphys[length(tphys)], 0.01)
rat <- Rfn(opt.res$par, tp=x-Thalf, m=pion.effmass$opt.res$par[1])
lines(x=x, y=rat)

cat("deltaE = ", opt.res$par[2], "+-", sd(opt.tsboot[,2]), "\n")
cat("chi^2 = ", opt.res$value, "\n")
cat("dof = ", length(tt)-2, "\n")
cat("chi^2/dof = ", opt.res$value/(length(tt)-2), "\n")
cat("Qval = ", 1-pchisq(opt.res$value, length(tt)-2), "\n")

cat("mpi*a_0 = ", -opt.res$par[2]*pion.effmass$opt.res$par[1]^2*Thalf^3/4/pi, "+-",
    sd(opt.tsboot[,2]*pion.effmass$massfit.tsboot[,1]^2*Thalf^3/4./pi), "\n")
