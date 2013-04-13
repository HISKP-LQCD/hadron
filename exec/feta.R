feta <- function(amps, mpi, mK, ml, ms) {
  ## separate light and strange angles
  pl <- atan(amps$a3/amps$a1)
  phil <- pl/2/pi*360
  ps <-  - atan(amps$a2/amps$a4)
  phis <- ps/2/pi*360
  ## mean angle
  p <- atan(sqrt(-amps$a3/amps$a1*amps$a2/amps$a4))
  phi <- p/2/pi*360
  ## decay constants using separate angles
  fl <- sqrt(2)*2*ml*amps$a1/cos(pl)/mpi^2
  fs <- sqrt(2)*2*ms*amps$a2/(-sin(ps))/(2*mK^2-mpi^2)
  ##fl <- 2*ml*amps$a3/sin(pl)/mpi^2
  ##fs <- 2*ms*amps$a4/(cos(ps))/(2*mK^2-mpi^2)
  flovfs <- amps$a1/amps$a4/cos(pl)*cos(ps)
  ## decay constants using mean angle
  flp <- 2*ml*amps$a1/cos(p)/mpi^2
  fsp <- 2*ms*amps$a2/(-sin(p))/(2*mK^2-mpi^2)
  return(invisible(c(phil, phis, phi, fl, fs, flp, fsp, flovfs)))
}

compfs <- function(ml, ms, path, debug=FALSE) {
  ##Ameta <- read.table(paste(path, "eta_trick/L_2x2.blocking.amplitudes.cosh.state.0.all", sep=""))
  ##Ametap <- read.table(paste(path, "eta_trick/L_2x2.blocking.amplitudes.cosh.state.1.all", sep=""))
  ##am <- data.frame(a1 = Ameta$V2, a2 = Ameta$V3, a3 = Ametap$V2, a4 = Ametap$V3)
  am <- read.table(paste(path,"eta_trick/L_2x2.blocking.amplitudes.all", sep=""),
                   col.names=c("i", "a1", "da1", "a2", "da2", "a3", "da3", "a4", "da4"))
  ##meta <- read.table(paste(path,"eta_trick/factorize.L_2x2.blocking.masses.all", sep=""))
  meta <- read.table(paste(path,"eta_trick/L_2x2.blocking.masses.cosh.state.0.all", sep=""))
  metap <- read.table(paste(path,"eta_trick/L_2x2.blocking.masses.cosh.state.1.all", sep=""))
  mpi <- read.table(paste(path,"pion+/LF_4x4.blocking.masses.cosh.state.0.all", sep=""))
  if(file.exists(paste(path,"kaon/LF_4x4.blocking.masses.cosh.state.0.all", sep=""))) {
    mK <- read.table(paste(path,"kaon/LF_4x4.blocking.masses.cosh.state.0.all", sep=""))
  }
  if(file.exists(paste(path,"kaon/LF_8x8.masses.cosh.state.0.all", sep=""))) {
    mK <- read.table(paste(path,"kaon/LF_8x8.masses.cosh.state.0.all", sep=""))
  }
  else if(file.exists(paste(path,"kaon/resampling.dat", sep=""))) {
    mK <- read.table(paste(path,"kaon/resampling.dat", sep=""))
  }
  xbar <- feta(amps=am[1,], mpi$V2[1], mK$V2[1], ml, ms)
  N <- length(am$a2)
  xb <- array(0, dim=c(N, length(xbar)+2))
  xb[1,c(1:8)] <- xbar
  xb[1,9] <- mpi$V2[1]
  xb[1,10] <- meta$V2[1]^2 + metap$V2[1]^2 - 2*mK$V2[1]^2
  for(i in 2:N) {
    xb[i,c(1:8)] <- feta(amps=am[i,], mpi$V2[i], mK$V2[i], ml, ms)
    xb[i,9] <- mpi$V2[i]
    xb[i,10] <- meta$V2[i]^2 + metap$V2[i]^2 - 2*mK$V2[i]^2
  }
  if(debug) {
    ii <- c(2:N)
    ##print(xb[ii,3])
    cat("phi_l    =\t", xbar[1], "+-", sd(xb[ii,1]), "\n");
    cat("phi_s    =\t", xbar[2], "+-", sd(xb[ii,2]), "\n");
    cat("phi      =\t", xbar[3], "+-", sd(xb[ii,3]), "\n");
    cat("f_l      =\t", xbar[4], "+-", sd(xb[ii,4]), "\n");
    cat("f_s      =\t", xbar[5], "+-", sd(xb[ii,5]), "\n");
    cat("f_l(phi) =\t", xbar[6], "+-", sd(xb[ii,6]), "\n");
    cat("f_s(phi) =\t", xbar[7], "+-", sd(xb[ii,7]), "\n");
  }
  return(invisible(xb))
}

ensembles <- dir("/dsk/lattice02-0/ottnad/data/",
                 pattern = "^[ABD]",full.names=FALSE, ignore.case = FALSE)
Nens <- length(ensembles)
mul <- c(0.010, 0.010, 0.0030, 0.0040, 0.0040, 0.0040, 0.0050, 0.0060, 0.0080, 0.0080,
         0.0025, 0.0035, 0.0035, 0.0055, 0.0075, 0.0085,
         0.0015, 0.0020, 0.0030, 0.0045)
mus <- c(0.022643, 0.0179509, rep(0.022643, times=7), 0.0179509,
         rep(0.018397, times=6),
         rep(0.01622195, times=3), 0.01300039)
r0 <- c(rep(5.231, times=10),
        rep(5.710, times=6),
        rep(7.538, times=5))
dr0 <- c(rep(0.038, times=10),
        rep(0.041, times=6),
        rep(0.058, times=5))

bres <- array(0, dim=c(Nens, 1000, 12))
dres <- array(0, dim=c(Nens, 12))
for(i in 1:Nens) {
  if(!(i %in% c(4,7,13,18,21))) {
    path <- paste("/dsk/lattice02-0/ottnad/data/", ensembles[i], "/analyse/", sep="")
    cat(ensembles[i], i, "\n")
    bres[i,,c(1:10)] <- compfs(mul[i], mus[i], path)
    bres[i,,11] <- bres[i,,1] - bres[i,,2]
    bres[i,,12] <- bres[i,,4]/bres[i,,5]
    dres[i,] <- apply(bres[i,,], 2, sd)
  }
  else {
    bres[i,,] <- NA
    dres[i,] <- NA
  }
}
res <- bres[,1,]
result <- data.frame(ensemble=ensembles, r0=r0, dr0=dr0,
                     phi_l = res[,1], dphi_l = dres[,1],
                     phi_s = res[,2], dphi_s = dres[,2],
                     phi = res[,3], dphi = dres[,3],
                     f_l = res[,4], df_l = dres[,4],
                     f_s = res[,5], df_s = dres[,5],
                     mpi = res[,9], dmpi = dres[,9],
                     dphi = res[,11], ddphi = dres[,11],
                     flovfs = res[,8], dflovfs = dres[,8],
                     mU = res[,10], dmU = dres[,10],
                     flovfs2 = res[,12], dflovfs2 = dres[,12])



write.table(result, file="feta_res.dat", col.names=FALSE, row.names=FALSE, quote=FALSE)
save(result, file="feta_res.Rdata")
bfetares <- bres
save(bfetares, file="feta_bootstrap.Rdata")
rm(bres)
