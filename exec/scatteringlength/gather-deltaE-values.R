source("analyse_pipi.R")
source("summary.R")

ens <- c("A40.32", "A40.24", "A40.20")
ensdir <- c("A40.32cont", "A40.24cont", "A40.20cont")
filename <- c("Energies-contR2.Rdata")
tablefile <- c("dEres-all-contR2.dat")

##ensdir <- c("A40.32lat", "A40.24lat", "A40.20lat")
##filename <- c("Energies-lat.Rdata")
##tablefile <- c("dEres-all-lat.dat")

Llist <- c(32, 24, 20)
irreplist <- c("A1", "E", "T2")
tplist <- list(c("TP0", "TP1", "TP2", "TP3"), c("TP0"), c("TP0"))
idlist <- list(c(2,2,2,2), c(2), c(1))
boot.R <- 1500
R2 <- TRUE

boot.data <- c()
bdata <- array(0., dim=c(boot.R+1,2))

reslist <- c()
Epi <- data.frame()
Epipi <- data.frame()

for(d in c(1:length(ens))) {
  path <- paste(ensdir[d], "/", sep="")
  L <- Llist[d]
  for(i in c(1:length(irreplist))) {
    irrep <- irreplist[i]
    for(t in c(1:length(tplist[[i]]))) {
      tp <- tplist[[i]][t]
      for(pc.id in c(1:idlist[[i]][t])) {
        mergeEnergies(irrep, ens[d], pc.id, tp, path=path, R2=R2)
        rl <- summariseEnergies(irrep, ens[d], pc.id, tp, path=path, R2=R2)
        Epipi <- rbind(Epipi, data.frame(ens=ens[d], irrep=i, frame=t, level=pc.id, L=L, E=rl$E[1], dE=rl$E[2], dmE=rl$E[3], dpE=rl$E[4], irrepname=irrep, tpname=tp))
        Epi <- rbind(Epi, data.frame(ens=ens[d], irrep=i, frame=t, level=pc.id, L=L, E=rl$Mpi[1], dE=rl$Mpi[2], dmE=rl$Mpi[3], dpE=rl$Mpi[4], irrepname=irrep, tpname=tp))
        bdata[,1] <- rl$Eboots
        bdata[,2] <- rl$Mpiboots
        boot.data <- cbind(boot.data, bdata)
        cat(boot.data[1,], "\n")
        reslist <- rbind(reslist, c(ens[d], irrep, tp, pc.id, L, rl$deltaE, rl$E))
      }
    }
  }
}

print(Epipi)
print(Epi)
## covariance matrix for Epipi                                                                                                                                                                                                                
save(Epipi, Epi, boot.data, file=filename)

write.table(reslist, file=tablefile, quote=FALSE, sep=" ", row.names=FALSE, col.names=FALSE)

