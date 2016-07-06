source("/hiskp2/urbach/head/hadron/exec/rho-phaseshift/phaseshift.rho.R")
source("/hiskp2/urbach/head/hadron/exec/rho-phaseshift/summarise.R")
source("parameters.R")

res.pc1 <- summarise.rho(ens, frame, PC="pc1")
res.pc2 <- summarise.rho(ens, frame, PC="pc2")

## correct for negative delta where needed
##res.pc1 <- shift.delta(res.pc1)
##res.pc2 <- shift.delta(res.pc2)

res.all.pc1 <- compute.error.rho(res.pc1, PC="pc1")
res.all.pc2 <- compute.error.rho(res.pc2, PC="pc2")

res.boot.pc1 <- array(0, dim=c(boot.R+1, 3))
res.boot.pc2 <- array(0, dim=c(boot.R+1, 3))
for(i in c(1:3)) {
  res.boot.pc1[,i] <- compute.boots(res.pc1, index=i)
  res.boot.pc2[,i] <- compute.boots(res.pc2, index=i)
}
cat("Ecm1:", res.all.pc1$Ecm, "\n")
cat("delta1:", res.all.pc1$delta, "\n")
cat("Ecm2:", res.all.pc2$Ecm, "\n")
cat("delta2", res.all.pc2$delta, "\n")

save(res.pc1, res.pc2, res.all.pc1, res.all.pc2, res.boot.pc1, res.boot.pc2, file=paste("res.", ens, frame, ".Rdata", sep=""))

