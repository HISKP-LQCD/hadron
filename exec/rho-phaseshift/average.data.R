source("/hiskp2/urbach/head/hadron/exec/rho-phaseshift/phaseshift.rho.R")
source("/hiskp2/urbach/head/hadron/exec/rho-phaseshift/summarise.R")
source("../../ens.R")
hint <- rep("no", times=5)
source("parameters.R")
source("../../../detect_irrep_frame.R")

pdf(file=paste("histograms", ens, frame, irrep, "pdf", sep="."))
res <- list()
res.all <- list()
res.boot <- list()
for(i in c(1:min(2, N))) {
  if(i == 1) PC <- "pc1"
  if(i == 2) PC <- "pc2"
  if(i == 3) PC <- "pc3"
  cat("hint", i, hint[i], "\n")
  res[[i]] <- summarise.rho(ens=ens, frame=frame, irrep=irrep, PC=PC, hint=hint[i])
  res.all[[i]] <- compute.error.rho(res[[i]], PC=PC)
  res.boot[[i]] <- array(0, dim=c(boot.R+1, 3))
  for(j in c(1:3)) {
    res.boot[[i]][,j] <- compute.boots(res[[i]], index=j)
  }
  cat("Ecm:", i, res.all[[i]]$Ecm, "\n")
  cat("delta:", i, res.all[[i]]$delta, "\n")
}
dev.off()
rm(i,j)

save.image(file=paste("res", ens, frame, irrep, "Rdata", sep="."))

