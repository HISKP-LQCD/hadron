source("/hiskp2/werner/pipi_I1/2017-08-03_analyse/analyse-A40.32.R")
ens <-  args$ens
boot.R <- args$boot.R

hint <- rep("no", times=5)

source(paste(args$path.to.hadron, "phaseshift.rho.R", sep="/"))
source(paste(args$path.to.hadron, "summarise.R", sep="/"))

## extracts irrep and frame from directory name
## also defines N for th ematrix size
## and path
source(paste(args$path.to.hadron, "/detect_irrep_frame.R", sep="/"))

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

