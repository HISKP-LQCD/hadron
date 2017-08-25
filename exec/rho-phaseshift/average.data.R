args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)!=1) {
    stop("Please specify the input file!", call.=FALSE)
} else  {
    # default output file
    input.file = args[1]
}

source(input.file)

ens <-  args$ens
boot.R <- args$boot.R
dirs <- args$dirs

hint <- rep("no", times=5)

source(paste(args$path.to.hadron, "/exec/rho-phaseshift/phaseshift.rho.R", sep="/"))
source(paste(args$path.to.hadron, "/exec/rho-phaseshift/summarise.R", sep="/"))

for(dir in dirs){

  setwd(paste(args$output.path, "5_fit-data", dir, sep="/"))
  sink("average.data.log", append=FALSE)

  ## extracts irrep and frame from directory name
  ## also defines N for the matrix size
  ## and path
  source(paste(args$path.to.hadron, "/exec/rho-phaseshift/detect_irrep_frame.R", sep="/"))
  
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
}  

sink()
