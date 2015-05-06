## needs also ratio.R and analyse_pipi.R
source("summary.R")

boot.R <- 1500
boot.l <- 1

ens <- c("cA2.09.48")

if(file.exists("parameters.R")) {
  source("parameters.R")
}


if(!file.exists(paste("summary-res.averx.", ens, ".Rdata", sep=""))) {
  res <- compile.averxdata(ens=ens)
  
  save(res, file=paste("summary-res.averx.", ens, ".Rdata", sep=""))
}
load(paste("summary-res.averx.", ens, ".Rdata", sep=""))

if(!interactive()) pdf(file=paste("whist.", ens, ".pdf", sep=""))

x1 <- estimate.error(res, index=1, main=c("<x>_1 weighted histogram"), Qval=res[1,,c(4,5)])
x1boot <- compute.boots(res, index=1, Qval=res[1,,c(4,5)])

x2 <- estimate.error(res, index=2, main=c("<x>_2 weighted histogram"), Qval=res[1,,c(4,5)])
x2boot <- compute.boots(res, index=2, Qval=res[1,,c(4,5)])

Mpi <- estimate.error(res, index=3, main=c("Mpi weighted histogram"), Qval=res[1,,c(4,5)], piononly=TRUE)
Mpiboot <- compute.boots(res, index=3, Qval=res[1,,c(4,5)], piononly=TRUE)

save(x1, x2, Mpi, file=paste("final-res.averx.", ens, ".Rdata", sep=""))

if(!interactive()) dev.off()
