source("summary.R")
source("fit_finite_range.R")

reread <- FALSE
ensembles <- c("A30.32", "A40.32", "A40.24", "A40.20", "A60.24", "A80.24", "A100.24", "B55.32", "D45.32", "B35.32", "B85.24")
##ensembles <- c("A30.32", "A40.32", "A40.24", "A40.20", "A60.24", "A80.24", "A100.24", "B55.32", "D45.32")

## now for the effective mass data

tres <- array(0, dim=c(length(ensembles),12))

for(i in c(1:length(ensembles))) {
  cat("processing ensemble", ensembles[i], "using effective masses\n")
  if(file.exists(paste("liuming/", ensembles[i], "/", "res-efm.", ensembles[i], ".Rdata", sep="")) && !reread) {
    load(paste("liuming/", ensembles[i], "/", "res-efm.", ensembles[i], ".Rdata", sep=""))
  }
  else {
    load(paste("liuming/", ensembles[i], "/", "res.pc1.TP0.Rdata", sep=""))
    source(paste("liuming/", ensembles[i], "/", "parameters.R", sep=""))
    res <- compile.efm.sldata(ens=ensembles[i], path=paste("liuming/", ensembles[i], "/", sep=""), data=res, L=L)
  }
  if(!interactive()) pdf(onefile=TRUE, file=paste("whists-efm.", ensembles[i], ".pdf", sep=""))
  tmp <- compute.error.piL(res)
  tres[i,] <- c(tmp$deltaE, tmp$a0, tmp$mpia0)
  if(!interactive()) dev.off()
}
res <- tres
rm(tres)

save(ensembles, res, file="res-efm-allens.Rdata")
write.table(cbind(ensembles, format(res[,c(1:4)], digits=3, scientific=FALSE)), sep=" & ", quote=FALSE, row.names=FALSE, col.names=c("ens", '$\\delta E$', '$d\\delta E$', '$d^-\\delta E$', '$d^+\\delta E$'), file="res-deltaE-efm-allens.dat", eol=" \\\\ \n")
write.table(cbind(ensembles, format(res[,c(9:12)], digits=3, scientific=FALSE)), sep=" & ", quote=FALSE, row.names=FALSE, col.names=c("ens", '$M_\\pi a_0$', '$dM_\\pi a_0$', '$d^-M_\\pi a_0$', '$d^+M_\\pi a_0$'), file="res-mpia0-efm-allens.dat", eol=" \\\\ \n")

rm(res)

res.finite.range.fit <- fit.finite.range(path=paste("liuming/", sep=""), type="-efm")

save(res.finite.range.fit, file="res-finite-range-fit-efm.Rdata")

rm(res.finite.range.fit)
