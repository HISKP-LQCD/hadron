source("summary.R")
source("fit_finite_range.R")

ensembles <- c("A30.32", "A40.32", "A40.24", "A40.20", "A60.24", "A80.24", "A100.24", "B55.32", "D45.32", "B35.32", "B85.24", "D15.48")

## compute the errors for ratio data

res <- array(0, dim=c(length(ensembles),12))

for(i in c(1:length(ensembles))) {
  tmp <- compile.ratio.sldata(ens=ensembles[i], path=paste(ensembles[i], "/", sep=""))
  if(!interactive()) pdf(onefile=TRUE, file=paste("whists-ratio.", ensembles[i], ".pdf", sep=""))
  tmp <- compute.error.piL(tmp)
  res[i,] <- c(tmp$deltaE, tmp$a0, tmp$mpia0)
  if(!interactive()) dev.off()
}

save(ensembles, res, file="res-ratio-allens.Rdata")
write.table(cbind(ensembles, format(res[,c(1:4)], digits=3, scientific=FALSE)), sep=" & ", quote=FALSE, row.names=FALSE, col.names=c("ens", '$\\delta E$', '$d\\delta E$', '$d^-\\delta E$', '$d^+\\delta E$'), file="res-deltaE-ratio-allens.dat", eol=" \\\\ \n")
write.table(cbind(ensembles, format(res[,c(9:12)], digits=3, scientific=FALSE)), sep=" & ", quote=FALSE, row.names=FALSE, col.names=c("ens", '$M_\\pi a_0$', '$dM_\\pi a_0$', '$d^-M_\\pi a_0$', '$d^+M_\\pi a_0$'), file="res-mpia0-ratio-allens.dat", eol=" \\\\ \n")

rm(res)

res.finite.range.fit <- fit.finite.range(typ="-ratio")

save(res.finite.range.fit, file="res-finite-range-fit-ratio.Rdata")

rm(res.finite.range.fit)

