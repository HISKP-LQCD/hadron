source("summary.R")
source("fit_finite_range.R")

ensembles <- c("A30.32", "A40.32", "A40.24", "A40.20", "A60.24", "A80.24", "A100.24", "B55.32", "D45.32", "B35.32", "B85.24")

## compute the errors for ratio and efm data
## requires to run get_summary-efm and get_summary-ratio first

tres <- array(0, dim=c(length(ensembles),12))

for(i in c(1:length(ensembles))) {
  cat("processing ensemble", ensembles[i], "for combined error analysis\n")
  load(paste("liuming/", ensembles[i], "/", "res-efm.", ensembles[i], ".Rdata", sep=""))
  res.efm <- res
  dim.efm <- dim(res.efm)
  load(paste(ensembles[i], "/", "res.", ensembles[i], ".Rdata", sep=""))
  res.ratio <- res
  dim.ratio <- dim(res.ratio)
  n1 <- min(dim.efm[1], dim.ratio[1])
  tmp <- array(0, dim=c(n1, dim.efm[2]+dim.ratio[2], dim.efm[3]))
  tmp[,c(1:dim.efm[2]),] <- res.efm[c(1:n1),,]
  tmp[,c((dim.efm[2]+1):(dim.efm[2]+dim.ratio[2])),] <- res.ratio[c(1:n1),,]

  if(!interactive()) pdf(onefile=TRUE, file=paste("whists.", ensembles[i], ".pdf", sep=""))
  tmp <- compute.error.piL(tmp)
  tres[i,] <- c(tmp$deltaE, tmp$a0, tmp$mpia0)
  if(!interactive()) dev.off()
}

res <- tres
save(ensembles, res, file="res-allens.Rdata")
write.table(cbind(ensembles, format(res[,c(1:4)], digits=3, scientific=FALSE)), sep=" & ", quote=FALSE, row.names=FALSE, col.names=c("ens", '$\\delta E$', '$d\\delta E$', '$d^-\\delta E$', '$d^+\\delta E$'), file="res-deltaE-allens.dat", eol=" \\\\ \n")
write.table(cbind(ensembles, format(res[,c(9:12)], digits=3, scientific=FALSE)), sep=" & ", quote=FALSE, row.names=FALSE, col.names=c("ens", '$M_\\pi a_0$', '$dM_\\pi a_0$', '$d^-M_\\pi a_0$', '$d^+M_\\pi a_0$'), file="res-mpia0-allens.dat", eol=" \\\\ \n")

rm(res)

res.finite.range.fit <- fit.finite.range()

save(res.finite.range.fit, file="res-finite-range-fit.Rdata")

rm(res.finite.range.fit)

