source("summary.R")

ensembles <- c("A30.32", "A40.32", "A40.24", "A40.20", "A60.24", "A80.24", "A100.24", "B55.32", "D45.32", "B35.32", "B85.24")

## compute the errors for ratio data

res <- array(0, dim=c(length(ensembles),24))
bootres <- array()
                        
for(i in c(1:length(ensembles))) {
  tmpres <- compile.ratio.sldata(ens=ensembles[i], path=paste(ensembles[i], "/", sep=""))
  if(i == 1) {
    R <- length(tmpres[,1,1])
    bootres <- array(0, dim=c(R, length(ensembles), 6))
  }
  else if(R != length(tmpres[,1,1])) {
    stop(paste("number of bootstrap samples for ensemble ", ensembles[i], " does not match\n"), sep="")
  }
  if(!interactive()) pdf(onefile=TRUE, file=paste("whists-ratio.", ensembles[i], ".pdf", sep=""))
  tmp <- compute.error.piL(tmpres)
  res[i,] <- c(tmp$deltaE, tmp$a0, tmp$mpia0, tmp$qcotdelta, tmp$Ecm, tmp$qsq)
  k <- 1
  for(j in c(1,2,3,7,8,9)) {
    bootres[,i,k] <- compute.boots(res=tmpres, index=j)
    k <- k+1
  }
  if(!interactive()) dev.off()
}

save(ensembles, res, bootres, file="res-ratio-allens.Rdata")
write.table(cbind(ensembles, format(res[,c(1:4)], digits=3, scientific=FALSE)), sep=" & ", quote=FALSE, row.names=FALSE, col.names=c("ens", '$\\delta E$', '$d\\delta E$', '$d^-\\delta E$', '$d^+\\delta E$'), file="res-deltaE-ratio-allens.dat", eol=" \\\\ \n")
write.table(cbind(ensembles, format(res[,c(9:12)], digits=3, scientific=FALSE)), sep=" & ", quote=FALSE, row.names=FALSE, col.names=c("ens", '$M_\\pi a_0$', '$dM_\\pi a_0$', '$d^-M_\\pi a_0$', '$d^+M_\\pi a_0$'), file="res-mpia0-ratio-allens.dat", eol=" \\\\ \n")
write.table(cbind(ensembles, format(res[,c(13:16)], digits=3, scientific=FALSE)), sep=" & ", quote=FALSE, row.names=FALSE, col.names=c("ens", '$q\\cot\\delta_0$', '$dq\\cot\\delta_0$', '$d^-q\\cot\\delta_0$', '$d^+q\\cot\\delta_0$'), file="res-qcotdelta-ratio-allens.dat", eol=" \\\\ \n")
write.table(cbind(ensembles, format(res[,c(21:24)], digits=3, scientific=FALSE)), sep=" & ", quote=FALSE, row.names=FALSE, col.names=c("ens", '$q^2$', '$dq^2$', '$d^-q^2$', '$d^+q^2$'), file="res-qsq-ratio-allens.dat", eol=" \\\\ \n")
rm(res, bootres, tmp, tmpres)

