source("../analyse_pipi.R")
redo <- FALSE
source("parameters.R")

cargs <- commandArgs(trailingOnly = TRUE)
if(length(cargs) > 0) {
  tr1 <- as.integer(cargs[1])
}

if(length(cargs) > 1) {
  t1 <- as.integer(cargs[2])
}

if(length(cargs) > 2) {
  tr2 <- as.integer(cargs[3])
}

if(length(cargs) > 3) {
  t2 <- as.integer(cargs[4])
}


cat("tr1 =", tr1, ", t1 =", t1, "\n")

data <- read.pipidata(path=srcpath, T=T, Tformat=TRUE)
if(!interactive()) pdf(onefile=TRUE, file=paste("L-analysis", t1, ".", t2, ".", tr1, ".", tr2, ".pdf", sep=""))

data$pion.cor <- bootstrap.cf(data$pion.cor, boot.R=boot.R, boot.l=boot.l)
data$pion.effmass <- bootstrap.effectivemass(data$pion.cor, boot.R=boot.R, boot.l=boot.l, type="solve")
data$pion.effmass <- fit.effectivemass(data$pion.effmass, t1=t1, t2=t2, useCov=useCov)


for(ta in c(tr1:tr2)) {
  if(tr2 - ta < 5) break
  if(file.exists(paste("data.", ens, ".", t1, ".", t2, ".tr1", ta, ".tr2", tr2, ".Rdata", sep="")) && !redo) {
    load(paste("data.", ens, ".", t1, ".", t2, ".tr1", ta, ".tr2", tr2, ".Rdata", sep=""))
  }
  else {
    pipi.data <- try(run.pipi.analysis.ratio(data, t1=t1, t2=t2, tr1=ta, tr2=tr2, ens=ens, L=L, boot.R=boot.R, boot.l=boot.l))
  }

  if(!inherits(pipi.data, "try-error")) {
    meta.data <- list(ens=ens, T=T, L=L, srcpath=srcpath, t1=t1, t2=t2, tr1=ta, tr2=tr2)
    save(meta.data, pipi.data, file=paste("data.", ens, ".", t1, ".", t2, ".tr1", ta, ".tr2", tr2, ".Rdata", sep=""))
  
    summary(pipi.data)
    plot(pipi.data)
  }
  else {
    cat("Skiping combination", ta, "-", tr2, "due to non converging fit\n")
  }
}

if(!interactive()) dev.off()
