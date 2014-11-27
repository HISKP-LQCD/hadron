## read pion data from file
if(!file.exists(paste(path,"pidata.Rdata", sep=""))) {
  pidata <- readtextcf(paste(srcpath, "pi_corr_p0.dat", sep=""), T=T, check.t=1, path=path)
  pidata <- bootstrap.cf(pidata, boot.R=boot.R, boot.l=boot.l, seed=seed)
  save(pidata, file=paste(path,"pidata.Rdata", sep=""))
}
load(file=paste(path, "pidata.Rdata", sep=""))

if(redofit) {
  cat("...determining pion mass from effective mass and matrix fit\n")
  t1 <- 9
  t2 <- T/2
  for(ta in c(t1:t2)) {
    if(t2-ta <= 5) break
    pion.effectivemass <- fit.effectivemass(bootstrap.effectivemass(cf=pidata, type="solve"), t1=ta, t2=t2-1, useCov=TRUE)
    summary(pion.effectivemass)
    if(pion.effectivemass$Qval > 0.01 && pion.effectivemass$Qval < 0.99) save(pion.effectivemass, file=paste(path, "pion.p0.effectivemass.", ta, ".", t2, ".Rdata", sep=""))
    pion.matrixfit <- matrixfit(pidata, t1=ta, t2=t2, useCov=TRUE)
    if(pion.matrixfit$Qval > 0.01 && pion.matrixfit$Qval < 0.99) save(pion.matrixfit, file=paste(path, "pion.p0.matrixfit.", ta, ".", t2, ".Rdata", sep=""))
  }
}
