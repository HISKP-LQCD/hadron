files <- getorderedfilelist(basename="disc.0.13872.0.003.k0v.", last.digits=3)
vdata <- readcmiloopfiles(files)
save(vdata, file="vdata.Rdata")
cat("read done\n")
pi0loops <- extract.loop(vdata, obs=9)
save(pi0loops, file="pi0loops.Rdata")
cat("extract done\n")
## local-local
pi0disc.ll <- computeDisc(pi0loops, smeared=FALSE, real=TRUE, subtract.vev=TRUE)
## local-fuzzed
pi0disc.lf <- computeDisc(cf = pi0loops, cf2 = pi0loops, smeared=FALSE, smeared2=TRUE, real=TRUE, subtract.vev=TRUE)
## pi0disc.fl per construction equal to pi0disc.lf
## fuzzed-fuzzed
pi0disc.ff <- computeDisc(pi0loops, smeared=TRUE, real=TRUE, subtract.vev=TRUE)
save(pi0disc.ll, file="pi0disc.ll.Rdata")
save(pi0disc.lf, file="pi0disc.lf.Rdata")
save(pi0disc.ff, file="pi0disc.ff.Rdata")
#load("pi0disc.ll.Rdata")
#load("pi0disc.lf.Rdata")
#load("pi0disc.fl.Rdata")
#load("pi0disc.ff.Rdata")

files <- getorderedfilelist(basename="outprcvn.")
cmicor <- readcmidatafiles(files)
save(cmicor, file="neutral.Rdata")
#load("neutral.Rdata")
## get the connected bit...
pion0.con <- extract.obs(cmicor,  vec.obs=c(5))
## get the disconnected bit as a concatenation...
pion0.disc <- c(pi0disc.ll, pi0disc.lf, pi0disc.lf, pi0disc.ff)
## now add connected and disconnected bits...
## factors depend on your normalisation
pion0.cor <- add.cf(pion0.con, pion0.disc, a=-0.5, b=2.)
## bootstrap the correlators
pion0.cor <- bootstrap.cf(pion0.cor, boot.R=400, boot.l=1)
