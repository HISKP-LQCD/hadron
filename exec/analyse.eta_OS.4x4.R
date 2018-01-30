## light connected
files <- getorderedfilelist(path="/hiskp2/fzimmerm/D45/", basename="neutral", last.digits=4)
conf.index <- getorderedconfignumbers(path="/hiskp2/fzimmerm/D45/", basename="neutral", last.digits=4)
cmicor <- readcmidatafiles(files)

con.l <- extract.obs(cmicor,  vec.obs=c(5))
con.l <- addConfIndex2cf(con.l, conf.index)
save(con.l, file="con.l.Rdata")

## strange connected
files <- getorderedfilelist(path="/hiskp2/fzimmerm/D45/etaS/", basename="neutral", last.digits=4)
conf.index <- getorderedconfignumbers(path="/hiskp2/fzimmerm/D45/etaS", basename="neutral", last.digits=4)
cmicor <- readcmidatafiles(files, verbose=TRUE)

con.s <- extract.obs(cmicor,  vec.obs=c(5))
con.s <- addConfIndex2cf(con.s, conf.index)
save(con.s, file="con.s.Rdata")

## strange connected Kp matching
files <- getorderedfilelist(path="/hiskp2/fzimmerm/D45/kp/", basename="neutral", last.digits=4)
conf.index <- getorderedconfignumbers(path="/hiskp2/fzimmerm/D45/kp", basename="neutral", last.digits=4)
cmicor <- readcmidatafiles(files, verbose=TRUE)

con.s.kp <- extract.obs(cmicor,  vec.obs=c(5))
con.s.kp <- addConfIndex2cf(con.s.kp, conf.index)
save(con.s.kp, file="con.s.kp.Rdata")


## light loops
files <- getorderedfilelist(path="/hiskp2/fzimmerm/D45/", basename="disc.0.156315.0.0045.k0vv.", last.digits=4)
vvloops <- readcmiloopfiles(files, obs=9, verbose=TRUE)

loops.l <- extract.loop(vvloops, obs=9)
save(loops.l, file="loops.l.Rdata")

## strange loops
files <- getorderedfilelist(path="/hiskp2/fzimmerm/D45/etaS/", basename="disc.0.156315.0.0118.k0vv.", last.digits=4)
vvloops <- readcmiloopfiles(files, obs=9, verbose=TRUE)

loops.s <- extract.loop(vvloops, obs=9)
save(loops.s, file="loops.s.Rdata")

## strange loops Kp matching
files <- getorderedfilelist(path="/hiskp2/fzimmerm/D45/kp/", basename="disc.0.156315.0.01488.k0vv.", last.digits=4)
vvloops <- readcmiloopfiles(files, obs=9, verbose=TRUE)

loops.s.kp <- extract.loop(vvloops, obs=9)
save(loops.s.kp, file="loops.s.kp.Rdata")


## light local-local
disc.l.ll <- computeDisc(loops.l, smeared=FALSE, real=FALSE, subtract.vev=FALSE)
## local-fuzzed
disc.l.lf <- computeDisc(cf = loops.l, cf2 = loops.l, smeared=FALSE, smeared2=TRUE, real=FALSE, real2=FALSE,
                       subtract.vev=FALSE, subtract.vev2=FALSE)
## disc.l.fl per construction equal to disc.l.lf
## fuzzed-fuzzed
disc.l.ff <- computeDisc(loops.l, smeared=TRUE, real=FALSE, subtract.vev=FALSE)
save(disc.l.ll, file="disc.l.ll.Rdata")
save(disc.l.lf, file="disc.l.lf.Rdata")
save(disc.l.ff, file="disc.l.ff.Rdata")

## strange local-local
disc.s.ll <- computeDisc(loops.s, smeared=FALSE, real=FALSE, subtract.vev=FALSE)
## local-fuzzed
disc.s.lf <- computeDisc(cf = loops.s, cf2 = loops.s, smeared=FALSE, smeared2=TRUE, real=FALSE, real2=FALSE,
                       subtract.vev=FALSE, subtract.vev2=FALSE)
## disc.s.fl per construction equal to disc.s.lf
## fuzzed-fuzzed
disc.s.ff <- computeDisc(loops.s, smeared=TRUE, real=FALSE, subtract.vev=FALSE)
save(disc.s.ll, file="disc.s.ll.Rdata")
save(disc.s.lf, file="disc.s.lf.Rdata")
save(disc.s.ff, file="disc.s.ff.Rdata")

use.samples <- min(loops.l$nrSamples, loops.s$nrSamples)

## light-strange (strange-light is the same...)
disc.ls.ll <- computeDisc(cf = loops.l, cf2 = loops.s, smeared=FALSE, smeared2=FALSE, real=FALSE, real2=FALSE,
                          subtract.vev=FALSE, subtract.vev2=FALSE, subtract.equal=FALSE)
disc.ls.lf <- computeDisc(cf = loops.l, cf2 = loops.s, smeared=FALSE, smeared2=TRUE, real=FALSE, real2=FALSE,
                          subtract.vev=FALSE, subtract.vev2=FALSE, subtract.equal=FALSE)
disc.ls.fl <- computeDisc(cf = loops.l, cf2 = loops.s, smeared=TRUE, smeared2=FALSE, real=FALSE, real2=FALSE,
                          subtract.vev=FALSE, subtract.vev2=FALSE, subtract.equal=FALSE)
disc.ls.ff <- computeDisc(cf = loops.l, cf2 = loops.s, smeared=TRUE, smeared2=TRUE, real=FALSE, real2=FALSE,
                          subtract.vev=FALSE, subtract.vev2=FALSE, subtract.equal=FALSE)

save(disc.ls.ll, file="disc.ls.ll.Rdata")
save(disc.ls.lf, file="disc.ls.lf.Rdata")
save(disc.ls.fl, file="disc.ls.fl.Rdata")
save(disc.ls.ff, file="disc.ls.ff.Rdata")

eta.l.cor <- add.cf(con.l, c(disc.l.ll, disc.l.lf, disc.l.lf, disc.l.ff), a=-0.5, b=-2.)
eta.s.cor <- add.cf(con.s, c(disc.s.ll, disc.s.lf, disc.s.lf, disc.s.ff), a=-0.25, b=-0.5)

eta.ls.cor <- mul.cf(c(disc.ls.ll, disc.ls.lf, disc.ls.fl, disc.ls.ff), a=-1.)
eta.sl.cor <- mul.cf(c(disc.ls.ll, disc.ls.fl, disc.ls.lf, disc.ls.ff), a=-1.)

Cmatrix <- c(eta.l.cor, eta.ls.cor, eta.sl.cor, eta.s.cor)
Cmatrix <- bootstrap.cf(Cmatrix, boot.R=999, boot.l=10)
save(Cmatrix, file="Cmatrix.Rdata")

## we use element.order to bring the matrix into the right order
Cmatrix.bootstrap.gevp <- bootstrap.gevp(Cmatrix, matrix.size=4, t0=1,
                                         element.order=c(
                                            1, 2, 5, 6,   
                                            3, 4, 7, 8,  
                                            9,11,13,14,
                                           10,12,15,16 ) )

## solve the GEVP
## eta
eta <- gevp2cf(Cmatrix.bootstrap.gevp, id=1)
eta.effectivemass <- bootstrap.effectivemass(cf=eta, type="solve")
eta.effectivemass <- fit.effectivemass(eta.effectivemass, t1=10, t2=23, useCov=FALSE)
plot(eta.effectivemass, ylim=c(0.1,0.3))
summary(eta.effectivemass)

## eta prime
etap <- gevp2cf(Cmatrix.bootstrap.gevp, id=2)
etap.effectivemass <- bootstrap.effectivemass(cf=etap, type="solve")
etap.effectivemass <- fit.effectivemass(etap.effectivemass, t1=10, t2=23, useCov=FALSE)
X11()
plot(etap.effectivemass, ylim=c(0.1,0.6))
summary(etap.effectivemass)


## now with excited state subtraction
## need to carefully chose the fit-range to get good results in the subtraction
con.l.ll <- extractSingleCor.cf(con.l, id=1)
con.l.ll <- mul.cf(con.l.ll, -1)
con.l.ll <- bootstrap.cf(con.l.ll, boot.R=999, boot.l=10)
con.l.ll.matrixfit <- matrixfit(con.l.ll, t1=15, t2=32, symmetrise=TRUE, useCov=FALSE)
summary(con.l.ll.matrixfit)

con.s.ll <- extractSingleCor.cf(con.s, id=1)
con.s.ll <- mul.cf(con.s.ll, -1)
con.s.ll <- bootstrap.cf(con.s.ll, boot.R=999, boot.l=10)
con.s.ll.matrixfit <- matrixfit(con.s.ll, t1=15, t2=32, symmetrise=TRUE, useCov=FALSE)
summary(con.s.ll.matrixfit)

con.l.subtracted <- subtract.excitedstates(con.l.ll, con.l.ll.matrixfit, from.samples=FALSE)
con.s.subtracted <- subtract.excitedstates(con.s.ll, con.s.ll.matrixfit, from.samples=FALSE)

## note the different sign of a because of the mul.cf with -1 above
eta.l.cor.subtracted <- add.cf(con.l.subtracted, disc.l.ll, a=0.5, b=-2.)
eta.s.cor.subtracted <- add.cf(con.s.subtracted, disc.s.ll, a=0.25, b=-0.5)

Cmatrix.subtracted <- c(eta.l.cor.subtracted, mul.cf(disc.ls.ll, -1), mul.cf(disc.ls.ll, -1), eta.s.cor.subtracted)
Cmatrix.subtracted <- bootstrap.cf(Cmatrix.subtracted, boot.R=999, boot.l=10)

Cmatrix.subtracted.bootstrap.gevp <- bootstrap.gevp(Cmatrix.subtracted, matrix.size=2, t0=1,
                                                    element.order=c(1,2,3,4))

eta.sub <- gevp2cf(Cmatrix.subtracted.bootstrap.gevp, id=1)
eta.sub.effectivemass <- bootstrap.effectivemass(cf=eta.sub, type="solve")
eta.sub.effectivemass <- fit.effectivemass(eta.sub.effectivemass, t1=2, t2=8, useCov=FALSE)
plot(eta.sub.effectivemass, ylim=c(0.1,0.3))
summary(eta.sub.effectivemass)

etap.sub <- gevp2cf(Cmatrix.subtracted.bootstrap.gevp, id=2)
etap.sub.effectivemass <- bootstrap.effectivemass(cf=etap.sub, type="solve")
etap.sub.effectivemass <- fit.effectivemass(etap.sub.effectivemass, t1=2, t2=8, useCov=FALSE)
plot(etap.sub.effectivemass, ylim=c(0.1,0.6))
summary(etap.sub.effectivemass)


