## if the data is concatenated into a single file "pion.dat" already, use this
cmicor <- readcmicor("pion.dat")
## otherwise you need this here for reading the single files
### files <- getorderedfilelist(basename="outprcv.")
### cmicor <- readcmidatafiles(files)

## now extract the gamma matrix combination of interest
pion.cor <- extract.obs(cmicor,  vec.obs=c(1))
## which will extract gamma5 from the file for all smearings available
pion.cor <- bootstrap.cf(pion.cor, boot.R=400, boot.l=1)

## now we can attempt a constrained fit to the matrix
pion.cor.matrixfit <- matrixfit(pion.cor, t1=10, t2=23, symmetrise=TRUE, useCov=FALSE)
## compute the ps decay constant for the twisted mass case, need mu and kappa
pion.cor.matrixfit <- computefps(pion.cor.matrixfit, mu1=0.003, kappa=0.13782)
X11()
plot(pion.cor.matrixfit, xlab=c("t/a"), ylab=c("C(t)"))
summary(pion.cor.matrixfit)

## or an effective mass analysis of the correlators
pion.cor.effectivemass <- bootstrap.effectivemass(pion.cor, type="acosh")
pion.cor.effectivemass <- fit.effectivemass(pion.cor.effectivemass, t1=10, t2=23, useCov=TRUE)
X11()
plot(pion.cor.effectivemass, xlab=c("t/a"), ylab=c("aM"))
summary(pion.cor.effectivemass)

## apply a GEVP analysis
pion.cor.gevp <- bootstrap.gevp(pion.cor, t0=1)
## extract the first principal correlator
pion.pc1 <- gevp2cf(pion.cor.gevp, id=1)
## which can now be treated like a bootstrapped correlation function
pion.pc1.effectivemass <- bootstrap.effectivemass(cf=pion.pc1, type="acosh")
pion.pc1.effectivemass <- fit.effectivemass(pion.pc1.effectivemass, t1=10, t2=23, useCov=TRUE)
summary(pion.pc1.effectivemass)
X11()
plot(pion.pc1.effectivemass, xlab=c("t/a"), ylab=c("aM"))
