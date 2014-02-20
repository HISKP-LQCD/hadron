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
plot(pion.cor.matrixfit, xlab=c("t/a"), ylab=c("C(t)"))
summary(pion.cor.matrixfit)

## or an effective mass analysis of the correlators
pion.cor.effectivemass <- bootstrap.effectivemass(pion.cor, boot.R=400, boot.l=1, type="acosh")
pion.cor.effectivemass <- fit.effectivemass(pion.cor.effectivemass, t1=10, t2=23, useCov=TRUE)
X11()
plot(pion.cor.effectivemass, xlab=c("t/a"), ylab=c("aM"))
summary(pion.cor.effectivemass)

