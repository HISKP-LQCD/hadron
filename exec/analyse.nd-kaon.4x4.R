
## we analyse 8x8 matrix for the kaon
## we need to read four files corresponding to heavy propagators as
file.strings <- c("ss", "sc", "cs", "cc")

if(FALSE) {
  for(i in c(1:4)) {
    files <- getorderedfilelist(basename=paste("outprcv.l", file.strings[i], ".", sep=""))
    ## set skip to 1 for libcvcpp
    assign(file.strings[i], readcmidatafiles(files, skip=1, verbose=TRUE))
  }
}
load("ss.Rdata")
load("sc.Rdata")
load("cs.Rdata")
load("cc.Rdata")

flavour.strings <- array(c("cc", "cs", "cc", "cs", 
                           "sc", "ss", "sc", "ss", 
                           "cc", "cs", "cc", "cs", 
                           "sc", "ss", "sc", "ss"),
                         dim=c(4,4))


## the following have to be chosen for the heavyheavy code
flavour.factors.hh <- 0.5*array(c(-1, 1, 1, 1,
                                  1, -1, -1, -1,
                                  1, -1, 1, 1,
                                  1, -1, 1, 1), dim=c(4,4))

gamma.indices.hh <- array(c(4, 4, 3, 3,
                            4, 4, 3, 3,
                            2, 2, 1, 1,
                            2, 2, 1, 1),
                          dim=c(4,4))

## and this here for libcvcpp (might change in the future...)
flavour.factors <- 0.5*array(c(-1, 1, -1, -1,
                               1, -1, 1, 1,
                               -1, 1, 1, 1,
                               -1, 1, 1, 1), dim=c(4,4))

gamma.indices <- array(c(5, 5, 7, 7,
                         5, 5, 7, 7,
                         6, 6, 1, 1,
                         6, 6, 1, 1),
                       dim=c(4,4))


Cmatrix <- cf()

for(i in c(3:4)) {
  for(j in c(3:4)) {
    tmp <- extract.obs(eval(as.name(flavour.strings[i,j])), vec.obs=c(gamma.indices[i,j]))
    Cmatrix <- c(Cmatrix, mul.cf(tmp, a=flavour.factors[i,j]))
  }
}

## we bootstrap the matrix and save
Cmatrix <- bootstrap.cf(Cmatrix, boot.R=400, boot.l=2)
save(Cmatrix, file="Cmatrix.Rdata")

##load("Cmatrix.Rdata")

## we use element.order to bring the matrix into the right order
##Cmatrix.bootstrap.gevp <- bootstrap.gevp(Cmatrix, matrix.size=8,
##                                         element.order=c(
##                                           1, 2, 5, 6,   9,10,13,14,
##                                           3, 4, 7, 8,  11,12,15,16,
##                                           17,18,21,22, 25,26,29,30,
##                                           19,20,23,24, 27,28,31,32,
##                                           33,35,49,51, 41,42,45,46,
##                                           34,36,50,52, 43,44,47,48,
##                                           37,39,53,55, 57,58,61,62,
##                                           38,40,54,56, 59,60,63,64))

Cmatrix.bootstrap.gevp <- bootstrap.gevp(Cmatrix, matrix.size=4,
                                         element.order=c(
                                           1,2,5,6,
                                           3,4,7,8,
                                           9,11,13,14,
                                           10,12,15,16
                                           ))

## solve the GEVP
etass.pc1 <- gevp2cf(Cmatrix.bootstrap.gevp, id=1)
etass.pc1.effectivemass <- bootstrap.effectivemass(cf=etass.pc1, type="acosh")
etass.pc1.effectivemass <- fit.effectivemass(etass.pc1.effectivemass, t1=12, t2=23, useCov=TRUE)
plot(etass.pc1.effectivemass, ylim=c(0.1,0.3))
summary(etass.pc1.effectivemass)
