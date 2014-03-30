
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
flavour.factors <- 0.5*array(c(-1, 1, 1, 1,
                                  1, -1, -1, -1,
                                  1, -1, 1, 1,
                                  1, -1, 1, 1), dim=c(4,4))

gamma.indices.hh <- array(c(4, 4, 3, 3,
                            4, 4, 3, 3,
                            2, 2, 1, 1,
                            2, 2, 1, 1),
                          dim=c(4,4))

## and this here for libcvcpp (might change in the future...)
flavour.factors.cvc <- 0.5*array(c(-1, 1, -1, -1,
                               1, -1, 1, 1,
                               -1, 1, 1, 1,
                               -1, 1, 1, 1), dim=c(4,4))

gamma.indices <- array(c(5, 5, 7, 7,
                         5, 5, 7, 7,
                         6, 6, 1, 1,
                         6, 6, 1, 1),
                       dim=c(4,4))

if(!file.exists("Cmatrix.Rdata")) {
  Cmatrix <- cf()
  
  for(i in c(1:4)) {
    for(j in c(1:4)) {
      tmp <- extract.obs(eval(as.name(flavour.strings[j,i])), vec.obs=c(gamma.indices[j,i]))
      Cmatrix <- c(Cmatrix, mul.cf(tmp, a=flavour.factors[j,i]))
    }
  }

  ## we bootstrap the matrix and save
  Cmatrix <- bootstrap.cf(Cmatrix, boot.R=400, boot.l=2)
  save(Cmatrix, file="Cmatrix.Rdata")
}
load("Cmatrix.Rdata")

## we use element.order to bring the matrix into the right order
Cmatrix.bootstrap.gevp <- bootstrap.gevp(Cmatrix, matrix.size=8,
                                         element.order=c(
                                           1, 2, 5, 6,   9,10,13,14,
                                           3, 4, 7, 8,  11,12,15,16,
                                           17,19,21,22, 25,26,29,30,
                                           18,20,23,24, 27,28,31,32,
                                           33,35,37,39, 41,42,45,46,
                                           34,36,38,40, 43,44,47,48,
                                           49,51,53,55, 57,59,61,62,
                                           50,52,54,56, 58,60,63,64))

## solve the GEVP
kaon.pc1 <- gevp2cf(Cmatrix.bootstrap.gevp, id=1)
kaon.pc1.effectivemass <- bootstrap.effectivemass(cf=kaon.pc1, type="solve")
kaon.pc1.effectivemass <- fit.effectivemass(kaon.pc1.effectivemass, t1=10, t2=23, useCov=FALSE)
plot(kaon.pc1.effectivemass, ylim=c(0.1,0.3))
summary(kaon.pc1.effectivemass)
