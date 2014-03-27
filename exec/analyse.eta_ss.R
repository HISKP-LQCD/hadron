## we analyse the connected only heavy-heavy contributions to
## the eta correlation matrix
## we need the following flavour combinations

flavour.strings <- array(c("cscs","cssc","sccs","scsc", "cscc","csss","sccc","scss",
                           "cccs","ccsc","sscs","sssc", "cccc","ccss","sscc","ssss"),
                         dim=c(4,2,2))

elements.strings <- array(c("ss", "sp", "ps", "pp"), dim=c(2,2))

## mapping needed to create the filenames of the input files
## mapping is needed because of gamma_5 trick etc.
flavour.mapping <- function(s) {
  return(paste(substr(s,1,1), substr(s,4,4), substr(s,2,2), substr(s,3,3), sep=""))
}

flavour.factors <- array(c(+.25,-.25,-.25,+.25, +.25,+.25,-.25,-.25,
                           -.25,+.25,-.25,+.25, +.25,+.25,+.25,+.25),
                         dim=c(4,2,2))
## the following have to be chosen for the heavyheavy code
gamma.indices <- array(c(4,4,4,4, 3,3,3,3,
                         2,2,2,2, 1,1,1,1),
                       dim=c(4,2,2))
## and this here for libcvcpp (might change in the future...)
##gamma.indices <- array(c(5,5,5,5, 7,7,7,7,
##                         6,6,6,6, 1,1,1,1),
##                       dim=c(4,2,2))


# set reread = TRUE when you want to read the data again
reread <- FALSE
if(!file.exists("Cmatrix.Rdata") || reread) {
  for(i in c(1:2)) {
    for(j in c(1:2)) {
      for(k in c(1:4)) {
        files <- getorderedfilelist(basename=paste("outprcvn.", flavour.mapping(flavour.strings[k,i,j]), ".", sep=""))
        ## set skip to 1 for libcvcpp
        cmicor <- readcmidatafiles(files, skip=0, verbose=TRUE)
        assign(flavour.strings[k,i,j], extract.obs(cmicor,  vec.obs=c(gamma.indices[k,i,j])))
        if(k == 1) {
          assign("tmp", eval(as.name(flavour.strings[k,i,j])))
          tmp <- mul.cf(tmp, flavour.factors[k,i,j])
        }
        else {
          tmp <- add.cf(tmp, eval(as.name(flavour.strings[k,i,j])), a=1., b=flavour.factors[k,i,j])
        }
      }
      assign(elements.strings[i,j], tmp)
      rm(tmp)
    }
  }

  ## now we coerce to obtain the full matrix
  ## note that here we have smearing as fastest index
  Cmatrix <- c(eval(as.name(elements.strings[1,1])), eval(as.name(elements.strings[1,2])),
               eval(as.name(elements.strings[2,1])), eval(as.name(elements.strings[2,2])))
  
  ## we bootstrap the matrix and save
  Cmatrix <- bootstrap.cf(Cmatrix, boot.R=400, boot.l=2)
  save(Cmatrix, file="Cmatrix.Rdata")
}
load("Cmatrix.Rdata")

## we use element.order to bring the matrix into the right order
Cmatrix.bootstrap.gevp <- bootstrap.gevp(Cmatrix, matrix.size=4,
                                         element.order=c(
                                           1,2,5,6,
                                           3,4,7,8,
                                           9,11,13,14,
                                           10,12,15,16))

## solve the GEVP
etass.pc1 <- gevp2cf(Cmatrix.bootstrap.gevp, id=1)
etass.pc1.effectivemass <- bootstrap.effectivemass(cf=etass.pc1, type="acosh")
etass.pc1.effectivemass <- fit.effectivemass(etass.pc1.effectivemass, t1=12, t2=23, useCov=TRUE)
plot(etass.pc1.effectivemass, ylim=c(0.2,0.4))
summary(etass.pc1.effectivemass)
