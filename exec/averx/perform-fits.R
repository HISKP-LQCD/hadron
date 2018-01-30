boot.R <- 1500
boot.l <- 1
seed <- 123476

useCov <- TRUE
ens <- c("cA2.09.48")

tlower <- seq(4,20,2)
tupper <- seq(28,44,2)
piontlower <- seq(8,20,2)
piontupper <- c(44,46,48)
redo <- FALSE

if(file.exists("parameters.R")) {
  source("parameters.R")
}

filename <- paste("averx-data-", ens, ".Rdata", sep="")
if(!file.exists(filename)) {
  data2pt <- bootstrap.cf(convert2cf(read.table("ppcorlikepoint.dat")), boot.R=boot.R, boot.l=boot.l, seed=seed)
  data3pt <- bootstrap.cf(mul.cf(convert2cf(read.table("momf2elikepoint.dat")), -1.), boot.R=boot.R, boot.l=boot.l, seed=seed)
  
  save(data2pt, data3pt, file=filename)
}
load(filename)
for(piont1 in piontlower) {
  for(piont2 in piontupper) {  
    pionfit <- matrixfit(data2pt, t1=piont1, t2=piont2, symmetrise=TRUE, useCov=useCov,
                         matrix.size=1, parlist=array(c(1,1), dim=c(2,1)))
    for(t1 in tlower) {
      for(t2 in tupper) {
        filename <- paste("res.averx.t1", t1, ".t2", t2, ".piont1", piont1, ".piont2", piont2, ".", ens, ".Rdata", sep="")
        if(!file.exists(filename) && !redo) {
          res.averx <- averx(data3pt, data2pt, pionfit, boot.R=boot.R, boot.l=boot.l, piont1=piont1, piont2=piont2, t1=t1, t2=t2, useCov=useCov)
          
          summary(res.averx)
          plot(res.averx)
          
          
          save(res.averx, file=filename)
        }
      }
    }
  }
}

