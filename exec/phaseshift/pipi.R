## start with zero total momentum
## and analyse the matrix
## matrix size
N <- 5
## read data into Cmatrix
if(!file.exists(paste(path, "Cmatrix.Rdata", sep=""))) {
  Cmatrix <- cf()
  
  for(i in c(1:N)) {
    for(j in c(1:N)) {
      filename <- paste(srcpath, "pipi_pipi_A1_corr_TP0_", i-1, j-1, ".dat", sep="")
      tmp <- readtextcf(filename, T=T, check.t=1, path=path)
      Cmatrix <- c(Cmatrix, tmp)
    }
  }
  
  ## we bootstrap the matrix and save
  Cmatrix <- bootstrap.cf(Cmatrix, boot.R=boot.R, boot.l=boot.l, seed=seed)
  save(Cmatrix, file=paste(path,"Cmatrix.Rdata", sep=""))
}
load(file=paste(path,"Cmatrix.Rdata", sep=""))

if(redofit) {
  cat("...solving the GEVP\n")
  ## we use element.order to bring the matrix into the right order
  Cmatrix.bootstrap.gevp <- bootstrap.gevp(Cmatrix, matrix.size=N,
                                           element.order=c(1:(N^2)))
  
  ## extract the principal correlators
  ## determine masses from effective mass plus temporal state
  cat("...determining masses for pi pi states\n")
  t1 <- c(10, 8, 5)
  t2 <- c(31, 15, 13)
  for(i in c(1:3)) {
    pc <- gevp2cf(Cmatrix.bootstrap.gevp, id=i)
    for(ta in c(t1[i]:t2[i])) {
      if(t2[i] - ta <= 5) break
      pc.effectivemass <- fit.effectivemass(bootstrap.effectivemass(cf=pc, type="temporal"), t1=ta, t2=t2[i], useCov=TRUE)
      summary(pc.effectivemass)
      if(pc.effectivemass$Qval > 0.01 && pc.effectivemass$Qval < 0.99) {
        save(pc.effectivemass, file=paste(path, "pc", i, ".effectivemass.", ta, ".", t2[i], ".Rdata", sep=""))
      }
    }
  }
}
