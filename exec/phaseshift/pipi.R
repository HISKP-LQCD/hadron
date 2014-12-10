## matrix size: N
## tp: reference frame
energies.pipi <- function(N=5, tp="TP0", basename="pipi_pipi_A1_corr_", redofit=TRUE, N.ids=3,
                          t1 = c(10, 8, 5), t2 = c(31, 15, 13), seed=12345,
                          boot.R=400, boot.l=1, T=64,
                          srcpath="./", path="./") {
  
  
  ## read data into Cmatrix
  if(!file.exists(paste(path, "Cmatrix.", tp, ".Rdata", sep=""))) {
    Cmatrix <- cf()
    
    for(i in c(1:N)) {
      for(j in c(1:N)) {
        filename <- paste(basename, tp, "_", i-1, j-1, ".dat", sep="")
        tmp <- readtextcf(filename, T=T, check.t=1, path=srcpath)
        Cmatrix <- c(Cmatrix, tmp)
      }
    }
    
    ## we bootstrap the matrix and save
    Cmatrix <- bootstrap.cf(Cmatrix, boot.R=boot.R, boot.l=boot.l, seed=seed)
    save(Cmatrix, file=paste(path,"Cmatrix.", tp, ".Rdata", sep=""))
  }
  load(file=paste(path,"Cmatrix.", tp, ".Rdata", sep=""))
  
  if(redofit) {
    cat("...solving the GEVP\n")
    ## we use element.order to bring the matrix into the right order
    Cmatrix.bootstrap.gevp <- bootstrap.gevp(Cmatrix, matrix.size=N,
                                             element.order=c(1:(N^2)))

    save(Cmatrix.bootstrap.gevp, file=paste(path, "Cmatrix.bootstrap.gevp.", tp, ".Rdata", sep=""))
    ## extract the principal correlators
    ## determine masses from effective mass plus temporal state
    cat("...determining masses for pi pi states\n")
    for(i in c(1:N.ids)) {
      pc <- gevp2cf(Cmatrix.bootstrap.gevp, id=i)
      for(ta in c(t1[i]:t2[i])) {
        if(t2[i] - ta < 5) break
        pc.effectivemass <- fit.effectivemass(bootstrap.effectivemass(cf=pc, type="temporal"), t1=ta, t2=t2[i], useCov=TRUE)
        summary(pc.effectivemass)
        ## lets save all of them
        save(pc.effectivemass, file=paste(path, "pc", i, ".", tp, ".effectivemass.", ta, ".", t2[i], ".Rdata", sep=""))
      }
    }
  }
  return(invisible(Cmatrix))
}
