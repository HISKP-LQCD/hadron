## matrix size: N
## tp: reference frame
energies.pipi <- function(N=5, tp="TP0", irrep="A1", basename="pipi_pipi_", redofit=TRUE, N.ids=3,
                          t1 = c(10, 8, 5), t2 = c(31, 15, 13), seed=12345, t0=1,
                          boot.R=400, boot.l=1, T=64,
                          srcpath="./", path="./", ens, ind.vector=c(2,3)) {
  
  
  ## read data into Cmatrix
  if(!file.exists(paste(path, "Cmatrix.", irrep, ".", tp, ".Rdata", sep=""))) {
    Cmatrix <- cf()
    
    for(i in c(1:N)) {
      for(j in c(1:N)) {
        filename <- paste(basename, irrep, "_corr_", tp, "_", i-1, j-1, ".dat", sep="")
        tmp <- readtextcf(filename, T=T, check.t=1, path=srcpath, ind.vector=ind.vector)
        Cmatrix <- c(Cmatrix, tmp)
      }
    }
    
    ## we bootstrap the matrix and save
    Cmatrix <- bootstrap.cf(Cmatrix, boot.R=boot.R, boot.l=boot.l, seed=seed)
    save(Cmatrix, file=paste(path,"Cmatrix.", irrep, ".", tp, ".Rdata", sep=""))
  }
  load(file=paste(path,"Cmatrix.", irrep, ".", tp, ".Rdata", sep=""))
  
  if(redofit) {
    cat("...solving the GEVP\n")
    ## we use element.order to bring the matrix into the right order
    Cmatrix.bootstrap.gevp <- bootstrap.gevp(Cmatrix, matrix.size=N, t0=t0,
                                             element.order=c(1:(N^2)))

    save(Cmatrix.bootstrap.gevp, file=paste(path, "Cmatrix.bootstrap.gevp.", irrep, ".", tp, ".Rdata", sep=""))
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
        meta.data <- list(ens=ens, T=T, srcpath=srcpath, t1=ta, t2=t2[i], tp=tp, pcid=i, irrep=irrep)
        save(pc.effectivemass, meta.data, file=paste(path, "pc", i, ".", irrep, ".", tp, ".effectivemass.", ta, ".", t2[i], ".Rdata", sep=""))

	plot(pc.effectivemass, xlab=c("t/a"), ylab=c("Meff"), main=paste("pc", i, ".", irrep, ".", tp, ".effectivemass.", ta, ".", t2[i], "p=", pc.effectivemass$Qval, sep=""))
      }
    }
  }
  return(invisible(Cmatrix))
}

