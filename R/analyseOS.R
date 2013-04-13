analyse.os <- function(N, kappa, t1, t2, boot.R, boot.l,
                       muvalues, ename, starti=0) {
  if(length(t1) == 1) {
    t1 <- rep(t1[1], times=((N+1)*(N+2)/2))
  }
  else if(length(t1) != ((N+1)*(N+2)/2)) {
    stop("Error! t1 must be of lenght 1 or (N+1)*(N+2)/2!")
  }
  if(length(t2) == 1) {
    t2 <- rep(t2[1], times=((N+1)*(N+2)/2))
  }
  else if(length(t2) != ((N+1)*(N+2)/2)) {
    stop("Error! t2 must be of lenght 1 or (N+1)*(N+2)/2!")
  }
  
  result <- data.frame(mu1 = numeric((N+1)*(N+2)/2), mu2 = numeric((N+1)*(N+2)/2),
                       m = numeric((N+1)*(N+2)/2), dm = numeric((N+1)*(N+2)/2),
                       f = numeric((N+1)*(N+2)/2), df = numeric((N+1)*(N+2)/2),
                       fsinh = numeric((N+1)*(N+2)/2), L = numeric((N+1)*(N+2)/2), 
                       chisqr = numeric((N+1)*(N+2)/2), dof = numeric((N+1)*(N+2)/2),
                       t1 = t1, t2 = t2)
  
  c <- 1
  bootsamples <- array(0., dim=c(boot.R, 2, (N+1)*(N+2)/2))

  for (i in starti:N) {
    for (j in i:N) {
      mu1 <- muvalues$V1[i+1]
      mu2 <- muvalues$V1[j+1]
      filename=paste("ppcorrel.",sprintf("%.02d", i),".",sprintf("%.02d", j),".dat",sep="")
      cat(i, j, filename, "\n")
      cmicor <- read.table(filename, header=F, colClasse=c("integer","integer","integer","numeric","numeric"));
      res <- smearedpion(cmicor, t1=t1[c], t2=t2[c], debug=FALSE, method="all", mu1=mu1, mu2=mu2, kappa=kappa, boot.R=boot.R, boot.l=boot.l)
      result$L <- max(cmicor[3])
      result$mu1[c] <- mu1
      result$mu2[c] <- mu2
      result$m[c] <- abs(res$fitresult$par[3])
      result$dm[c] <- res$uwerrresultmps$dvalue
      result$f[c] <- 2*kappa*2*(mu1+mu2)/2/sqrt(2)*abs(res$fitresult$par[1])/sqrt(abs(result$m[c])^3)
      result$fsinh[c] <- 2*kappa*2*(mu1+mu2)/2/sqrt(2)*abs(res$fitresult$par[1])/sqrt(abs(result$m[c]))/sinh(abs(result$m[c]))
      result$df[c] <- res$uwerrresultfps$dvalue*2*kappa*2*(mu1+mu2)/2./sqrt(2)
      result$chisqr[c] <- res$fitresult$value
      result$dof[c] <- res$dof
      bootsamples[,1,c] <- res$tsboot$t[,1]
      bootsamples[,2,c] <- res$tsboot$t[,2]
      c <- c+1
      save(result, file=paste("result", ename, ".Rdata", sep=""))
      save(bootsamples, file=paste("bootsamples", ename, ".Rdata", sep=""))
    }
  }
  write.table(result, file=paste("result", ename, ".dat", sep=""), quote=FALSE, sep="\t")
}
