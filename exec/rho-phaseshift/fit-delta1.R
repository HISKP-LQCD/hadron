clargs = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(clargs)!=1) {
    stop("Please specify the input file!", call.=FALSE)
} else  {
    # default output file
    input.file = clargs[1]
}

source(input.file)

ensembles <- args$ens
#ensembles <- c("A40.24")
frames <- c("cmf", "mf1", "mf2", "mf3")
momenta <- c("p0", "p1", "p2", "p3")
dirs <- args$dirs
path.to.data <- args$path.to.data
output.path <- args$output.path

setwd(output.path)

## Starting values for fit parameters
## g_{\rho\pi\pi} and M_\rho
par <- c(6,0.4)

labels <- apply( cbind(unlist(lapply(strsplit(dirs, "/"), function(x){x[1]})), unlist(lapply(strsplit(dirs, "/"), function(x){x[2]}))), 1, paste, collapse=" ")

load( paste(args$path.to.data, "/pion.Rdata", sep="") )

boot.R <- pion.cor$boot.R

MKdata <- read.table( paste(args$output.path, "mk.dat", sep="") )

Mk <- MKdata$V2[which(MKdata$V1 == ensembles[1])]
Mpi <- pion.matrixfit$opt.res$par[1]
Mpiboot <- c(Mpi, pion.matrixfit$opt.tsboot[1,])

cat("Mk", Mk, "\n")
Ndirs <- length(dirs)

alldata <- array(NA, dim=c(boot.R+1,4*Ndirs))

rr <- c(1:boot.R)+1

for(i in c(1:Ndirs)) {
  splitwd <- strsplit(dirs[i], "/")[[1]]
  momentum <- splitwd[1]
  frameid <- which(momenta == momentum)
  frame <- frames[frameid]
  Ensemble <- ensembles[1]
  for(j in c(1:length(ensembles))){
    cat(length(ensembles), ensembles, ensembles[j], dirs[i], "\n")
    if(grepl(ensembles[j], dirs[i])) {
      Ensemble <- ensembles[j]
      break
    }
  }
  
  file <- paste(path.to.data, dirs[i], "/res.", Ensemble, ".", frame, ".", splitwd[2], ".Rdata", sep="")

  cat("loading file", file, "\n")
  load(file)
  j <- 2*(i-1)+1
  ## write first the Ecm values and then delta
  ## alldata[,c(1:N)] -> Ecm
  ## alldata[,c((N+1):(2*N))] -> delta
  n <- length(res.boot)

  alldata[,c(j, 2*Ndirs+j)] <- res.boot[[1]][,c(1:2)]
  alldata[,c(j+1, 2*Ndirs+j+1)] <- res.boot[[2]][,c(1:2)]
}

ii <- which(alldata[1,c(1:(2*Ndirs))] < 2*Mk)
data <- alldata[,c(ii, ii+2*Ndirs)]
N <- length(ii)

M <- invertCovMatrix(data, boot.samples=TRUE)

L <- chol(M)

par <- c(par, data[1,c(1:length(ii))])

deltaovEcm <- function(par, Ecm, Mpi) {
  x <- atan(par[1]^2/6/pi*sqrt(Ecm^2/4.-Mpi^2)^3/Ecm/(par[2]^2-Ecm^2))
  x[which(x < 0)] <- x[which(x < 0)] + pi
  return(x)
}


Chi <- function(par, y, L, Mpi) {
  N <- (length(par)-2)
  z <- deltaovEcm(par[c(1:2)], Ecm=par[c(1:N)+2], Mpi=Mpi)
  i <- which(z < 0)
  ##print(L)
  ##cat(z, "\n")
  ##cat(par[c(1:N)+2], "\n")
  ##cat (c(par[c(1:N)+2]-y[c(1:N)], z-y[c(1:N)+N]), "\n")
  ##cat (dim(L), "\n")
  if(any(is.na(z))) return(rep(10000, times=2*N))
  return(L %*% c(par[c(1:N)+2]-y[c(1:N)], z-y[c(1:N)+N]))
}

lm.avail <- require(minpack.lm)


npar <- length(par)
Mrho.res <- array(NA, dim=c(boot.R+1, npar+1))

##if(FALSE) {
filename <- paste("Mrho-res", ensembles[1], ".", ensembles[2], ".Rdata", sep="")
if(file.exists(filename) && file_test("-nt", filename, "parameters-A40.R")) {
  load(filename)
}
if(!(file.exists(filename) && file_test("-nt", filename, "parameters-A40.R"))) {
  for(i in c(1:(boot.R+1))) {
    opt.res <- nls.lm(par, fn=Chi, y=data[i,], L=L, Mpi=Mpiboot[i], control = nls.lm.control(ftol=1.e-8, ptol=1.e-8, maxiter=500))
    Mrho.res[i, c(1:npar)] <- opt.res$par
    Mrho.res[i, npar+1] <- opt.res$rsstrace[length(opt.res$rsstrace)]
  }
  
  save(Mrho.res, data, file=filename)
}

dof <- npar-4
cat("\n************************************\n")
cat(ensembles, "\n\n")
cat("Mrho", Mrho.res[1,2], "+-", sd(Mrho.res[,2]), "\n")
cat("g", Mrho.res[1,1], "+-", sd(Mrho.res[,1]), "\n")
cat("chi^2", Mrho.res[1,dim(Mrho.res)[2]], "dof", dof, "\n")
cat("\n************************************\n")

x <- seq(2*Mpi, 2*Mk, 0.005)
y <- deltaovEcm(par=Mrho.res[1, c(1:npar)], x, Mpi)
y[which(y<0)] <- y[which(y<0)] + pi

dy <- apply(apply(Mrho.res[, c(1:npar)], 1, deltaovEcm, x, Mpi), 1, sd)
##}


tikzfiles <- tikz.init(basename="fit-delta1", width=4.5, height=5.)
##output.path, 
n <- 2*length(dirs)
cols <- rep(c("red", "blue", "darkgreen", "magenta", "brown", "black", "cyan", "green", "orange"), times=3)
# must match the number of plotted irreps for different volumes to have the same color

pchs <- rep(c(21:25,21:25,21:25), times=3)
plot(NA, xlim=c(2, 4), ylim=c(-0.4,pi+.4), xlab=c("$E_\\mathrm{CM}/M_\\pi$"), ylab=c("$\\delta_1$"))
polygon(x=c(x,rev(x))/Mpi, y=c(y+dy, rev(y-dy)), col="gray", lty=0, lwd=0.001, border="gray")
for(i in seq(1,n-1,2)) {
  plotwitherror(x=alldata[1,c(i,i+1)]/Mpi, y=alldata[1,c(i,i+1)+n], dx=apply(alldata[,c(i,i+1)], 2, sd)/Mpi, dy=apply(alldata[,c(i,i+1)+n], 2, sd),
                pch=pchs[i/2+1], col=cols[i/2+1], rep=TRUE)
}
lines(x/Mpi,y)
abline(v=2*Mk/Mpi, lty=c(2))
text(x=2*Mk/Mpi+0.25, y=0.4, "$2M_K/M_\\pi$")
legend("topleft", legend=labels, bty="n", col=cols[1:(n/2)], pch=pchs[1:(n/2)])
tikz.finalize(tikzfiles=tikzfiles)

if(FALSE) {
  tikzfiles <- tikz.init(basename=paste("delta-fitqq", ensembles[1], sep=""), width=5., height=4.)
  s <- seq(0,1,1./boot.R)
  x <- qchisq(p=s, df=dof)
  qqplot(x=x, y=Mrho.res[,dim(Mrho.res)[2]], xlab="Theoretical Quantiles", ylab="Sample Quantiles", main="QQ-Plot chisqr Values")
  
  
  tikz.finalize(tikzfiles=tikzfiles)
}
