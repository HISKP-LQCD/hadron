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
all.dirs <- args$all.dirs
max.pcs <- args$pcs
data.paths <- args$data.paths
output.path <- paste(args$output.path, "6_phaseshift", sep="/")

# Paramters for fit layout
# must match the number of plotted irreps for different volumes to have the same color
pchs <- c(21:25) # vector with filled plot symbols
cols <- c("red", "blue", "darkgreen", "magenta", "brown", "black", "cyan", "green", "orange")
unique.dirs <- unique(unlist(all.dirs))
cols <- array(cols[1:length(unique.dirs)], dimnames=list(unique.dirs))

dir.create(output.path, showWarnings=FALSE, recursive=TRUE)
setwd(output.path)
## WARNING: Name of sink may not have the same basename as tikz uses below
sink("fitresults.log", append=FALSE, split=TRUE)

## Starting values for fit parameters
## g_{\rho\pi\pi} and M_\rho
par <- c(6,0.4)

load( paste(data.paths[length(data.paths)], "/pion.Rdata", sep="") )

## TODO: That is dangerous, espacially when equallity of boot.R is not checked
boot.R <- pion.cor$boot.R

MKdata <- read.table( paste(data.paths[length(data.paths)], "mk.dat", sep="") )

Mk <- MKdata$V2[which(MKdata$V1 == ensembles[length(data.paths)])]
cat("Mk", Mk, "\n")

Mpi <- pion.matrixfit$opt.res$par[1]
Mpiboot <- c(Mpi, pion.matrixfit$opt.tsboot[1,])


E.cm <- array(dim=c(boot.R+1, 0))
delta <- array(dim=c(boot.R+1, 0))

labels <- c()
pch <- c()
col <- c()

all.pch <- c()
all.col <- c()

N.ensembles <- length(ensembles)
for( e in c(1:N.ensembles) ){

  Ensemble <- ensembles[e]
  data.path <- data.paths[e]
 
  for( d in c(1:length(all.dirs[[e]])) ){

    dir <- all.dirs[[e]][d]

    splitwd <- strsplit(dir, "/")[[1]]
    momentum <- splitwd[1]
    frameid <- which(momenta == momentum)
    frame <- frames[frameid]

    labels <- append(labels, paste(Ensemble, dir))
    pch <- append(pch, pchs[e])
    col <- append(col, cols[dir])

    for( pc.counter in c(1:max.pcs[[e]][d]) ){

      file <- paste(data.path, dir, "/res.", Ensemble, ".", frame, ".", splitwd[2], ".Rdata", sep="")
      cat("loading file", file, "\n")
      # dir is overwritten by version without / as last character. Workaround: No / in infile
      load(file)
  
      E.cm <- cbind(E.cm, res.boot[[pc.counter]][,1])
      delta <- cbind(delta, res.boot[[pc.counter]][,2])

      all.pch <- append(all.pch, pchs[e])
      all.col <- append(all.col, cols[dir])

    }
  }
}

ii <- which(E.cm[1,] < 2*Mk)

E.cm.below.KKbar <- E.cm[,ii]
delta.below.KKbar <- delta[,ii]
data <- cbind(E.cm.below.KKbar, delta.below.KKbar)

M <- invertCovMatrix(data, boot.samples=TRUE)

L <- chol(M)

par <- c(par, data[1,c(1:length(ii))])
npar <- length(par)

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

Mrho.res <- array(NA, dim=c(boot.R+1, npar+1))
filename <- paste("Mrho-res", ensembles[1], ".", ensembles[2], ".Rdata", sep="")
for(i in c(1:(boot.R+1))) {
  opt.res <- nls.lm(par, fn=Chi, y=data[i,], L=L, Mpi=Mpiboot[i], control = nls.lm.control(ftol=1.e-8, ptol=1.e-8, maxiter=500))
  Mrho.res[i, c(1:npar)] <- opt.res$par
  Mrho.res[i, npar+1] <- opt.res$rsstrace[length(opt.res$rsstrace)]
}

save(Mrho.res, data, file=filename)

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

## Create plot
tikzfiles <- tikz.init(basename="fit-delta1", width=4.5, height=5.)
plot(NA, xlim=c(2, 4), ylim=c(-0.4,pi+.4), xlab=c("$E_\\mathrm{CM}/M_\\pi$"), ylab=c("$\\delta_1$"))
abline(v=2*Mk/Mpi, lty=c(2))
text(x=2*Mk/Mpi+0.25, y=0.4, "$2M_K/M_\\pi$")

## Plot fit to data
polygon(x=c(x,rev(x))/Mpi, y=c(y+dy, rev(y-dy)), col="gray", lty=0, lwd=0.001, border="gray")
lines(x/Mpi,y)

## Plot all data
for( i in c(1:length(E.cm[1,])) ){
  plotwitherror(x=E.cm[1,i]/Mpi, y=delta[1,i], dx=sd(E.cm[,i])/Mpi, dy=sd(delta[,i]),
                pch=all.pch[i], col=all.col[i], rep=TRUE)
}

## Create legend
## Warning: Will not work if colors are symbols are used repeatedly
legend("topleft", legend=labels, bty="n", col=col, pch=pch)

tikz.finalize(tikzfiles=tikzfiles)

#if(FALSE) {
#  tikzfiles <- tikz.init(basename=paste("delta-fitqq", ensembles[1], sep=""), width=5., height=4.)
#  s <- seq(0,1,1./boot.R)
#  x <- qchisq(p=s, df=dof)
#  qqplot(x=x, y=Mrho.res[,dim(Mrho.res)[2]], xlab="Theoretical Quantiles", ylab="Sample Quantiles", main="QQ-Plot chisqr Values")
#  
#  tikz.finalize(tikzfiles=tikzfiles)
#}

sink()
