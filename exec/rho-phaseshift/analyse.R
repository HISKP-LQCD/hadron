## "Rscript analyse.R infile-analyse.R"
library(hadron)

## see phaseshift.rho.R for the implemented irreps

clargs = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(clargs)!=1) {
  stop("Please specify the input file!", call.=FALSE)
} else {
  # default output file
  input.file = clargs[1]
}

## Standard parameters

## t0 for the GEVP
t0 <- 2
## minimal number of timeslices in fit range
dof <- rep(4, times=6)
## momentum of the moving frame
p <- c(0,0,0)
## maximal number of principal correlators to analyse
maxpcs <- 2

## read input file and use namespace args
source(input.file)
ens <-  args$ens
disp <- args$disp
maxpcs <- args$maxpcs
hint <- rep("no", times=maxpcs)
L <- args$L
T <- args$T
## the left boundary of the fit interval runs from t10 to t11 in steps of 1
## the right from t11 to t21
## all these can be set specifically in parameters.R
t10 <- args$t10
t11 <- args$t11
t21 <- args$t21
t.step <- args$t.step
## force a data reread
reread <- args$reread
boot.R <- args$boot.R
boot.l <- args$boot.l
seed <- args$seed
dirs <- args$dirs
output.path <- paste(args$output.path, "5_fit-data", sep="/")

path.to.hadron <- getwd()

source(paste(path.to.hadron, "/phaseshift.rho.R", sep="/"))
source(paste(path.to.hadron, "/summarise.R", sep="/"))

## With R >= 3.2 dir.exists(output.path) can replace showWarnings=FALSE
dir.create(output.path, showWarnings=FALSE)
setwd(output.path)

## pion analysis
pion.filename <- paste(output.path, "/pion.Rdata", sep="")
if(reread || !file.exists(pion.filename)) {
  pion.cor <- bootstrap.cf(readtextcf("pi_p0.dat", T=T, check.t=0, path=args$path.to.data, ind.vector=c(2,2)), boot.R=boot.R, boot.l=boot.l, seed=seed)
  pion.matrixfit <- matrixfit(pion.cor, t1=16, t2=24, useCov=TRUE, parlist=array(c(1,1), dim=c(2,1)), sym.vec=c("cosh"), fit.method="lm")
  
  save(pion.cor, pion.matrixfit, file=pion.filename)
} else {
  load(pion.filename)
}

for(dir in dirs){

  dir.create(paste(output.path, dir, sep="/"), showWarnings=FALSE, recursive=TRUE)
  setwd(paste(output.path, dir, sep="/"))

  ## extracts irrep and frame from directory name
  ## also defines N for th ematrix size
  ## and path
  source(paste(path.to.hadron, "/detect_irrep_frame.R", sep="/"))

  ## full momentum is given by p*2*pi/L
  ## n is needed for function phaseshift.rho
  ## while p is needed for weighting and shifting
  n <- 1
  if(momentum == "p1") p <- c(0,0,1)
  if(momentum == "p2") p <- c(0,1,1)
  if(momentum == "p3") p <- c(1,1,1)
  if(momentum == "p4") {
    p <- c(0,0,2)
    n <- 2
  }
  
  model <- "single"
  sym.vec <- c("exp")
  if(momentum == "p0") {
    model <- "shifted"
    sym.vec<- c("cosh")
  }
  type="subtracted"
  if(frame != "cmf") type="weighted"

  cat("\nRunning for", ens, momentum, "with irrep", irrep, "size", N, "with momentum", p, "\n")

  if(file.exists("parameters.R")){
    print("Found local file parameters.R: Overwriting default")
    source("parameters.R")
    cat(readLines("parameters.R"), sep="\n")
  }

  sink("analyse.log", append=FALSE)
 
  ## Create correlation matrix 
  if(reread || !file.exists(paste("Cmatrix", ens, frame, irrep, "Rdata", sep="."))) {
    cat("reading raw data\n")
    Cmatrix <- cf()
    for(i in c(0:(N-1))) {
      for(j in c(0:(N-1))) {
        tmp <- readtextcf(paste("rho", ".", i, ".", j, ".dat", sep=""), T=T, check.t=1, path=path, ind.vector=c(2,2))
        Cmatrix <- c(Cmatrix, tmp)
      }
    }
    cat("reading done, bootstrapping now...\n")
  
    Cmatrix <- bootstrap.cf(Cmatrix, boot.R=boot.R, boot.l=boot.l, seed=seed)
    Cmatrix.rt <- removeTemporal.cf(Cmatrix, p1=p, p2=c(0,0,0), single.cf1=pion.matrixfit, L=L, lat.disp=TRUE)
  
    Cmatrix.bootstrap.gevp.n <- bootstrap.gevp(Cmatrix, matrix.size=N, element.order=c(1:N^2), t0=t0)
    Cmatrix.bootstrap.gevp <- bootstrap.gevp(Cmatrix.rt, matrix.size=N, element.order=c(1:N^2), t0=t0)
    cat("...done\n")
    save(Cmatrix, Cmatrix.rt, Cmatrix.bootstrap.gevp, Cmatrix.bootstrap.gevp.n, irrep, frame, N, file=paste("Cmatrix", ens, frame, irrep, "Rdata", sep="."))

  ## If it is created anew, solve gevp and plot eigenvalues and their effective mass 
#  if(!file.exists(paste(ens, ".", frame, ".", irrep, "-efm.pdf", sep=""))) {
  pdf(file=paste(ens, ".", frame, ".", irrep, "-efm.pdf", sep=""))
  plot(NA, xlim=c(1,T/2), ylim=c(0,1), xlab="t/a", ylab="Meff")
  clrs <- c("red", "blue", "darkgreen", "black", "purple", "orange")
  pchs <- c(21, 22, 23, 24, 25, 26)
  pcs <- c("pc1", "pc2", "pc3", "pc4", "pc5")
  for(id in c(1:N)) {
    pc <- gevp2cf(Cmatrix.bootstrap.gevp, id=id)
    effmass <- bootstrap.effectivemass(pc, boot.R=boot.R, boot.l=boot.l, seed=seed, type=type)
    plot(effmass, rep=TRUE, col=clrs[id], pch=pchs[id])
  }
  legend("topright", legend=pcs[1:N], col=clrs[1:N], bty="n", pch=pchs[1:N])
  
  plot(NA, xlim=c(1,T/2), ylim=c(0.000001,1), xlab="t/a", ylab="C(t)", log=c("y"))
  for(id in c(1:N)) {
    pc <- gevp2cf(Cmatrix.bootstrap.gevp, id=id)
    plot(pc, rep=TRUE, col=clrs[id], pch=pchs[id])
  }
  legend("topright", legend=pcs[1:N], col=clrs[1:N], bty="n", pch=pchs[1:N])
  dev.off()
#  }
  } else {
    cat("loading Rdata\n")
    load(paste("Cmatrix", ens, frame, irrep, "Rdata", sep="."))
    cat("...done\n")
  }
  
  for(id in c(1:min(N, maxpcs))) {
    if(id == 1) PC="pc1"
    if(id == 2) PC="pc2"
    if(id == 3) PC="pc3"
    if(id == 4) PC="pc4"
    if(id == 5) PC="pc5"

    pc <- gevp2cf(Cmatrix.bootstrap.gevp, id=id)
    if(file.exists(paste(PC, ".R", sep=""))) {
      source(paste(PC, ".R", sep=""))
    }

    cat(PC, "\n")
    cat("Mpi = ", pion.matrixfit$opt.res$par[1], sd(pion.matrixfit$opt.tsboot[1,]), "\n")
    pdf(file=paste(ens, PC, frame, "-fits.pdf", sep=""))
  
    for(t1 in seq(t10[id], t11[id], t.step)) {
      for(t2 in seq(t11[id], t21[id], t.step)) {
        if(t2-t1 > dof[id]) {
          cat(t1, t2, "\n")

          filename <- paste("rhoana", PC, t1, t2, ens, frame, irrep, "Rdata", sep=".")
          if(reread || !file.exists(filename)) {
            pc.matrixfit <- matrixfit(pc, t1=t1, t2=t2, useCov=TRUE, parlist=array(c(1,1), dim=c(2,1)), sym.vec=sym.vec, fit.method="lm", model=model)
          } else {
            load(filename)
          }

          plot(pc.matrixfit, main=paste(PC, "t1", t1, "t2", t2, "Qval", format(pc.matrixfit$Qval, digits=3), "E=", format(pc.matrixfit$opt.res$par[1], digits=3)))
  
          gs <- phaseshift.rho(pcfit =pc.matrixfit, L=L, Mpi=pion.matrixfit$opt.res$par[1], Mpiboot=pion.matrixfit$opt.tsboot[1,], frame=frame, irrep=irrep, disp=disp, n=n)
  
          save(pc.matrixfit, type, L, T, pion.matrixfit, gs, frame, ens, disp, PC, file=filename)
        
          cat(PC, t1, t2, pc.matrixfit$opt.res$par[1], sd(pc.matrixfit$opt.tsboot[1,]), gs$Ecm, sd(gs$Ecmboot, na.rm=TRUE),
          gs$delta, sd(gs$deltaboot, na.rm=TRUE), sin(gs$delta)^2, sd(sin(gs$deltaboot)^2, na.rm=TRUE),
          gs$tandelta, sd(gs$tandeltaboot, na.rm=TRUE), "Qval=", pc.matrixfit$Qval, "\n")
        }
      }
    }
    dev.off()
  }

  sink()

  ## average bootstrapsamples and write result to average.data.log
  source(paste(path.to.hadron, "/average.data.R", sep="/"))
  sink()

}
