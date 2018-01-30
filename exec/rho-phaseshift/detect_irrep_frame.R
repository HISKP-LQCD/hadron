## detects the irrep and frame from the directory structure
## its supposed to be called from a directory with name
## .../p<N>/<irrep>
## with <N> being the squared momentum in units of 2pi/L
## and <irrep> the corresponding irreducible representation

## PATH must be set by the calling script to a directory such that
## the data can be found in PATH/p<N>/<irrep>/

wd <- getwd()
splitwd <- strsplit(wd, "/")[[1]]
momentum <- splitwd[length(splitwd)-1]
irrep <- splitwd[length(splitwd)]
rm(wd, splitwd)

momenta <- c("p0", "p1", "p2", "p3", "p4")

frames <- c("cmf", "mf1", "mf2", "mf3", "mf1")
frameid <- which(momenta == momentum)
frame <- frames[frameid]

irreps <- list(c("T1"), c("A1", "E2"), c("A1", "B1", "B2"), c("A1"), c("A1"))
irrepid <- which(irreps[[frameid]] == irrep)
sizes <- list(c(3), c(5, 3), c(4, 4, 3), c(4), c(2))

path <- paste(PATH, "/", momentum, "/", irrep, "/", sep="")

N <- 1
for(N in c(0:10)) {
  if(!file.exists(paste(path, "/", "Rho_Gevp_", momentum, "_", irrep, ".", N, ".", N, ".dat", sep=""))) {
    break
  }
}
