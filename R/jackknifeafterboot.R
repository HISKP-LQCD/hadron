find.duplicates <- function(xstar, x) {
  duplicated(c(xstar, x))[(length(x) + 1):(2 * length(x))]
}

jack.boot <- function(indices, xstar, f) {
  if(is.null(dim(xstar)))   apply(xstar[!indices], MARGIN=2L, FUN=f)
  else apply(xstar[!indices, ], MARGIN=2L, FUN=f)
}

jackknifeafterboot <- function(cf, m=1) {
  jab(cf=cf, m=m)
}

jab.cf <- function(cf, m=1) {
  if(!any(class(cf) == "cf")) {
    stop("bootstrap.cf requires an object of class cf as input! Aborting!\n")
  }
  if(!cf$boot.sample) {
    stop("cf must be bootstrapped already using bootstrap.cf! Aborting!\n")
  }
  if(cf$cf.tsboot$sim != "fixed") {
    stop("JAB only implemented for 'sim=fixed' at the moment")
  }

  ## save random number generator state
  if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
    temp <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  else temp <- NULL
  ## we set the seed as used for tsboot
  set.seed(cf$seed)

  m <- 1
  ## the resampling block indices
  blockind <- boot:::ts.array(n=cf$cf.tsboot$n, n.sim=cf$cf.tsboot$n.sim,
                              R=cf$boot.R, l=cf$boot.l, sim=cf$sim, endcorr=cf$cf.tsboot$endcorr)

  duplicates <- t(apply(blockind$starts, MARGIN=1L, FUN=find.duplicates, c(1:ncol(blockind$starts))))
  jack.boot.values <- apply(duplicates, MARGIN=2L, FUN=jack.boot, xstar=cf$cf.tsboot$t, f=sd)

  ## total number of blocks
  N <- nrow(jack.boot.values)
  M <- N - m + 1
  phitilde <- (N*cf$tsboot.se - (N-m)*jack.boot.values)/m - cf$tsboot.se
  jack.boot.se <- sqrt(m/(N-m)/M * apply(phitilde, MARGIN=1L,
                                           FUN=function(x) {sum(x^2)}))

  
  ## restore random number generator state
  if (!is.null(temp))
    assign(".Random.seed", temp, envir = .GlobalEnv)
  else rm(.Random.seed, pos = 1)
  cf$jack.boot.se <- jack.boot.se
  return(invisible(cf))
}
