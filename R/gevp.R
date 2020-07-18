## generates all permutations of c(1:n)
permutations <- function(n){
  if(n==1) {
    return(matrix(1))
  }
  
  sp <- permutations(n-1)
  p <- nrow(sp)
  A <- matrix(nrow=n*p,ncol=n)
  for(i in 1:n){
    A[(i-1)*p+1:p,] <- cbind(i,sp+(sp>=i))
  }
  return(A)
}

## t0 starts counting at 0!
##                           ( 1 , 2 )
## order c(1,2,3,4) goes to 
##                           ( 3 , 4 )



#' solve GEVP for correlator matrix
#' 
#' solve GEVP for a real, symmetric correlator matrix
#' 
#' The generalised eigenvalue problem\cr \eqn{ }{ C(t) v(t,t0) =
#' C(t0)lambda(t,t0) v(t,t0)}\eqn{ C(t) v(t,t_0) = C(t_0) \lambda(t,t_0)
#' v(t,t_0) }{ C(t) v(t,t0) = C(t0)lambda(t,t0) v(t,t0)}\eqn{ }{ C(t) v(t,t0) =
#' C(t0)lambda(t,t0) v(t,t0)}\cr is solved by performing a Cholesky
#' decomposition of \eqn{C(t_0)=L^t }{C(t0)=t(L) L}\eqn{ L}{C(t0)=t(L) L} and
#' transforming the GEVP into a standard eigenvalue problem for all values of
#' \eqn{t}. The matrices \eqn{C} are symmetrised for all \eqn{t}. So we solve
#' for \eqn{\lambda}{lambda}\cr \eqn{(L^t)^{-1} C(t) L^{-1} w = \lambda
#' w}{solve(t(L)) C(t) solve(L) w = lambda w}\cr with\cr \eqn{w = L v} or the
#' wanted \eqn{v = L^{-1} w}.
#' 
#' The amplitudes can be computed from\cr \eqn{ A_i^{(n)}(t) =
#' \sum_{j}C_{ij}(t) v_j^{(n)}(t,t_0)/(\sqrt{(v^{(n)}, Cv^{(n)})(\exp(-mt)\pm
#' \exp(-m(t-t)))}) } and this is what the code returns up to the factor\cr
#' \eqn{ 1/\sqrt{\exp(-mt)\pm \exp(-m(t-t))} } The states are sorted by their
#' eigenvalues when "values" is chosen. If "vectors" is chosen, we take \eqn{
#' \max( \sum_i \langle v(t_0,i), v(t, j)\rangle) } with \eqn{v} the
#' eigenvectors. For sort type "det" we compute \eqn{ \max(...)  }
#' 
#' @param cf correlation matrix preferably obtained with a call to
#' \code{extrac.obs} (or at leas with the same structure) or an already
#' averaged one.
#' 
#' cf is supposed to be an array of \code{dim=c(N, n*(Time/2+1))}, where
#' \code{N} is the number of observations and \code{n} is the number of single
#' correlators in the matrix. E.g. for a 2x2 matrix \code{n} would be 4.
#' @param Time time extent of the lattice.
#' @param t0 initial time value of the GEVP, must be in between 0 and
#' \code{Time/2-2}. Default is 1.
#' @param element.order specifies how to fit the \code{n} linearly ordered
#' single correlators into the correlator matrix.
#' \code{element.order=c(1,2,3,4)} leads to a matrix
#' \code{matrix(cf[element.order], nrow=2)}.
#' @param sort.type Sort the eigenvalues either in descending order, or by
#' using the scalar product of the eigenvectors with the eigenvectors at
#' \eqn{t=t_0+1}{t=t0+1}. Possible values are "values", "vectors" or "det".
#' @param for.tsboot for internal use of \code{\link{bootstrap.gevp}}. Alters
#' the returned values, see details.
#' @param sort.t0 if true (default), sort with respect to data at t0, otherwise
#' with respect to t-1.
#' @return Returns a list with the sorted eigenvalues, sorted eigenvectors and
#' sorted (reduced) amplitudes for all t > t0.
#' 
#' In case \code{for.tsboot=TRUE} the same is returned as one long vector with
#' first all eigenvalues concatenated, then all eigenvectors and then all
#' (reduced) amplitudes concatenated.
#' @author Carsten Urbach, \email{curbach@@gmx.de}
#' @seealso \code{boostrap.gevp}, \code{extract.obs}
#' @references Michael, Christopher and Teasdale, I., Nucl.Phys.B215 (1983)
#' 433, DOI: 10.1016/0550-3213(83)90674-0\cr Blossier, B. et al., JHEP 0904
#' (2009) 094, DOI: 10.1088/1126-6708/2009/04/094, arXiv:0902.1265
#' @keywords GEVP
#' @export gevp
gevp <- function(cf, Time, t0 = 1, element.order = 1:cf$nrObs,
                 for.tsboot=TRUE, sort.type="vectors", sort.t0=TRUE) {
  if(t0 < 0 || t0 > (Time/2-2)) {
    stop("t0 must be in between 0 and T/2-2. Aborting ...\n")
  }
  if(!any(c("values", "vectors", "det") == sort.type)) {
    stop("possible values for sort.ype are values or vectors. Aborting\n")
  }
  Thalf <- Time/2
  if(length(dim(cf)) == 2) {
    Cor <- apply(cf, 2, mean)
  }
  else {
    Cor <- cf
  }

  ## need to check consistency of cf here!
  ## can only operate on a square matrix

  ## number of correlators in cf
  Ncor <- length(Cor)/(Thalf+1)
  matrix.size <- as.integer(round(sqrt(length(element.order))))
  if(length(element.order) != matrix.size^2) {
    stop("gevp can only operate on square matrices, please adjust element.order! Aborting!\n")
  }
  if(max(element.order) > Ncor) {
    stop("element.order tries to index beyond the available correlators in cf! Aborting...\n")
  }
  ## index array for indexing the linear data
  ii <- c()
  for(i in c(1:Ncor)) {
    ii <- c(ii, (i-1)*(Thalf+1)+1)
  }
  ## re-order as to match the input order
  ii <- ii[element.order]
  
  evalues <-  array(NA, dim=c(Thalf+1, matrix.size))
  evectors <- array(NA, dim=c(Thalf+1, matrix.size, matrix.size))
  amplitudes <- array(NA, dim=c(Thalf+1, matrix.size, matrix.size))
  ## matrix at t=t0 (ii takes care of the indices starting at 1 and not 0)
  ## and symmetrise
  cM <- 0.5*matrix(Cor[ii+t0], nrow=matrix.size, ncol=matrix.size)
  cM <- cM + t(cM)
  ## check for positive definiteness
  ev.cM <- eigen(cM, symmetric=TRUE, only.values = TRUE)
  if(any(ev.cM$values < 0)) {
    stop("gevp: matrix at t0 is not positive definite. Aborting...\n")
  }
  ## compute Cholesky factorisation
  L <- chol(cM)
  invL <- solve(L)
  
  ## now the time dependence for t != t0
  ## we need to multiply from the left with t(invL) and from the right with invL
  for(t in c((t0 + 1):(Thalf), (t0-1):0)) {
    ## the +2 comes from the fact that evectors are stored at evectors[t+1,,]
    ## C to R index convention
    t.sort <- t0+2
    ## if wanted sort by the previous t, i.e. t.sort <- t + 1 - 1 for t > t0
    if((t > t0+1) && !sort.t0) t.sort <- t
    ## and t.sort <- t + 1 + 1 for t < t0
    if((t < t0-1) && !sort.t0) t.sort <- t + 2
    ## for t=t0-1 t.sort = t0+2, the default
    ## matrix at t and symmetrise
    cM <- 0.5*matrix(Cor[ii+t], nrow=matrix.size, ncol=matrix.size)
    cM <- cM + t(cM)
    if(any(is.na(cM))) {
      evalues[t+1,] <- NA
      evectors[t+1,,] <- NA
      amplitudes[t+1,,] <- NA
    }
    else {
      ## determine eigenvalues and vectors
      
      variational.solve <- eigen(t(invL) %*% cM %*% invL,
                                 symmetric=TRUE, only.values = FALSE, EISPACK=FALSE)
      ## sort depending on input by values or vectors
      sortindex <- integer(matrix.size)
      decreasing <- (t >= t0)
      if(sort.type == "values" || t == t0+1) {
        sortindex <- order(variational.solve$values, decreasing=decreasing)
      }
      else if(sort.type == "vectors") {
        ## compute the scalar product of eigenvectors with those at t.sort
        ## for each column apply order
        idx <- apply(abs( t(variational.solve$vectors) %*% evectors[t.sort,,] ), 1, order, decreasing=TRUE)
        ## obtain the indices of the maxima (depending on variable decreasing) per row
        sortindex <- idx[1,]
        ## if that fails, i.e. sortindex not a permutation, use t0+2 (i.e. physical t0+1) as reference time slice
        if(anyDuplicated(sortindex) && !sort.t0) {
          idx <- apply(abs( t(variational.solve$vectors) %*% evectors[t0+2,,] ), 1, order, decreasing=TRUE)
          sortindex <- idx[1,]
        }
        ## Fallback is simply sort by eigenvalues, i.e. sort.type="values"
        if(anyDuplicated(sortindex)) {
          sortindex <- order(variational.solve$values, decreasing=decreasing)
        }
      }
      else {
        Perms <- permutations(matrix.size)
        NPerms <- dim(Perms)[1]
        DMp <- integer(NPerms)
        Mp <- matrix(0, nrow=matrix.size, ncol=matrix.size)
        for(p in c(1:NPerms)) {
          for(i in c(1:matrix.size)) {
            ij <- Perms[p,c(1:matrix.size)[-i]]
            if(i == 1) Mp <- matrix(c(variational.solve$vectors[,i], as.vector(evectors[t.sort,,ij])), nrow=matrix.size, ncol=matrix.size)
            else {
              Mp <- Mp %*% matrix(c(variational.solve$vectors[,i], as.vector(evectors[t.sort,,ij])), nrow=matrix.size, ncol=matrix.size)
            }
          }
          DMp[p] <- determinant(Mp, logarithm=FALSE)$modulus
        }
        sortindex <- Perms[which.max(DMp),]
        ##if(!decreasing) sortindex <- rev(sortindex)
      }
      evalues[t+1,] <- variational.solve$values[sortindex]
      evectors[t+1,,] <- variational.solve$vectors[, sortindex]
    }
  }
  for(t in c((0:(t0-1)), (t0 + 1):(Thalf))) {
    evectors[t+1,,] <- invL %*% evectors[t+1,,]
    tmp <- matrix(Cor[ii+t], nrow=matrix.size, ncol=matrix.size) %*% evectors[t+1,,]
    ## t(evectors[t+1,,]) %*% tmp must be proportional to delta_ij
    ## these are the amplitudes up to a factor sqrt(exp(-mt) \pm exp(-m(T-t)))
    ## diag(t(evectors[t+1,,]) %*% tmp) might get negative due to fluctuations
    ## we set them to NA first
    d <- diag(t(evectors[t+1,,]) %*% tmp)
    d[d < 0] <- NA
    amplitudes[t+1,,] <- t(t(tmp)/sqrt(d))
    rm(tmp)
  }
  ## at t0 we set to 1
  evalues[t0+1,] <- 1.
  ## in case of bootstrapping everything (eigenvalues and eigenvectors)
  ## is concatenated into a single vector
  if(for.tsboot) {
    return(c(as.vector(evalues), as.vector(amplitudes), as.vector(evectors)))
  }
  
  return(invisible(list(evalues=evalues, evectors=evectors, amplitudes=amplitudes)))
}




#' perform a bootstrap analysis of a GEVP
#' 
#' perform a bootstrap analysis of a GEVP for a real, symmetric correlator
#' matrix
#' 
#' Say something on "det" sorting method.
#' 
#' @param cf correlation matrix obtained with a call to \code{extrac.obs}.
#' @param t0 initial time value of the GEVP, must be in between 0 and
#' \code{Time/2-2}. Default is 1.
#' @param element.order specifies how to fit the \code{n} linearly ordered
#' single correlators into the correlator matrix.
#' \code{element.order=c(1,2,3,4)} leads to a matrix
#' \code{matrix(cf[element.order], nrow=2)}.  Double indexing is allowed.
#' @param sort.type Sort the eigenvalues either in descending order, or by
#' using the scalar product of the eigenvectors with the eigenvectors at
#' \eqn{t=t_0+1}{t=t0+1}. Possible values are "values", "vectors" and "det".
#' The last one represents a time consuming, but in principle better version of
#' sorting by vectors.
#' @param sort.t0 for \code{sort.type} "vectors" use \eqn{t_0}{t0} as reference
#' or \eqn{t-1}{t-1}.
#' @return Returns an object of class \code{gevp} with member objects:
#' 
#' \code{cf}:\cr The input data, if needed bootstrapped with
#' \code{\link{bootstrap.cf}}.
#' 
#' \code{res.gevp}:\cr The object returned from the call to \code{\link{gevp}}.
#' For the format see \code{\link{gevp}}.
#' 
#' \code{gevp.tsboot}:\cr The bootstrap samples of the GEVP. For the format see
#' \code{\link{gevp}}.
#' @author Carsten Urbach, \email{curbach@@gmx.de}
#' @seealso \code{gevp}, \code{extract.obs}, \code{bootstrap.cf}
#' @references Michael, Christopher and Teasdale, I., Nucl.Phys.B215 (1983)
#' 433, DOI: 10.1016/0550-3213(83)90674-0\cr Blossier, B. et al., JHEP 0904
#' (2009) 094, DOI: 10.1088/1126-6708/2009/04/094, arXiv:0902.1265
#' @keywords GEVP
#' @examples
#' 
#' data(correlatormatrix)
#' ## bootstrap the correlator matrix
#' correlatormatrix <- bootstrap.cf(correlatormatrix, boot.R=99, boot.l=1, seed=132435)
#' ## solve the GEVP
#' t0 <- 4
#' correlatormatrix.gevp <- bootstrap.gevp(cf=correlatormatrix, t0=t0, element.order=c(1,2,3,4))
#' ## extract the ground state and plot
#' pc1 <- gevp2cf(gevp=correlatormatrix.gevp, id=1)
#' plot(pc1, log="y")
#' ## determine the corresponding effective masses
#' pc1.effectivemass <- bootstrap.effectivemass(cf=pc1)
#' pc1.effectivemass <- fit.effectivemass(cf=pc1.effectivemass, t1=5, t2=20)
#' ## summary and plot
#' summary(pc1.effectivemass)
#' plot(pc1.effectivemass)
#' 
#' ## we can also use matrixfit with a special model for a principal
#' ## correlators
#' pc1.matrixfit <- matrixfit(pc1, t1=2, t2=24, fit.method="lm", model="pc", useCov=FALSE,
#'                            parlist=array(c(1,1), dim=c(2,1)), sym.vec=c("cosh"), neg.vec=c(1))
#' summary(pc1.matrixfit)
#' plot(pc1.matrixfit)
#' 
#' ## the same can be achieved using bootstrap.nlsfit
#' model <- function(par, x, t0, ...) {
#'   return(exp(-par[1]*(x-t0))*(par[3]+(1-par[3])*exp(-par[2]*(x-t0))))
#' }
#' ii <- c(2:4, 6:25)
#' fitres <- parametric.nlsfit(fn=model, par.guess=c(0.5, 1, .9),
#'                             y=pc1$cf0[ii], dy=pc1$tsboot.se[ii],
#'                             x=ii-1, boot.R=pc1$boot.R, t0=t0)
#' summary(fitres)
#' plot(fitres, log="y")
#' 
#' @export bootstrap.gevp
bootstrap.gevp <- function(cf, t0 = 1, element.order = 1:cf$nrObs,
                           sort.type = "vectors", sort.t0 = TRUE) {
  stopifnot(inherits(cf, 'cf_meta'))
  stopifnot(inherits(cf, 'cf_boot'))
  stopifnot(inherits(cf, 'cf_orig'))
  stopifnot(sort.type %in% c("values", "vectors", "det"))

  N <- length(cf$cf[,1])
  seed <- cf$seed
  boot.R <- cf$boot.R
  boot.l <- cf$boot.l
  matrix.size <- as.integer(round(sqrt(length(element.order))))
  res <- gevp(cf$cf0, Time=cf$Time, t0=t0, element.order=element.order, for.tsboot=FALSE, sort.type=sort.type, sort.t0=sort.t0)

  gevp.tsboot <- t(apply(cf$cf.tsboot$t, 1, gevp, Time=cf$Time, t0=t0,
                         element.order=element.order,
                         for.tsboot=TRUE, sort.type=sort.type, sort.t0=sort.t0))

  ## gevp.tsboot contains first the N*(Thalf+1) eigenvalues
  ## and the the N*N*(Thalf+1) eigenvectors
  
  ret <- list(cf=cf,
              res.gevp=res,
              gevp.tsboot=gevp.tsboot,
              boot.R=boot.R,
              boot.l=boot.l,
              seed=seed,
              matrix.size=matrix.size,
              sort.type=sort.type,
              t0=t0,
              sort.t0=sort.t0)
  class(ret) <- c("gevp", class(ret))
  return(invisible(ret))
}



#' Extracts a principle correlator from a GEVEP
#' 
#' Extracts a principle correlator from a GEVP and converts it into an object
#' of class \code{cf}
#' 
#' 
#' @param gevp An object returned by \code{\link{bootstrap.gevp}}.
#' @param id The index of the principal correlator to extract.
#' @return An object of class \code{cf}, which contains bootstrap samples
#' already. So a call to \code{bootstrap.cf} is neither needed nor possible. It
#' can be treated further by \code{\link{bootstrap.effectivemass}} or
#' \code{\link{matrixfit}} to extract a mass value.
#' @author Carsten Urbach, \email{curbach@@gmx.de}
#' @seealso \code{\link{gevp}}, \code{\link{matrixfit}},
#' \code{\link{bootstrap.effectivemass}}
#' @keywords GEVP
#' @examples
#'
#' data(correlatormatrix)
#' ## bootstrap the correlator matrix
#' correlatormatrix <- bootstrap.cf(correlatormatrix, boot.R=99, boot.l=1, seed=132435)
#' ## solve the GEVP
#' t0 <- 4
#' correlatormatrix.gevp <- bootstrap.gevp(cf=correlatormatrix, t0=t0, element.order=c(1,2,3,4))
#' ## extract the ground state and plot
#' pc1 <- gevp2cf(gevp=correlatormatrix.gevp, id=1)
#' plot(pc1, log="y")
#'
#' @export gevp2cf
gevp2cf <- function(gevp, id=1) {
  stopifnot(inherits(gevp, "gevp"))
  stopifnot(!(id > gevp$matrix.size || id < 1))

  # Base `cf` properties.
  cf <- cf_meta(nrObs = 1,
                Time = gevp$cf$Time,
                nrStypes = 1,
                symmetrised = gevp$cf$symmetrised)

  # Add the `cf_boot` mixin.
  tt <- (id-1)*(cf$Time/2+1)+seq(1, cf$Time/2+1)
  cf.tsboot <- list(t = gevp$gevp.tsboot[,tt],
                    t0 = gevp$res.gevp$evalues[,id])

  cf <- cf_boot(cf,
                boot.R = gevp$boot.R,
                boot.l = gevp$boot.l,
                seed = gevp$seed,
                sim = gevp$cf$sim,
                endcorr = gevp$cf$endcorr,
                cf.tsboot = cf.tsboot,
                resampling_method = gevp$cf$resampling_method)

  cf <- cf_principal_correlator(cf,
                                id = id,
                                gevp_reference_time = gevp$t0)

  if (inherits(gevp$cf, 'cf_shifted')) {
    cf <- cf_shifted(cf,
                     deltat = gevp$cf$deltat,
                     forwardshift = gevp$cf$forwardshift)
  }

  if (inherits(gevp$cf, 'cf_weighted')) {
    cf <- cf_weighted(cf,
                      weight.factor = gevp$cf$weight.factor,
                      weight.cosh = gevp$cf$weight.cosh)
  }

  # Add some other stuff
  cf$N <- length(gevp$cf$cf[,1])

  return (invisible(cf))
}



#' Extracts physical amplitudes from a GEVP
#' 
#' Given a GEVP generated with \code{bootstrap.gevp} and masses determined from
#' the principle correlator with given \code{id}, the physical amplitudes are
#' extracted and bootstraped. The man amplitude is determined from a constant
#' fit to the data in the specified time range.
#' 
#' 
#' @param gevp An object of class \code{gevp} as generated with a call to
#' \code{bootstrap.gevp}.
#' @param mass Optimally, this is an object either of class
#' \code{effectivemassfit} generated using \code{\link{fit.effectivemass}} or
#' of class \code{matrixfit} generated with \code{\link{matrixfit}} to the
#' principal correlator extracted using \code{\link{gevp2cf}} applied to
#' \code{gevp} using the same value of \code{id}.
#' 
#' It can also be given as a numerical vector with the bootstrap samples as
#' entries. The mean will then be computed as the bootstrap mean over this
#' vector. The number of samples must agree with the number of bootstrap
#' samples in \code{gevp}.
#' @param id The index of the principal correlator to extract, i.e. the
#' physical state to extract.
#' @param op.id The index of the operator for which to extract the amplitude.
#' @param type The symmetry of the pricipal correlator in time, can be either
#' "cosh" or "sinh".
#' @param t1,t2 The time range in which to fit the amplitude starting with 0.
#' If not given it will be tried to infer these from the \code{mass} object.
#' @param useCov Use the covariance matrix for fitting the constant to the
#' amplitude data.
#' @param fit perform a fit to the data.
#' @author Carsten Urbach, \email{curbach@@gmx.de}
#' @seealso \code{\link{matrixfit}}, \code{\link{fit.effectivemass}},
#' \code{\link{gevp}}, \code{\link{gevp2cf}}, \code{\link{computefps}}
#' @keywords GEVP
#'
#' @return
#' Returns an object of S3 class `gevp.amplitude`, generated as a list with named
#' elements `amplitude` the numeric vector of amplitudes, `amplitude.tsboot`
#' the corresponding bootstrap samples, `damplitude` the estimates for the
#' standard errors, `fit` the object returned by the fit routine,
#' `meanAmplitude` and `meanAmplitude.tsboot` mean amplitude and its
#' bootstrap samples, `chisqr` the residual sum of squares, `dof` the numberi of
#' degrees of freedom, `t1` and `t2` the fit range, and then all the input
#' objects.
#' 
#' @examples
#'
#' data(correlatormatrix)
#' ## bootstrap the correlator matrix
#' correlatormatrix <- bootstrap.cf(correlatormatrix, boot.R=99, boot.l=1, seed=132435)
#' ## solve the GEVP
#' t0 <- 4
#' correlatormatrix.gevp <- bootstrap.gevp(cf=correlatormatrix, t0=t0, element.order=c(1,2,3,4))
#' ## extract the ground state and plot
#' pion.pc1 <- gevp2cf(gevp=correlatormatrix.gevp, id=1)
#' pion.pc1.effectivemass <- bootstrap.effectivemass(cf=pion.pc1, type="solve")
#' pion.pc1.effectivemass <- fit.effectivemass(pion.pc1.effectivemass, t1=8, t2=23,
#'                                             useCov=FALSE)
#' ## now determine the amplitude
#' pion.pc1.amplitude <- gevp2amplitude(correlatormatrix.gevp, pion.pc1.effectivemass,
#'                                      useCov=FALSE, t1=8, t2=14)
#' plot(pion.pc1.amplitude)
#' summary(pion.pc1.amplitude)
#'
#' @export gevp2amplitude
gevp2amplitude <- function(gevp, mass, id=1, op.id=1, type="cosh", t1, t2, useCov=TRUE, fit=TRUE) {
  if(id > gevp$matrix.size || id < 1 || op.id > gevp$matrix.size || op.id < 1) {
    stop("gevp2cf: id and op.id must be <= matrix.size and > 0. Aborting...\n")
  }
  if(missing(t1) || missing(t2)) {
    if(inherits(mass, "effectivemassfit") || inherits(mass, "matrixfit")) {
      t1 <- mass$t1
      t2 <- mass$t2
    }
    else {
      stop("gevp2amplitude: t1 and/or t2 missing... Aborting\n")
    }
  }
  if((t2 <= t1) || (t1 < 0) || (t2 > (gevp$cf$Time/2-1))) {
    stop("gevp2amplitude: t1 < t2 and both in 0...Time/2-1 is required. Aborting...\n")
  }
  
  sign <- +1
  if(type != "cosh") sign <- -1.
  if(!inherits(gevp, "gevp")) {
    stop("gevp2amplitude requires an element of class gevp as input. Aborting...\n")
  }
  if(is.numeric(mass)){
    if(length(mass) != gevp$boot.R) {
      stop(paste("gevp2amplitude: gevp and mass differ in number of bootstrap samples:", gevp$boot.R, "and", length(mass), ". Aborting...\n"))
    }
    m0 <- mean(mass)
    m <-  mass
  }
  else if(inherits(mass, "effectivemassfit") || inherits(mass, "matrixfit")) {
    if(gevp$boot.R != mass$boot.R) {
      stop(paste("gevp2amplitude: gevp and mass differ in number of bootstrap samples:", gevp$boot.R, "and", mass$boot.R, ". Aborting...\n"))
    }
    m0 <- mass$opt.res$par[1]
    m  <- mass$massfit.tsboot[,1]
  }
  else {
    stop("gevp2amplitude requires a numeric vector or an object either of type effectivemassfit or matrixfit as input. Abortgin...\n")
  }
  Time <- gevp$cf$Time
  t <- c(0:(Time/2))
  amplitude <- abs(gevp$res.gevp$amplitudes[,id,op.id])/sqrt(.5*(exp(-m0*t)+ sign*exp(-m0*(Time-t))))
  tt <- gevp$matrix.size*(Time/2+1) + ((id-1)*gevp$matrix.size+(op.id-1))*(Time/2+1) + seq(1, Time/2+1)
  amplitude.tsboot <- array(NA, dim=c(gevp$boot.R, Time/2+1))
  for(i in c(1:gevp$boot.R)) {
    amplitude.tsboot[i,] <- abs(gevp$gevp.tsboot[i,tt])/sqrt(.5*(exp(-m[i]*t)+ sign*exp(-m[i]*(Time-t))))
  }
  damplitude <- apply(amplitude.tsboot, 2, gevp$cf$error_fn)
  
  ## now we perform a constant fit
  ii <- c((t1+1):(t2+1))
  if(any(is.na(amplitude[ii])) || any(is.na(damplitude[ii]))) {
    stop("At least one amplitude or its error take the value NA, change t1 and t2! Aborting!")
  }
  M <- diag(1/damplitude[ii]^2)
  if(useCov) {
    ## compute correlation matrix and compute the correctly normalised inverse
    M <- try(invertCovMatrix(amplitude.tsboot[,ii], boot.samples=TRUE), silent=TRUE)
    if(inherits(M, "try-error")) {
      warning("inversion of variance covariance matrix failed in gevp2amplitude, continuing with uncorrelated chi^2\n")
      M <- diag(1/damplitude[ii]^2)
      useCov <- FALSE
    }
  }
  meanAmplitude <- NA
  opt.res <- list()
  opt.res$value <- NA
  if(fit) {
    ## the chisqr function
    fn <- function(par, y, M) { sum((y-par[1]) %*% M %*% (y-par[1]))}
    
    par <- c(amplitude[t1+1])
    opt.res <- optim(par, fn = fn,
                     method="BFGS", M=M, y = amplitude[ii])
    opt.res <- optim(opt.res$par, fn = fn,
                     control=list(parscale=1/opt.res$par),
                     method="BFGS", M=M, y = amplitude[ii])
    meanAmplitude <- par[1]
    par <- opt.res$par
  }
  meanAmplitude.tsboot <- array(0, dim=c(gevp$boot.R, 2))
  if(fit) {
    for(i in 1:gevp$boot.R) {
      opt <- optim(par, fn = fn,
                   control=list(parscale=1/par),
                   method="BFGS", M=M, y = amplitude.tsboot[i,ii])
      meanAmplitude.tsboot[i, 1] <- opt$par[1]
      meanAmplitude.tsboot[i, 2] <- opt$value
    }
  }
  res <- list(amplitude=amplitude,
              amplitude.tsboot=amplitude.tsboot,
              damplitude=damplitude, fit=fit,
              meanAmplitude=meanAmplitude,
              meanAmplitude.tsboot=meanAmplitude.tsboot,
              chisqr = opt.res$value,
              dof=t2-t1, t1=t1, t2=t2,
              mass=mass, gevp=gevp,
              boot.R=gevp$boot.R, boot.l=gevp$boot.l, seed=gevp$seed,
              id=id, op.id=op.id,
              Time=Time, m0=m0, m0.tsboot=m, useCov=useCov,
              Qval=1-pchisq(opt.res$value, t2-t1),
              error_fn = gevp$cf$error_fn)
  attr(res, "class") <- c("gevp.amplitude", class(res))
  return(invisible(res))
}

#' summary.gevp.amplitude
#'
#' @param object Object of type `gevp.amplitude`.
#' @param ... Generic Parameters to be passed on.
#'
#' @return
#' No return values.
#' 
#' @export
summary.gevp.amplitude <- function (object, ...) {
  amp <- object
  cat("\n ** Result of a GEVP analysis for the amplitude **\n\n")
  cat("time range from", amp$t1, " to ", amp$t2, "\n")
  cat("mass:\n")
  cat("m \t=\t", amp$m0, "\n")
  cat("dm\t=\t", amp$error_fn(amp$m0.tsboot), "\n")
  cat("\nAmplitude:\n")
  cat("operator id:", amp$op.id, "\n")
  cat("state id   :", amp$id, "\n")
  cat(" P[", amp$id, ",", amp$op.id, "] = ", amp$meanAmplitude, "\n")
  cat("dP[", amp$id, ",", amp$op.id, "] = ", amp$error_fn(amp$meanAmplitude.tsboot[,1]), "\n")
  cat("\n")
  cat("boot.R\t=\t", amp$gevp$boot.R, "\n")
  cat("boot.l\t=\t", amp$gevp$boot.l, "\n")
  cat("useCov\t=\t", amp$useCov, "\n")
  cat("chisqr\t=\t", amp$chisqr, "\ndof\t=\t", amp$dof, "\nchisqr/dof=\t",
      amp$chisqr/amp$dof, "\n")
  cat("Quality of the fit (p-value):", amp$Qval, "\n")
  if(any(names(amp) == "fps")) {
    cat("\nDecay Constant (derived quantity):\n")
    cat("mu1 \t=\t", amp$mu1, "\n")
    cat("mu2 \t=\t", amp$mu2, "\n")
    if(amp$normalisation == "cmi") cat("kappa\t=\t", amp$kappa,"\n")
    cat("fps \t=\t", amp$fps, "\n")
    cat("dfps\t=\t", amp$error_fn(amp$fps.tsboot), "\n")
  }
}

#' plot.gevp.amplitude
#'
#' @param x Object of type `gevp.amplitude`.
#' @param xlab x axis label
#' @param ylab y axis label
#' @param ... Graphical parameters to be passed on.
#'
#' @return
#' No return value.
#' 
#' @export
plot.gevp.amplitude <- function (x, xlab="t",
                                 ylab=paste0("P[,",  x$id, ",", x$op.id, "]"),
                                 ...) {
  amp <- x
  plotwitherror(c(0:(amp$Time/2)), amp$amplitude, amp$damplitude,
                xlab=xlab, ylab=ylab, ...)
  if(amp$fit) {
    arrows(x0=amp$t1, y0=amp$meanAmplitude,
           x1=amp$t2, y1=amp$meanAmplitude, col=c("red"), length=0)
    arrows(x0=amp$t1, y0=amp$meanAmplitude+amp$error_fn(amp$meanAmplitude.tsboot[,1]),
           x1=amp$t2, y1=amp$meanAmplitude+amp$error_fn(amp$meanAmplitude.tsboot[,1]),
           col=c("red"), length=0, lwd=c(1))
    arrows(x0=amp$t1, y0=amp$meanAmplitude-amp$error_fn(amp$meanAmplitude.tsboot[,1]),
           x1=amp$t2, y1=amp$meanAmplitude-amp$error_fn(amp$meanAmplitude.tsboot[,1]),
           col=c("red"), length=0, lwd=c(1))
  }
}
