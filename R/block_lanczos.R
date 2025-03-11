#' @title Block-Lanczos method for LQCD correlators
#' 
#' @description
#'   Taking a correlation matrix as input, the method
#'   determines the first Neig energy levels plus its bootstrap uncertainty.
#'
#' @param cf object of type \link{cf}, optimally returned by
#'   \code{\link{bootstrap.cf}}
#' @param N Integer. Maximal time index in correlation function to be used in
#'                   Lanczos analysis
#' @param bias_correction boolean. If set to 'TRUE', the median of the bootstrap
#'   distribution is used as estimator for the energy values.
#'   This will be set to TRUE for errortyp is equal 'dbboot'
#' @param errortype string. Determines the treatment of the bootstrap
#'   histograms to determine the statistical error on eigenvalues. Can
#'   be: 1. 'outlier-removal' for which outliers are removed according to
#'   the 0.25 and 0.75 quantiles and the inter-quantile-range,
#'   i.e. only values are kept which are in the interval
#'   \eqn{[Q_25-1.5IQR, Q_75+1.5IQR]}
#'   and the error is computed from the standard deviation of the bootstrap distribution.
#'   2. 'quantiles' for which the error is estimated from the difference
#'   between the 0.16 and 0.84 quantile of the original bootstrap distribution
#'   3. 'dbboot' which works only, if the 'cf' is double bootstrapped. It will
#'   estimate the error from the true error of the median
#' @param probs numeric. Vector of probabilities for the error estimation method
#'   'quantiles'.
#' @seealso \code{\link{plot.effectivemass}}, \code{\link{bootstrap.effectivemass}}
#' @references M. Wagman, 'Lanczos, the transfer matrix, and the signal-to-noise problem',
#'   arXiv:2406.20009 , D.C. Hackett, M. Wagman, 'Block Lanczos for lattice QCD spectroscopy and matrix elements'
#'   arXiv:2412.04444
#' @return
#'   Returns an object of S3 class `effectivemass`.
#' 
#' @family lanczos
#' @export
#' @examples
#' data(pscor.sample)
#' newcf <- cf_orig(cf=t(array(pscor.sample[,2], dim=c(48, 316))))
#' newcf <- cf_meta(newcf, nrObs=1, Time=48, symmetrised=FALSE)
#' newcf.boot <- bootstrap.cf(newcf)
#' ncf.boot <- symmetrise.cf(newcf.boot)
#' ncf.effmass <- bootstrap.effectivemass(ncf.boot)
#' plot(ncf.effmass, ylim=c(0.1,0.2))
#' res <- bootstrap.lanczos(newcf.boot, N=newcf$Time)
#' plot(res, rep=TRUE, col="red", pch=22, xshift=0.2)
bootstrap.block_lanczos <- function(cf, N = (cf$Time/2+1), Neig=1, bias_correction=FALSE,
                              errortype="outlier-removal", probs=c(0.16,0.84)) {
  ## wrapper function, not yet bootstrapping...
  stopifnot(inherits(cf, 'cf_meta'))
  stopifnot(inherits(cf, 'cf_boot'))
  stopifnot(errortype %in% c("outlier-removal", "quantiles", "dbboot"))

  dbboot <- inherits(cf, 'cf_dbboot')
  if(errortype == "dbboot" & (!dbboot)) {
    if(!dbboot) cat("errortype dbboot needs a doubly bootstrapped cf\n")
    stopifnot(dbboot)
  }
  if(errortype == "dbboot") bias_correction = TRUE
  else(dbboot = FALSE)
  
  seed <- cf$seed
  boot.R <- cf$boot.R
  boot.l <- cf$boot.l
  res <- block_lanczos.solve(cf=cf$cf0,Time=cf$Time, Neig=Neig,element.order=1:cf$nrObs)
  effMass <- -log(res$eigvalues)
  lanczos.tsboot.orig <- t(apply(cf$cf.tsboot$t, 1, block_lanczos.solve, Time=cf$Time, Neig=Neig, element.order=1:cf$nrObs))
  lanczos.tsboot <- lanczos.tsboot.orig
  lanczos.dbboot <- array()
  deffMass <- rep(NA, length(effMass))
  if(errortype=="outlier-removal") {
    remove_outliers <- function(x, probs=c(0.25,0.75)) {
      for (line in 1:Neig)
      Q <- quantile(x, probs=probs, na.rm=TRUE)
      iqr <- Q[2]-Q[1]
      x[x<(Q[1]-1.5*iqr) | x > (Q[2] + 1.5*iqr)] <- NA
      return(invisible(x))
    }
    deffMass <- NULL
    lanczos.tsboot <- NULL
    for (lineeig in 1:Neig){
      bs_samples <- NULL
      for (linebs in 1:length(lanczos.tsboot.orig)){
        tmp <- filter(lanczos.tsboot.orig[[linebs]],index==lineeig)$eigvalues
        bs_samples <- c(bs_samples, tmp)
      }
      ncol=length(bs_samples)/length(lanczos.tsboot.orig)
      nrow=length(lanczos.tsboot.orig)
      bsm_samples <- matrix(0,nrow=nrow,ncol=ncol)
      for (line1 in 1:nrow){
	for (line2 in 1:ncol){
	  bsm_samples[line1,line2] <- bs_samples[(line1-1)*ncol+line2]
        }
      }
      bsm_samplessecond <- apply(bsm_samples, 2, remove_outliers)
      lanczos.tsboot <- c(lanczos.tsboot, bsm_samplessecond)
      deffMass <- c(deffMass, apply(-log(bsm_samplessecond), 2L, cf$error_fn, na.rm=TRUE))
    }
  }
  else if(errortype == "quantiles") {
    error_fn <- function(x, probs=c(0.16, 0.84)) {
      str(x)
      Q <- quantile(x, probs=probs)
      return(Q[2]-Q[1])
    }
    
    deffMass <- NULL

    for (lineeig in 1:Neig){
      bs_samples <- NULL
      for (linebs in 1:length(lanczos.tsboot.orig)){
	tmp <- filter(lanczos.tsboot.orig[[linebs]],index==lineeig)$eigvalues
        bs_samples <- c(bs_samples, tmp)
      }
      ncol=length(bs_samples)/length(lanczos.tsboot.orig)
      nrow=length(lanczos.tsboot.orig)
      bsm_samples <- matrix(0,nrow=nrow,ncol=ncol)
      for (line1 in 1:nrow){
        for (line2 in 1:ncol){
          bsm_samples[line1,line2] <- bs_samples[(line1-1)*ncol+line2]
        }
      }
      deffMass <- c(deffMass, apply(-log(bsm_samples), 2L, error_fn, probs=probs))
    }
  }
  ret <- list(t.idx=c((res$m)*2-1), cf=cf, res.lanczos=res, #bias=bias,
              lanczos.tsboot.orig=lanczos.tsboot.orig, lanczos.tsboot=lanczos.tsboot,
              effMass=effMass, deffMass=deffMass, effMass.tsboot=-log(lanczos.tsboot),
              effMass.dbboot=lanczos.dbboot,
              opt.res=NULL, t1=NULL, t2=NULL, type="log", useCov=NULL, CovMatrix=NULL, invCovMatrix=NULL,
              boot.R = boot.R, boot.l = boot.l, seed = seed,
              massfit.tsboot=NULL, Time=cf$Time, nrObs=1, dof=NULL,
              chisqr=NULL, Qval=NULL
             )
  ret$t0 <- effMass
  ret$t <- ret$effMass.tsboot
  ret$se <- deffMass
  attr(ret, "class") <- c("effectivemass", "lanczos", class(ret))
  return(invisible(ret))
}

#' @title Block-Lanczos solver
#' 
#' @description
#'   blub ...
#'
#' @param cf Numeric vector (this will generally be a correlation function or a bootstrap sample thereof).
#' @param Time time extent of the lattice. 
#' @param Neig Integer. The number of eigenvectors to be used in the Lanczos analysis
#' @param element.order specifies how to fit the \code{n} linearly ordered
#' single correlators into the correlator matrix.
#' \code{element.order=c(1,2,3,4)} leads to a matrix
#' \code{matrix(cf[element.order], nrow=2)}.
#' @return
#'   tbw
#' 
#' @family lanczos
block_lanczos.solve <- function(cf, Time, Neig=1, element.order) {
  Cor <- cf
  Thalf <- Time/2

  N <- Time/2+1

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


  ## container for the eigenvalues and overlaps per m
  evs <- data.frame()

  ## container storing the approximation to the transfer matrix T
  M <- c()


  #Taking the real symmetric part of the correlator
  cM00 <- 0.5*matrix(Cor[ii], nrow=matrix.size, ncol=matrix.size)
  cM00 <- cM00 + t(cM00)

  #Normalize the correlator
  for (i in 0:(N-1)){
    cM0 <- 0.5*matrix(Cor[ii+i], nrow=matrix.size, ncol=matrix.size)
    cM0 <- cM0 + t(cM0)
    cM1 <- matrix(0,nrow=matrix.size, ncol=matrix.size)
    for (linea in 1:matrix.size){
      for (lineb in 1:matrix.size){
        cM1[linea, lineb] <- cM0[linea,lineb]/sqrt( cM00[linea,linea]* cM00[lineb,lineb])
      }
    }
    Cor[ii+i]<- cM1
  }
	   
  #Using the symmetric convention 
  #TBD: check other conventions
  gamma1 <- sqrtm(matrix(Cor[ii], nrow=matrix.size, ncol=matrix.size))
  beta1  <- gamma1

  for(m in c(1:(N/2))) {

    Aj <- vector("list",N-1)

    for (i in 1:(N-1)){
      cM <- matrix(Cor[ii+i], nrow=matrix.size, ncol=matrix.size)
      Aj[[i]] <- solve(beta1)%*%cM%*%solve(gamma1)
    }

    Aj <- do.call(c, Aj)


    Bj <- rep(0, times=length(Aj))
    Bjp1 <- Bj
    Gj <- Bj
    Gjp1 <- Gj
    Ajm1 <- Bj
    Ajp1 <- Bj

    alpha <- rep(NA, times=matrix.size*matrix.size*N)
    beta <- alpha
    gamma <- alpha

    alpha[1:Ncor] <- Aj[1:Ncor]
    beta[1:Ncor]  <- rep(0,Ncor)
    gamma[1:Ncor] <- rep(0,Ncor)
    if(m > 1) {
      for(j in c(1:(m-1))) {
        ## eq.(53) suppl. mat.
        Aj_2 <- matrix(Aj[(1+Ncor):(2*Ncor)],ncol=matrix.size,nrow=matrix.size)
        alpha_j <- matrix(alpha[((j-1)*Ncor+1):(j*Ncor)],nrow=matrix.size,ncol=matrix.size)
        alphaj_squared <- alpha_j%*%alpha_j
 
        gamma_j <- matrix(gamma[((j-1)*Ncor+1):(j*Ncor)],nrow=matrix.size,ncol=matrix.size)

        beta_j <- matrix(beta[((j-1)*Ncor+1):(j*Ncor)],nrow=matrix.size,ncol=matrix.size)

        delta_jp1 <- Aj_2-alphaj_squared-gamma_j%*%beta_j
        delta_jp1 <- Re(delta_jp1)
	gamma_jp1 <- sqrtm(delta_jp1)
	beta_jp1  <- gamma_jp1

	#Using symmetric convention
        gamma[(j*Ncor+1):((j+1)*Ncor)] <- gamma_jp1

        beta[(j*Ncor+1):((j+1)*Ncor)] <- beta_jp1 #diag(1,ncol=matrix.size,nrow=matrix.size)


        ## below eq.(28) suppl. mat.
        kmax <- 2*(m-j) + 1
        ## B_{j+1}^k and G_{j+1}^k

        Gjp1 <- vector("list", kmax)  # Preallocate list to store results

        Bjp1 <- vector("list", kmax)  # Preallocate list to store results

        for (k in c(1:kmax)){

          idx_next     <- (k * Ncor + 1):((k + 1) * Ncor)
          idx_curr     <- ((k - 1) * Ncor + 1):(k * Ncor)


          Aj_tp1 <- matrix(Aj[idx_next],ncol=matrix.size,nrow=matrix.size)

          Aj_t <- matrix(Aj[idx_curr],ncol=matrix.size,nrow=matrix.size)
          Bj_t <- matrix(Bj[idx_curr],ncol=matrix.size,nrow=matrix.size)
          Gj_t <- matrix(Gj[idx_curr],ncol=matrix.size,nrow=matrix.size)

          Gjp1[[k]] <- solve(beta_jp1) %*%( Aj_tp1-alpha_j%*%Aj_t-gamma_j%*%Bj_t)
          Bjp1[[k]] <- (Aj_tp1-Aj_t%*%alpha_j-Gj_t%*%beta_j) %*% solve(gamma_jp1)

        }

        Gjp1 <- do.call(c, Gjp1)
        Bjp1 <- do.call(c, Bjp1)


	Ajp1 <- vector("list", kmax)  # Preallocate list to store results

        for (k in seq_len(kmax)) {
          idx_next     <- (k * Ncor + 1):((k + 1) * Ncor)
          idx_nextnext <- ((k + 1) * Ncor + 1):((k + 2) * Ncor)
          idx_curr     <- ((k - 1) * Ncor + 1):(k * Ncor)

          # Extract matrices efficiently
          Aj_tp1 <- matrix(Aj[idx_next], ncol = matrix.size, nrow = matrix.size)
          Bj_tp1 <- matrix(Bj[idx_next], ncol = matrix.size, nrow = matrix.size)
          Gj_tp1 <- matrix(Gj[idx_next], ncol = matrix.size, nrow = matrix.size)

          Aj_tp2 <- matrix(Aj[idx_nextnext], ncol = matrix.size, nrow = matrix.size)

          Aj_t   <- matrix(Aj[idx_curr], ncol = matrix.size, nrow = matrix.size)
          Ajm1_t <- matrix(Ajm1[idx_curr], ncol = matrix.size, nrow = matrix.size)
          Bj_t   <- matrix(Bj[idx_curr], ncol = matrix.size, nrow = matrix.size)
          Gj_t   <- matrix(Gj[idx_curr], ncol = matrix.size, nrow = matrix.size)

          # Compute the result using matrix operations
          Ajp1[[k]] <- solve(beta_jp1) %*% ( Aj_tp2 - (alpha_j %*% Aj_tp1 + Aj_tp1 %*% alpha_j) +
		  alpha_j %*% Aj_t %*% alpha_j +
	      	  gamma_j %*% Ajm1_t %*% beta_j -
	      	  (gamma_j %*% Bj_tp1 + Gj_tp1 %*% beta_j) +
	      	  gamma_j %*% Bj_t %*% alpha_j +
	      	  alpha_j %*% Gj_t %*% beta_j ) %*% solve(gamma_jp1)
	}
	
	Ajp1 <- do.call(c, Ajp1)

        alpha[(j*Ncor+1):((j+1)*Ncor)] <- Ajp1[1:Ncor]

        ## don't do the copying if not needed
        if(j == m-1) break
        Ajm1 <- Aj
        Aj <- Ajp1
        Bj <- Bjp1
        Gj <- Gjp1
      }
    }
    if(m == 1) {
      ## ndimxndim case
      M <- matrix(alpha[1:Ncor],ncol=matrix.size,nrow=matrix.size)
      eigens <- try(eigen(M, symmetric=FALSE, only.values = FALSE, EISPACK = FALSE), TRUE)
      if(inherits(eigens, "try-error")) {
        warning("eigen failed in lanczos.solvesasaaaaaaaaa\n")
      }
      else {
	 overlap_ZL <- NULL
         overlap_ZR <- NULL
	 checkeigenvalue <- NULL
	 delta_ZCV <- NULL
	 for (line in 1:length(eigens$values)){
		 delta_ZCV <- c(delta_ZCV, sum(abs(eigens$vector[,line]*solve(eigens$vector)[line,])))
	 }
         evstmp <- data.frame(eigvalues=eigens$values,
			      delta_ZCV=delta_ZCV,
                              m=rep(1, length(delta_ZCV)))
	 evstmp <- filter(evstmp, Im(eigvalues)==0)
	 evstmp <- filter(evstmp, eigvalues<1,eigvalues>0)
         newevstmp <-  evstmp[order(evstmp$eigvalues,decreasing=TRUE),]
         newevstmp$index <- seq(1,length(delta_ZCV))
	 newevstmp <- filter(newevstmp, index<=Neig)
         evs <- rbind(evs, newevstmp)
      }
    }
    else{
      Mnew <-matrix(0, ncol=m*matrix.size,nrow=m*matrix.size)

      Mnew[1:((m-1)*matrix.size),1:((m-1)*matrix.size)] <- M
      Mnew[((m-1)*matrix.size+1):(m*matrix.size),((m-1)*matrix.size+1):(m*matrix.size)] <- matrix(alpha[((m-1)*Ncor+1):(m*Ncor)],ncol=matrix.size,nrow=matrix.size)
      Mnew[((m-2)*matrix.size+1):((m-1)*matrix.size),((m-1)*matrix.size+1):(m*matrix.size)] <- matrix(gamma[((m-1)*Ncor+1):(m*Ncor)],ncol=matrix.size,nrow=matrix.size)
      Mnew[((m-1)*matrix.size+1):(m*matrix.size),((m-2)*matrix.size+1):((m-1)*matrix.size)] <- matrix(beta[((m-1)*Ncor+1):(m*Ncor)],ncol=matrix.size,nrow=matrix.size)

      ## eigensolve and extract the lowest eigenvalues
      eigens <- try(eigen(Mnew, symmetric=FALSE, only.values = FALSE, EISPACK = FALSE), TRUE)
      if(inherits(eigens, "try-error")) {
        warning("eigen failed in lanczos.solve\n")
      }
      else {
         eigen_inv <- solve(eigens$vectors)
	 delta_ZCV <- Re(colSums(eigens$vectors[1:matrix.size, ] * t(eigen_inv[,1:matrix.size ])))
	
	 evstmp <- data.frame(
			      eigvalues=eigens$values,
                              delta_ZCV=delta_ZCV,
                              m=rep(1, length(delta_ZCV)))
         evstmp <- filter(evstmp, abs(Im(eigvalues))<1e-4)
         evstmp <- filter(evstmp, Re(eigvalues)<1,Re(eigvalues)>0)
	 evstmp$eigvalues <- Re(evstmp$eigvalues)

         newevstmp <-  evstmp[order(evstmp$eigvalues,decreasing=TRUE),]

         #Filtering out eigenvectors with significant overlap with the first lanzos vector
         newevstmp <- filter(newevstmp, delta_ZCV>0.1)
         newevstmp$index <- seq(1,nrow(newevstmp))
         newevstmp <- filter(newevstmp, index<=Neig)
         evs <- rbind(evs, newevstmp)
      }
      M <- Mnew
    }
  }
  return(evs)
}
