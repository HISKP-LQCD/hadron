## sum_i (a_i cosh(m_i*t))
## a_i are the amplitudes, m_i the masses and t is a vector of times
sum.cosh <- function(masses, amplitudes, t){
	colSums(amplitudes*cosh(masses%o%t))
}

## extract the effective mass from the sum of coshs
cosh.to.effmass <- function(masses, amplitudes, t, Thalf, type){
	if(type == "solve"){
		sapply(t, FUN=function(x) {invcosh(sum.cosh(masses, amplitudes, x-Thalf)/sum.cosh(masses, amplitudes, x+1-Thalf), timeextent=2*Thalf, t=x+1)})
	}else if(type == "acosh"){
		acosh((sum.cosh(masses, amplitudes, t-Thalf-1)+sum.cosh(masses, amplitudes, t-Thalf+1))/(2*sum.cosh(masses, amplitudes, t-Thalf)))
	}
	## the other types might be included as well...
}

## this fit-function works completely analogous to fit.effecitvemass
## but it fits the correlator directly to a sum of cosh-terms
fit.cosh <- function(effMass, t1, t2, useCov=FALSE, replace.na=TRUE, m.init=0.01, par, n.cosh=2, every, ...){
	stopifnot(inherits(effMass, 'effectivemass'))

	Thalf = effMass$Time/2
	if(effMass$cf$symmetrised) {
		tmax <- Thalf
	}else{
		tmax <- effMass$Time-1
	}

	stopifnot(t1 < t2)
	stopifnot(0 <= t1)
	stopifnot(t2 <= tmax)
	
	if(missing(par)){
		## some initial values, probably not optimal...
		masses = c(0:(n.cosh-1))
		masses = masses+m.init
		amplitudes = effMass$cf$cf0[t1+1]/n.cosh/cosh(masses*(t1-Thalf))
	}else{
		masses = par[1:n.cosh]
		amplitudes = par[(n.cosh+1):(2*n.cosh)]
	}

	cf.save <- effMass$cf$cf.tsboot$t
	cf0.save <- effMass$cf$cf0

	## extract the relevant time slices
	ii <- c((t1+1):(t2+1))
	if(!missing(every)){
		ii <- ii[which(ii%%every == 0)]
	}
	if(any(is.na(cf0.save[ii]))){
		ii.na <- which(is.na(cf0.save[ii]))
		ii <- ii[-ii.na]
	}
	ii.remove <- c()
	if(replace.na && any(is.na(cf.save))) {
		for(k in ii) {
		  if(any(is.na(cf.save[,k]))) {
			## indices of NA elements
			jj <- which(is.na(effMass$t[,k]))
			## if there are more NA elements than non-NA elements, the replacement
			## below will fail
			if( length(cf.save[-jj, k]) > length(jj) ) {
			  ## random indices in the non-NA elements of t
			  rj <- sample.int(n=length(cf.save[-jj, k]), size=length(jj), replace=FALSE)
			  ## replace
			  cf.save[jj, k] <- cf.save[-jj, k][rj]
			}
			else {
			  ## so we remove this column from the analysis below
			  ii.remove <- c( ii.remove, which( ii == k ) )
			}
		  }
		}
	}
	if( length( ii.remove ) > 0 ) {
		## remove the columns that should be excluded from the fit below
		ii <- ii[ -ii.remove ]
		cat("Due to NAs we have removed the time slices", ii.remove-1, " from the fit\n")
	}
	tt = ii-1-Thalf

	## function to fit
	sum.cosh.fit <- function(par, x) {
		sum.cosh(par[1:n.cosh], abs(par[(n.cosh+1):(2*n.cosh)]), x)
	}
	## corresponding Jacobian
	sum.cosh.jac <- function(par, x, ...) {
		df.dm <- x * t(abs(par[(n.cosh+1):(2*n.cosh)])* sinh(par[1:n.cosh]%o%x)) 
		df.da <- t(sign(par[(n.cosh+1):(2*n.cosh)]) * cosh(par[1:n.cosh]%o%x))
		return(cbind(df.dm, df.da))
	}

	dy <- effMass$cf$tsboot.se[ii]
	fit.res <- bootstrap.nlsfit(fn=sum.cosh.fit, gr=sum.cosh.jac, par.guess=c(masses, amplitudes), sim="direct", boot.R=effMass$cf$boot.R, y=cf0.save[ii], dy=dy, x=tt, bsamples=cf.save[,ii], useCov=useCov, use.minpack.lm=FALSE, ...)

	opt.res <- fit.res$t0[1:(2*n.cosh)]
	opt.res[1:n.cosh] <- abs(opt.res[1:n.cosh])
	opt.res[(n.cosh+1):(2*n.cosh)] <- abs(opt.res[(n.cosh+1):(2*n.cosh)][order(opt.res[1:n.cosh])])
	opt.res[1:n.cosh] <- opt.res[1:n.cosh][order(opt.res[1:n.cosh])]

    ## now we bootstrap the fit
	massfit.tsboot <- abs(fit.res$t[,1:(2*n.cosh)])
	if(n.cosh >= 2){
		massfit.tsboot[,(n.cosh+1):(2*n.cosh)] <- t(apply(massfit.tsboot, FUN=function(res) {res[(n.cosh+1):(2*n.cosh)][order(res[1:n.cosh])]}, MARGIN=1))
		massfit.tsboot[,1:n.cosh] <- t(apply(massfit.tsboot[,1:n.cosh], FUN=function(res) {res[order(res)]}, MARGIN=1))
	}

	effMass$coshfit <- list()
	effMass$coshfit$t0 <- opt.res
    effMass$coshfit$t <- massfit.tsboot
	effMass$coshfit$se <- apply(massfit.tsboot, 2, sd)
	effMass$coshfit$n.cosh <- n.cosh
	effMass$t1 <- t1
	effMass$coshfit$t1 <- t1
	effMass$t2 <- t2
	effMass$coshfit$t2 <- t2
	effMass$coshfit$ii <- ii
	effMass$coshfit$useCov <- useCov
	effMass$coshfit$invCovMatrix <- fit.res$invCovMatrix
	effMass$coshfit$replace.na <- replace.na
	effMass$dof <- length(ii)-2*n.cosh
	effMass$chisqr <- fit.res$chisqr
	effMass$coshfit$dof <- length(ii)-2*n.cosh
	effMass$coshfit$chisqr <- fit.res$chisqr
	effMass$coshfit$Qval <- fit.res$Qval
	attr(effMass, "class") <- c("coshfit", class(effMass))
	return(invisible(effMass))
}

plot.coshfit <- function(effMass, col.fitline="black", plot.mass=TRUE, plot.corr=FALSE, ...) {
	t.all <- effMass$t.idx
	t <- c(effMass$t1:effMass$t2)
	Thalf <- effMass$Time/2
	n.cosh <- effMass$coshfit$n.cosh

	pcol <- col2rgb(col.fitline,alpha=TRUE)/255
	pcol[4] <- 0.65
	pcol <- rgb(red=pcol[1],green=pcol[2],blue=pcol[3],alpha=pcol[4])

	if(plot.mass){
		op <- options()
		options(warn=-1)
		plotwitherror(x=t.all-1, y=effMass$effMass[t.all], dy=effMass$deffMass[t.all], ...)
		options(op)

		if(!is.null(effMass$coshfit)){
			Y <- cosh.to.effmass(effMass$coshfit$t0[1:n.cosh], effMass$coshfit$t0[(n.cosh+1):(2*n.cosh)], t, Thalf, type=effMass$type)
			Y.boot <- apply(effMass$coshfit$t, FUN=function(x) {cosh.to.effmass(x[1:n.cosh], x[(n.cosh+1):(2*n.cosh)], t, Thalf, type=effMass$type)}, MARGIN=1)
			se <- apply(Y.boot, MARGIN=1, FUN=sd, na.rm=TRUE)

			## plot it
			polyval <- c(Y+se, rev(Y-se))
			polygon(x=c(t, rev(t)), y=polyval, col=pcol, lty=0, lwd=0.001, border=pcol)

			## plot the fitted curve on top
			lines(x=t, y=Y, col=col.fitline, lwd=1.3)
		}
	}

	if(plot.corr){
		plot(effMass$cf, ...)

		if(!is.null(effMass$coshfit)){
			Y <- sum.cosh(effMass$coshfit$t0[1:n.cosh], effMass$coshfit$t0[(n.cosh+1):(2*n.cosh)], t-Thalf)
			Y.boot <- apply(effMass$coshfit$t, FUN=function(x) {sum.cosh(x[1:n.cosh], x[(n.cosh+1):(2*n.cosh)], t-Thalf)}, MARGIN=2)
			se <- apply(Y.boot, MARGIN=1, FUN=sd, na.rm=TRUE)

			## plot it
			polyval <- c(Y+se, rev(Y-se))
			polygon(x=c(t, rev(t)), y=polyval, col=pcol, lty=0, lwd=0.001, border=pcol)

			## plot the fitted curve on top
			lines(x=t, y=Y, col=col.fitline, lwd=1.3)
		}
	}
}

summary.coshfit <- function(effMass, verbose=FALSE) {
  cat("\n ** Result of", effMass$coshfit$n.cosh, "-fold cosh-fit **\n\n")
  cat("no. measurements\t=\t", effMass$N, "\n")
  cat("boot.R\t=\t", effMass$boot.R, "\n")
  cat("boot.l\t=\t", effMass$boot.l, "\n")
  cat("Time extend\t=\t", effMass$Time, "\n")
  cat("NA count in fitted bootstrap samples:\t", length(which(is.na(effMass$cf$cf.tsboot$t[,effMass$coshfit$ii]))),
      "(",100*length(which(is.na(effMass$cf$cf.tsboot$t[,effMass$coshfit$ii])))/ length(effMass$cf$cf.tsboot$t[,effMass$coshfit$ii]), "%)\n")
  cat("time range from", effMass$coshfit$t1, " to ", effMass$coshfit$t2, "\n")
  if(verbose) {
    cat("\nmasses:\n")
    print(data.frame(m = effMass$coshfit$t0[1:effMass$coshfit$n.cosh], dm = effMass$coshfit$se[1:effMass$coshfit$n.cosh]))
    cat("\namplitudes:\n")
    print(data.frame(a = effMass$coshfit$t0[(effMass$coshfit$n.cosh+1):(2*effMass$coshfit$n.cosh)], da = effMass$coshfit$se[(effMass$coshfit$n.cosh+1):(2*effMass$coshfit$n.cosh)]))
	cat("\ncorrelations:\n")
	fit.cor <- cor(effMass$coshfit$t[,1:(2*effMass$coshfit$n.cosh)], effMass$coshfit$t[,1:(2*effMass$coshfit$n.cosh)])
	print(fit.cor)
	cat("\n")
  }else{
	  cat("m\t=\t", effMass$coshfit$t0[1], "\n")
	  cat("dm\t=\t", effMass$coshfit$se[1], "\n")
  }
  cat("chisqr\t=\t", effMass$coshfit$chisqr, "\n")
  cat("dof\t=\t", effMass$coshfit$dof, "\n")
  cat("chisqr/dof=\t",
      effMass$coshfit$chisqr/effMass$coshfit$dof, "\n")
  cat("Quality of the fit (p-value):",   effMass$coshfit$Qval, "\n")

}
