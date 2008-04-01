print.chiralfit <- function(fit, ...) {
  summary(fit, ...)
}

summary.chiralfit <- function(fit, show.input=FALSE) {
  N <- length(fit$data)
  npar <- length(fit$par)
  if(!is.null(fit$boot.result)) {
    cat("Errors computed using", fit$boot.R, "bootsamples \n")
  }
  cat("chisqr       = ", fit$result$chisqr, "\n")
  cat("dof          = ", fit$result$dof, "\n")
  cat("red. chisqr  = ", fit$result$chisqr/fit$result$dof, "\n")
  if(!is.null(fit$boot.result)) {
    if(fit$fit.l12) {
      cat("l1           = ", fit$result$l1, "+-", sd(fit$boots[,(8+2*N)], na.rm=TRUE), "\n")
      cat("l2           = ", fit$result$l2, "+-", sd(fit$boots[,(9+2*N)], na.rm=TRUE), "\n")
    }
    cat("l3           = ", fit$result$l3, "+-", sd(fit$boots[,(1+2*N)], na.rm=TRUE), "\n")
    cat("l4           = ", fit$result$l4, "+-", sd(fit$boots[,(2+2*N)], na.rm=TRUE), "\n")
    cat("F0           = ", fit$result$F, "+-", sd(fit$boots[,(5+2*N)], na.rm=TRUE), "MeV \n")
    cat("B0           = ", fit$result$B0, "+-", sd(fit$boots[,(6+2*N)], na.rm=TRUE), "MeV \n")
    cat("mN0          = ", fit$result$mN0, "+-", sd(fit$boots[,(7+2*N)], na.rm=TRUE), "MeV \n")
    cat("mN           = ", fit$result$mN, "+-", sd(fit$boots[,(8+2*N)], na.rm=TRUE), "MeV \n")
    cat("c1           = ", fit$result$c1, "+-", sd(fit$boots[,(8+2*N+6+2*N)], na.rm=TRUE), "\n")
    cat("gA           = ", fit$result$gA, "+-", sd(fit$boots[,(8+2*N+7+2*N)], na.rm=TRUE), "\n")
    cat("mu_phys      = ", fit$result$mu.phys[1], "+-", sd(fit$boots[,1], na.rm=TRUE), "MeV \n\n")
    for(i in 1:N) {
      cat("lattice spacing", i, ":\n")
      cat("lattice spacing at r0/a = ",fit$r0data$r0[i], ": a = ", fit$result$a[i], "+-",
          sd(fit$boots[,(N+i)], na.rm=TRUE),"fm \n")
      cat("            fitted r0/a = ", fit$par[4+i], "\n")
      cat("            fitted ZP   = ", fit$par[4+N+i], "\n")
      if(show.input) {
        cat("Raw data used:\n")
        print(fit$data[[i]])
        cat("Raw r0 data used:\n")
        cat(fit$r0data$r0[i], fit$r0data$dr0[i], "\n")
        cat("Raw ZP data used:\n")
        cat(fit$ZPdata$ZP[i], fit$ZPdata$dZP[i], "\n")
      }
      cat("\n")
    }
  }
  else {
    if(fit$fit.l12) {
      cat("l1           = ", fit$result$l1, "\n")
      cat("l2           = ", fit$result$l2, "\n")
    }
    cat("l3           = ", fit$result$l3, "\n")
    cat("l4           = ", fit$result$l4, "\n")
    cat("F0           = ", fit$result$F, "MeV \n")
    cat("B0           = ", fit$result$B0, "MeV \n")
    cat("mN0          = ", fit$result$mN0, "MeV \n")
    cat("mN           = ", fit$result$mN, "MeV \n")
    cat("c1           = ", fit$result$c1, "\n")
    cat("gA           = ", fit$result$gA, "\n")
    cat("mu_phys      = ", fit$result$mu.phys[1], "MeV \n")
    for(i in 1:N) {
      cat("lattice spacing", i, ":\n")
      cat("lattice spacing at r0/a = ",fit$r0data$r0[i], ": a = ", fit$result$a[i], "fm \n")
      cat("            fitted r0/a = ", fit$par[4+i], "\n")
      cat("            fitted ZP   = ", fit$par[4+N+i], "\n")
      if(show.input) {
        cat("Raw data used:\n")
        print(fit$data[[i]])
        cat("Raw r0 data used:\n")
        cat(fit$r0data$r0[i], fit$r0data$dr0[i], "\n")
        cat("Raw ZP data used:\n")
        cat(fit$ZPdata$ZP[i], fit$ZPdata$dZP[i], "\n")
      }
      cat("\n")
    }
  }
}
