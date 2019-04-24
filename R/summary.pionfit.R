#' pionfit print
#'
#' @param x class pionfit. object to print
#' @param ... additional parameters to be passed on
print.pionfit <- function (x, ...) {
  summary(x, ...)
}

#' pionfit summary
#'
#' @param object class pionfit. object to summarise
#' @param ... additional parameters to be passed on
summary.pionfit <- function (object, ...) {
  fit <- object
  kappa <- fit$kappa
  mu <- fit$mu
  t1 <- fit$t1
  t2 <- fit$t2
  fit.mass <- abs(fit$fitresult$par[fit$matrix.size+1])
  fit.fpi <- 2*kappa*2*mu/sqrt(2)*abs(fit$fitresult$par[1])/sqrt(fit.mass^3)
  fit.chisqr <- fit$fitresult$value
  fit.dof <- length(fit$fitdata$t)-length(fit$fitresult$par)
  
  cat("mu     = ", mu, "\n")
  cat("kappa  = ", kappa, "\n")
  cat("Nr of measurements = ", fit$N, "\n")
  cat("No of replica = ", length(fit$nrep), "\n")
  cat("no or measurements per replicum: ", fit$nrep, "\n")
  cat("fitrange = ", t1, "-", t2, "\n")
  cat("chi^2    = ", fit.chisqr, "\n")
  cat("dof    = ", fit.dof, "\n")
  cat("chi^2/dof = ", fit.chisqr/fit.dof, "\n")
  
  cat("\nmpi    = ", fit.mass, "\n")
  cat("fpi    = ", fit.fpi, "\n")
  if(fit$matrix.size > 2) {
    cat("mpcac  =", fit.mass*fit$fitresult$par[3]/fit$fitresult$par[1]/2., "\n")
  }
  if(fit$matrix.size > 4) {
    cat("Z_V    =", 2.*mu/fit.mass*fit$fitresult$par[1]/fit$fitresult$par[5], "\n")
  }
  ii <- seq(1, fit$no.masses*(fit$matrix.size+1), by=fit$matrix.size+1)
  cat("\nstate  =", seq(1,fit$no.masses), "\n", sep="\t")
  cat("masses =", abs(fit$fitresult$par[ii+fit$matrix.size]), "\n", sep="\t")
  cat("P_L    =", fit$fitresult$par[ii], "\n", sep="\t")
  cat("P_F    =", fit$fitresult$par[ii+1], "\n", sep="\t")
  if(fit$matrix.size >2) {
    cat("A_L    =", fit$fitresult$par[ii+2], "\n", sep="\t")
    cat("A_F    =", fit$fitresult$par[ii+3], "\n", sep="\t")
  }
  if(fit$matrix.size >4) {
    cat("4_L    =", fit$fitresult$par[ii+4], "\n", sep="\t")
    cat("4_F    =", fit$fitresult$par[ii+5], "\n", sep="\t")
  }

  
  if(!is.null(fit$uwerrresultmps)) {
    cat("\n--- Autocorrelation analysis for m_ps ---\n")
    cat("\nS        = ", fit$uwerrresultmps$S, "\n")
    cat("mps      = ", fit$uwerrresultmps$value[1], "\n")
    cat("dmps     = ", fit$uwerrresultmps$dvalue[1], "\n")
    cat("ddmps    = ", fit$uwerrresultmps$ddvalue[1], "\n")
    cat("tauint   = ", fit$uwerrresultmps$tauint[1], "\n")
    cat("dtauint  = ", fit$uwerrresultmps$dtauint[1], "\n")
    cat("Wopt     = ", fit$uwerrresultmps$Wopt[[1]], "\n")
    if(fit$uwerrresultmps$R>1) {
      cat("Qval     =", fit$uwerrresultmps$Qval[1], "\n")
    }
  }
  if(!is.null(fit$uwerrresultmps2)) {
    cat("\n--- Autocorrelation analysis for m_ps ---\n")
    cat("\nS        = ", fit$uwerrresultmps2$S, "\n")
    cat("mps2     = ", fit$uwerrresultmps2$value[1], "\n")
    cat("dmps     = ", fit$uwerrresultmps2$dvalue[1], "\n")
    cat("ddmps    = ", fit$uwerrresultmps2$ddvalue[1], "\n")
    cat("tauint   = ", fit$uwerrresultmps2$tauint[1], "\n")
    cat("dtauint  = ", fit$uwerrresultmps2$dtauint[1], "\n")
    cat("Wopt     = ", fit$uwerrresultmps2$Wopt[[1]], "\n")
    if(fit$uwerrresultmps2$R>1) {
      cat("Qval     =", fit$uwerrresultmps2$Qval[1], "\n")
    }
  }
  if(!is.null(fit$uwerrresultmps3)) {
    cat("\n--- Autocorrelation analysis for m_ps ---\n")
    cat("\nS        = ", fit$uwerrresultmps3$S, "\n")
    cat("mps3     = ", fit$uwerrresultmps3$value[1], "\n")
    cat("dmps     = ", fit$uwerrresultmps3$dvalue[1], "\n")
    cat("ddmps    = ", fit$uwerrresultmps3$ddvalue[1], "\n")
    cat("tauint   = ", fit$uwerrresultmps3$tauint[1], "\n")
    cat("dtauint  = ", fit$uwerrresultmps3$dtauint[1], "\n")
    cat("Wopt     = ", fit$uwerrresultmps3$Wopt[[1]], "\n")
    if(fit$uwerrresultmps3$R>1) {
      cat("Qval     =", fit$uwerrresultmps3$Qval[1], "\n")
    }
  }
  if(!is.null(fit$uwerrresultfps)) {
    cat("\n--- Autocorrelation analysis for f_ps ---\n")    
    cat("\nS        = ", fit$uwerrresultfps$S, "\n")
    cat("fps      = ", fit$uwerrresultfps$value[1]*2*kappa*2*mu/sqrt(2), "\n")
    cat("dfps     = ", fit$uwerrresultfps$dvalue[1]*2*kappa*2*mu/sqrt(2), "\n")
    cat("ddfps    = ", fit$uwerrresultfps$ddvalue[1]*2*kappa*2*mu/sqrt(2), "\n")
    cat("tauint   = ", fit$uwerrresultfps$tauint[1], "\n")
    cat("dtauint  = ", fit$uwerrresultfps$dtauint[1], "\n")
    cat("Wopt     = ", fit$uwerrresultfps$Wopt[[1]], "\n")
    if(fit$uwerrresultfps$R>1) {
      cat("Qval     =", fit$uwerrresultfps$Qval[1], "\n")
    }
  }
  if(!is.null(fit$uwerrresultmpcac)) {
    cat("\n--- Autocorrelation analysis for m_pcac ---\n")
    cat("\nS        = ", fit$uwerrresultmpcac$S, "\n")
    cat("mpcac    = ", fit$uwerrresultmpcac$value[1], "\n")
    cat("dmpcac   = ", fit$uwerrresultmpcac$dvalue[1], "\n")
    cat("ddmpcac  = ", fit$uwerrresultmpcac$ddvalue[1], "\n")
    cat("tauint   = ", fit$uwerrresultmpcac$tauint[1], "\n")
    cat("dtauint  = ", fit$uwerrresultmpcac$dtauint[1], "\n")
    cat("Wopt     = ", fit$uwerrresultmpcac$Wopt[[1]], "\n")
    if(fit$uwerrresultmpcac$R>1) {
      cat("Qval     =", fit$uwerrresultmpcac$Qval[1], "\n")
    }
  }
  if(!is.null(fit$uwerrresultzv)) {
    cat("\n--- Autocorrelation analysis for Z_V ---\n")
    cat("\nS        = ", fit$uwerrresultzv$S, "\n")
    cat("Z_V      = ", fit$uwerrresultzv$value, "\n")
    cat("dzv      = ", fit$uwerrresultzv$dvalue, "\n")
    cat("ddzv     = ", fit$uwerrresultzv$ddvalue, "\n")
    cat("tauint   = ", fit$uwerrresultzv$tauint, "\n")
    cat("dtauint  = ", fit$uwerrresultzv$dtauint, "\n")
    cat("Wopt     = ", fit$uwerrresultzv$Wopt, "\n")
    if(fit$uwerrresultmpcac$R>1) {
      cat("Qval     =", fit$uwerrresultmpcac$Qval, "\n")
    }
  }
  
  if(!is.null(fit$boot)) {
    cat("--- Bootstrap analysis  ---\n")
    cat("---", fit$boot$R, "samples  ---\n")
    cat("          mean        -err           +err            stderr        bias\n")
    fit$boot.ci <- boot::boot.ci(fit$boot, type = c("norm"), index=1)
    cat("mpi    = ", fit$boot$t0[1], "(", (fit$boot.ci$normal[1,2]-fit$boot$t0[1])/1.96
	, ",", -(fit$boot$t0[1]-fit$boot.ci$normal[1,3])/1.96, ")", sd(fit$boot$t[,1]),
	mean(fit$boot$t[,1])-fit$boot$t0[1],"\n")
    
    fit$boot.ci <- boot::boot.ci(fit$boot, type = c("norm"), index=2)
    cat("fpi    = ", fit$boot$t0[2], "(", (fit$boot.ci$normal[1,2]-fit$boot$t0[2])/1.96
	, ",", -(fit$boot$t0[2]-fit$boot.ci$normal[1,3])/1.96, ")", sd(fit$boot$t[,2]),
	mean(fit$boot$t[,2])-fit$boot$t0[2], "\n")

    if(fit$matrix.size > 2) {
      fit$boot.ci <- boot::boot.ci(fit$boot, type = c("norm"), index=fit$matrix.size+3)
      cat("mpcac  = ", fit$boot$t0[fit$matrix.size+3], "(", (fit$boot.ci$normal[1,2]-fit$boot$t0[fit$matrix.size+3])/1.96
          , ",", -(fit$boot$t0[fit$matrix.size+3]-fit$boot.ci$normal[1,3])/1.96, ")", sd(fit$boot$t[,(fit$matrix.size+3)]),
          mean(fit$boot$t[,(fit$matrix.size+3)])-fit$boot$t0[fit$matrix.size+3], "\n")
    }

    if(fit$matrix.size > 4) {
      fit$boot.ci <- boot::boot.ci(fit$boot, type = c("norm"), index=fit$matrix.size+4)
      cat("Z_V    = ", fit$boot$t0[fit$matrix.size+4], "(", (fit$boot.ci$normal[1,2]-fit$boot$t0[fit$matrix.size+4])/1.96
          , ",", -(fit$boot$t0[fit$matrix.size+4]-fit$boot.ci$normal[1,3])/1.96, ")", sd(fit$boot$t[,(fit$matrix.size+4)]),
          mean(fit$boot$t[,(fit$matrix.size+4)])-fit$boot$t0[fit$matrix.size+4], "\n")      
    }
    
    fit$boot.ci <- boot::boot.ci(fit$boot, type = c("norm"), index=4)
    cat("P_L    = ", fit$boot$t0[4], "(", (fit$boot.ci$normal[1,2]-fit$boot$t0[4])/1.96
	, ",", -(fit$boot$t0[4]-fit$boot.ci$normal[1,3])/1.96, ")", sd(fit$boot$t[,4]),
	mean(fit$boot$t[,4])-fit$boot$t0[4], "\n")
    
    fit$boot.ci <- boot::boot.ci(fit$boot, type = c("norm"), index=5)
    cat("P_F    = ", fit$boot$t0[5], "(", (fit$boot.ci$normal[1,2]-fit$boot$t0[5])/1.96
	, ",", -(fit$boot$t0[5]-fit$boot.ci$normal[1,3])/1.96, ")", sd(fit$boot$t[,5]),
	mean(fit$boot$t[,5])-fit$boot$t0[5],"\n")
  }
  if(!is.null(fit$tsboot)) {
    cat("\n--- Bootstrap analysis with blocking ---\n")
    cat("---", fit$boot$R, "samples  ---\n")
    cat("--- block size", fit$tsboot$l, "---\n")
    fit$tsboot.ci <- boot::boot.ci(fit$tsboot, type = c("norm"), index=1)
    cat("mpi    = ", fit$tsboot$t0[1], "(", (fit$tsboot.ci$normal[1,2]-fit$tsboot$t0[1])/1.96
	, ",", -(fit$tsboot$t0[1]-fit$tsboot.ci$normal[1,3])/1.96, ")", sd(fit$tsboot$t[,1]),
	mean(fit$tsboot$t[,1])-fit$tsboot$t0[1], "\n")
    
    fit$tsboot.ci <- boot::boot.ci(fit$tsboot, type = c("norm"), index=2)
    cat("fpi    = ", fit$tsboot$t0[2], "(", (fit$tsboot.ci$normal[1,2]-fit$tsboot$t0[2])/1.96
	, ",", -(fit$tsboot$t0[2]-fit$tsboot.ci$normal[1,3])/1.96, ")", sd(fit$tsboot$t[,2]),
	mean(fit$tsboot$t[,2])-fit$tsboot$t0[2], "\n")

    if(fit$matrix.size > 2) {
      fit$tsboot.ci <- boot::boot.ci(fit$tsboot, type = c("norm"), index=fit$matrix.size+3)
      cat("mpcac  = ", fit$tsboot$t0[fit$matrix.size+3], "(", (fit$tsboot.ci$normal[1,2]-fit$tsboot$t0[fit$matrix.size+3])/1.96
          , ",", -(fit$tsboot$t0[fit$matrix.size+3]-fit$tsboot.ci$normal[1,3])/1.96, ")", sd(fit$tsboot$t[,(fit$matrix.size+3)]),
          mean(fit$tsboot$t[,(fit$matrix.size+3)])-fit$tsboot$t0[fit$matrix.size+3], "\n")
    }

    if(fit$matrix.size > 4) {
      fit$tsboot.ci <- boot::boot.ci(fit$tsboot, type = c("norm"), index=fit$matrix.size+4)
      cat("Z_V    = ", fit$tsboot$t0[fit$matrix.size+4], "(", (fit$tsboot.ci$normal[1,2]-fit$tsboot$t0[fit$matrix.size+4])/1.96
          , ",", -(fit$tsboot$t0[fit$matrix.size+4]-fit$tsboot.ci$normal[1,3])/1.96, ")", sd(fit$tsboot$t[,(fit$matrix.size+4)]),
          mean(fit$tsboot$t[,(fit$matrix.size+4)])-fit$tsboot$t0[fit$matrix.size+4], "\n")      
    }
    
    fit$tsboot.ci <- boot::boot.ci(fit$tsboot, type = c("norm"), index=4)
    cat("P_L    = ", fit$tsboot$t0[4], "(", (fit$tsboot.ci$normal[1,2]-fit$tsboot$t0[4])/1.96
	, ",", -(fit$tsboot$t0[4]-fit$tsboot.ci$normal[1,3])/1.96, ")", sd(fit$tsboot$t[,4]),
	mean(fit$tsboot$t[,4])-fit$tsboot$t0[4], "\n")
    
    fit$tsboot.ci <- boot::boot.ci(fit$tsboot, type = c("norm"), index=5)
    cat("P_F    = ", fit$tsboot$t0[5], "(", (fit$tsboot.ci$normal[1,2]-fit$tsboot$t0[5])/1.96
	, ",", -(fit$tsboot$t0[5]-fit$tsboot.ci$normal[1,3])/1.96, ")", sd(fit$tsboot$t[,5]),
	mean(fit$tsboot$t[,5])-fit$tsboot$t0[5], "\n")
    
  }

  if(!is.null(fit$variational.masses)) {
    cat("\n--- Variational analysis ---\n")
    cat("masses:", fit$variational.masses, "\n")
  }
  
}
