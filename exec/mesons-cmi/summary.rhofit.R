print.rhofit <- function(fit) {
  summary(fit)
}

summary.rhofit <- function(fit) {
  kappa <- fit$kappa
  mu <- fit$mu
  t1 <- fit$t1
  t2 <- fit$t2
  ij <- seq(1, fit$no.masses*(fit$matrix.size+1), by=fit$matrix.size+1)
  sortindex <- order(abs(fit$fitresult$par[ij+fit$matrix.size]))
  ii <- ij[sortindex]

  fit.mass <- abs(fit$fitresult$par[ii[1] + fit$matrix.size])
  fit.chisqr <- fit$fitresult$value
  fit.dof <- length(fit$fitdata$t)-length(fit$fitresult$par)
  
  cat("mu     = ", mu, "\n")
  cat("kappa  = ", kappa, "\n")
  cat("No of measurements = ", fit$N, "\n")
  cat("No of replica = ", length(fit$nrep), "\n")
  cat("no or measurements per replicum: ", fit$nrep, "\n")
  cat("fitrange = ", t1, "-", t2, "\n")
  cat("chi^2    = ", fit.chisqr, "\n")
  cat("dof    = ", fit.dof, "\n")
  cat("chi^2/dof = ", fit.chisqr/fit.dof, "\n")
  
  cat("mv     = ", fit.mass, "\n")
  ii <- seq(1, fit$no.masses*(fit$matrix.size+1), by=fit$matrix.size+1)

  cat("\nstate  =", seq(1,fit$no.masses), "\n", sep="\t")
  cat("masses =", abs(fit$fitresult$par[ii+fit$matrix.size]), "\n", sep="\t")
  cat("4_L    =", fit$fitresult$par[ii], "\n", sep="\t")
  cat("4_F    =", fit$fitresult$par[ii+1], "\n", sep="\t")
  if(fit$matrix.size >2) {
    cat("A_L    =", fit$fitresult$par[ii+2], "\n", sep="\t")
    cat("A_F    =", fit$fitresult$par[ii+3], "\n", sep="\t")
  }
  if(fit$matrix.size >4) {
    cat("V_L    =", fit$fitresult$par[ii+4], "\n", sep="\t")
    cat("V_F    =", fit$fitresult$par[ii+5], "\n", sep="\t")
  }

  if(!is.null(fit$uwerrresultmv)) {
    cat("\n--- Autocorrelation analysis for m_v ---\n")
    cat("\nS        = ", fit$uwerrresultmv$S, "\n")
    cat("mv       = ", fit$uwerrresultmv$value, "\n")
    cat("dmv      = ", fit$uwerrresultmv$dvalue, "\n")
    cat("ddmv     = ", fit$uwerrresultmv$ddvalue, "\n")
    cat("tauint   = ", fit$uwerrresultmv$tauint, "\n")
    cat("dtauint  = ", fit$uwerrresultmv$dtauint, "\n")
    cat("Wopt     = ", fit$uwerrresultmv$Wopt, "\n")
    if(fit$uwerrresultmv$R>1) {
      cat("Qval     =", fit$uwerrresultmv$Qval, "\n")
    }
    if(fit$no.masses > 1) {
      cat("\n--- Autocorrelation analysis for m_v2 ---\n")
      cat("\nS        = ", fit$uwerrresultmv2$S, "\n")
      cat("mv2      = ", fit$uwerrresultmv2$value, "\n")
      cat("dmv2     = ", fit$uwerrresultmv2$dvalue, "\n")
      cat("ddmv2    = ", fit$uwerrresultmv2$ddvalue, "\n")
      cat("tauint2  = ", fit$uwerrresultmv2$tauint, "\n")
      cat("dtauint2 = ", fit$uwerrresultmv2$dtauint, "\n")
      cat("Wopt2    = ", fit$uwerrresultmv2$Wopt, "\n")
      if(fit$uwerrresultmv2$R>1) {
        cat("Qval     =", fit$uwerrresultmv2$Qval, "\n")
      }
    }
  }

  if(!is.null(fit$mv.boot)) {
    cat("--- Bootstrap analysis  ---\n")
    cat("---", fit$mv.boot$R, "samples  ---\n")
    cat("          mean        -err           +err            stderr        bias\n")
    for(no in 1:fit$no.masses) {
      index <- no*(fit$matrix.size+1)
      b.ci <- boot.ci(fit$mv.boot, type = c("norm"), index=index)
      cat("mv[",no,"] = ", abs(fit$mv.boot$t0[index]), "(", (b.ci$normal[1,2]-fit$mv.boot$t0[index])/1.96
          , ",", -(fit$mv.boot$t0[index]-b.ci$normal[1,3])/1.96, ")", sd(fit$mv.boot$t[,index]),
          mean(fit$mv.boot$t[,index])-fit$mv.boot$t0[index],"\n")
    }
  }
  if(!is.null(fit$mv.tsboot)) {
    cat("\n--- Bootstrap analysis  with blocking ---\n")
    cat("---", fit$mv.tsboot$R, "samples  ---\n")
    cat("--- block size", fit$mv.tsboot$l, "---\n")
    for(no in 1:fit$no.masses) {
      index <- no*(fit$matrix.size+1)
      tsb.ci <- boot.ci(fit$mv.tsboot, type = c("norm"), index=index)
      cat("mv[",no,"] = ", fit$mv.tsboot$t0[index], "(", (tsb.ci$normal[1,2]-fit$mv.tsboot$t0[index])/1.96
          , ",", -(fit$mv.tsboot$t0[index]-tsb.ci$normal[1,3])/1.96, ")", sd(fit$mv.tsboot$t[,index]),
          mean(fit$mv.tsboot$t[,index])-fit$mv.tsboot$t0[index], "\n")
    }
  }
  if(!is.null(fit$variational.masses)) {
    cat("\n--- Variational analysis ---\n")
    cat("masses:", fit$variational.masses, "\n")
  }
  
}
