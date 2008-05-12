print.b1fit <- function(fit) {
  summary(fit)
}

summary.b1fit <- function(fit) {
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

  if(fit$no.masses == 1) {
    cat("mb1    = ", fit.mass, "\n")
  }

  cat("\nstate  =", seq(1,fit$no.masses), "\n", sep="\t")
  cat("masses =", abs(fit$fitresult$par[ii+fit$matrix.size]), "\n", sep="\t")
  cat("B_L    =", fit$fitresult$par[ii], "\n", sep="\t")
  cat("B_F    =", fit$fitresult$par[ii+1], "\n", sep="\t")

  if(!is.null(fit$uwerrresultmb1)) {
    cat("\n--- Autocorrelation analysis for mb1 ---\n")
    cat("\nS        = ", fit$uwerrresultmb1$S, "\n")
    cat("mb1      = ", fit$uwerrresultmb1$value, "\n")
    cat("dmb1     = ", fit$uwerrresultmb1$dvalue, "\n")
    cat("ddmb1    = ", fit$uwerrresultmb1$ddvalue, "\n")
    cat("tauint   = ", fit$uwerrresultmb1$tauint, "\n")
    cat("dtauint  = ", fit$uwerrresultmb1$dtauint, "\n")
    cat("Wopt     = ", fit$uwerrresultmb1$Wopt, "\n")
    if(fit$uwerrresultmb1$R>1) {
      cat("Qval     =", fit$uwerrresultmb1$Qval, "\n")
    }
    if(fit$no.masses > 1) {
      cat("\n--- Autocorrelation analysis for mb1 ---\n")
      cat("\nS        = ", fit$uwerrresultmb12$S, "\n")
      cat("mb12     = ", fit$uwerrresultmb12$value, "\n")
      cat("dmb12    = ", fit$uwerrresultmb12$dvalue, "\n")
      cat("ddmb12   = ", fit$uwerrresultmb12$ddvalue, "\n")
      cat("tauint2  = ", fit$uwerrresultmb12$tauint, "\n")
      cat("dtauint2 = ", fit$uwerrresultmb12$dtauint, "\n")
      cat("Wopt2    = ", fit$uwerrresultmb12$Wopt, "\n")
      if(fit$uwerrresultmb12$R>1) {
        cat("Qval     =", fit$uwerrresultmb12$Qval, "\n")
      }
    }
  }

  if(!is.null(fit$boot)) {
    cat("--- Bootstrap analysis  ---\n")
    cat("---", fit$boot$R, "samples  ---\n")
    cat("          mean        -err           +err            stderr        bias\n")
    for(no in 1:fit$no.masses) {
      index <- (no-1)*(fit$matrix.size+1)+1
      b.ci <- boot.ci(fit$boot, type = c("norm"), index=index)
      cat("mb1[",no,"] = ", abs(fit$boot$t0[index]), "(", (b.ci$normal[1,2]-fit$boot$t0[index])/1.96
          , ",", -(fit$boot$t0[index]-b.ci$normal[1,3])/1.96, ")", sd(fit$boot$t[,index]),
          mean(fit$boot$t[,index])-fit$boot$t0[index],"\n")
    }
  }
  if(!is.null(fit$tsboot)) {
    cat("\n--- Bootstrap analysis  with blocking ---\n")
    cat("---", fit$tsboot$R, "samples  ---\n")
    cat("--- block size", fit$tsboot$l, "---\n")
    for(no in 1:fit$no.masses) {
      index <- (no-1)*(fit$matrix.size+1)+1
      tsb.ci <- boot.ci(fit$tsboot, type = c("norm"), index=index)
      cat("mb1[",no,"] = ", fit$tsboot$t0[index], "(", (tsb.ci$normal[1,2]-fit$tsboot$t0[index])/1.96
          , ",", -(fit$tsboot$t0[index]-tsb.ci$normal[1,3])/1.96, ")", sd(fit$tsboot$t[,index]),
          mean(fit$tsboot$t[,index])-fit$tsboot$t0[index], "\n")
    }
  }
  if(!is.null(fit$variational.masses)) {
    cat("\n--- Variational analysis ---\n")
    cat("masses:", fit$variational.masses, "\n")
  }
  
}
