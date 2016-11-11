print.ofit <- function(fit) {
  summary(fit)
}

summary.ofit <- function(fit) {
  kappa <- fit$kappa
  mu <- fit$mu
  t1 <- fit$t1
  t2 <- fit$t2
  fit.mass <- abs(fit$fitresult$par[3])
  fit.fpi <- 2*kappa*2*mu/sqrt(2)*abs(fit$fitresult$par[1])/sqrt(fit.mass^3)
  fit.chisqr <- fit$fitresult$value
  fit.dof <- length(fit$fitdata$t)-length(fit$fitresult$par)
  
  cat("mu        = ", mu, "\n")
  cat("kappa     = ", kappa, "\n")
  cat("Nr of measurements = ", fit$N, "\n")
  cat("No of replica = ", length(fit$nrep), "\n")
  cat("no or measurements per replicum: ", fit$nrep, "\n")
  cat("fitrange  = ", t1, "-", t2, "\n")
  cat("chi^2     = ", fit.chisqr, "\n")
  cat("dof       = ", fit.dof, "\n")
  cat("chi^2/dof = ", fit.chisqr/fit.dof, "\n")
  
  cat("\nmass     = ", fit.mass, "\n")
  cat("mpcac    = ", fit.mass*fit$fitresult$par[2]/fit$fitresult$par[1]/2., "\n")
  cat("\nP_L      = ", fit$fitresult$par[1], "\n")
  cat("A_L      = ", fit$fitresult$par[2], "\n")
  
  if(!is.null(fit$uwerrresultmps)) {
    cat("\n--- Autocorrelation analysis for m_ps ---\n")
    cat("\nS        = ", fit$uwerrresultmps$S, "\n")
    cat("mass     = ", fit$uwerrresultmps$value, "\n")
    cat("dmass    = ", fit$uwerrresultmps$dvalue, "\n")
    cat("ddmass   = ", fit$uwerrresultmps$ddvalue, "\n")
    cat("tauint   = ", fit$uwerrresultmps$tauint, "\n")
    cat("dtauint  = ", fit$uwerrresultmps$dtauint, "\n")
    cat("Wopt     = ", fit$uwerrresultmps$Wopt, "\n")
    if(fit$uwerrresultmps$R>1) {
      cat("Qval     =", fit$uwerrresultmps$Qval, "\n")
    }    
  }
  if(!is.null(fit$uwerrresultfps)) {
    cat("\n--- Autocorrelation analysis for f_ps ---\n")    
    cat("\nS        = ", fit$uwerrresultfps$S, "\n")
    cat("fps      = ", fit$uwerrresultfps$value*2*kappa*2*mu/sqrt(2), "\n")
    cat("dfps     = ", fit$uwerrresultfps$dvalue*2*kappa*2*mu/sqrt(2), "\n")
    cat("ddfps    = ", fit$uwerrresultfps$ddvalue*2*kappa*2*mu/sqrt(2), "\n")
    cat("tauint   = ", fit$uwerrresultfps$tauint, "\n")
    cat("dtauint  = ", fit$uwerrresultfps$dtauint, "\n")
    cat("Wopt     = ", fit$uwerrresultfps$Wopt, "\n")
    if(fit$uwerrresultfps$R>1) {
      cat("Qval     =", fit$uwerrresultfps$Qval, "\n")
    }
  }
  if(!is.null(fit$uwerrresultmpcac)) {
    cat("\n--- Autocorrelation analysis for m_pcac ---\n")
    cat("\nS        = ", fit$uwerrresultmpcac$S, "\n")
    cat("mpcac    = ", fit$uwerrresultmpcac$value, "\n")
    cat("dmpcac   = ", fit$uwerrresultmpcac$dvalue, "\n")
    cat("ddmpcac  = ", fit$uwerrresultmpcac$ddvalue, "\n")
    cat("tauint   = ", fit$uwerrresultmpcac$tauint, "\n")
    cat("dtauint  = ", fit$uwerrresultmpcac$dtauint, "\n")
    cat("Wopt     = ", fit$uwerrresultmpcac$Wopt, "\n")
    if(fit$uwerrresultmpcac$R>1) {
      cat("Qval     =", fit$uwerrresultmpcac$Qval, "\n")
    }
  }
  if(!is.null(fit$boot)) {
    cat("--- Bootstrap analysis  ---\n")
    cat("---", fit$boot$R, "samples  ---\n")
    cat("          mean        -err           +err            stderr        bias\n")
    fit$boot.ci <- boot.ci(fit$boot, type = c("norm"), index=1)
    cat("mpi    = ", fit$boot$t0[1], "(", (fit$boot.ci$normal[1,2]-fit$boot$t0[1])/1.96
	, ",", -(fit$boot$t0[1]-fit$boot.ci$normal[1,3])/1.96, ")", sd(fit$boot$t[,1]),
	mean(fit$boot$t[,1])-fit$boot$t0[1],"\n")
    
#    fit$boot.ci <- boot.ci(fit$boot, type = c("norm"), index=2)
#    cat("fpi    = ", fit$boot$t0[2], "(", (fit$boot.ci$normal[1,2]-fit$boot$t0[2])/1.96
#	, ",", -(fit$boot$t0[2]-fit$boot.ci$normal[1,3])/1.96, ")", sd(fit$boot$t[,2]),
#	mean(fit$boot$t[,2])-fit$boot$t0[2], "\n")

    fit$boot.ci <- boot.ci(fit$boot, type = c("norm"), index=2)
    cat("mpcac  = ", fit$boot$t0[2], "(", (fit$boot.ci$normal[1,2]-fit$boot$t0[2])/1.96
        , ",", -(fit$boot$t0[2]-fit$boot.ci$normal[1,3])/1.96, ")", sd(fit$boot$t[,2]),
        mean(fit$boot$t[,2])-fit$boot$t0[2], "\n")
    
    fit$boot.ci <- boot.ci(fit$boot, type = c("norm"), index=3)
    cat("P_L    = ", fit$boot$t0[3], "(", (fit$boot.ci$normal[1,2]-fit$boot$t0[3])/1.96
	, ",", -(fit$boot$t0[3]-fit$boot.ci$normal[1,3])/1.96, ")", sd(fit$boot$t[,3]),
	mean(fit$boot$t[,3])-fit$boot$t0[3], "\n")
    
    fit$boot.ci <- boot.ci(fit$boot, type = c("norm"), index=4)
    cat("A_L    = ", fit$boot$t0[4], "(", (fit$boot.ci$normal[1,2]-fit$boot$t0[4])/1.96
	, ",", -(fit$boot$t0[4]-fit$boot.ci$normal[1,3])/1.96, ")", sd(fit$boot$t[,4]),
	mean(fit$boot$t[,4])-fit$boot$t0[4],"\n")
  }
  if(!is.null(fit$tsboot)) {
    cat("\n--- Bootstrap analysis  with blocking ---\n")
    cat("---", fit$boot$R, "samples  ---\n")
    cat("--- block size", fit$tsboot$l, "---\n")
    fit$tsboot.ci <- boot.ci(fit$tsboot, type = c("norm"), index=1)
    cat("mpi    = ", fit$tsboot$t0[1], "(", (fit$tsboot.ci$normal[1,2]-fit$tsboot$t0[1])/1.96
	, ",", -(fit$tsboot$t0[1]-fit$tsboot.ci$normal[1,3])/1.96, ")", sd(fit$tsboot$t[,1]),
	mean(fit$tsboot$t[,1])-fit$tsboot$t0[1], "\n")
    
                                        #    fit$tsboot.ci <- boot.ci(fit$tsboot, type = c("norm"), index=2)
                                        #    cat("fpi    = ", fit$tsboot$t0[2], "(", (fit$tsboot.ci$normal[1,2]-fit$tsboot$t0[2])/1.96
                                        #	, ",", -(fit$tsboot$t0[2]-fit$tsboot.ci$normal[1,3])/1.96, ")", sd(fit$tsboot$t[,2]),
                                        #	mean(fit$tsboot$t[,2])-fit$tsboot$t0[2], "\n")

    fit$tsboot.ci <- boot.ci(fit$tsboot, type = c("norm"), index=2)
    cat("mpcac  = ", fit$tsboot$t0[2], "(", (fit$tsboot.ci$normal[1,2]-fit$tsboot$t0[2])/1.96
        , ",", -(fit$tsboot$t0[2]-fit$tsboot.ci$normal[1,3])/1.96, ")", sd(fit$tsboot$t[,2]),
        mean(fit$tsboot$t[,2])-fit$tsboot$t0[2], "\n")
    
    fit$tsboot.ci <- boot.ci(fit$tsboot, type = c("norm"), index=3)
    cat("P_L    = ", fit$tsboot$t0[3], "(", (fit$tsboot.ci$normal[1,2]-fit$tsboot$t0[3])/1.96
	, ",", -(fit$tsboot$t0[3]-fit$tsboot.ci$normal[1,3])/1.96, ")", sd(fit$tsboot$t[,3]),
	mean(fit$tsboot$t[,3])-fit$tsboot$t0[3], "\n")
    
    fit$tsboot.ci <- boot.ci(fit$tsboot, type = c("norm"), index=4)
    cat("A_L    = ", fit$tsboot$t0[4], "(", (fit$tsboot.ci$normal[1,2]-fit$tsboot$t0[4])/1.96
	, ",", -(fit$tsboot$t0[4]-fit$tsboot.ci$normal[1,3])/1.96, ")", sd(fit$tsboot$t[,4]),
	mean(fit$tsboot$t[,4])-fit$tsboot$t0[4], "\n")
    
  }  
}
