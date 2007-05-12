print.cfit <- function(fit) {
  summary(fit)
}

summary.cfit <- function(fit) {
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
  cat("fitrange = ", t1, "-", t2, "\n")
  cat("chi^2    = ", fit.chisqr, "\n")
  cat("dof    = ", fit.dof, "\n")
  cat("chi^2/dof = ", fit.chisqr/fit.dof, "\n")
  
  if(is.null(fit$uwerrresultmv) && is.null(fit$mv.boot)) {
    cat("\nmpi    = ", fit.mass, "\n")
    cat("fpi    = ", fit.fpi, "\n")
    if(fit$matrix.size > 2) {
      cat("mpcac  = ", fit.mass*fit$fitresult$par[3]/fit$fitresult$par[1]/2., "\n")
    }
    cat("\nP_L    = ", fit$fitresult$par[1], "\n")
    cat("P_F    = ", fit$fitresult$par[2], "\n")
    if(fit$matrix.size >2) {
      cat("A_L    = ", fit$fitresult$par[3], "\n")
      cat("A_F    = ", fit$fitresult$par[4], "\n")
    }
    if(fit$matrix.size >4) {
      cat("4_L    = ", fit$fitresult$par[5], "\n")
      cat("4_F    = ", fit$fitresult$par[6], "\n")
    }
  }
  if(!is.null(fit$uwerrresultmv)) {
    cat("mv     = ", abs(fit$fitresult$par[fit$matrix.size+1]), "\n")
    cat("4_L    = ", fit$fitresult$par[1], "\n")
    cat("4_F    = ", fit$fitresult$par[2], "\n")
    if(fit$matrix.size == 2 && fit$no.masses == 2) {
      cat("mv2    = ", abs(fit$fitresult$par[2*(fit$matrix.size+1)]), "\n")
      cat("4_L2   = ", fit$fitresult$par[(fit$matrix.size+2)], "\n")
      cat("4_F2   = ", fit$fitresult$par[(fit$matrix.size+3)], "\n")
    }
    if(fit$matrix.size == 4) {
      cat("A_L    = ", fit$fitresult$par[3], "\n")
      cat("A_F    = ", fit$fitresult$par[4], "\n")
    }
    if(fit$matrix.size == 6) {
      cat("V_L    = ", fit$fitresult$par[5], "\n")
      cat("V_F    = ", fit$fitresult$par[6], "\n")
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
    
  }
  if(!is.null(fit$uwerrresultmps)) {
    cat("\n--- Autocorrelation analysis for m_ps ---\n")
    cat("\nS        = ", fit$uwerrresultmps$S, "\n")
    cat("mps      = ", fit$uwerrresultmps$value, "\n")
    cat("dmps     = ", fit$uwerrresultmps$dvalue, "\n")
    cat("ddmps    = ", fit$uwerrresultmps$ddvalue, "\n")
    cat("tauint   = ", fit$uwerrresultmps$tauint, "\n")
    cat("dtauint  = ", fit$uwerrresultmps$dtauint, "\n")
    cat("Wopt     = ", fit$uwerrresultmps$Wopt, "\n")
    
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
    if(fit$no.masses > 1) {
      cat("\n--- Autocorrelation analysis for m_v2 ---\n")
      cat("\nS        = ", fit$uwerrresultmv2$S, "\n")
      cat("mv2      = ", fit$uwerrresultmv2$value, "\n")
      cat("dmv2     = ", fit$uwerrresultmv2$dvalue, "\n")
      cat("ddmv2    = ", fit$uwerrresultmv2$ddvalue, "\n")
      cat("tauint2  = ", fit$uwerrresultmv2$tauint, "\n")
      cat("dtauint2 = ", fit$uwerrresultmv2$dtauint, "\n")
      cat("Wopt2    = ", fit$uwerrresultmv2$Wopt, "\n")
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
  }
  if(!is.null(fit$boot)) {
    cat("--- Bootstrap analysis  ---\n")
    cat("---", fit$boot$R, "samples  ---\n")
    cat("          mean        -err           +err            stderr        bias\n")
    fit$boot.ci <- boot.ci(fit$boot, type = c("norm"), index=1)
    cat("mpi    = ", fit$boot$t0[1], "(", (fit$boot.ci$normal[1,2]-fit$boot$t0[1])/1.96
	, ",", -(fit$boot$t0[1]-fit$boot.ci$normal[1,3])/1.96, ")", sd(fit$boot$t[,1]),
	mean(fit$boot$t[,1])-fit$boot$t0[1],"\n")
    
    fit$boot.ci <- boot.ci(fit$boot, type = c("norm"), index=2)
    cat("fpi    = ", fit$boot$t0[2], "(", (fit$boot.ci$normal[1,2]-fit$boot$t0[2])/1.96
	, ",", -(fit$boot$t0[2]-fit$boot.ci$normal[1,3])/1.96, ")", sd(fit$boot$t[,2]),
	mean(fit$boot$t[,2])-fit$boot$t0[2], "\n")

    if(fit$matrix.size > 2) {
      fit$boot.ci <- boot.ci(fit$boot, type = c("norm"), index=fit$matrix.size+3)
      cat("mpcac  = ", fit$boot$t0[fit$matrix.size+3], "(", (fit$boot.ci$normal[1,2]-fit$boot$t0[fit$matrix.size+3])/1.96
          , ",", -(fit$boot$t0[fit$matrix.size+3]-fit$boot.ci$normal[1,3])/1.96, ")", sd(fit$boot$t[,(fit$matrix.size+3)]),
          mean(fit$boot$t[,(fit$matrix.size+3)])-fit$boot$t0[fit$matrix.size+3], "\n")
    }
    
    fit$boot.ci <- boot.ci(fit$boot, type = c("norm"), index=4)
    cat("P_L    = ", fit$boot$t0[4], "(", (fit$boot.ci$normal[1,2]-fit$boot$t0[4])/1.96
	, ",", -(fit$boot$t0[4]-fit$boot.ci$normal[1,3])/1.96, ")", sd(fit$boot$t[,4]),
	mean(fit$boot$t[,4])-fit$boot$t0[4], "\n")
    
    fit$boot.ci <- boot.ci(fit$boot, type = c("norm"), index=5)
    cat("P_F    = ", fit$boot$t0[5], "(", (fit$boot.ci$normal[1,2]-fit$boot$t0[5])/1.96
	, ",", -(fit$boot$t0[5]-fit$boot.ci$normal[1,3])/1.96, ")", sd(fit$boot$t[,5]),
	mean(fit$boot$t[,5])-fit$boot$t0[5],"\n")
  }
  if(!is.null(fit$tsboot)) {
    cat("\n--- Bootstrap analysis  with blocking ---\n")
    cat("---", fit$boot$R, "samples  ---\n")
    cat("--- block size", fit$tsboot$l, "---\n")
    fit$tsboot.ci <- boot.ci(fit$tsboot, type = c("norm"), index=1)
    cat("mpi    = ", fit$tsboot$t0[1], "(", (fit$tsboot.ci$normal[1,2]-fit$tsboot$t0[1])/1.96
	, ",", -(fit$tsboot$t0[1]-fit$tsboot.ci$normal[1,3])/1.96, ")", sd(fit$tsboot$t[,1]),
	mean(fit$tsboot$t[,1])-fit$tsboot$t0[1], "\n")
    
    fit$tsboot.ci <- boot.ci(fit$tsboot, type = c("norm"), index=2)
    cat("fpi    = ", fit$tsboot$t0[2], "(", (fit$tsboot.ci$normal[1,2]-fit$tsboot$t0[2])/1.96
	, ",", -(fit$tsboot$t0[2]-fit$tsboot.ci$normal[1,3])/1.96, ")", sd(fit$tsboot$t[,2]),
	mean(fit$tsboot$t[,2])-fit$tsboot$t0[2], "\n")

    if(fit$matrix.size > 2) {
      fit$tsboot.ci <- boot.ci(fit$tsboot, type = c("norm"), index=fit$matrix.size+3)
      cat("mpcac  = ", fit$tsboot$t0[fit$matrix.size+3], "(", (fit$tsboot.ci$normal[1,2]-fit$tsboot$t0[fit$matrix.size+3])/1.96
          , ",", -(fit$tsboot$t0[fit$matrix.size+3]-fit$tsboot.ci$normal[1,3])/1.96, ")", sd(fit$tsboot$t[,(fit$matrix.size+3)]),
          mean(fit$tsboot$t[,(fit$matrix.size+3)])-fit$tsboot$t0[fit$matrix.size+3], "\n")
    }

    
    fit$tsboot.ci <- boot.ci(fit$tsboot, type = c("norm"), index=4)
    cat("P_L    = ", fit$tsboot$t0[4], "(", (fit$tsboot.ci$normal[1,2]-fit$tsboot$t0[4])/1.96
	, ",", -(fit$tsboot$t0[4]-fit$tsboot.ci$normal[1,3])/1.96, ")", sd(fit$tsboot$t[,4]),
	mean(fit$tsboot$t[,4])-fit$tsboot$t0[4], "\n")
    
    fit$tsboot.ci <- boot.ci(fit$tsboot, type = c("norm"), index=5)
    cat("P_F    = ", fit$tsboot$t0[5], "(", (fit$tsboot.ci$normal[1,2]-fit$tsboot$t0[5])/1.96
	, ",", -(fit$tsboot$t0[5]-fit$tsboot.ci$normal[1,3])/1.96, ")", sd(fit$tsboot$t[,5]),
	mean(fit$tsboot$t[,5])-fit$tsboot$t0[5], "\n")
    
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
