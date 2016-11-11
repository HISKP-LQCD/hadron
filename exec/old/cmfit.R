cmfit <- function(fitpar, optim.func=NULL, optim.method="BFGS", optim.control=list(trace=0),
                  gsl.prec=c(1.e-10,1.e-3), no.masses=1,
                  Thalf, x, y, err, tr, N, fit.routine="gsl") {

  if(missing(fitpar)) {
    stop("Error, parameter list must be given!")
  }


  if(fit.routine == "gsl") {
    fit <- gsl_fit_correlator_matrix(par=fitpar, Thalf=Thalf, x=x, y=y,
                                       err=err, tr=tr, N=N, no_masses=no.masses,
                                       prec=c(1.e-10,1.e-3))
  }
  else if(fit.routine == "optim") {
    fit <- optim(par=fitpar, fn=optim.func, method=optim.method, control=optim.control, Thalf=Thalf,
                 x=x, y=y, err=err, tr=tr, N=N)

  }
  else {
    npar <- length(fitpar)
    parsave <- numeric(npar)
    for(i in 1:npar) parsave[i] <- fitpar[i]
    
    state <- .Call("multifit_cor", fitpar, Thalf, x, y, err, tr, prec, N, 500, no_masses)
    if(state[5] >= 0) {
      return(invisible(list(par=state[6:(npar+5)], value=state[1],
                            convergence=state[5], counts = state[3], dof = state[4])))
    }
    else if(state[5] < 0) {
      state <- .Call("multimin_cor", parsave, Thalf, x, y, err, tr, prec, N, 500, no_masses)
      if(state[5] >= 0) {
        return(invisible(list(par=state[6:(npar+5)], value=state[1],
                              convergence=state[5], counts = state[3], dof = state[4])))
      }
      if(state[5] < 0) {
        fit <- optim(par=parsave, fn=optim.func, method=optim.method, control=optim.control, Thalf=Thalf,
                     x=x, y=y, err=err, tr=tr, N=N)
      }
    }
  }

  return(invisible(fit))
}
