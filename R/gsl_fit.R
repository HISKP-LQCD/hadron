# ...


gsl_fit_correlator_matrix <- function(par, Thalf, x, y, err, tr, N, no_masses=1, prec=c(1.e-10,1.e-4)) {
  if(missing(par)) {
    stop("Error, parameter list must be given!")
  }
                                        #etc

  npar <- length(par)
  parsave <- numeric(npar)
  
  state <- .Call("multifit_cor", par, Thalf, x, y, err, tr, prec, N, 500, no_masses)
  if(state[5] < 0) {
    warning("gsl_multifit did not converge ", state[5])
  }

  return(invisible(list(par=state[6:(npar+5)], value=state[1], convergence=state[5], counts = state[3], dof = state[4])))
}
