# ...


gsl_fit_correlator_matrix <- function(par, Thalf, x, y, err, tr, N, no_masses=1, prec=c(1.e-10,1.e-4)) {
  if(missing(par)) {
    stop("Error, parameter list must be given!")
  }
                                        #etc

  value = 0.

  state <- .Call("multifit_cor", par, Thalf, x, y, err, tr, value, prec, N, 500, no_masses)
  convergence <- state[5]
  return(invisible(list(par=par, value=value, convergence=convergence, counts = state[3], dof = state[4])))
}
