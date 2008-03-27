# ...


gsl_fit_correlator_matrix <- function(par, Thalf, x, y, err, tr, N, prec=c(1.e-4,1.e-4)) {
  if(missing(par)) {
    stop("Error, parameter list must be given!")
  }
                                        #etc

  value = 0.

  state <- .Call("multifit_S1", par, Thalf, x, y, err, tr, value, prec, N, 100)
  return(invisible(list(par=par, value=value, convergence=state[5], counts = state[3], dof = state[4])))
}
