# ...

gsl_fit_smeared_correlator <- function(par, Thalf, x, y, err, tr, prec=c(1.e-10,1.e-4), N=1) {
  if(missing(par)) {
    stop("Error, parameter list must be given!")
  }
                                        #etc
  
  npar <- length(par)
  parsave <- numeric(npar)
  parsave <- par

  state <- .Call("multifit_smearedcor", par, Thalf, x, y, err, tr, prec, N, 500, 1)
  if(state[5] >= 0) {
    return(invisible(list(par=state[6:(npar+5)], value=state[1],
                          convergence=state[5], counts = state[3], dof = state[4])))
  }
  else if(state[5] < 0) {
    ##state <- .Call("multimin_cor", parsave, Thalf, x, y, err, tr, prec, N, 500, no_masses)
    ##if(state[5] < 0) {
      warning("gsl_multifit did not converge ", state[5], " chisqr ",
              state[1], "\npars ", state[6:(npar+5)],
              "\nResults may be unreliable!\nPlease try to vary initial guesses!");
    ##}
  }

  return(invisible(list(par=state[6:(npar+5)], value=state[1],
                        convergence=state[5], counts = state[3], dof = state[4])))
}

gsl_fit_correlator_matrix <- function(par, Thalf, x, y, err, tr, N, no_masses=1, prec=c(1.e-10,1.e-4)) {
  if(missing(par)) {
    stop("Error, parameter list must be given!")
  }
                                        #etc

  npar <- length(par)
  parsave <- numeric(npar)
  parsave <- par
  ##  for(i in 1:npar) parsave[i] <- par[i]

  state <- .Call("multifit_cor", par, Thalf, x, y, err, tr, prec, N, 500, no_masses)
  if(state[5] >= 0) {
    return(invisible(list(par=state[6:(npar+5)], value=state[1],
                          convergence=state[5], counts = state[3], dof = state[4])))
  }
  else if(state[5] < 0) {
    state <- .Call("multimin_cor", parsave, Thalf, x, y, err, tr, prec, N, 500, no_masses)
    if(state[5] < 0) {
      warning("gsl_multifit and multimin did not converge ", state[5], " chisqr ",
              state[1], "\npars ", state[6:(npar+5)],
              "\nResults may be unreliable!\nPlease try to vary initial guesses!");
    }
  }

  return(invisible(list(par=state[6:(npar+5)], value=state[1],
                        convergence=state[5], counts = state[3], dof = state[4])))
}

gsl_min_correlator_matrix <- function(par, Thalf, x, y, err, tr, N, no_masses=1, prec=c(1.e-10,1.e-4)) {
  if(missing(par)) {
    stop("Error, parameter list must be given!")
  }
                                        #etc

  npar <- length(par)
  parsave <- numeric(npar)
  
  state <- .Call("multimin_cor", par, Thalf, x, y, err, tr, prec, N, 500, no_masses)
  if(state[5] < 0) {
    warning("gsl_multimin did not converge ", state[5], " chisqr ", state[1], "\npars ", state[6:(npar+5)])
    state <- .Call("multifit_cor", parsave, Thalf, x, y, err, tr, prec, N, 500, no_masses)
    if(state[5] < 0) warning("failed again...")
  }

  return(invisible(list(par=state[6:(npar+5)], value=state[1],
                        convergence=state[5], counts = state[3], dof = state[4])))
}
