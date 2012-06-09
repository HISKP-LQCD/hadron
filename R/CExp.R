CExp <- function(m, Time, x, sign=1.) {
  return(0.5*(exp(-m*(Time-x)) + sign*exp(-m*x)))
}

dCExpdm <- function(m, Time, x, sign=1.) {
  return(0.5*(-(Time-x)*exp(-m*(Time-x)) -x* sign*exp(-m*x)))
}

