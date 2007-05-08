CExp <- function(m, Time, x, sign=1.) {
  return(0.5*(exp(-m*x) + sign*exp(-m*(Time-x))))
}
