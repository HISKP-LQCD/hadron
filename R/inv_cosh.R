

invcosh <- function(ratio, timeextent, t, eps=1.e-9, maxiterations=1000) {

  if(ratio < 1 || is.na(ratio)) {
    return(NA);
#    stop("Error: ratio is smaller than 1 in invcosh!")
  }

  return(.Call("invcosh", ratio, timeextent, t, eps, maxiterations))
}
