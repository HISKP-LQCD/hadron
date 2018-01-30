

invcosh <- function(ratio, timeextent, t, eps=1.e-9, maxiterations=1000) {

  if(ratio < 1 || is.na(ratio)) {
    return(NA);
#    stop("Error: ratio is smaller than 1 in invcosh!")
  }
  mass <- log(ratio)
  newmass <- log(ratio)
  i <-  0
  repeat {
    mass = newmass
    i=i+1
    if(t <= timeextent/2) {
      r <-  (1+exp(-mass*(timeextent - 2*(t-1) - 2.)))/(1+exp(-mass*(timeextent - 2.*(t-1))))
    }
    else {
      r <-  (1+exp(-mass*(timeextent - 2*(t-1) + 2.)))/(1+exp(-mass*(timeextent - 2.*(t-1))))
    }
    newmass = log(ratio*r)
#    cat("iteration = ",i, maxiterations,"\n")
#    cat("mass = ",mass, "newmass = ", newmass, "\n")
    if((abs(mass-newmass) < eps) || (i > maxiterations)) {
      return(newmass)
    }
  }

  return(mass)
}
