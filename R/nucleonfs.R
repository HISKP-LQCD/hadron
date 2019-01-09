# this is only the asymptotic formula for the FS corrections
# of the nucleon mass as found in
# hep-lat/0403015
nucleonfs <- function(amN, gA, gDelta = 1.4, Delta, L, ampi, aF0) {
  # eq(17)
  deltaNA <- (9.*gA^2*ampi/(8*pi*aF0^2) + 4.*gDelta^2*ampi^{5/2}/((2*pi)^(3/2)*aF0^2*Delta*sqrt(L))) *
    exp(-ampi*L)/L
  return(invisible(list(amNV = amN - deltaNA)))
}

