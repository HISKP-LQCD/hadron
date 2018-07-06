## The functions in this file date back to the very beginnings of hadron
## and are still used in a number of legacy routines, so they are kept
## around


# this was originally in pp.R
effmass <- function(data, timeextent, t) {
  mass <- invcosh((data[1]+data[4])/(data[2]+data[3]), timeextent=timeextent, t=t)
  return(invisible(mass))
}

effmass2 <- function(data, timeextent, t) {
  mass <- invcosh(ratio=(data[1])/(data[2]), timeextent=timeextent, t=t)
  return(invisible(mass))
}

effectivemass <- function(from, to, Time, Z, pl=TRUE, S,...) {
  L <- (to-from+1)
  i <- 1
  result <- data.frame(t = array(0.,dim=c(L)), mass = array(0.,dim=c(L)), dmass = array(0.,dim=c(L)),
                       ddmass = array(0.,dim=c(L)), tauint = array(0.,dim=c(L)), dtauint = array(0.,dim=c(L)))
  for(t in from:to) {
    try(mass <- uwerrderived(f=effmass2, data=t(Z[t:(t+1),]), S=S, pl=F, timeextent=Time, t=t, ...))

   cat(sprintf("Time t=%d\n", t)) 
       print(str(mass))
       cat("\n\n")

    result$t[i] <- t-1
    if( length(mass$res$value) > 0 ){
      result$mass[i] <- mass$res$value[1]
    }
    if( length(mass$res$dvalue) > 0 ){
      result$dmass[i] <- mass$res$dvalue[1]
    }
    if( length(mass$res$ddvalue) > 0 ){
      result$ddmass[i] <- mass$res$ddvalue[1]
    }
    if( length(mass$res$tauint) > 0 ){
      result$tauint[i] <- mass$res$tauint[1]
    }
    if( length(mass$res$dtauint) > 0 ){  
      result$dtauint[i] <- mass$res$dtauint[1]
    }
    i = i+1
  }
  rm(mass)
  rm(Z)
  attr(result, "class") <- c("massfit", "data.frame")
  if(pl == T) {
    if(interactive() && (grepl(pattern="X11", x=names(dev.cur()), ignore.case=TRUE) || grepl(pattern="null", x=names(dev.cur()), ignore.case=TRUE))) {
      X11()
    }
    plot(result)
  }
  return(invisible(result))
}

