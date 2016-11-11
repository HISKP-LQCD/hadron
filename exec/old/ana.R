anaoutput <- function(data, from, to, npsf = 2, S=1.5) {

  plaq <- uwerrprimary(data$V1[from:to], plot=FALSE, S=S)
  dH <- uwerrprimary(data$V2[from:to]*data$V2[from:to], plot=FALSE, S=S)
  edH <- uwerrprimary(data$V3[from:to], plot=FALSE, S=S)
  iter0 <- uwerrprimary(data$V4[from:to]+data$V5[from:to]+data$V6[from:to], plot=FALSE, S=S)
  if(npsf == 1) {
    iter <- iter0
    acc <- uwerrprimary(data$V7[from:to], plot=FALSE, S=S);
    cat(plaq$value, plaq$dvalue, plaq$ddvalue, plaq$tauint, plaq$dtauint, iter$value, iter$dvalue, acc$value, acc$dvalue, edH$value, edH$dvalue, sqrt(dH$value), sqrt((dH$dvalue)^2/4./dh$value), "\n")
  }

  if(npsf == 2) {
    iter1 <- uwerrprimary(data$V7[from:to]+data$V8[from:to]+data$V9[from:to], plot=FALSE, S=S)
    iter <-uwerrprimary(data$V4[from:to]+data$V5[from:to]+data$V6[from:to]+data$V7[from:to]+data$V8[from:to]+data$V9[from:to], plot=FALSE, S=S)
    acc <- uwerrprimary(data$V10[from:to], plot=FALSE, S=S);
    cat(plaq$value, plaq$dvalue, plaq$ddvalue, plaq$tauint, plaq$dtauint, iter$value, iter$dvalue, acc$value, acc$dvalue, edH$value, sqrt(dH$value), sqrt((dH$dvalue)^2/4./dh$value), "\n")
  }
  if (npsf == 3) {
    iter1 <- uwerrprimary(data$V7[from:to]+data$V8[from:to]+data$V9[from:to], plot=FALSE, S=S)
    iter2 <- uwerrprimary(data$V10[from:to]+data$V11[from:to]+data$V12[from:to], plot=FALSE, S=S)
    iter <-uwerrprimary(data$V4[from:to]+data$V5[from:to]+data$V6[from:to]+data$V7[from:to]+data$V8[from:to]+data$V9[from:to]+data$V10[from:to]+data$V11[from:to]+data$V12[from:to],
                        plot=FALSE, S=S)
    acc <- uwerrprimary(data$V13[from:to], plot=FALSE, S=S);
#    cat(plaq$value, plaq$dvalue, dH$value, dH$dvalue, acc$value, acc$dvalue, ddH$value, ddH$dvalue, ddU$value, ddU$Dvalue, iter$value, iter$dvalue, iter0$value, iter0$dvalue, iter1$value, iter1$dvalue, iter2$value, iter2$dvalue, "\n", file="res.dat")
    cat(plaq$value, plaq$dvalue, plaq$ddvalue, plaq$tauint, plaq$dtauint, iter$value, iter$dvalue, acc$value, acc$dvalue, edH$value, sqrt(dH$value), sqrt((dH$dvalue)^2/4./dh$value), "\n")
  }

}

anaoutputaver <- function(file = "output.data", npsf = 2) {

  data <- read.table(file);
  ret <- read.table("return_check.data")
  plaq <- average(data$V1)
  dH <- average(data$V2*data$V2)
  iter0 <- average(data$V4+data$V5+data$V6)
  ddH <- average(sqrt(ret$V3*ret$V3))
  ddU <- average(ret$V5)
  if(npsf == 1) {
    iter <- iter0
    acc <- average(data$V7)
    cat(plaq$value, plaq$dvalue, sqrt(dH$value), sqrt((dH$dvalue)^2/4/dH$value), acc$value, acc$dvalue, ddH$value, ddH$dvalue, ddU$value, ddU$dvalue, iter0$value, iter0$dvalue, "\n", file="res.dat")
  }

  if(npsf == 2) {
    iter1 <- average(data$V7+data$V8+data$V9)
    iter <-average(data$V4+data$V5+data$V6+data$V7+data$V8+data$V9)
    acc <- average(data$V10)
    cat(plaq$value, plaq$dvalue, sqrt(dH$value), sqrt((dH$dvalue)^2/4/dH$value), acc$value, acc$dvalue, ddH$value, ddH$dvalue, ddU$value, ddU$dvalue, iter$value, iter$dvalue, iter0$value, iter0$dvalue, iter1$value, iter1$dvalue, "\n", file="res.dat")
  }
  if (npsf == 3) {
    iter1 <- average(data$V7+data$V8+data$V9)
    iter2 <- average(data$V10+data$V11+data$V12)
    iter <-average(data$V4+data$V5+data$V6+data$V7+data$V8+data$V9+data$V10+data$V11+data$V12)
    acc <- average(data$V13)
    cat(plaq$value, plaq$dvalue, sqrt(dH$value), sqrt((dH$dvalue)^2/4/dH$value), acc$value, acc$dvalue, ddH$value, ddH$dvalue, ddU$value, ddU$Dvalue, iter$value, iter$dvalue, iter0$value, iter0$dvalue, iter1$value, iter1$dvalue, iter2$value, iter2$dvalue, "\n", file="res.dat")
  }

}
