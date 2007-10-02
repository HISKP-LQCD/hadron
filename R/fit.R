fitcoshnlm <- function(data, T2, cutoff, A=0.01, m=0.01, debug=FALSE) {
  S <- cutoff
  E <- T2+2-cutoff
  H <- (T2/2+1)

  x=c(S:E)
  y=data
  fnct <- function(p) sum((y - p[1]*cosh(p[2]*(x-T2)))^2)
  cfit <- nlm(fnct, p = c(A,m), hessian=FALSE)
  if(debug == TRUE) {
    plot(x, y,log="y")
    xfit <- c(S:E)
    yfit <- cfit$estimate[1]*cosh(cfit$estimate[2]*(xfit-H))
    lines(spline(xfit,yfit), col = "red")
    readline("Please type <RETURN> to continue")
  }
                                        #  detach(df)
                                        #  rm(df)
  rm(fnct)
  
  res <- list(A=cfit$estimate[1], m=cfit$estimate[2])
  return(invisible(res))
}

# We know the derivative, so let's use it.
coshfit <- deriv(~A*cosh(m*(x-T2)),
                 c('A', 'm'),
                 function(A,m,T2,x){}
                 )

doublecoshfit <- deriv(~(A*cosh(m*(x-T2)) + B*cosh(m2*(x-T2))),
                       c('A', 'm', 'B', 'm2'),
                       function(A, m, B, m2, T2, x){}
                       )

fitdoublecoshnls <- function(data, T2, cutoff, A=0.01, B=1., m=0.1, m2=.6, debug=FALSE) {
  S <- cutoff
  E <- T2+2-cutoff
  H <- (T2/2+1)
  df <- data.frame(x=c(S:(H-15)), y=data, T2=H)
  ##   nls.control(maxiter = 500)
  cfit <- nls( y ~ doublecoshfit(A, m, B, m2, T2[1], x), data = df, start = list(A = A, m = m, B = B, m2 = m2), trace=T, control=nls.control(maxiter = 500000, minFactor = 1/10000000000), algorithm="port")
##  cfit <- nls( y ~ (A*cosh(m*(x-T2)) + B*cosh(m2*(x-T2))), data = df, start = list(A = A, m = m, B = B, m2 = m2), trace=F, control=nls.control(maxiter = 500000, minFactor = 1/10000000000))
  if(debug == TRUE) {
    print(summary(cfit))
    print(cfit)
    plot(df$x, df$y,log="y", ylim=c(0.001,0.11), xlim=c((S-1),(H+3)))
    xfit <- c(S:(H+3))
    yfit <- coef(cfit)[[1]]*cosh(coef(cfit)[[2]]*(xfit-H)) + coef(cfit)[[3]]*cosh(coef(cfit)[[4]]*(xfit-H))
##    yfit <- A*cosh(m*(xfit-H)) + B*cosh(m2*(xfit-H))
    lines(spline(xfit,yfit), col = "red")
    readline("Please type <RETURN> to continue")
  }
  rm(df)
  res <- list(A=coef(cfit)[[1]], m=coef(cfit)[[2]], B=coef(cfit)[[3]], m2=coef(cfit)[[4]])
                                        #  options(error = NULL)
  return(invisible(res))
}

fitcoshnls <- function(data, T2, cutoff, A=0.01, m=0.01, debug=FALSE) {
  S <- cutoff
  E <- T2+2-cutoff
  H <- (T2/2+1)
  df <- data.frame(x=c(S:E), y=data, T2=H)
                                        #  cfit <- nls(y ~ A*cosh(m*(x-T2[1])), data = df, start = list(A = 0.0009, m = 0.04), trace=TRUE)
  # We have to catch the error
                                        #  try(cfit <- nls(y ~ coshfit(A,m,T2[1],x), data = df, start = list(A = A, m = m), trace=F), silent=F)
  cfit <- nls(y ~ coshfit(A,m,T2[1],x), data = df, start = list(A = A, m = m), trace=F)
  if(debug == TRUE) {
    print(summary(cfit))
    plot(df$x, df$y,log="y")
    xfit <- c(S:E)
    yfit <- coef(cfit)[[1]]*cosh(coef(cfit)[[2]]*(xfit-H))
    lines(spline(xfit,yfit), col = "red")
    readline("Please type <RETURN> to continue")
  }
  rm(df)
  # catch the error here again, but the expression does not make sense?!?
                                        #  options(error = quote(res <- list(A = 0.01, m = 0.01)))
                                        #  try(res <- list(A=coef(cfit)[[1]], m=coef(cfit)[[2]]),silent=T)
  res <- list(A=coef(cfit)[[1]], m=coef(cfit)[[2]])
                                        #  options(error = NULL)
  return(invisible(res))
}

fitcoshnls2 <- function(data, T2, from, to, A=0.01, m=0.01, debug=FALSE) {

#  H <- (T2/2+1)
  H <- (T2/2)
  df <- data.frame(x=c((from-1):(to-1)), y=data, T2=H)
#  plot(df$x, df$y,log="y")
  cfit <- nls(y ~ coshfit(A,m,T2[1],x), data = df, start = list(A = A, m = m), trace=F, algorithm="port")
  if(debug == TRUE) {
    print(summary(cfit))
    print(cfit)
    chisq <- 0.
    for(i in 1:length(data)) {
      chisq = chisq + ((fitted(cfit)[[i]]-data[i])^2)
    }
#    print(chisq)

  }
  rm(df)
  res <- list(A=coef(cfit)[[1]], m=coef(cfit)[[2]])
                                        #  options(error = NULL)
  return(invisible(res))
}

fitconst <- function(data) {

  df <- data.frame(y=data)
  cfit <- lm(y ~ 1, df)
  rm(df)
  res <- list(c = coef(cfit)[[1]])
  return(invisible(res))
}

getmass <- function(data, T2, cutoff, A=0.01, m=0.01, debug=T) {
  fit <- fitcoshnls(data=data, T2=T2, cutoff=cutoff, A=1., m=0.1, debug=debug)
  return(fit$m)
}

getamp <- function(data, T2, cutoff, A=0.01, m=0.01, debug=FALSE) {
  fit <- fitcoshnls(data=data, T2=T2, cutoff=cutoff, A=1., m=0.1, debug=debug)
  return(fit$A)
}

getfps <- function(data, T2, cutoff, A=0.01, m=0.01, debug=FALSE, mq) {
  fit <- fitcoshnls(data=data, T2=T2, cutoff=cutoff, A=1., m=0.1, debug=debug)
  return(2.*mq*sqrt(fit$A/abs(fit$m)^3)*exp(abs(fit$m)*T2/4.))
}

getmass2 <- function(data, T2, from, to, A=0.01, m=0.01, debug=T, mu=0.) {
  fit <- fitcoshnls2(data=data, T2=T2, from=from, to=to, A=1., m=0.1, debug=debug)
  return(fit$m)
}
getamp2 <- function(data, T2, from, to, A=0.01, m=0.01, debug=T, mu=0.) {
  fit <- fitcoshnls2(data=data, T2=T2, from=from, to=to, A=1., m=0.1, debug=debug)
  return(fit$A)
}
getfpi2 <- function(data, T2, from, to, A=0.01, m=0.01, debug=T, mu) {
  fit <- fitcoshnls2(data=data, T2=T2, from=from, to=to, A=1., m=0.1, debug=debug)
  Fpi=2.*mu*sqrt(fit$A/abs(fit$m)^3)*exp(abs(fit$m)*T2/4.)
  return(Fpi)
}
getsigma2 <- function(data, T2, from, to, A=0.01, m=0.01, debug=T, mu) {
  fit <- fitcoshnls2(data=data, T2=T2, from=from, to=to, A=1., m=0.1, debug=debug)
  Fpi=2.*mu*sqrt(fit$A/abs(fit$m)^3)*exp(abs(fit$m)*T2/4.)
  sigma = Fpi^2 * fit$m^2 / 4. / mu
  return(sigma)
}

getSigmaInter <- function(data, from1, to1, from2, to2,
                         T2, A, m, mu1, mu2, r0, r0mpssq) {
  range1 <- to1-from1+1
  range2 <- to2-from2+1
  fit1 <- fitcoshnls2(data=data[1:range1], T2=T2, from=from1, to=to1,
                      A=1., m=0.1, debug=FALSE)
  fit2 <- fitcoshnls2(data=data[(range1+1):(range1+range2)], T2=T2, from=from2, to=to2,
                      A=1., m=0.1, debug=FALSE)
  r0mpssq1 = (r0*fit1$m)^2
  r0mpssq2 = (r0*fit2$m)^2
  Fpi1=2.*mu1*sqrt(fit1$A/abs(fit1$m)^3)*exp(abs(fit1$m)*T2/4.)
  Sigma1 = Fpi1^2 * fit1$m^2 / 4. / mu1
  Fpi2=2.*mu2*sqrt(fit2$A/abs(fit2$m)^3)*exp(abs(fit2$m)*T2/4.)
  Sigma2 = Fpi2^2 * fit2$m^2 / 4. / mu2
  a <- (Sigma2-Sigma1)/(r0mpssq2-r0mpssq1)
  b = Sigma1-a*r0mpssq1
#  fitdata <- data.frame(y = c(Sigma1, Sigma2),
#                        x = c(r0mpssq1, r0mpssq2))
#  fit <- glm(y ~ 1 + x, data=fitdata)
#  return(coef(fit)[[2]]*r0mpssq+coef(fit)[[1]])
  return(a*r0mpssq+b)
}

getMuIntermps <- function(data, from1, to1, from2, to2,
                         T2, A, m, mu1, mu2, r0, r0mpssq) {
  range1 <- to1-from1+1
  range2 <- to2-from2+1
  fit1 <- fitcoshnls2(data=data[1:range1], T2=T2, from=from1, to=to1,
                      A=1., m=0.1, debug=FALSE)
  fit2 <- fitcoshnls2(data=data[(range1+1):(range1+range2)], T2=T2, from=from2, to=to2,
                      A=1., m=0.1, debug=FALSE)
  r0mpssq1 = (r0*fit1$m)^2
  r0mpssq2 = (r0*fit2$m)^2
  fitdata <- data.frame(y = c(mu1, mu2),
                        x = c(r0mpssq1, r0mpssq2))
  fit <- glm(y ~ 1 + x, data=fitdata)
  return(coef(fit)[[2]]*r0mpssq+coef(fit)[[1]])
}

getMuInterfps <- function(data, from1, to1, from2, to2,
                         T2, A, m, mu1, mu2, r0, r0fps) {
  range1 <- to1-from1+1
  range2 <- to2-from2+1
  fit1 <- fitcoshnls2(data=data[1:range1], T2=T2, from=from1, to=to1,
                      A=1., m=0.1, debug=FALSE)
  fit2 <- fitcoshnls2(data=data[(range1+1):(range1+range2)], T2=T2, from=from2, to=to2,
                      A=1., m=0.1, debug=FALSE)
  r0Fpi1=abs(r0*2.*mu1*sqrt(fit1$A/abs(fit1$m)^3)*exp(abs(fit1$m)*T2/4.))
  r0Fpi2=abs(r0*2.*mu2*sqrt(fit2$A/abs(fit2$m)^3)*exp(abs(fit2$m)*T2/4.))
  fitdata <- data.frame(y = c(mu1, mu2),
                        x = c(r0Fpi1, r0Fpi2))
  fit <- glm(y ~ 1 + x, data=fitdata)
  return(coef(fit)[[2]]*r0fps+coef(fit)[[1]])
}


# stderr <- function(x) sqrt(var(x)/length(x))

tauint <- function(x,W,M=NULL) {
  gamma <- acf(x,lag.max=M)
  v <- gamma$acf
  ti <- 0.5*v[1] + sum(v[2:W])
  dti <- 2./length(x)*(2*W-3*ti+1+1/(4*ti))*ti^2
  res <- c(ti, dti)
  res
}

average <- function(x) {
  res <- list(value = mean(x), dvalue = sd(x)/sqrt(length(x)), variance = var(x))
  return(invisible(res))
}




#pp <- read.table(file="parplaq", col.names=c("plaq"), header=F, skip=0)
#Z <- array(x, dim=c(T,length(x)/T)
#psscar <- read.table(file="pssca_corr_1.dat", col.names=c("t","ps"), header=F, skip=0)
#T<-(max(t)-min(t)+1)
#Z <- array(ps, dim=c(T,length(ps)/T) 
