interpolateSigma <- function(filename1, filename2, skip=0, from1, to1, from2, to2,
                             S=1.5, A=0.01, m=0.01, plot=FALSE, debug=FALSE,
                             mu1=0.01, mu2=0.01, r0, r0mpssq, Zp, dZp) {
  if(!missing(filename1)) {
    psscar1 <- read.table(file=filename1, col.names=c("t","ps"), header=F)
  }
  else {
    stop("Error! Option filename1 is mandatory!")
  }
  if(!missing(filename2)) {
    psscar2 <- read.table(file=filename2, col.names=c("t","ps"), header=F)
  }
  else {
    stop("Error! Option filename2 is mandatory!")
  }
  T2 <- (max(psscar1$t)-min(psscar1$t)+1)
  Z <- array(psscar1$ps, dim=c(T2,length(psscar1$ps)/T2))
  W1 <- array(0, dim=c((T2/2+1), length(psscar1$ps)/T2))
  for(i in 1:(T2/2+1)) {
    if(i==1 || i == (T2/2+1)) {
      W1[i,] <- Z[i,]
    }
    else {
      W1[i,] = 0.5*(Z[i,]+Z[(T2+2-i),])
    }
  }
  rm(Z)
  Nalpha <- length(W1[(from1):(to1),1])
  N <- length(W1[1,skip:(length(psscar2$ps)/T2)])
  abb <- c(1:Nalpha)
  for (i in 1:Nalpha) {
    abb[i] <- mean(W1[(i-1+from1),])
  }
  cat("Performing a first fit for first correlator:\n\n")
  print(abb)
  fitresult <- fitcoshnls2(abb, T2=T2, from=from1, to=to1, A=A, m=m, debug=T)
  Z <- array(psscar2$ps, dim=c(T2,length(psscar1$ps)/T2))
  W2 <- array(0, dim=c((T2/2+1), length(psscar1$ps)/T2))
  for(i in 1:(T2/2+1)) {
    if(i==1 || i == (T2/2+1)) {
      W2[i,] <- Z[i,]
    }
    else {
      W2[i,] = 0.5*(Z[i,]+Z[(T2+2-i),])
    }
  }
  Nalpha <- length(W2[(from2):(to2),1])
  N <- length(W2[1,skip:(length(psscar2$ps)/T2)])
  abb <- c(1:Nalpha)
  for (i in 1:Nalpha) {
    abb[i] <- mean(W2[(i-1+from2),])
  }
  cat("Performing a first fit for second correlator:\n\n")
  print(abb)
  fitresult <- fitcoshnls2(abb, T2=T2, from=from2, to=to2, A=A, m=m, debug=T)
  W <- rbind(W1[from1:to1,], W2[from2:to2,])
  rm(W1, W2, Z)
  result <- data.frame(r0mpssq = r0mpssq, r0=r0, Zp=Zp, dZp=dZp, sigma=0., dsigma=0., ddsigma=0.,
                       sigmatauint=0., sigmadtauint=0., sigmaZpr0c=0., dsigmaZpr0c=0.,
                       mu=0., dmu=0., ddmu=0., mutauint=0., mudtauint=0., r0muZp =0.,
                       dr0muZp=0.)
  cat("Found", length(psscar1$ps)/T2, "measurements, skipping", skip, " \n")
  sigma <- uwerrderived(getSigmaInter, t(W[,skip:length(W[1,(skip):(length(psscar1$ps)/T2)])]), S, plot=F, from1, to1,
                        from2, to2, T2, A, m, mu1, mu2, r0, r0mpssq)
  mu <- uwerrderived(getMuIntermps, t(W[,skip:length(W[(skip):(length(psscar1$ps)/T2)])]), S, plot=F, from1, to1,
                     from2, to2, T2, A, m, mu1, mu2, r0, r0mpssq)
  
  result$sigma <- sigma$res$value[1]
  result$dsigma <- sigma$res$dvalue[1]
  result$ddsigma <- sigma$res$ddvalue[1]
  result$sigmatauint <- sigma$res$tauint[1]
  result$sigmadtauint <- sigma$res$dtauint[1]
  result$sigmaZpr0c <- sigma$res$value[1]*Zp*r0^3/2.
  result$dsigmaZpr0c <- r0^3*(sqrt(sigma$res$value[1]^2*dZp^2 + sigma$res$dvalue[1]^2*Zp^2))/2.
  result$mu <- mu$res$value
  result$dmu <- mu$res$dvalue
  result$ddmu <- mu$res$ddvalue
  result$mutauint <- mu$res$tauint
  result$mudtauint <- mu$res$dtauint
  result$r0muZp <- r0*mu$res$value/Zp
  result$dr0muZp <- r0*sqrt((mu$res$dvalue/Zp)^2 + (mu$res$value/Zp^2*dZp)^2)
  return(result)
}

interpolateMu <- function(filename1, filename2, skip=0, from1, to1, from2, to2,
                          S=1.5, A=0.01, m=0.01, plot=FALSE, debug=FALSE,
                          mu1=0.01, mu2=0.01, r0, r0mpssq, r0fps, Zp, dZp) {
  if(!missing(filename1)) {
    psscar1 <- read.table(file=filename1, col.names=c("t","ps"), header=F)
  }
  else {
    stop("Error! Option filename1 is mandatory!")
  }
  if(!missing(filename2)) {
    psscar2 <- read.table(file=filename2, col.names=c("t","ps"), header=F)
  }
  else {
    stop("Error! Option filename2 is mandatory!")
  }
  T2 <- (max(psscar1$t)-min(psscar1$t)+1)
  Z <- array(psscar1$ps, dim=c(T2,length(psscar1$ps)/T2))
  W1 <- array(0, dim=c((T2/2+1), length(psscar1$ps)/T2))
  for(i in 1:(T2/2+1)) {
    if(i==1 || i == (T2/2+1)) {
      W1[i,] <- Z[i,]
    }
    else {
      W1[i,] = 0.5*(Z[i,]+Z[(T2+2-i),])
    }
  }
  rm(Z)
  Z <- array(psscar2$ps, dim=c(T2,length(psscar1$ps)/T2))
  W2 <- array(0, dim=c((T2/2+1), length(psscar1$ps)/T2))
  for(i in 1:(T2/2+1)) {
    if(i==1 || i == (T2/2+1)) {
      W2[i,] <- Z[i,]
    }
    else {
      W2[i,] = 0.5*(Z[i,]+Z[(T2+2-i),])
    }
  }
  W <- rbind(W1[from1:to1,], W2[from2:to2,])
  rm(W1, W2, Z)
  result <- data.frame(r0mpssq = r0mpssq, r0fps = r0fps, r0=r0, Zp=Zp, dZp=dZp,
                       mumps = 0., dmumps = 0., ddmumps = 0., mumpstauint=0., mumpsdtauint=0., r0mumpsZp=0., dr0mumpsZp=0.,
                       mufps = 0., dmufps = 0., ddmufps = 0., mufpstauint=0., mufpsdtauint=0., r0mufpsZp=0., dr0mufpsZp=0.)
  cat("Found", length(psscar1$ps)/T2, "measurements, skipping", skip, " \n")
  mumps <- uwerrderived(getMuIntermps, t(W[,skip:length(W[1,])]), S, plot=F, from1, to1,
                        from2, to2, T2, A, m, mu1, mu2, r0, r0mpssq)
  mufps <- uwerrderived(getMuInterfps, t(W[,skip:length(W[1,])]), S, plot=F, from1, to1,
                        from2, to2, T2, A, m, mu1, mu2, r0, r0fps)
  
  result$mumps <- mumps$res$value[1]
  result$dmumps <- mumps$res$dvalue[1]
  result$ddmumps <- mumps$res$ddvalue[1]
  result$mumpstauint <- mumps$res$tauint[1]
  result$mumpsdtauint <- mumps$res$dtauint[1]
  result$r0mumpsZp <- r0*mumps$res$value/Zp
  result$dr0mumpsZp <- r0*sqrt((mumps$res$dvalue/Zp)^2 + (mumps$res$value/Zp^2*dZp)^2)
  result$mufps <- mufps$res$value
  result$dmufps <- mufps$res$dvalue
  result$ddmufps <- mufps$res$ddvalue
  result$mufpstauint <- mufps$res$tauint
  result$mufpsdtauint <- mufps$res$dtauint
  result$r0mufpsZp <- r0*mufps$res$value/Zp
  result$dr0mufpsZp <- r0*sqrt((mufps$res$dvalue/Zp)^2 + (mufps$res$value/Zp^2*dZp)^2)
  return(result)
}


ppfit <- function(filename, skip=0, from, to, S=1.5, A=0.01, m=0.01, plot=FALSE, debug=FALSE, mu=0.01) {
  if(!missing(filename)) {
    psscar <- read.table(file=filename, col.names=c("t","ps"), header=F)
  }
  else {
    stop("Error! Option psfilename is mandatory!")
  }
  T2 <- (max(psscar$t)-min(psscar$t)+1)
  if(missing(from)) {
    from <- 1
  }
  if(missing(to) || to > T2/2+1) {
    to <- T2/2
  }
  Z <- array(psscar$ps, dim=c(T2,length(psscar$ps)/T2))
  W <- array(0, dim=c((T2/2+1), length(psscar$ps)/T2))
                                        #Fold the data
  for(i in 1:(T2/2+1)) {
    if(i==1 || i == (T2/2+1)) {
      W[i,] <- Z[i,]
    }
    else {
      W[i,] = 0.5*(Z[i,]+Z[(T2+2-i),])
    }
  }
  Nalpha <- length(W[(from):(to),1])
  N <- length(W[1,skip:(length(psscar$ps)/T2)])
  abb <- c(1:Nalpha)
  for (i in 1:Nalpha) {
    abb[i] <- mean(W[(i-1+from),])
  }
  cat("Performing a first fit:\n\n")
  fitresult <- fitcoshnls2(abb, T2=T2, from=from, to=to, A=A, m=m, debug=T)
  
  result <- data.frame(from = from, to = to, mass = 0., dmass = 0.,
                       amp = 0., damp = 0., fpi = 0., dfpi = 0.,
                       sigma = 0., dsigma = 0.,
                       ddmass = 0., masstauint = 0.,
                       massdtauint = 0.,
                       ddamp = 0., amptauint = 0., ampdtauint = 0.,
                       ddfpi = 0., fpitauint = 0., fpidtauint = 0.,
                       ddsigma = 0., sigmatauint = 0., sigmadtauint = 0.)

  cat("\nComputing the error...\n")
  cat("Found", length(psscar$ps)/T2, "measurements, skipping", skip, " \n")
  cat("fitting from timeslice", from, " to timeslice ", to, " !\n")

  try(mass <- uwerrderived(getmass2, data=t(W[(from):to,skip:(length(psscar$ps)/T2)]),
                           S, plot=debug, T2, from, to, A, m, debug=debug))

  try(amp <- uwerrderived(getamp2, data=t(W[(from):to,skip:(length(psscar$ps)/T2)]),
                           S, plot=debug, T2, from, to, A, m, debug=debug))

  try(fpi <- uwerrderived(getfpi2, data=t(W[(from):to,skip:(length(psscar$ps)/T2)]),
                           S, plot=debug, T2, from, to, A, m, debug=debug, mu=mu))

  try(sigma <- uwerrderived(getsigma2, data=t(W[(from):to,skip:(length(psscar$ps)/T2)]),
                          S, plot=debug, T2, from, to, A, m, debug=debug, mu=mu))
  
  result$mass <- mass$res$value[1]
  result$dmass <- mass$res$dvalue[1]
  result$amp <- amp$res$value[1]
  result$damp <- amp$res$dvalue[1]
  result$ddmass <- mass$res$ddvalue[1]
  result$masstauint <- mass$res$tauint[1]
  result$massdtauint <- mass$res$dtauint[1]
  result$ddamp <- amp$res$ddvalue[1]
  result$amptauint <- amp$res$tauint[1]
  result$ampdtauint <- amp$res$dtauint[1]
  result$fpi <- fpi$res$value[1]
  result$dfpi <- fpi$res$dvalue[1]
  result$ddfpi <- fpi$res$ddvalue[1]
  result$fpitauint <- fpi$res$tauint[1]
  result$fpidtauint <- fpi$res$dtauint[1]
  result$sigma <- sigma$res$value[1]
  result$dsigma <- sigma$res$dvalue[1]
  result$ddsigma <- sigma$res$ddvalue[1]
  result$sigmatauint <- sigma$res$tauint[1]
  result$sigmadtauint <- sigma$res$dtauint[1]

  attr(result, "class") <- c("data.frame")
  if(plot == TRUE) {
    cat("Plotting the function...\n")
    correl <- ppcorr(filename, skip=skip, S=S, plotit=T)
    xfit <- c(0:T2)
    yfit <- result$amp*cosh(result$mass*(xfit-T2/2))
    lines(spline(xfit,yfit), col = "red")

  }
  rm(psscar)
  rm(Z)
  rm(W)
  return(invisible(result))  
  
}

