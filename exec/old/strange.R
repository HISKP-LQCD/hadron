mucont <- function(j=4) {
  sigma6.20 <- read.table("b6.20_L24T48_new/sigmainterpolated.dat", header=TRUE)
  sigma6.00 <- read.table("b6.00_L16T32_new/sigmainterpolated6.0.dat", header=TRUE)
  sigma5.85 <- read.table("b5.85_L16T32_new/sigmainterpolated5.85.dat", header=TRUE)
  sigma5.7 <- read.table("b5.70_L12T32_new/sigmainterpolated5.7.dat", header=TRUE)
  df <- data.frame(x = c(1/sigma6.20$r0[(j)], 1/sigma6.00$r0[(j)], 1/sigma5.85$r0[(j)]),
                   y= c(sigma6.20$r0muZp[(j)], sigma6.00$r0muZp[(j)], sigma5.85$r0muZp[(j)]))
  fit <- glm(y ~ 1 , data=df, weights=c(1/sigma6.20$dr0muZp[(j)]^2, 1/sigma6.00$dr0muZp[(j)]^2, 1/sigma5.85$dr0muZp[(j)]^2))
  fit2 <- glm(y ~ 1 +x, data=df, weights=c(1/sigma6.20$dr0muZp[(j)]^2, 1/sigma6.00$dr0muZp[(j)]^2, 1/sigma5.85$dr0muZp[(j)]^2))
  df2 <- data.frame(x = c(seq(0.,0.25,0.01)))
  prefit <- predict(fit, newdata=df2, se.fit=TRUE)
  prefit2 <- predict(fit2, newdata=df2, se.fit=TRUE)
  plotwitherror(c(1/sigma6.20$r0[(j)], 1/sigma6.00$r0[(j)], 1/sigma5.85$r0[(j)], 1/sigma5.7$r0[(j)]), c(sigma6.20$r0muZp[(j)], sigma6.00$r0muZp[(j)], sigma5.85$r0muZp[(j)], sigma5.7$r0muZp[(j)]), c(sigma6.20$dr0muZp[(j)], sigma6.00$dr0muZp[(j)], sigma5.85$dr0muZp[(j)], sigma5.7$dr0muZp[(j)]), xlim=c(0.,0.35), xlab="a/r0", ylab="r0*mu_s/Zp", tck=0.01, mgp=c(3,1.,0), las=1, ylim=c(0.1,0.16))
#  lines(spline(df2$x, prefit$fit+prefit$se.fit), col="red")
#  lines(spline(df2$x, prefit$fit-prefit$se.fit), col="red")
#  lines(spline(df2$x, prefit$fit), col="red")
  lines(spline(df2$x, prefit2$fit+prefit2$se.fit), col="blue", lty=2)
  lines(spline(df2$x, prefit2$fit-prefit2$se.fit), col="blue", lty=2)
  lines(spline(df2$x, prefit2$fit), col="blue", lty=2)
  dev.print(device=postscript, file="ms_contextr.eps", horizontal=FALSE, onefile=FALSE, width=8, height=8)
  return(data.frame(ms=prefit2$fit[[1]], dms=prefit2$se.fit[[1]]))
}

mucontfps <- function() {
  mustrange <- read.table(file="mustrange.dat", header=T)
  df <- data.frame(x = 1/mustrange$r0[1:3], y=mustrange$r0mufpsZp[1:3])
  fit2 <- glm(y ~ 1 +x, data=df, weights=1/mustrange$dr0mufpsZp[1:3]^2)
  df2 <- data.frame(x = c(seq(0.,0.25,0.01)))
  prefit2 <- predict(fit2, newdata=df2, se.fit=TRUE)
  plotwitherror(1/mustrange$r0, mustrange$r0mufpsZp, mustrange$dr0mufpsZp, xlim=c(0,0.35), xlab="a/r0", ylab="r0*mu_s/Zp", tck=0.01,
                mgp=c(3,1.,0), las=1, ylim=c(0.06,0.2))
  lines(spline(df2$x, prefit2$fit+prefit2$se.fit), col="blue", lty=2)
  lines(spline(df2$x, prefit2$fit-prefit2$se.fit), col="blue", lty=2)
  lines(spline(df2$x, prefit2$fit), col="blue", lty=2)
  dev.print(device=postscript, file="msfps_contextr.eps", horizontal=FALSE, onefile=FALSE, width=8, height=8)
  return(data.frame(ms=prefit2$fit[[1]], dms=prefit2$se.fit[[1]]))
}
sigmacont <- function(j=1) {
  sigma6.20 <- read.table("b6.20_L24T48_new/sigmainterpolated.dat", header=TRUE)
  sigma6.00 <- read.table("b6.00_L16T32_new/sigmainterpolated6.0.dat", header=TRUE)
  sigma5.85 <- read.table("b5.85_L16T32_new/sigmainterpolated5.85.dat", header=TRUE)
  sigma5.7 <- read.table("b5.70_L12T32_new/sigmainterpolated5.7.dat", header=TRUE)
  df <- data.frame(x = c(1/sigma6.20$r0[(j)], 1/sigma6.00$r0[(j)], 1/sigma5.85$r0[(j)]),
                   y= c(sigma6.20$sigmaZpr0c[(j)], sigma6.00$sigmaZpr0c[(j)], sigma5.85$sigmaZpr0c[(j)]))
  fit <- glm(y ~ 1 , data=df, weights=c(1/sigma6.20$dsigmaZpr0c[(j)]^2, 1/sigma6.00$dsigmaZpr0c[(j)]^2, 1/sigma5.85$dsigmaZpr0c[(j)]^2))
  fit2 <- glm(y ~ 1 +x, data=df, weights=c(1/sigma6.20$dsigmaZpr0c[(j)]^2, 1/sigma6.00$dsigmaZpr0c[(j)]^2, 1/sigma5.85$dsigmaZpr0c[(j)]^2))
  df2 <- data.frame(x = c(seq(0.,0.25,0.01)))
  prefit <- predict(fit, newdata=df2, se.fit=TRUE)
  prefit2 <- predict(fit2, newdata=df2, se.fit=TRUE)
  plotwitherror(c(1/sigma6.20$r0[(j)], 1/sigma6.00$r0[(j)], 1/sigma5.85$r0[(j)], 1/sigma5.7$r0[(j)]), c(sigma6.20$sigmaZpr0c[(j)], sigma6.00$sigmaZpr0c[(j)], sigma5.85$sigmaZpr0c[(j)], sigma5.7$sigmaZpr0c[(j)]), c(sigma6.20$dsigmaZpr0c[(j)], sigma6.00$dsigmaZpr0c[(j)], sigma5.85$dsigmaZpr0c[(j)], sigma5.7$dsigmaZpr0c[(j)]), xlim=c(0.,0.35), xlab="a/r0", ylab="r0*mu_s/Zp", tck=0.01, mgp=c(3,1.,0), las=1, ylim=c(0.25,0.55))
  lines(spline(df2$x, prefit$fit+prefit$se.fit), col="red")
  lines(spline(df2$x, prefit$fit-prefit$se.fit), col="red")
  lines(spline(df2$x, prefit$fit), col="red")
  lines(spline(df2$x, prefit2$fit+prefit2$se.fit), col="blue", lty=2)
  lines(spline(df2$x, prefit2$fit-prefit2$se.fit), col="blue", lty=2)
  lines(spline(df2$x, prefit2$fit), col="blue", lty=1)
  dev.print(device=postscript, file="sigmamatmu_contextr.eps", horizontal=FALSE, onefile=FALSE, width=8, height=8)
  return(data.frame(r0mpssq = sigma6.20$r0mpssq[(j)], sigma = fit2$coeff[[1]], dsigma=prefit2$se.fit[[1]]))
}
