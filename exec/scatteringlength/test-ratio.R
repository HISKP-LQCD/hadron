T <- 64
L <- 32
type <- "log"

E0 <- 0.1415
E1 <- 0.2420
deltaE <- 0.0050
Epipi <- E0+E1+deltaE
ddE <- E1-E0

t <- c(0:(T/2))

Cpipi <- cosh(ddE*(t-T/2)) + cosh(Epipi*(t-T/2))
Cpi0 <- cosh(E0*(t-T/2))
Cpi1 <- cosh(E1*(t-T/2))

t <- c(0:(T/2-1))
tt <- t+1
R1 <- ((Cpipi[tt]*exp(ddE*t)-Cpipi[tt+1]*exp(ddE*(t+1)))*exp(-ddE*t) )/( Cpi0[tt]^2 - Cpi0[tt+1]^2 )
R2 <- ((Cpipi[tt]*exp(ddE*t)-Cpipi[tt+1]*exp(ddE*(t+1)))*exp(-ddE*t) )/((Cpi0[tt]*Cpi1[tt]*exp(ddE*t)-Cpi0[tt+1]*Cpi1[tt+1]*exp(ddE*(t+1))) )
R3 <- ((Cpipi[tt]*exp(ddE*t)-Cpipi[tt+1]*exp(ddE*(t+1)))*exp(-ddE*t) )/( Cpi0[tt]*Cpi1[tt] - Cpi0[tt+1]*Cpi1[tt+1] )
R4 <- ((Cpipi[tt]*exp(ddE*t)-Cpipi[tt+1]*exp(ddE*(t+1)))*exp(-ddE*t) )/( (Cpi0[tt]^2 - Cpi0[tt+1]^2)*exp(-ddE) )
R1[T/2+1] <- NA
R2[T/2+1] <- NA
R3[T/2+1] <- NA
R4[T/2+1] <- NA

mR1 <- effectivemass.cf(cf=R1, Thalf=T/2, type=type)
mR2 <- effectivemass.cf(cf=R2, Thalf=T/2, type=type)
mR3 <- effectivemass.cf(cf=R3, Thalf=T/2, type=type)
mR4 <- effectivemass.cf(cf=R4, Thalf=T/2, type=type)

#tt <- c(0:(T/2))
#plot(tt, R1, log="y")
#points(tt, R2, col="red")

tikzfiles <- tikz.init(paste("testdata", sep=""),width=6,height=5)
plot(NA, ylim=c(0.9,1.1), xlim=c(1,26), xlab=c("$t/a$"), ylab=c("$\\delta E_\\mathrm{ratio}/\\delta E_\\mathrm{true}$"))
points((mR1-E1+E0)/deltaE, pch=21, col="black", bg="black")
points(c(0:(T/2-1)), (mR2-E1+E0)/deltaE, col="red", bg="red", pch=22)
points(c(0:(T/2-1)), mR3/deltaE, col="blue", bg="blue", pch=23)
#points(c(0:(T/2-1)), (mR4-E1+E0)/deltaE, col="darkgreen")
abline(h=1)
legend("topleft", legend=c("$R_1$", "$R_2$", "$R_3$"), col=c("black", "red", "blue"), pt.bg=c("black", "red", "blue"), bty="n", pch=c(21,22,23))
              
tikz.finalize(tikzfiles=tikzfiles,clean=TRUE, crop=FALSE)

