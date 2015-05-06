

tikzfiles <- tikz.init(paste("compare-em", sep=""),width=5,height=5)

plot(NA, ylim=c(0.05,0.15), xlim=c(5,26.5), main="", xlab=c("$t/a$"), ylab=c("$aE_\\mathrm{eff}$"))

plotwitherror(x=c(0:30), effmass[1,], apply(effmass, 2, sd), col="black", rep=TRUE)
plotwitherror(x=c(0:30)+0.2, effmassR2[1,], apply(effmassR2, 2, sd), rep=TRUE, col="blue")
plotwitherror(x=c(0:30)-0.2, effmass3[1,], apply(effmass3, 2, sd), rep=TRUE, col="red")
plotwitherror(x=c(0:30)+0.4, effmassR3[1,], apply(effmassR3, 2, sd), rep=TRUE, col="darkgreen")
legend("bottomleft", legend=c("$R_1, t_0=1$", "$R_2, t_0=1$", "$R_1, t_0=8$", "$R_2, t_0=8$"), bty="n", col=c("black", "blue", "red", "darkgreen"), pch=c(21,21,21,21))
legend("topleft", legend=c("$\\delta E$ from $R_1, t_0=1$"), lty=c(1), col=c("black"), bty="n")
abline(h=dEres$opt.tsboot[1,2], col="black")

tikz.finalize(tikzfiles=tikzfiles,clean=TRUE, crop=FALSE)

tikzfiles <- tikz.init(paste("compare-ratio", sep=""),width=5,height=5)

tt <- c((t0+1):(T/2))
shift <- 0.5
if(mf) shift <- 0.

plot(NA, xlab="$t/a$", ylab="$R_i$", log="y", ylim=c(0.05,1.1), xlim=c(8,32))
plotwitherror(x=tt-shift, y=dEres$Rpipi.tsboot[1,tt]/dEres$Rpipi.tsboot[1,9], dy=dEres$dRpipi[tt]/dEres$Rpipi.tsboot[1,9], rep=TRUE, col="black")
plotwitherror(x=tt-shift+0.2, y=dEresR2$Rpipi.tsboot[1,tt]/dEresR2$Rpipi.tsboot[1,9], dy=dEresR2$dRpipi[tt]/dEresR2$Rpipi.tsboot[1,9], rep=TRUE, col="blue")
plotwitherror(x=tt-shift-0.2, y=dEres3$Rpipi.tsboot[1,tt]/dEres3$Rpipi.tsboot[1,9], dy=dEres3$dRpipi[tt]/dEres3$Rpipi.tsboot[1,9], rep=TRUE, col="red")
plotwitherror(x=tt-shift+0.4, y=dEresR3$Rpipi.tsboot[1,tt]/dEresR3$Rpipi.tsboot[1,9], dy=dEresR3$dRpipi[tt]/dEresR3$Rpipi.tsboot[1,9], rep=TRUE, col="darkgreen")

legend("bottomleft", legend=c("$R_1, t_0=1$", "$R_2, t_0=1$", "$R_1, t_0=8$", "$R_2, t_0=8$"), col=c("black", "blue", "red", "darkgreen"), pch=c(21,21, 21, 21), bty="n")

tikz.finalize(tikzfiles=tikzfiles,clean=TRUE, crop=FALSE)
