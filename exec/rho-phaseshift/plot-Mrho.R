ens <- c("A40.32", "A30.32", "A60.24", "A80.24", "B55.32")

pch <- c(21, 22, 23, 24, 25)
col <- c("red", "darkblue", "darkgreen", "purple", "orange", "black")

r0data <- read.table("r0.dat")

tikzfiles <- tikz.init(basename=paste("MrhoMpisq", sep=""), width=5., height=4.)

plot(NA, xlim=c(0,1.2), ylim=c(1,2.5), xlab=c("$(r_0M_\\pi)^2$"), ylab=c("$r_0 M_\\rho$"))

for(i in c(1:length(ens))) {
  cat(paste(ens[i], "/pion.Rdata", sep=""), "\n")
  load(paste(ens[i], "/pion.Rdata", sep=""))
  cat(paste(ens[i], "/Mrho-res", ens[i], ".Rdata", sep=""), "\n")
  load(paste(ens[i], "/Mrho-res", ens[i], ".Rdata", sep=""))

  k <- 1
  if(grepl("B", ens[i])) k <- 2
  if(grepl("D", ens[i])) k <- 3

  r0 <- r0data$V1[k]
  dr0 <- r0data$V2[k]
  boot.R <- 5000
  r0.boot <- rnorm(n=boot.R, mean=r0, sd=dr0)

  cat(r0*Mrho.res[1,2], "\n")
  cat((r0*pion.matrixfit$opt.res$par[1])^2, "\n")
  plotwitherror(x = (r0*pion.matrixfit$opt.res$par[1])^2,
                dx = sd((r0.boot*pion.matrixfit$opt.tsboot[1,])^2),
                y = r0*Mrho.res[1,2],
                dy = sd(r0.boot*Mrho.res[c(2:(boot.R+1)),2]),
                rep=TRUE, col=col[k], pch=pch[k], bg=col[k])
}

legend("bottomright", legend=c("A ensembles", "B ensemble"), pch=pch[c(1:2)], col=col[c(1:2)], pt.bg=col[c(1:2)], bty="n")

tikz.finalize(tikzfiles=tikzfiles)


tikzfiles <- tikz.init(basename=paste("gMpisq", sep=""), width=5., height=4.)

plot(NA, xlim=c(0,1.2), ylim=c(4,8), xlab=c("$(r_0M_\\pi)^2$"), ylab=c("$r_0 M_\\rho$"))

for(i in c(1:length(ens))) {
  cat(paste(ens[i], "/pion.Rdata", sep=""), "\n")
  load(paste(ens[i], "/pion.Rdata", sep=""))
  cat(paste(ens[i], "/Mrho-res", ens[i], ".Rdata", sep=""), "\n")
  load(paste(ens[i], "/Mrho-res", ens[i], ".Rdata", sep=""))

  k <- 1
  if(grepl("B", ens[i])) k <- 2
  if(grepl("D", ens[i])) k <- 3

  r0 <- r0data$V1[k]
  dr0 <- r0data$V2[k]
  boot.R <- 5000
  r0.boot <- rnorm(n=boot.R, mean=r0, sd=dr0)

  cat(Mrho.res[1,1], "\n")
  cat((r0*pion.matrixfit$opt.res$par[1])^2, "\n")
  plotwitherror(x = (r0*pion.matrixfit$opt.res$par[1])^2,
                dx = sd((r0.boot*pion.matrixfit$opt.tsboot[1,])^2),
                y = Mrho.res[1,1],
                dy = sd(Mrho.res[c(2:(boot.R+1)),1]),
                rep=TRUE, col=col[k], pch=pch[k], bg=col[k])
}

legend("bottomright", legend=c("A ensembles", "B ensemble"), pch=pch[c(1:2)], col=col[c(1:2)], pt.bg=col[c(1:2)], bty="n")
tikz.finalize(tikzfiles=tikzfiles)
