source("parameters.R")

source("singlepi.R")

singlepi(t1=14, t2=T/2, T=T, p="p0", srcpath=srcpath)
singlepi(t1=14, t2=T/2, T=T, p="p1", srcpath=srcpath)
singlepi(t1=11, t2=22, T=T, p="p2", srcpath=srcpath)
singlepi(t1=10, t2=15, T=T, p="p3", srcpath=srcpath)
singlepi(t1=10, t2=15, T=T, p="p4", srcpath=srcpath)

source("pipi.R")

energies.pipi(boot.R=boot.R, boot.l=boot.l, redofit=redofit, tp="TP0", N=5, seed=seed, T=T, N.ids=2, t1=c(9,8,5), t2=c(26, 15, 13), srcpath=srcpath)
energies.pipi(boot.R=boot.R, boot.l=boot.l, redofit=redofit, tp="TP1", N=4, seed=seed, T=T, N.ids=1, t1=c(9,8,5), t2=c(24, 15, 13), srcpath=srcpath)
energies.pipi(boot.R=boot.R, boot.l=boot.l, redofit=redofit, tp="TP2", N=5, seed=seed, T=T, N.ids=1, t1=c(8,6,5), t2=c(17, 15, 13), srcpath=srcpath)
energies.pipi(boot.R=boot.R, boot.l=boot.l, redofit=redofit, tp="TP3", N=3, seed=seed, T=T, N.ids=1, t1=c(5,6,6), t2=c(16, 15, 13), srcpath=srcpath)

