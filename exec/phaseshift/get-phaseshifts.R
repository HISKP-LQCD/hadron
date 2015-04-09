source("parameters.R")

source("../phaseshift.pipiswave.R")
phaseshift.pipi.swave(PC="pc1", tp="TP0", boot.R=boot.R, boot.l=boot.l, L=L, T=T, dvec=c(0,0,0), debug=TRUE, p1=c(0,0,0), p2=c(0,0,0))
phaseshift.pipi.swave(PC="pc2", tp="TP0", boot.R=boot.R, boot.l=boot.l, L=L, T=T, dvec=c(0,0,0), debug=TRUE, p1=c(0,0,0), p2=c(0,0,0))
phaseshift.pipi.swave(PC="pc1", tp="TP1", boot.R=boot.R, boot.l=boot.l, L=L, T=T, dvec=c(0,0,1), debug=TRUE, p1=c(0,0,1), p2=c(0,0,0))
phaseshift.pipi.swave(PC="pc1", tp="TP2", boot.R=boot.R, boot.l=boot.l, L=L, T=T, dvec=c(0,1,1), debug=TRUE, p1=c(0,1,1), p2=c(0,0,0))
phaseshift.pipi.swave(PC="pc1", tp="TP3", boot.R=boot.R, boot.l=boot.l, L=L, T=T, dvec=c(1,1,1), debug=TRUE, p1=c(1,1,1), p2=c(0,0,0))

