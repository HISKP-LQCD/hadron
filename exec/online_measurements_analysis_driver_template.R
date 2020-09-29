skip <- 0
stat_skip <- 0

L <- 48
Time <- 96
beta <- 1.778
type <- "iwa"
kappa <- 0.1394267
mul <- 0.0025
csw <- 1.69
musigma <- 0.1246864
mudelta <- 0.1315052

evals <- 5
cg_col <- 20

boot.l <- 2
boot.R <- 1000

t1 <- 15
t2 <- 42

for( path in c("cB211a.25.48","cB211b.25.48") ){
  analysis_online(type=type, beta=beta, L=L, Time=Time, kappa=kappa, mul=mul,
               t1=t1, t2=t2, csw=csw, musigma=musigma, mudelta=mudelta,
               skip=skip,
               stat_skip=stat_skip,
               addon="", title=FALSE,
               evals=evals,
               cg_col=cg_col, 
               plotsize=4.5,
               rundir=path,
               boot.l=boot.l, boot.R=boot.R, method="all",
               acc=TRUE)
  try(analysis_gradient_flow(path=path,
                         basename="gradflow",
                         outputbasename=path,
                         pl=TRUE,
                         read.data=TRUE,
                         skip=skip/4,
                         scale=4,
                         dbg=FALSE))
}

