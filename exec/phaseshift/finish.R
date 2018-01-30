source("parameters.R")
pdf(file="finish.pdf")

pc <- "pc1"
TP <- "TP0"
source("../summary.R")
sr <- array(c(qsqovmpisq, qcotdeltaovmpi, delta, qsq, qcotdelta, Epi, Epipi, q), dim=c(1,32))

save(sr, file="sr.Rdata")

pc <- "pc2"
TP <- "TP0"
source("../summary.R")
sr <- rbind(sr, c(qsqovmpisq, qcotdeltaovmpi, delta, qsq, qcotdelta, Epi, Epipi, q))

pc <- "pc1"
TP <- "TP1"
source("../summary.R")
sr <- rbind(sr, c(qsqovmpisq, qcotdeltaovmpi, delta, qsq, qcotdelta, Epi, Epipi, q))

pc <- "pc1"
TP <- "TP2"
source("../summary.R")
sr <- rbind(sr, c(qsqovmpisq, qcotdeltaovmpi, delta, qsq, qcotdelta, Epi, Epipi, q))

pc <- "pc1"
TP <- "TP3"
source("../summary.R")
sr <- rbind(sr, c(qsqovmpisq, qcotdeltaovmpi, delta, qsq, qcotdelta, Epi, Epipi, q))

sr

save(sr, file="sr.Rdata")
dev.off()

