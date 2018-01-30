source("fit_finite_range.R")

res.finite.range.fit <- fit.finite.range(typ="-ratio")

save(res.finite.range.fit, file="res-finite-range-fit-ratio.Rdata")

rm(res.finite.range.fit)


res.finite.range.fit.qcotdelta <- fit.finite.range.qcotdelta(typ="-ratio")

save(res.finite.range.fit.qcotdelta, file="res-finite-range-fit-qcotdelta-ratio.Rdata")
