\name{samplecf}
\alias{samplecf}
\docType{data}
\title{
Sample cf data
}
\description{
Sample data for a correlation function for a 24 cube times 48 lattice
QCD simulation representing a pion propagation. It is stored in form of
an object of class \code{cf}, which is derived from \code{list}.
}
\usage{data("samplecf")}
\format{
  The format is:
List of 15
 $ cf               : num [1:1018, 1:25] 521 533 532 531 561 ...
 $ icf              : num [1:1018, 1:25] 521 533 532 531 561 ...
 $ Time             : num 48
 $ nrStypes         : num 1
 $ nrObs            : num 1
 $ boot.samples     : logi TRUE
 $ jackknife.samples: logi FALSE
 $ symmetrised      : logi TRUE
 $ boot.R           : num 1500
 $ boot.l           : num 2
 $ seed             : num 1442556
 $ sim              : chr "geom"
 $ cf0              : num [1:25] 519 375 274 221 185 ...
 $ cf.tsboot        :List of 11
  ..$ t0       : num [1:25] 519 375 274 221 185 ...
  ..$ t        : num [1:1500, 1:25] 521 518 520 519 519 ...
  ..$ R        : num 1500
  ..$ data     : num [1:1018, 1:25] 521 533 532 531 561 ...
  ..$ seed     : int [1:626] 403 624 -867935848 1692432057 -1535150298 -1438296209 912697060 1838233749 1438572626 999279531 ...
  ..$ statistic:function (x)  
  ..$ sim      : chr "geom"
  ..$ n.sim    : int 1018
  ..$ call     : language tsboot(tseries = cf$cf, statistic = function(x) {     return(apply(x, MARGIN = 2L, FUN = mean)) ...})
  ..$ l        : num 2
  ..$ endcorr  : logi TRUE
  ..- attr(*, "class")= chr "boot"
  ..- attr(*, "boot_type")= chr "tsboot"
 $ tsboot.se        : num [1:25] 1.001 0.615 0.572 0.537 0.499 ...
 - attr(*, "class")= chr [1:2] "cf" "list"
}
%%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%%}
%%\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
%%}
%%\references{
%%  ~~ possibly secondary sources and usages ~~
%%}
\examples{
data(samplecf)
plot(samplecf)
}
\keyword{datasets}
