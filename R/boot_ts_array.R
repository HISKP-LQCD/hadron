#' Copy of boot:::ts.array
#'
#' The function `ts.array` of the `boot` package is not exported, yet we depend
#' on it. It is bad to depend on private functions of a library, so we copy it
#' here in the current version such that the `boot` package is free to change
#' it.
#'
#' @param n integer. Length of original data
#' @param n.sim The length of the simulated time series.  Typically this will
#'          be equal to the length of the original time series but there
#'          are situations when it will be larger.  One obvious situation
#'          is if prediction is required.  Another situation in which
#'          ‘n.sim’ is larger than the original length is if ‘tseries’ is
#'          a residual time series from fitting some model to the
#'          original time series. In this case, ‘n.sim’ would usually be
#'          the length of the original time series.
#' @param R A positive integer giving the number of bootstrap replicates
#'          required.
#' @param l If ‘sim’ is ‘"fixed"’ then ‘l’ is the fixed block length used
#'          in generating the replicate time series.  If ‘sim’ is
#'          ‘"geom"’ then ‘l’ is the mean of the geometric distribution
#'          used to generate the block lengths. ‘l’ should be a positive
#'          integer less than the length of ‘tseries’.  This argument is
#'          not required when ‘sim’ is ‘"model"’ but it is required for
#'          all other simulation types.
#' @param sim The type of simulation required to generate the replicate
#'          time series.  The possible input values are ‘"model"’ (model
#'          based resampling), ‘"fixed"’ (block resampling with fixed
#'          block lengths of ‘l’), ‘"geom"’ (block resampling with block
#'          lengths having a geometric distribution with mean ‘l’) or
#'          ‘"scramble"’ (phase scrambling).
#' @param endcorr boolean. whether or not to apply end correction
boot_ts_array <- function (n, n.sim, R, l, sim, endcorr) 
{
    endpt <- if (endcorr) 
        n
    else n - l + 1
    cont <- TRUE
    if (sim == "geom") {
        len.tot <- rep(0, R)
        lens <- NULL
        while (cont) {
            temp <- 1 + rgeom(R, 1/l)
            temp <- pmin(temp, n.sim - len.tot)
            lens <- cbind(lens, temp)
            len.tot <- len.tot + temp
            cont <- any(len.tot < n.sim)
        }
        dimnames(lens) <- NULL
        nn <- ncol(lens)
        st <- matrix(sample.int(endpt, nn * R, replace = TRUE), 
            R)
    }
    else {
        nn <- ceiling(n.sim/l)
        lens <- c(rep(l, nn - 1), 1 + (n.sim - 1)%%l)
        st <- matrix(sample.int(endpt, nn * R, replace = TRUE), 
            R)
    }
    list(starts = st, lengths = lens)
}
