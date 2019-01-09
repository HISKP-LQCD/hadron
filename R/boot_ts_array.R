#' Copy of boot:::ts.array
#'
#' The function `ts.array` of the `boot` package is not exported, yet we depend
#' on it. It is bad to depend on private functions of a library, so we copy it
#' here in the current version such that the `boot` package is free to change
#' it.
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
