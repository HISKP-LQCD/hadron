#' Reweighting factor container
#'
#' This function `rw()` creates containers for reweighting factors
#' of class `rw`. This class is particularly designed to store reweighting 
#' factors for each gauge configuration, that can be applied on 
#' correlation functions emerging in statistical and quantum field theory
#' simulations. Arithmetic operations are defined for this class in
#' several ways, as well as concatenation and \link{is.rw}.
#'
#' @details
#'
#' And last but not least, these are the fields that are used somewhere in the library but we have not figured out which mixin these should belong to:
#'
#' @family rw constructors
#'
#' @export
rw_orig <- function (.rw = rw(), rw, max_value) {
  stopifnot(inherits(.rw, 'rw'))

  .rw$rw <- rw
  .rw$max_value <- max_value

  class(.rw) <- append(class(.rw), 'rw_orig')
  return (.rw)
}
#' @export
rw <- function () {
  rw <- list()
  class(rw) <- append(class(rw), 'rw')
  return (rw)
}

#' rw metadata mixin constructor
#'
#' @param .rw `rw` object to extend.
#' @param conf.index list of Integers, containing the index of gauge configurations
#'
#' @family cf constructors
#'
#' @export
rw_meta <- function (.rw = rw(), conf.index = NA) {
   stopifnot(inherits(.rw, 'rw'))

   .rw$conf.index = conf.index

   class(.rw) <- append(class(.rw), 'rw_meta')
   return (.rw)
}

#' Checks whether the cf object contains no data
#'
#' @param .rw `rw` object.
#'
#' @examples
#' # The empty cf object must be empty:
#' is_empty.cf(cf())
#'
#' # The sample cf must not be empty:
#' is_empty.cf(samplecf)
is_empty.rw <- function (.rw) {
  setequal(class(.rw), class(rw())) &&
    is.null(names(.rw))
}


#' Arithmetically multiplies two reweighting factors functions
#'
#' @param rw1,rw2 `rw_orig` object.
#' @param nf1,nf2 Integer. Factors that determines the number of flavours you 
#'        have for reweighting factor 1 and 2
#'
#' @return
#' The value is
#' \deqn{exp(C_1*nf1)*exp(C_2*nf2) \,.}
#'
#' @export
multiply.rw <- function(rw1, rw2, nf1 = 1, nf2 = 1) {
  stopifnot(inherits(rw1, 'rw'))
  stopifnot(inherits(rw2, 'rw'))
  stopifnot(inherits(rw1, 'rw_orig'))
  stopifnot(inherits(rw2, 'rw_orig'))
  stopifnot(inherits(rw1, 'rw_meta'))
  stopifnot(inherits(rw2, 'rw_meta'))
  stopifnot(all(dim(rw$rw) == dim(rw2$rw)))
  stopifnot(rw1$conf.index == rw2$conf.index)

  rw <- rw1
  rw$rw <- exp(nf1*rw1$rw + nf2*rw2$rw)
  return(rw)
}

#' Arithmetically multiply reweighting factors
#'
#' @param rw1,rw2 `rw_orig` objects.
#'
#' @export
'*.cf' <- function (rw1, rw2) {
  multiply.rw(cf1, cf2, nf1 = 1, nf2 = 1)
}


#' Checks whether an object is a rw
#'
#' @param x Object, possibly of class `rw`.
#'
#' @export
is.rw <- function (x) {
  inherits(x, "rw")
}


#' Plot a reweighting factor function
#'
#' @param x `rw` object
#' @param neg.vec Numeric vector of length `rw$rw0`. This allows switching the
#' sign for certain gauge configurations such that displaying in
#' log-scale is sensible.
#' @param rep See \code{\link{plotwitherror}}.
#'
#' @inheritParams plotwitherror
#'
#' @export
plot.rw <- function(x, neg.vec = rep(1, times = length(rw$rw)), rep = FALSE, ...) {
  rw <- x
  stopifnot(inherits(rw, 'rw'))

  val <- rw$rw/mean(rw$rw)
  err <- rep(x=0,length(rw$rw))

  df <- data.frame(t = c(1:length(val)),
                   CF = val,
                   Err = err)

  plotwitherror(x = c(1:length(val)), y = neg.vec * df$CF, dy = df$Err, rep = rep, ...)

  return(invisible(df))
}

addStat.rw <- function(rw1, rw2,reverse1=FALSE, reverse2=FALSE) {
  stopifnot(inherits(rw1, 'rw'))
  stopifnot(inherits(rw2, 'rw'))

  if (is_empty.rw(rw1)) {
    return (invisible(rw2))
  }
  if (is_empty.cf(rw2)) {
    return (invisible(rw1))
  }

  stopifnot(inherits(rw1, 'rw_meta'))
  stopifnot(inherits(rw2, 'rw_meta'))

  #stopifnot(identical(rw1$conf.index == rw2$conf.index))
  #stopifnot(dim(rw1$rw)[2] == dim(rw2$rw)[2])

  rw    <- rw1

  if (reverse1 == TRUE){
    rw_1   <- rev(rw1$rw)
    conf_1 <- rev(rw1$conf.index)
  }
  else{
    rw_1   <- rw1$rw
    conf_1 <- rw1$conf.index 
  }
  if (reverse2 == TRUE){
    rw_2   <- rev(rw2$rw)
    conf_2 <- rev(rw2$conf.index)
  }
  else{
    rw_2   <- rw2$rw
    conf_2 <- rw2$conf.index 
  }
  rw_1 <- rw_1*rw$max_value
  rw_2 <- rw_2*rw$max_value

  rw$rw <- rbind(rw_1, rw_2)
  max_value <- max(rw$rw)
  rw$rw <- rw$rw/max_value
  rw$conf.index <- rbind(conf_1, conf_2) 
  rw$max_value  <- max_value 

  return (invisible(rw))
}
