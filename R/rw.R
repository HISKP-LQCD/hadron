#' Reweighting factor container
#'
#' This function `rw()` creates containers for reweighting factors
#' of class `rw`. This class is particularly designed to store reweighting 
#' factors for each gauge configuration, that can be applied on 
#' correlation functions emerging in statistical and quantum field theory
#' simulations. Multiplication operation is defined for this class,
#' as well as for increasing statistics and \link{is.rw}.
#'
#' @details
#'#'
#' @family rw constructors
#'
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
#' @family rw constructors
#'
#' @export
rw_meta <- function (.rw = rw(), conf.index ) {
   stopifnot(inherits(.rw, 'rw'))

   .rw$conf.index <- conf.index

   class(.rw) <- append(class(.rw), 'rw_meta')
   return (.rw)
}


#' @export
rw_orig <- function (.rw = rw(), rw, conf.index, max_value) {
  stopifnot(inherits(.rw, 'rw'))

  .rw$rw <- rw
  .rw$max_value <- max_value
  .rw$conf.index <- conf.index

  class(.rw) <- append(class(.rw), 'rw_orig')
  class(.rw) <- append(class(.rw), 'rw_meta')

  ## Norm is needed in order to be able to extend statistics
  ## If the object does not have `rw_norm` it cannot be extended
  class(.rw) <- append(class(.rw), 'rw_norm')
  return (.rw)
}

#' @export
rw_unit <- function (.rw = rw(), conf.index) {
  stopifnot(inherits(.rw, 'rw'))

  .rw$conf.index <- conf.index
  .rw$rw <- rep(1.0,length(conf.index))
  .rw$max_value <- 1.0

  class(.rw) <- append(class(.rw), 'rw_orig')
  class(.rw) <- append(class(.rw), 'rw_meta')
  class(.rw) <- append(class(.rw), 'rw_norm')
  return (.rw)
}

#' Checks whether the cf object contains no data
#'
#' @param .rw `rw` object.
#'
#' @examples
#' # The empty rw object must be empty:
#' is_empty.rw(rw())
#'
is_empty.rw <- function (.rw) {
  setequal(class(.rw), class(rw())) &&
    is.null(names(.rw))
}


#' Arithmetically multiplies two reweighting factors 
#'
#' @param rw1,rw2 `rw_orig` objects to be muplitplied
#' @param nf1,nf2 Integer. Factors that determines the number of flavours we
#'        have for reweighting factor 1 and 2, the default value is 1, because
#'        usually we compute the reweighting factors for Q: the product of 
#'        up and down or strange and charm determinant, there are no additional 
#'        terms in the sea determinant that have to be reweighted with the same
#'        rw factor. In case we compute the rw1 factor only for the up quark, we 
#'        for example, we have to set nf1=2 to obtain the rw factor
#'        for the light determinant.  
#'
#' @return
#' The value is
#' \deqn{rw_1*nf1*rw_2*nf2) \,.}
#'
#' @export
multiply.rw <- function(rw1, rw2, nf1 = 1, nf2 = 1) {

  stopifnot(inherits(rw1, 'rw'))
  stopifnot(inherits(rw2, 'rw'))
 
  stopifnot(inherits(rw1, 'rw_orig'))
  stopifnot(inherits(rw2, 'rw_orig'))

  stopifnot(inherits(rw1, 'rw_meta'))
  stopifnot(inherits(rw2, 'rw_meta'))
 
  stopifnot(identical( rw1$conf.index, rw2$conf.index ))

  rw <- rw()
  rw$rw <- nf1*rw1$rw * nf2*rw2$rw
  rw$max_value <- NA
  rw$conf.index <- rw1$conf.index 

  class(rw) <- append(class(rw), 'rw_orig')
  class(rw) <- append(class(rw), 'rw_meta')
  
  return(rw)
}

#' Arithmetically multiply reweighting factors
#'
#' @param rw1,rw2 `rw_orig` objects.
#'
#' @export
'*.rw' <- function (rw1, rw2) {
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

meanval <- function(d, w) {return(mean(d$rw[w]))}

#' Plot a reweighting factor function
#'
#' @param x `rw` object
#' @param rep See \code{\link{plotwitherror}}.
#'
#' @inheritParams plotwitherror
#'
#' @export
plot.rw <- function(x, rep = FALSE, ...) {
  rw <- x
  stopifnot(inherits(rw, 'rw'))
  negs <- which(rw$rw < 0)
  neg.vec <- rep(1,times=length(rw$rw))
  neg.vec[negs] <- -1
 
  colneg <- rep("black",times=length(rw$rw))
  colneg[negs] <- "red"

# Estimating the statistical error on the reweighting factor

  rw_boot <- boot(rw, meanval, R=1500)
  rw_error <- sd(rw_boot)

  val <- rw$rw/mean(rw$rw)
  err <- rep(x=rw_error,length(rw$rw))

  df <- data.frame(t = c(1:length(val)),
                   CF = val,
                   Err = err)

  plotwitherror(x = c(1:length(val)), y = neg.vec * df$CF, dy = df$Err, col=colneg, rep = rep, ...)

  return(invisible(df))
}

#' Combine reweighting factors from different replicas
#' 
#' @param rw1 `rw` object: reweighting factor for replicum A
#' @param rw2 `cf` object: reweighting factor for replicum B
#' @param reverse1 `boolean` After the bifurcation point one of
#'                           the replicas (chain of reweighting 
#'                           factors in simulation time) has  
#'                           to be reversed.
#' @param reverse2 `boolean`
#'
#' @examples
#' # Suppose we have reweighting factors in replicum A from 0 to 500
#' # in steps of 4 and in replicum B from 4 to 500 in steps of 4.
#' # To combined the two replicas we have to use
#'
#' #addStat.rw(rw_replicumB, rw_replicumA, TRUE, FALSE)
#'
#' # which means
#' # combined=(rw500 from B, rw496 from B,...,rw004 from B, rw000 from A, ..
#' # rw500 from A) 
#' @export
addStat.rw <- function(rw1, rw2,reverse1=FALSE, reverse2=FALSE) {
  stopifnot(inherits(rw1, 'rw'))
  stopifnot(inherits(rw2, 'rw'))

  if (is_empty.rw(rw1)) {
    return (invisible(rw2))
  }
  if (is_empty.cf(rw2)) {
    return (invisible(rw1))
  }

  stopifnot(inherits(rw1, 'rw_orig'))
  stopifnot(inherits(rw2, 'rw_orig'))


  stopifnot(inherits(rw1, 'rw_meta'))
  stopifnot(inherits(rw2, 'rw_meta'))
 
  stopifnot(inherits(rw1, 'rw_norm'))
  stopifnot(inherits(rw2, 'rw_norm'))

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
  rw_1 <- rw_1*rw1$max_value
  rw_2 <- rw_2*rw2$max_value
  
  rw$rw <- c(rw_1, rw_2)
  max_value <- max(rw$rw)
  rw$rw <- rw$rw/max_value
  rw$conf.index <- c(conf_1, conf_2) 
  rw$max_value  <- max_value 

  return (invisible(rw))
}
