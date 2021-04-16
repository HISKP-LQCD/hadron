#' Reweighting factor container
#'
#' @details
#' This function `rw()` creates containers for reweighting factors
#' of class `rw`. This class is particularly designed to store reweighting 
#' factors for each gauge configuration, that can be applied on 
#' correlation functions emerging in statistical and quantum field theory
#' simulations. Note that the reweighting factors acts on bare cf objects
#' and not on bootstrapped correlation function. Multiplication operation 
#' is defined for this class, as well as increasing statistics and \link{is.rw}.
#'
#' @family rw constructors
#'
#' @examples
#' x <- rw()
#'
#' @return
#' returns an object of S3 class \code{rw} derived from a \code{list}
#'
#' @export
rw <- function() {
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
#'
#' @return 
#' returns the input object of class \code{rw} with the meta data mixin added
#'
#'
#' @example
#' newrw <- rw_meta(rw, conf.index=gauge_conf_list)
#'
#' 
#' @export
rw_meta <- function (.rw = rw(), conf.index ) {
   stopifnot(inherits(.rw, 'rw'))

   .rw$conf.index <- conf.index

   class(.rw) <- append(class(.rw), 'rw_meta')
   return (.rw)
}


#' rw metadata mixin constructor
#'
#' @param .rw `rw` object to extend.
#' @param rw vector of reweighting factors
#' @param conf.index list of Integers, containing the index of gauge configurations
#' @param max_value value used to be for normalization
#' @param stochastic_error error coming from the stochastic samples for each gauge conf
#' @family rw constructors
#'
#'
#' @return
#' returns the input object of class \code{rw} with the original data mixin added
#'
#'
#' @examples
#' rw_factor <- rw_orig(rw_factor, rw = rw_data, conf.index=gauge_conf_list, max_value = max(rw_data), stochastic_error=rw_data_error)
#'
#'
#' @export
rw_orig <- function (.rw = rw(), rw, conf.index, max_value, stochastic_error) {
  stopifnot(inherits(.rw, 'rw'))

  .rw$rw <- rw
  .rw$max_value <- max_value
  .rw$conf.index <- conf.index
  .rw$stochastic_error <- stochastic_error


  class(.rw) <- append(class(.rw), 'rw_orig')
  class(.rw) <- append(class(.rw), 'rw_meta')

  ## Norm is needed in order to be able to extend statistics
  ## If the object does not have `rw_norm` it cannot be extended
  class(.rw) <- append(class(.rw), 'rw_norm')
  return (.rw)
}

#' rw_unit
#'
#' generates an 'rw' object with all elements set to 1.
#' 
#' @param .rw  an object of class 'rw'
#' @param conf.index integer vector of configuration indices
#'
#' @return
#' an object of class 'rw' with all elements equal to unity
#'
#' @examples
#' x <- rw_unit(conf.index=c(1,2,3,4))
#' 
#' @export
rw_unit <- function (.rw = rw(), conf.index) {
  stopifnot(inherits(.rw, 'rw'))

  #we set the stochastic error to zero
  .rw <- rw_orig (.rw, rep(1.0, length(conf.index)), conf.index, 1.0, rep(0.0, length(conf.index)))

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
#' @return 
#' returns \code{FALSE} if \code{.rw} contains no data, \code{TRUE} otherwise
#'
#' @export
is_empty.rw <- function (.rw) {
  setequal(class(.rw), class(rw())) &&
    is.null(names(.rw))
}


#' Arithmetically multiplies two reweighting factors 
#'
#' @param rw1,rw2 `rw_orig` objects to be muplitplied
#' @param nf1,nf2 Integer. It is greater than one, when 
#'        we want to apply for more flavours than  we have
#'        actually computed. For example in case of nf=4 
#'        degenerate light quarks, we compute the reweighting
#'        factor for two using Q=D*D^dagger, and we actually want
#'        to apply it for all the four flavours, so we have to 
#'        choose in this case nf1=2,nf2=2
#' Note that in the process of multiplication we do not take
#' care of the normalization of the reweighting objects, since
#' it will drop out, when applying it to a correlation function.
#'
#'
#' @examples
#' data(samplerw)
#' ## the following is not useful, but
#' ## explains the usage
#' product_of_rew_factors <- multiply.rw(rw1=samplerw, rw2=samplerw)
#'
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
 
  #We can multiply reweighting factors, only if 
  #we have the same set of gauge configuration 
  #for both rw factors
  stopifnot(identical( rw1$conf.index, rw2$conf.index ))

  rw <- rw()
  rw$rw <- nf1*rw1$rw * nf2*rw2$rw
  rw$max_value <- NA
  rw$conf.index <- rw1$conf.index 

  #the stochastic error of the product of reweighting factors
  #this is approximately rw1*d(rw2)+rw2*d(rw1)
  rw$stochastic_error <- abs(nf1*rw1$rw*rw2$stochastic_error) + 
                         abs(nf2*rw2$rw*rw1$stochastic_error) 

  class(rw) <- append(class(rw), 'rw_orig')
  class(rw) <- append(class(rw), 'rw_meta')
  
  return(rw)
}

#' Arithmetically multiply reweighting factors
#'
#' @param rw1 `rw_orig` object
#' @param rw2 `rw_orig` object
#'
#' @return 
#' The value is
#' \deqn{rw1 * rw2 \,.}
#
#' @export
'*.rw' <- function (rw1, rw2) {
  multiply.rw(rw1, rw2, nf1 = 1, nf2 = 1)
}


#' Checks whether an object is a rw
#'
#' @param x Object, possibly of class `rw`.
#'
#' @examples
#' x <- c()
#' is.rw(x)
#'
#' @return
#' Returns TRUE if the input object is of class `rw`, FALSE otherwise.
#' 
#' @export
is.rw <- function (x) {
  inherits(x, "rw")
}


#' Plot a reweighting factor function
#'
#' @param x `rw` object
#' @param ... Generic graphical parameter to be passed on to \link{plotwitherror}
#' @param rep See \code{\link{plotwitherror}}.
#'
#' @return
#' Invisibly returns a data.frame with named columns `idx` with the indices of the gauge configurations , `rw` the mean values of the reweighting factor and `Err` its standard error calculated from the available stochastic samples.
#'
#' @export
plot.rw <- function(x, ..., rep = FALSE) {
  rw <- x
  stopifnot(inherits(rw, 'rw'))
  stopifnot(inherits(rw, 'rw_meta'))
  negs <- which(rw$rw < 0)
  neg.vec <- rep(1,times=length(rw$rw))
  neg.vec[negs] <- -1
 
  colneg <- rep("black",times=length(rw$rw))
  colneg[negs] <- "red"

  val <- rw$rw/mean(rw$rw)
  err <- rw$stochastic_error

  df <- data.frame(idx = rw$conf.index,
                   rw = val,
                   Err = err)

  plotwitherror(x = df$idx, y = neg.vec * df$rw, dy = df$Err, col=colneg, rep = rep, ...)

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
#' @param reverse2 `boolean`, see 'reverse1'
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
#'
#' @return an object of class \code{rw} with the statistics of the two input
#' \code{rw} objects combined
#'
#' @export
addStat.rw <- function(rw1, rw2, reverse1=FALSE, reverse2=FALSE) {
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
    stochastic_error_1 <- rev(rw1$stochastic_error)
  }
  else{
    rw_1   <- rw1$rw
    conf_1 <- rw1$conf.index 
    stochastic_error_1 <- rw1$stochastic_error
  }
  if (reverse2 == TRUE){
    rw_2   <- rev(rw2$rw)
    conf_2 <- rev(rw2$conf.index)
    stochastic_error_2 <- rev(rw2$stochastic_error)
  }
  else{
    rw_2   <- rw2$rw
    conf_2 <- rw2$conf.index 
    stochastic_error_2 <- rw2$stochastic_error
  }

  #Restore the proper normalization
  rw_1 <- rw_1*rw1$max_value
  rw_2 <- rw_2*rw2$max_value
  
  stochastic_error_1 <- stochastic_error_1*max_value
  stochastic_error_2 <- stochastic_error_2*max_value

  rw$rw <- c(rw_1, rw_2)
  rw$stochastic_error <- c(stochastic_error_1, stochastic_error_2)

  #New normalization
  max_value <- max(rw$rw)

  #Apply the new normalization
  rw$rw <- rw$rw/max_value
  rw$stochastic_error <- rw$stochastic_error/max_value

  rw$conf.index <- c(conf_1, conf_2) 
  rw$max_value  <- max_value 

  return (invisible(rw))
}
