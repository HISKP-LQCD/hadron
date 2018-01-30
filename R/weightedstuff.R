weighted.quantile <- function(x, w, probs = seq(0, 1, 0.25), na.rm = FALSE, names = TRUE) {
  ## n <- length(x)
  ## check:
  ## - if w exists (length(w)>0) => if not, then w <- rep(1,length(x))
  ## - if w is the weight or the inclusion probability
  ## + a test sum(w) > length(x)??????
  if(missing(x)) {
    stop("Usage: weighted.quantile(x,...) vector of values x required")    
  }
  if (missing(w))
    w <- rep(1, length(x))

  if (length(x) != length(w))
    stop("Weights and variable vectors are of unequal length!")
  if (na.rm) {
    valid.obs <- !is.na(x) & !is.na(w)
    x <- x[valid.obs]
    w <- w[valid.obs]
  }
  else if (any(is.na(x)) | any(is.na(w)))
    stop("Missing values and NaN's not allowed if `na.rm' is FALSE")
  if (any((p.ok <- !is.na(probs)) & (probs < 0 | probs > 1)))
    stop("probs outside [0,1]")
  if (na.p <- any(!p.ok)) {
    o.pr <- probs
    probs <- probs[p.ok]
  }
  np <- length(probs)

  ## check if weights are probs or unit correspondences
  ind <- order(x)
  cumw <- cumsum(w[ind])/sum(w)
  x <- x[ind]
  med <- numeric(0)
  for(i in 1:length(probs)) {
    if(probs[i] > 0 & probs[i] < 1) {
      op <- options()
      options(warn=-1)
      med[i] <- max(x[cumw<=probs[i]])
      options(op)
      if(is.infinite(med[i])) {
        warning("could not correctly estimate quantile ", probs[i])
        if(med[i] < 0) med[i] <- min(x)
        else med[i] <- max(x)
      }
    }
    else if (probs[i] == 0) med[i] <- min(x)
    else if (probs[i] == 1) med[i] <- max(x)
  }
  ## the treatment of names
  med
}

weighted.median1 <- function(x,w, na.rm = FALSE) {
  if (missing(w))
    w <- rep(1, length(x))
  if (na.rm) {
    w <- w[i <- !is.na(x)]
    x <- x[i]
    }
  weighted.quantile(x,w,.5)
}

############################################################################
# \name{weighted.median}
# \alias{weighted.median}
#
# \title{Weighted Median Value}
#
# \usage{weighted.median(x, w, na.rm=TRUE, ties=NULL)}
#
# \description{
#   Compute a weighted median of a numeric vector.
# }
#
# \arguments{
#   \item{x}{a numeric vector containing the values whose weighted median is
#            to be computed.}
#   \item{w}{a vector of weights the same length as \code{x} giving the weights
#            to use for each element of \code{x}. Default value is equal weight
#            to all values.}
#   \item{na.rm}{a logical value indicating whether \code{NA} values in
#            \code{x} should be stripped before the computation proceeds.}
#   \item{ties}{a character string specifying how to solve ties between two
#            \code{x}'s that are satisfying the weighted median criteria.
#            Note that at most two values can satisfy the criteria.
#            When \code{ties} is \code{"min"}, the smaller value of the two
#            is returned and when it is \code{"max"}, the larger value is
#            returned.
#            If \code{ties} is \code{"mean"}, the mean of the two values is
#            returned and if it is \code{"both"}, both values are returned.
#            Finally, if \code{ties} is \code{"weighted"} (or \code{NULL}) a
#            weighted average of the two are returned, where the weights are
#            weights of all values \code{x[i] <= x[k]} and \code{x[i] >= x[k]},
#            respectively. Default value is \code{NULL}.}
# }
#
# \value{
#   Returns the weighted median.
# }
#
# \details{
#  For the \code{n} elements \code{x = c(x[1], x[2], ..., x[n])} with positive
#  weights \code{w = c(w[1], w[2], ..., w[n])} such that \code{sum(w) = S},
#  the \emph{weighted median} is defined as the element \code{x[k]} for which
#initial  the total weight of all elements \code{x[i] < x[k]} is less or equal to
#  \code{S/2} and for which the total weight of all elements
#  \code{x[i] > x[k]} is less or equal to \code{S/2} (c.f. [1]).
#
#  If \code{w} is missing then all elements of \code{x} are given the same
#  positive weight. If all weights are zero, \code{NA} is returned.
#
#  When all the weights are the same and \code{ties} is \code{"weighted"} (or
#  \code{NULL}) \code{weighted.median} gives the same result as \code{median}.
#
#  If one or more weights are \code{Inf}, it is the same as these weights
#  have the same weight and the others has zero. This makes things easier for
#  cases where the weights are result of a division with zero.
#
#  The weighted median solves the following optimization problem:
#
#  \deqn{\alpha^* = \arg_\alpha \min \sum_{k=1}{K} w_k |x_k-\alpha|}
#  where \eqn{x=(x_1,x_2,\ldots,x_K)} are scalars and
#  \eqn{w=(w_1,w_2,\ldots,w_K)} are the corresponding "weights" for
#  each individual \eqn{x} value.
# }
#
# \examples{
#   x <- 1:10
#   n <- length(x)
#   median(x)                            # 5.5
#   weighted.median(x)                   # 5.5
#   w <- rep(1, n)
#   weighted.median(x, w)                # 5.5 (default)
#   weighted.median(x, ties="weighted")  # 5.5 (default)
#   weighted.median(x, ties="min")       # 5
#   weighted.median(x, ties="max")       # 6
#
#   # Pull the median towards zero
#   w[1] <- 5
#   weighted.median(x, w)                # 3.5
#   y <- c(rep(0,w[1]), x[-1])           # Only possible for integer weights
#   median(y)                            # 3.5
#
#   # Put even more weight on the zero
#   w[1] <- 8.5
#   weighted.median(x, w)                # 2
#
#   # All weight on the first value
#   w[1] <- Inf
#   weighted.median(x, w)                # 1
#
#   # All weight on the first value
#   w[1] <- 1
#   w[n] <- Inf
#   weighted.median(x, w)                # 10
#
#   # All weights set to zero
#   w <- rep(0, n)
#   weighted.median(x, w)                # NA
# }
#
# \seealso{
#   \code{\link[base]{median}}, \code{\link[base]{mean}} and
#   \code{\link[base]{weighted.mean}}
# }
#
# \reference{
#   [1]  T.H. Cormen, C.E. Leiserson, R.L. Rivest, Introduction to Algorithms,
#        The MIT Press, Massachusetts Institute of Technology, 1989.
# }
#
# \author{Henrik Bengtsson, hb at maths.lth.se with help from
#         Roger Koenker, roger at ysidro.econ.uiuc.edu}
#
#*/#########################################################################
###
weighted.median <- function(x, w, na.rm=TRUE, ties=NULL) {
  if (missing(w))
    w <- rep(1, length(x));

  # Remove values that are NA's
  if (na.rm == TRUE) {
    keep <- !(is.na(x) | is.na(w));
    x <- x[keep];
    w <- w[keep];
  } else if (any(is.na(x)))
    return(NA);

  # Assert that the weights are all non-negative.
  if (any(w < 0))
    stop("Some of the weights are negative; one can only have positive
weights.");

  # Remove values with weight zero. This will:
  #  1) take care of the case when all weights are zero,
  #  2) make sure that possible tied values are next to each others, and
  #  3) it will most likely speed up the sorting.
  n <- length(w);
  keep <- (w > 0);
  nkeep <- sum(keep);
  if (nkeep < n) {
    x <- x[keep];
    w <- w[keep];
    n <- nkeep;
  }

  # Are any weights Inf? Then treat them with equal weight and all others
  # with weight zero.
  wInfs <- is.infinite(w);
  if (any(wInfs)) {
    x <- x[wInfs];
    n <- length(x);
    w <- rep(1, n);
  }

  # Are there any values left to calculate the weighted me
  if (n == 0)
    return(NA);

  # Order the values and order the weights accordingly
  ord <- order(x);
  x <- x[ord];
  w <- w[ord];

  wcum <- cumsum(w);
  wsum <- wcum[n];
  wmid <- wsum / 2;

  # Find the position where the sum of the weights of the elements such that
  # x[i] < x[k] is less or equal than half the sum of all weights.
  # (these two lines could probably be optimized for speed).
  lows <- (wcum <= wmid);
  k  <- sum(lows);

  # Two special cases where all the weight are at the first or the
  # last value:
  if (k == 0) return(x[1]);
  if (k == n) return(x[n]);

  # At this point we know that:
  #  1) at most half the total weight is in the set x[1:k],
  #  2) that the set x[(k+2):n] contains less than half the total weight
  # The question is whether x[(k+1):n] contains *more* than
  # half the total weight (try x=c(1,2,3), w=c(1,1,1)). If it is then
  # we can be sure that x[k+1] is the weighted median we are looking
  # for, otherwise it is any function of x[k:(k+1)].

  wlow  <- wcum[k];    # the weight of x[1:k]
  whigh <- wsum - wlow;  # the weight of x[(k+1):n]
  if (whigh > wmid)
    return(x[k+1]);

  if (is.null(ties) || ties == "weighted") {  # Default!
    (wlow*x[k] + whigh*x[k+1]) / wsum;
  } else if (ties == "max") {
    x[k+1];
  } else if (ties == "min") {
    x[k];
  } else if (ties == "mean") {
    (x[k]+x[k+1])/2;
  } else if (ties == "both") {
    c(x[k], x[k+1]);
  }
}

# routine for computing the unbiased weighted sample variance
# following GSL implementation 
weighted.variance <- function(x, w, na.rm=FALSE) {
  mustar <- weighted.mean(x,w,na.rm)
  v2 <- sum(w^2)
  v1 <- sum(w)
  return ( v1/(v1^2-v2) * sum( w*(x-mustar)^2 ) )
}


get.breaks<-function(x,breaks) {
                                        # if a break computing function name is passed
  if(is.character(breaks)) 
    nbreaks<-do.call(paste("nclass",breaks,sep=".",collapse=""),list(x))
                                        # if breaks is numeric
  if(is.numeric(breaks)) {
                                        # if just the number of breaks is passed
    if(length(breaks) == 1) {
      nbreaks<-breaks
    }
                                        # otherwise assume that breaks specifies the breakpoints
    else return(breaks)
  }
  breakinc<-diff(range(x))/nbreaks
  breaks<-c(min(x),rep(breakinc,nbreaks))
  breaks<-cumsum(breaks)
  return(breaks)
}


## this function is taken from the plotrix package
## and extended to be able to deal with NA
weighted.hist<-function(x, w, breaks="Sturges",col=NULL,plot=TRUE,
                        freq=TRUE,ylim=NA,ylab=NULL,xaxis=TRUE,na.rm=FALSE,...) {
  
  if(missing(x))
    stop("Usage: weighted.hist(x,...) vector of values x required")
  if(missing(w)) w<-rep(1,length(x))
  if (length(x) != length(w))
    stop("Weights and variable vectors are of unequal length!")
  
  if(na.rm) {
    valid.obs <- !is.na(x) & !is.na(w)
    x <- x[valid.obs]
    w <- w[valid.obs]
  }
  else if (any(is.na(x)) | any(is.na(w)))
    stop("Missing values and NaN's not allowed if `na.rm' is FALSE")

  breaks<-get.breaks(x,breaks)
  width<-diff(breaks)
  diffx<-diff(range(x))
  equidist<-sum(width-width[1]) < diffx/1000
  nbreaks<-length(breaks)-1
                                        # make sure that the last break is greater than the maximum value
  lastbreak<-breaks[nbreaks+1]
  breaks[nbreaks+1]<-breaks[nbreaks+1]+diffx/1000
  if(diff(range(breaks)) < diffx)
    warning("Not all values will be included in the histogram")
  counts<-rep(0,nbreaks)
  for(bin in 1:nbreaks)
    counts[bin]<-sum(w[x >= breaks[bin] & x < breaks[bin+1]])
  density<-counts/sum(counts)
  if(freq) {
    if(is.null(ylab)) ylab<-"Frequency"
    heights<-counts
    if(!equidist) 
      warning("Areas will not relate to frequencies")
  }
  else {
    if(!equidist) {
      heights<-density*mean(width)/width
      heights<-heights/sum(heights)
    }
    else heights<-density/width
    if(is.null(ylab)) ylab<-"Density"
  }
  if(plot) {
    if(is.null(col)) col<-par("bg")
    if(is.na(ylim)) ylim<-c(0,1.1*max(heights,na.rm=TRUE))
    mids<-barplot(heights,width=width,col=col,space=0,ylim=ylim,ylab=ylab,...)
    tickpos<-c(mids-width/2,mids[length(mids)]+width[length(width)]/2)
    if(xaxis) axis(1,at=tickpos,labels=signif(c(breaks[1:nbreaks],lastbreak),3))
  }
  else mids<-breaks[-length(breaks)]+width/2
  invisible(list(breaks=breaks,counts=counts,density=density,
                 mids=mids,xname=deparse(substitute(x)),equidist=equidist))
}
