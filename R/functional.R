#' Folds the non-empty list with the binary function
#'
#' A right fold without the need for a neutral element. Does not work with
#' empty lists.
#'
#' @param f `function`. A binary function that takes two elements of the type
#'   contained in `xs` and returns another such element.
#' @param xs `list` or vector. Homogenious list or vector of elements.
#'
#' There is a `Reduce` function in base R that does left and right folds. It
#' always needs a starting element, which usually is the neutral element with
#' respect to the binary operation. We do not want to specify such a neutral
#' element for certain operations, like `+.cf`. Still a functional programming
#' style should be supported such that one can use maps and folds.
#'
#' @export
#'
#' @examples
#' # We generate some random numbers.
#' numbers <- rnorm(10)
#'
#' # The sum is easiest computed with the `sum` function:
#' sum(numbers)
#'
#' # If we wanted to implement `sum` ourselves, we can use a right fold to do
#' # so:
#' Reduce(`+`, numbers, 0.0)
#'
#' # With this new function we do not need a neutral element any more, but give
#' # up the possibility to fold empty lists.
#' foldr1(`+`, numbers)
foldr1 <- function (f, xs) {
  l <- length(xs)
  stopifnot(l > 0)

  if (l == 1)
    return (xs[[1]])
  else
    return (Reduce(f, xs[2:l], xs[[1]]))
}
