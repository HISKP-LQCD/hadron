#' Folds the non-empty list with the binary function
foldr1 <- function (f, xs) {
  l <- length(xs)
  stopifnot(l > 0)

  if (l == 1)
    return (xs[[1]])
  else
    return (Reduce(f, xs[2:l], xs[[1]]))
}
