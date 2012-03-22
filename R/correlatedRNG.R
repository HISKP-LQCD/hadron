corrnorm <- function(C, n) {
  U <- chol(C)
  nC <- dim(C)[1]
  N <- (n %% nC) + n
  result <- array(rnorm(N), dim=c(N/nC,nC)) %*% U
  return(invisible(result[1:n]))
}

corrnorm2 <- function(C, n) {
  U <- chol(C)
  nC <- dim(C)[1]
  N <- (n %% nC) + n
  result <- array(rnorm(N), dim=c(N/nC,nC)) %*% U
  return(invisible(result))
}
