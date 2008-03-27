readcmicor <- function(filename) {
  return(invisible(read.table(filename, header=F,
                              colClasse=c("integer","integer","integer","numeric","numeric","integer"))))
}

readhlcor <- function(filename) {
  return(invisible(read.table(filename, header=F,
                              colClasse=c("integer", "integer","integer","integer","numeric","numeric","numeric","numeric","integer"))))
}


