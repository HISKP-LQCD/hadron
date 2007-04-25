readcmicor <- function(filename) {
  return(invisible(read.table(filename, header=F,
                              colClasse=c("integer","integer","integer","numeric","numeric","numeric"))))
}


