if (!requireNamespace("tinytex", quietly = TRUE)) {
    if (!requireNamespace("devtools", quietly = TRUE)) {
        install.packages("jsonlite", repos = "https://cran.rstudio.com/", quiet = TRUE)
        install.packages("httr", repos = "https://cran.rstudio.com/", quiet = TRUE)
        install.packages("memoise", repos = "https://cran.rstudio.com/", quiet = TRUE)
        install.packages("devtools", repos = "https://cran.rstudio.com/", quiet = TRUE)
        cat('devtools installed\n')
    }
  devtools::install_github(c('yihui/tinytex'), quiet = TRUE)
  tinytex::install_tinytex()
}
