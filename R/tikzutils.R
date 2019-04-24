tikz.init <- function(basename, standAlone = TRUE, engine, ...) {
  havetikz <- requireNamespace("tikzDevice")
  if(!havetikz){
    stop("tikz.init: tikzDevice package was not found!")
  }
  if(missing(engine)) {
    engine <- getOption("tikzDefaultEngine")
  }

  temp <- sprintf("%s.%s",basename,c("tex","pdf","aux","log"))
  tikzfiles <- list(tex=temp[1], pdf=temp[2], aux=temp[3], log=temp[4], standAlone=standAlone)
  tikzDevice::tikz(tikzfiles$tex, standAlone = standAlone, engine=engine, ...)
  tikzfiles
}

tikz.finalize <- function(tikzfiles, crop=TRUE, margins=0, clean=TRUE) {
  dev.off()
  if(tikzfiles$standAlone) {
    tools::texi2dvi(tikzfiles$tex, pdf=T)
    if(crop){
      ## use pdfcrop tool to remove plot borders
      if( Sys.which("pdfcrop") != "" ){
        command <- sprintf("pdfcrop --margins=%s %s %s",margins,tikzfiles$pdf,tikzfiles$pdf)
        system(command)
      }
      else {
        warning(sprintf("tikz_finalize: crop requested for %s but 'pdfcrop' tool not found!",tikzfiles$pdf))
      }
    }
    if(clean){
      ## remove temporary files 
      command <- sprintf("rm %s %s %s", tikzfiles$tex, tikzfiles$log, tikzfiles$aux)
      system(command)
    }
  }
}

# from https://stackoverflow.com/questions/36338629/escaping-special-latex-characters-in-r
escapeLatexSpecials <- function(x) {
  x <- gsub("\\", "$\\backslash$", x, fixed = TRUE)
  x <- gsub("#", "\\\\#", x)
  x <- gsub("$", "\\\\$", x)
  x <- gsub("%", "\\\\%", x)
  x <- gsub("&", "\\\\&", x)
  x <- gsub("~", "\\\\~", x)
  x <- gsub("_", "\\\\_", x)
  x <- gsub("^", "\\\\^", x)
  x <- gsub("\\{", "\\\\{", x)
  x <- gsub("\\}", "\\\\}", x)
  x <- gsub(">", "$>$", x)
  x <- gsub("<", "$<$", x)
  return(x)
}
