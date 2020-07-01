#' tikz.init
#' 
#' Convenience Functions for \code{tikzDevice}
#'
#' @description
#' initialize and finalize a \code{tikzDevice} and carry out optional
#' post-processing
#' 
#' @param basename the base of the files which will be used by
#' \code{tikzDevice}, e.g. "basename" -> "basename.pdf", etc.
#' @param ...  optional arguments which are passed to \code{tikz}, see
#' \code{\link[tikzDevice:tikz]{tikzDevice::tikz}}
#' @param standAlone A logical value indicating whether the output file should
#' be suitable for direct processing by LaTeX. A value of \code{FALSE}
#' indicates that the file is intended for inclusion in a larger document.
#' @param engine used to specify the LaTex engine. If missing, the standard
#' engine of tikz is used.
#' @return \code{tikz.init} returns a list with character vector members, $pdf,
#' $tex, $aux $log containing the corresponding filenames
#' @author Bartosz Kostrzewa, \email{bartosz.kostrzewa@@desy.de}
#' @keywords file
#'
#' @family tikzutils
#' @examples
#' 
#' \donttest{tikzfiles <- tikz.init("plotname",width=3,height=4)}
#' \donttest{plot(x=c(1:3), y=c(1:3)^2, xlab="$x$", ylab="$y$")}
#' \donttest{tikz.finalize(tikzfiles=tikzfiles, clean=TRUE)}
#' \donttest{file.remove("plotname.pdf")}
#' 
#' @export
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

#' tikz.finalize
#' 
#' Convenience Functions for \code{tikzDevice}
#' 
#' @description
#' initialize and finalize a \code{tikzDevice} and carry out optional
#' post-processing
#' 
#' @param tikzfiles a list with members $pdf, $tex, $aux and $log, returned by
#' \code{tikz.init} which must be passed to \code{tikz.finalize}
#' @param crop boolean indicating whether \code{pdfcrop} should be called on
#' the resulting pdf ( existence of \code{pdfcrop} is checked before the
#' command is called ), default TRUE
#' @param margins margins argument for pdfcrop command, should be passed as a
#' string consisting of one or multiple numbers (e.g. "10" or "10.5 7.5 6.2
#' 10"), default 0
#' @param clean boolean indicating whether temporary files, e.g.
#' "basename.tex", "basename.aux" and "basename.log" should be deleted after
#' the pdf has been generated, default TRUE
#' @author Bartosz Kostrzewa, \email{bartosz.kostrzewa@@desy.de}
#' @keywords file
#'
#' @seealso \code{\link{tikz.init}} 
#' @family tikzutils
#'
#' @return
#' No return value, but the output PDF will be created and cropped.
#' 
#' @export
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


