require(tikzDevice)


tikz('deltaEovL.tex', standAlone = TRUE, width=5, height=5)



dev.off()
tools::texi2dvi('deltaEovL.tex',pdf=T)
