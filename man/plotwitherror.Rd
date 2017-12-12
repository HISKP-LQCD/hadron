\name{plotwitherror}
\alias{plotwitherror}
\title{Plot Command For XY Plots With Error Bars}
\description{
  Plot command for XY scatterplots based on plot and points which provides support for
  multiple, non-symmetric error bars. Error bars are drawn as vertical or horizontal
  lines originating from the point with narrow, perpendicular lines at the end of the
  error bar (end caps). When multiple errors are drawn, the width of the perpendicular line
  increases from the innermost error bar to the outermost one. Different summation
  methods for the individual errors are supported.
}
\usage{
plotwitherror(x, y, dy, ylim, dx, xlim, mdx, mdy, 
              errsum.method="linear.quadrature",
              rep=FALSE, col="black", ...)
}
\arguments{
  \item{x}{
    vector of x coordinates
  }
  \item{y}{
    vector of y coordinates
  }
  \item{dy}{
    one of:
    \itemize{
      \item{Vector of errors on y coordinates.}
      \item{Array, matrix or data frame if multiple error bars are to be drawn, such that each
            column refers to one error. The individual errors should be provided
            as is, because they are summed internally to draw the final error bars.
            A given column can also be provided with 0 entries, in which case the error
            bar will be drawn, but it will have zero length, such that only the end
            caps for this error will be visible.}
    }
  }
  \item{dx}{
    Same as \code{dy}, but for the x coordinate.
  }
  \item{mdx}{
    Support for non-symmetric error bars. Same as \code{dx}, but for errors in the
    negative x-direction. Errors should be provided as positive numbers, the correct
    sign will be added internally. If not provided, \code{dx} is used as a symmetric
    error.
  }
  \item{mdy}{
    Same as \code{mdx} but for the y coordinate.
  }
  \item{errsum.method}{
    Determines how the invidual errors should be summed for display purposes. Valid
    argument values are:
    \itemize{
      \item{"linear"}{
        \itemize{
          \item{Individual errors are summed linearly, such that the distance from the
                 point to the \eqn{i}'th error bar, \eqn{l_i}, is
                   \deqn{ l_i = \sum_{j=1}^i e_j }
                 Hence, the third error bar, for example, would be located at 
                   \deqn{ l_3 = e_1 + e_2 + e_3 }
                 while the second error bar is at
                   \deqn{ l_2 = e_1 + e_2 }
          }
        }
      }
      \item{"quadrature"}{
        \itemize{
          \item{Individual errors are summed in quadrature and error bars are drawn at
                the fractional position according to the following formula:
                \deqn{ l_{max} = \sqrt{ \sum_{j=1}^{max} e_j^2 } }
                \deqn{ l_i = \sum_{j=1}^i e_j^2 / l_{max} }
          }
        }
      }
      \item{"linear.quadrature"}{
        \itemize{
          \item{Errors are summed as for "linear", but the total error summed in
                quadrature is also indicated as an end cap of triple line width}
        }
      }
    }
  }
  \item{rep}{
    If set to \code{TRUE}, operate like "replot" in gnuplot. Allows
    adding points with error bars to the current plot. Switches the underlying
    plotting routine from \code{plot} to \code{points}.
  }
}
\value{
  a plot with error bars is drawn on the current device
}
\author{Carsten Urbach, \email{urbach@hiskp.uni-bonn.de} \cr
        Bartosz Kostrzewa, \email{bartosz.kostrzewa@desy.de}}

