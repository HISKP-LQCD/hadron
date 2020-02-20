#' apply parameter-free Gasser-Leutwyler FSE corrections
#' @details Apply FSE correction of a meson mass and decay constant
#' based on eqs. 8 and 9 of Colangeo, Duerr, Haefeli [Nucl.Phys.B 721 (2005)]
#' going back to Gasser and Leutwyler. Note that we employ the normalisation
#' in which f_pi = 130.4.
#' Note that \code{xi} and \code{lambda} below can be either of length one
#' or of the same length as the vectors \code{Mps_V} and \code{fps_V}. When
#' applying the correction to bootstrap samples, it might make sense
#' to pass the original data values for \code{xi} and \code{lambda}.
#' @param Mps_V Numeric vector, meson masses in any type of units. 
#' @param fps_V Numeric vector, meson decay constants in the same units as
#' \code{Mps_V}.
#' @param xi Numeric vector, the ratio \eqn{M_{ps}^2 / (4 \pi f_\pi)^2}. Either of length
#' one, in which case it is used for all \code{Mps_V} and \code{fps_V} or of the same
#' length as \code{Mps_V} and \code{fps_V}. 
#' @param lambda Numeric, the dimensionless product \eqn{M_{ps} \cdot L}. 
#' @return Data frame with two columns containing the correction factors \code{K_mps}
#' and \code{K_fps} for which we have:
#' \deqn{ M_{ps}(L) = M_{ps}(\infty) \cdot K_{M_{ps} + \mathcal{O}(M_{ps} \xi^2) }
#' \deqn{ f_{ps}(L) = f_{ps}(\infty) \cdot K_{f_{ps} + \mathcal{O}(f_{ps} \xi^2) }
#' @export
gl_fse_corr <- function(Mps_V, fps_V, xi, lambda){
  stopifnot( length(xi) == length(lambda) )
  stopifnot( length(Mps_V) == length(fps_V) )
  if( length(xi) > 1 ){
    stopifnot( length(xi) == length(Mps_V) )
  }

  g1_tilde <- g1(lambda)
  return(data.frame(K_Mps = (1+0.5*xi*g1_tilde), K_fps = (1-2*xi*g1_tilde)))
}

