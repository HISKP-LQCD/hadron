#' apply parameter-free Gasser-Leutwyler FSE corrections
#' @details Apply FSE correction of a meson mass and decay constant
#' based on eqs. 8 and 9 of Colangeo, Duerr, Haefeli [Nucl.Phys.B 721 (2005)]
#' going back to Gasser and Leutwyler. Note that we employ the normalisation
#' in which f_pi = 130.4.
#' @param xi Numeric vector, the ratio \eqn{M_{ps}^2 / (4 \pi f_\pi)^2}. 
#' @param lambda Numeric, the dimensionless product \eqn{M_{ps} \cdot L}. 
#' @return Data frame with two columns containing the correction factors \code{K_mps}
#' and \code{K_fps} for which we have:
#' \deqn{ M_{ps}(L) = M_{ps}(\infty) \cdot K_{M_{ps} + \mathcal{O}(M_{ps} \xi^2) }
#' \deqn{ f_{ps}(L) = f_{ps}(\infty) \cdot K_{f_{ps} + \mathcal{O}(f_{ps} \xi^2) }
#' @export
gl_fse_corr <- function(xi, lambda){
  stopifnot( length(xi) == length(lambda) )
  g1_tilde <- g1(lambda)
  return(data.frame(K_Mps = (1+0.5*xi*g1_tilde), K_fps = (1-2*xi*g1_tilde)))
}

