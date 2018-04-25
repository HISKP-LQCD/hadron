\name{LuescherZeta}
\alias{LuescherZeta}
\title{Computes Lueschers Zeta Function}
\description{
  Implementation of Lueschers Zeta Function for general \code{l},
  \code{m} and boost vectors \code{dvec}.
}
\usage{
  LuescherZeta(qsq, l=0, m=0, dvec=c(0,0,0), gamma=1., lambda=1.,
tol=0.000001, verbose=FALSE)
}
\arguments{
  \item{qsq}{
    a vector containing the squared scattering momenta to evaluate the
    Zeta function on
  }
  \item{l}{
    the angular momentum quantum number l
  }
  \item{m}{
    the angular momentum quantum number m
  }
  \item{dvec}{
    integer vector with three elements indicating the total momentum
  }
  \item{gamma}{
    boost factor
  }
  \item{lambda}{
    cut-off parameter for the integration
  }
  \item{tol}{
    error tolerance in numerical integrations
  }
  \item{verbose}{
    makes the function (very much) more verbose
  }
}
\value{
  a vector of complex values of the same length as \code{qsq} containing
  the value of the Zeta function
}
\details{
  The Zeta function is evaluated following arXiv:1107.5023. The
  evaluation is implemented in C using GSL routines.
}
\references{
  See arXiv:1107.5023
}
\examples{
  LuescherZeta(c(0.1, 0.2))
  ## next center of mass example
  Pcm <- c(0, 0, 0)
  q <- 0.1207*24/(2.*pi)
  gamma <- 1.0
  zeta <- Re(LuescherZeta(q^2, gamma = gamma, dvec = Pcm))
  delta <- atan2(pi^(3./2.)*q,zeta)*180./pi
  cat("delta = ", delta, "to compare with 136.6527\n")

  ## moving frame example
  Pcm <- c(0, 0, 1)
  L <- 32
  q <- 0.161*L/(2.*pi)
  E <- 0.440
  Ecm <- 0.396
  gamma <- E/Ecm
  Z00 = Re(LuescherZeta(q^2, gamma=gamma, dvec = Pcm))
  Z20 = Re(LuescherZeta(q^2, gamma=gamma, dvec = Pcm, l = 2))
  delta = atan2(gamma*pi**(3./2.) * q , (Z00 + (2./(q*q*sqrt(5)))*Z20))*180./pi
  cat("delta = ", delta, "to compare with 115.7653\n")

  ## moving frame 2
  Pcm <- c(1, 1, 0)
  L <- 32
  q <- 0.167*L/(2.*pi)
  E <- 0.490
  Ecm <- 0.407
  gamma <- E/Ecm
  Z00 <- Re(LuescherZeta(q^2, gamma=gamma, dvec = Pcm))
  Z20 <- Re(LuescherZeta(q^2, gamma=gamma, dvec = Pcm, l = 2))
  Z22  <- Im(LuescherZeta(q^2, gamma=gamma, dvec = Pcm, l = 2, m = 2))
  Z2_2 <- Im(LuescherZeta(q^2, gamma=gamma, dvec = Pcm, l = 2, m = -2))
  delta <- atan2(gamma*pi^(3./2.) * q , (Z00 - (1./(q^2*sqrt(5)))*Z20 + ((sqrt(3./10.)/(q^2))*(Z22-Z2_2))))*180./pi
  cat("delta = ", delta, "to compare with 127.9882\n")

}
\author{Carsten Urbach, \email{curbach@gmx.de}, Zahnlin Wang, \email{wangzhanlin@hiskp.uni-bonn.de}}
