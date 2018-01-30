## This file contains routines needed in the context of the Luescher method
## for investigating unstable particles in a finite box with
## Euclidean metric

## compute the lattice scattering momentum \tilde q

## default version with lattice dispersion relation
compute.qtildesq <- function(E, dvec, mpi, L) {
  ## center of mass energy
  ## cosh(Ecm) = cosh(E) - 2 sum sin^2(Pi/2)
  Ecm <- acosh(cosh(E) - 2*sum(sin(pi*dvec/L)^2))
  ## scattering momentum
  ## cosh(Ecm/2) = 2 sin^2(q*/2) + cosh(mpi)
  q = 2*asin(sqrt( (cosh(Ecm/2) - cosh(mpi))/2. ))
  ## qtsq = 2 pi q / L
  return(data.frame(gammaboost=E/Ecm, qtsq=(L*q/2./pi)^2, q=q, Ecm=Ecm))
}

## version with continuum dispersion relation
compute.qtildesq.contdisp <- function(E, dvec, mpi, L) {
  Ecmsq <- E^2 - sum((2*pi*dvec/L)^2)
  qsq <- Ecmsq/4.-mpi^2
  return(data.frame(gammaboost=E/sqrt(Ecmsq), qtsq=qsq*(L/2./pi)^2, q=sqrt(qsq), Ecm=sqrt(Ecmsq)))
}

