# Hadron

An R implementation of fitting routines used in lattice QCD. It provides useful
functions for extraction hadronic quantities and such like. 

The license is *GPL 3 or later*, even though the `DESCRIPTION` only shows
`GPL-3`.

[master branch](https://github.com/HISKP-LQCD/hadron): [![Build Status](https://travis-ci.org/HISKP-LQCD/hadron.svg?branch=master)](https://travis-ci.org/HISKP-LQCD/hadron)

# Installation

First clone `hadron` from github, e.g.

```{sh}
git clone https://github.com/HISKP-LQCD/hadron.git
```

Then change into the directory

```{sh}
cd hadron
```

`hadron` can be installed using the `install` script
provided. However, first the two packages `devtools` and `roxygen2`
need to be installed. Start `R` by typing

```{sh}
R
```

and then install the packages via

```{r}
install.packages(c("devtools", "roxygen2"), dependencies=TRUE)
```

following the instructions. `hadron` itself also depends on a few
libraries which can be installed as

```{r}
install.packages(c("Rcpp", "abind", "boot", "dplyr", "R6", "stringr"), dependencies=TRUE)
```

Thereafter, `hadron` can finally be installed from the linux command line 

```{sh}
./install
```

in the directory `hadron` was cloned into.

Alternatively, you may use the `install_github` function of the
`devtools` package to directly install from the github repository. 
