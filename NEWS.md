# Version 3.3.1

* added a double bootstrap routine for 'cf' called
  double_bootstrap.cf
* added a new interface for 'gevp.hankel' which evaluates as a function
  of 'n', the matrix size. The new interface is called
  'bootstrap.pgevm', which also acts on the double bootstrap data, if available
* 'pgevm2effectivemass' can also deal with the double bootstrap data
  to estimate the uncertainty for the media.
* added the oblique Lanczos method 'bootstrap.lanczos' for the
  analysis of Euclidean correlation functions
* 'gevp.hankel' received some more functionality
* 'plot.effectivemass' has a new optional parameter 'xshift' to shift
  data points in x-direction for legibility.
  
# Version 3.2.1

* fix problems in 'configure.ac'
  Now it should also work on Mac
  
# Version 3.2.0

* fix problems with roxygen docu in 'correlators_key_meson_3pt' and
  'cf_key_meson_3pt'
* faster implementation of 'plotwitherror'
* allow to plot a bootstrap nlsfit on top of existing plot
* Add some tests for 'computeDisc' function

