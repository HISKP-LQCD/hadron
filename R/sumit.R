# $Id$
# this does not help very much... :-(
# rather the contrary...
sumitinC <- function(y, par1, par2, v, err) {
#  .C("sumit",
#     as.double(y),
#     as.double(par1),
#     as.double(par2),
#     as.double(v),
#     as.double(err),
#     as.integer(length(y)),
#     res = double(1)
#     )$res
  .C("sumit",
     y,
     par1,
     par2,
     v,
     err,
     length(y),
     res = 0.,
     PACKAGE="hadron"
     )$res
}
