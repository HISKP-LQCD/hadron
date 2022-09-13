AC_INIT
AC_CONFIG_SRCDIR([src/alpha_s.c])

# This is from example from "writing R extensions"
# Now find the compiler and compiler flags to use
: ${R_HOME=`R RHOME`}
if test -z "${R_HOME}"; then
  echo "could not determine R_HOME"
  exit 1
fi
CC=`"${R_HOME}/bin/R" CMD config CC`
CXX=`"${R_HOME}/bin/R" CMD config CXX`
CFLAGS=`"${R_HOME}/bin/R" CMD config CFLAGS`
CPPFLAGS=`"${R_HOME}/bin/R" CMD config CPPFLAGS`

# Following lines kindly supplied by Dirk Eddelbuettel
## Use gsl-config to find arguments for compiler and linker flags
##
## Check for non-standard programs: gsl-config(1)
AC_PATH_PROG([GSL_CONFIG], [gsl-config])
## If gsl-config was found, let's use it
if test "${GSL_CONFIG}" != ""; then
   # Use gsl-config for header and linker arguments
   GSL_CFLAGS=`${GSL_CONFIG} --cflags`
   GSL_LIBS=`${GSL_CONFIG} --libs`
   CFLAGS="$CFLAGS $GSL_CFLAGS"  # kindly supplied by Ray Brownrigg
else
   AC_MSG_ERROR([gsl-config not found, is GSL installed?])
fi


#AC_CHECK_LIB(gsl, gsl_sf_airy_Ai_e,
#             [LIBS="${LIBS} -lgsl -lgslcblas"], 
#             [AC_MSG_ERROR([Cannot find Gnu Scientific Library.])],
#             [-lgslcblas])
#AC_CHECK_HEADER([gsl/gsl_version.h], ,
#             [AC_MSG_ERROR([Cannot find GSL headers.])])
AC_MSG_CHECKING([if GSL version >= 1.8])
AC_RUN_IFELSE([AC_LANG_SOURCE([[
#include <stdlib.h>
#include <string.h>
#include <gsl/gsl_version.h>
#include <stdio.h>
int main() {
#ifdef GSL_VERSION
   int major, minor;
   char *gslv = GSL_VERSION;
   if ((sscanf(gslv, "%d.%d", &major, &minor)) != 2) {
     exit (1);
   }
   exit( !( (major >= 2) || (major == 1 && minor >= 8 ) ) );
#else
  exit(1);
#endif
}
]])],
[gsl_version_ok=yes],
[gsl_version_ok=no],
[gsl_version_ok=yes])
if test "${gsl_version_ok}" = no; then
  AC_MSG_ERROR([Need GSL version >= 1.8])
else
  AC_MSG_RESULT([yes])
fi

AC_SUBST(GSL_CFLAGS) 
AC_SUBST(GSL_LIBS)

AC_CONFIG_FILES([src/Makevars])
AC_OUTPUT()

