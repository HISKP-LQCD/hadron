
AC_INIT([src/alpha_s.c])
AC_PROG_CXX
AC_CHECK_LIB(gsl, gsl_sf_airy_Ai_e,
             [LIBS="${LIBS} -lgsl -lgslcblas"], 
             [AC_MSG_ERROR([Cannot find Gnu Scientific Library.])],
             [-lgslcblas])
AC_CHECK_HEADER([gsl/gsl_version.h], ,
             [AC_MSG_ERROR([Cannot find GSL headers.])])
AC_MSG_CHECKING([if GSL version >= 1.8])
AC_RUN_IFELSE([AC_LANG_SOURCE([[
#include <stdlib.h>
#include <string.h>
#include <gsl/gsl_version.h>
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


AC_OUTPUT()

