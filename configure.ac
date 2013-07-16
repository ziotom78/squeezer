# Squeezer - compress LFI detector pointings and differenced data
# Copyright (C) 2013 Maurizio Tomasi (Planck collaboration)
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
# 02110-1301, USA.

AC_INIT([squeezer], [1.0], [tomasi@lambrate.inaf.it])

AM_INIT_AUTOMAKE([-Wall -Werror foreign])

######################################################################
# Check for a number of compiler's features

AC_PROG_CC
AC_PROG_CXX
AC_PROG_LIBTOOL

AX_CXX_COMPILE_STDCXX_11(,mandatory)

######################################################################
# Check for the presence of CFITSIO

AC_CHECK_LIB(cfitsio, ffopen,, AC_MSG_ERROR(Cannot find the CFITSIO library.), [-lz -lpthread])
LIBS="$LIBS -lz -lpthread"

######################################################################
# Check for the presence of GSL

AC_PATH_PROG(GSL_CONFIG, gsl-config, no)

if test "x$GSL_CONFIG" = "xno"; then
   AC_MSG_ERROR(Cannot configure GSL: gsl-config not found in path.)
else
   GSL_CFLAGS=`$GSL_CONFIG --cflags`
   GSL_LDFLAGS=`$GSL_CONFIG --libs`

   LDFLAGS="$GSL_LDFLAGS $LDFLAGS"
   AC_SUBST(GSL_CFLAGS)
   AC_SUBST(GSL_LDFLAGS)
fi

#AC_CHECK_LIB([m],[cos])
#AC_CHECK_LIB([gslcblas],[cblas_dgemm],,
#	[AC_MSG_ERROR(Cannot find the GSLCBLAS library.)],
#	[$GSL_LDFLAGS])
#AC_CHECK_LIB([gsl],[gsl_blas_dgemm],,
#        [AC_MSG_ERROR(Cannot find the GSL library)],
#	[$GSL_LDFLAGS])

######################################################################
# Check for the presence of CppUnit

AM_PATH_CPPUNIT(1.9.6, [cppunit=yes
		       AC_DEFINE(HAVE_CPPUNIT, 1, [Define to 1 if you have CppUnit available])],
		       [cppunit=no])

AM_CONDITIONAL([CPPUNIT_PRESENT], [test x$cppunit = xyes])

######################################################################
# Check for the presence of HPixLib

PKG_CHECK_MODULES([HPIXLIB], [libhpix],
	    [libhpix=yes AC_DEFINE(HAVE_HPIX, 1, [Define to 1 if you have HPixLib available])],
	    [libhpix=no])

AM_CONDITIONAL(HPIXLIB_PRESENT, [test x$libhpix == xyes])
AC_SUBST(HPIXLIB_CFLAGS)
AC_SUBST(HPIXLIB_LIBS)
AC_SUBST(HPIXLIB_LDFLAGS)

######################################################################
# Check for the presence of the TOODI library

AC_ARG_WITH([toodi-include-path],
	[AC_HELP_STRING([--with-toodi-include-path],
	                [location of the TOODI header files, defaults to $TOODI_JW])],
         [TOODI_CFLAGS="-I$withval"],
	 [TOODI_CFLAGS="-I$TOODI_JW"])

AC_ARG_WITH([toodi-lib-path],
	[AC_HELP_STRING([--with-toodi-lib-path],
	                [location of the TOODI library, defaults to $JW_LIB])],
         [TOODI_LDFLAGS="-L$withval -L$withval/client -L$withval/server -L$withval/native_threads -L$withval"],
	 [TOODI_LDFLAGS="-L$JW_LIB -L$JW_LIB/client -L$JW_LIB/server -L$JW_LIB/native_threads -L$JW_LIB"])

AC_MSG_CHECKING(for TOODI)
if test x$JW_LIB != x ; then :
   AC_MSG_RESULT(the JW_LIB variable is set, this looks promising)
   CPPFLAGS="$CPPFLAGS $TOODI_CFLAGS"
   LDFLAGS="$LDFLAGS $TOODI_LDFLAGS"
   LIBS="$LIBS -ljwrapper -ljvm"

   AC_CHECK_LIB(jwrapper, toodiInitializeLLIO,
	       [toodi=yes
	       AC_DEFINE(HAVE_TOODI, 1, [Define to 1 if you are compiling on Ironthrone])],
	       [toodi=no])
else
   AC_MSG_RESULT(no)
   toodi=no
fi

AM_CONDITIONAL([TOODI_PRESENT], [test x$toodi == xyes])

######################################################################
# Shut down and print a summary table

AC_CONFIG_HEADERS([src/config.hpp])
AC_CONFIG_FILES([
	Makefile
	src/Makefile
])
AC_OUTPUT

echo ""
echo "   Squeezer will be compiled with the following options:"
echo ""
echo "   HPixLib: $libhpix"
echo "   CppUnit: $cppunit"
echo "   TOODI: $toodi"
echo ""
