AC_PREREQ(2.59)

AC_INIT(sweepsims, [0.1.2], [https://github.com/molpopgen/sweepsims])
AC_CONFIG_SRCDIR([src/sweep.cc])
AM_INIT_AUTOMAKE
AC_CONFIG_HEADERS([config.h])

AC_CONFIG_MACRO_DIR([m4])

AC_PROG_CC
AC_C_CONST
AC_PROG_CXX
AX_CXX_COMPILE_STDCXX_11([noext],[mandatory])

AM_MAINTAINER_MODE([enable])

AC_LANG(C++)
AC_CONFIG_FILES([Makefile src/Makefile])

AC_ARG_ENABLE(debug,
[  --enable-debug    Turn on debugging],
[case "${enableval}" in
  yes) debug=true ;;
  no)  debug=false ;;
  *) AC_MSG_ERROR(bad value ${enableval} for --enable-debug) ;;
esac],[debug=false])
AM_CONDITIONAL(DEBUG, test x$debug = xtrue)

AC_ARG_ENABLE(profiling,
[  --enable-profiling    Turn on profiling],
[case "${enableval}" in
  yes) profiling=true ;;
  no)  profiling=false ;;
  *) AC_MSG_ERROR(bad value ${enableval} for --enable-profiling) ;;
esac],[profiling=false])
AM_CONDITIONAL(PROFILING, test x$profiling = xtrue)

dnl check for things that we need

dnl boost (only needed prior to 1.8.3)
dnl AC_CHECK_HEADER(boost/type_traits.hpp, LIBSEQUENCE_USE_BOOST=1 , [AC_MSG_ERROR([boost::type_traits required - either install boost (see http://www.boost.org) or write a replacement])])
dnl AC_CHECK_HEADER(boost/utility.hpp, LIBSEQUENCE_USE_BOOST=1 , [AC_MSG_ERROR([boost::noncopyable required (boost/utility.hpp) - either install boost (see http://www.boost.org) or write a replacement])])
dnl AC_CHECK_HEADER(boost/static_assert.hpp, LIBSEQUENCE_USE_BOOST=1 , [AC_MSG_ERROR([BOOST_STATIC_ASSERT required - either install boost (see http://www.boost.org) or write a replacement])])
dnl AC_CHECK_HEADER(boost/tuple/tuple.hpp,,[AC_MSG_ERROR([boot::tuple required - either install boost (see http://www.boost.org) or write a replacement])])

dnl gsl.  dependency removed in 1.8.3
dnl AC_CHECK_HEADER(gsl/gsl_rng.h, LIBSEQUENCE_USE_GSL=1 , [AC_MSG_ERROR([gsl/gsl_rng.h not found.  Please install GSL from http://www.gnu.org/software/gsl])])

AC_MSG_CHECKING(for namespace Sequence::coalsim)
AC_TRY_COMPILE(
[
#include <Sequence/Coalescent/Coalescent.hpp>

using namespace Sequence::coalsim;
],
[
Sequence::coalsim::arg x;
],
[AC_MSG_RESULT(yes)
AC_DEFINE([HAVE_SEQUENCE_COALSIM],[],"Have namespace Sequence::coalsim")],
[AC_MSG_RESULT(no)])

AC_CHECK_LIB(gslcblas,main,,[echo "GSL BLAS runtime library not found";exit 1])
AC_CHECK_LIB(gsl,main,,[echo "GSL runtime library not found";exit 1])

dnl zlib header
AC_CHECK_HEADER(zlib.h,,[AC_MSG_ERROR([zlib headers missing - cannot continue])])

dnl zlib runtime
AC_CHECK_LIB([z],gzungetc,,[echo "zlib run time library not found";exit 1])
AC_LANG_SAVE
          AC_LANG_CPLUSPLUS
          AC_REQUIRE_CPP
	  AC_CHECK_LIB(sequence, main,,[echo "libsequence runtime library does not seem to be here";exit 1])	
AC_LANG_RESTORE
dnl AM_CONFIG_HEADER(config.h)

AC_OUTPUT
