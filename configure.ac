#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([chenopis], [0.1], [fedor.baart@deltares.nl])
AC_CONFIG_SRCDIR([src/SwanCompUnstruc.F90])

AC_CONFIG_HEADER([config.h])
dnl Store generated files in config dir
AC_CONFIG_AUX_DIR([config])

dnl m4 macro's in m4 directory.
AC_CONFIG_MACRO_DIR([m4])
AM_INIT_AUTOMAKE

AC_PROG_FC

dnl Yes, we want free form Fortran!
AC_FC_FREEFORM
dnl Set unlimited line length (don't create lines longer than 132, because it breaks on some compilers,
dnl but because ifort and gfortran allow it let's not raise errors by default)
AC_FC_LINE_LENGTH([unlimited])

AC_PROG_F77
AC_FC_SRCEXT([F90])

# Checks for libraries.

# Checks for header files.

# Checks for typedefs, structures, and compiler characteristics.

# Checks for library functions.

# Check for features
AC_ARG_WITH(bmi,[AS_HELP_STRING([--with-bmi],[Build with bmi support @<:@default=no@:>@])], [],[with_bmi=yes])
AM_CONDITIONAL([HAVE_BMI], test "x$with_bmi" != "xno")

## ----------------------- ##
## Libtool initialisation. ##
## ----------------------- ##
LT_PREREQ(2.4)
LT_INIT



AC_CONFIG_FILES([
        Makefile
        src/Makefile
        tests/Makefile
        ])
AC_OUTPUT
