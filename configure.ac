#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.64])
AC_INIT([BSPT], [1.5.2-beta], [yanpc@bcc.ac.cn])
AC_CONFIG_SRCDIR([src/main.c])
AC_CONFIG_HEADERS([config.h])
AM_INIT_AUTOMAKE([foreign -Wall])

# Checks for programs.
AC_PROG_CC

# Checks for libraries.
AC_ARG_ENABLE(bzlib, [  --disable-bzlib         use this option to disable bzlib support], [
	AC_DEFINE(DISABLE_BZLIB,1,[Define 1 if we want to disable bzlib support.])], [
		AC_CHECK_HEADER([bzlib.h], [
			LIBS="${LIBS} -lbz2";AC_DEFINE(HAVE_BZLIB, 1, [Define to 1 if you have the <bzlib.h> header file.])
		],  AC_MSG_ERROR("could not find the bzlib library.  Please use --disable-bzlib if you wish to disable bzlib support.")
	)
])
AC_CHECK_LIB([z], [gzread], [ LIBS="${LIBS} -lz"; ], AC_MSG_ERROR("could not find the zlib libaray."))

#AC_CONDITIONAL([NO_BZ2], [test "$use_bz2" = yes])

# Checks for header files.
AC_CHECK_HEADERS([stdlib.h string.h sys/time.h unistd.h zlib.h])

# Checks for typedefs, structures, and compiler characteristics.

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_CHECK_FUNCS([bzero memset mkdir pow])

AC_CONFIG_FILES([Makefile src/Makefile])
AC_OUTPUT
