
AC_PREREQ([2.69])

AC_INIT([PPU-code],[2],[vlchaplin@gmail.com],[PPU-code])

AC_PREFIX_DEFAULT(/usr/local/PPU)

AC_CONFIG_SRCDIR([AnalyticalPPU.h])

AM_INIT_AUTOMAKE([foreign])

AC_ARG_VAR(HEADAS, [Location of include/fitsio.h and lib/libcfitsio.* (not needed if these are in the default header and & lib search paths, or if -I, -L flags are passed)])

AC_ARG_VAR(GSL, [Location of include/gsl/gsl*.h and lib/libgsl.* (not needed if [see above] )])

AC_ARG_VAR(HDF5, [Location of HDF5: i.e., the path containing 'include/H5Cpp.h' and 'lib/libhdf5_cpp.*' (not needed if [see above] )])
AC_ARG_VAR(HDF5_LDFLAGS, [Linker flags to succesfully compile HDF5-dependent code. If this variable isn't defined, a default is used: '-lz -lm -lhdf5 -lhdf5_cpp'. Consult your HDF5 library for the proper compression library (e.g., '-lsz' instead of '-lz'). ])


# Checks for programs.
AC_PROG_CXX
AC_PROG_CC

# Checks for library functions.
AC_HEADER_STDC
AC_TYPE_SIZE_T
AC_C_INLINE
AC_FUNC_MKTIME
AC_CHECK_HEADER_STDBOOL
AC_CHECK_HEADERS([stddef.h])
AC_CHECK_FUNCS([strpbrk])
AC_CHECK_FUNCS([floor])



if test x$HEADAS = x
then
AC_CHECK_HEADER( [fitsio.h], [have_fits=yes], [AC_MSG_WARN([fitsio.h not found in initial path])]  )
else
AC_CHECK_HEADER( [$HEADAS/include/fitsio.h], [have_fits=yes], [AC_MSG_WARN([fitsio.h not found in initial path])]  )
	if test x$have_fits = xyes
	then
		fits_inc=-I$HEADAS/include 
		fits_lib=-L$HEADAS/lib
	fi
fi

if test x$GSL = x
then
AC_CHECK_HEADER( [gsl/gsl_rng.h], [have_gsl=yes], [AC_MSG_WARN([gsl/gsl_rng.h not found in initial path])] )
else
AC_CHECK_HEADER( [$GSL/include/gsl/gsl_rng.h], [have_gsl=yes], [AC_MSG_WARN([gsl/gsl_rng.h not found in input path])] )
	if test x$have_gsl = xyes
	then
		gsl_inc=-I$GSL/include 
		gsl_lib=-L$GSL/lib
	fi
fi


if test x$HDF5 = x
then
AC_CHECK_HEADER( [H5version.h], [have_hdf5=yes], [AC_MSG_WARN([HDF5 headers not found in initial path])] )
else
AC_CHECK_HEADER( [$HDF5/include/H5version.h], [have_hdf5=yes], [AC_MSG_WARN([HDF5 headers not found in initial path])] )
	if test x$have_fits = xyes
	then
		h5_inc=-I$HDF5/include 
		h5_lib=-L$HDF5/lib
	fi
fi

if test "x$HDF5_LDFLAGS" = "x"
then
HDF5_LDFLAGS="-lz -lm -lhdf5 -lhdf5_cpp"
fi

AM_CONDITIONAL(ISDARWIN, test "x$is_darwin" = "xyes")

AC_SUBST([GSL])
AC_SUBST([HEADAS])
AC_SUBST([fits_inc])
AC_SUBST([fits_lib])
AC_SUBST([gsl_inc])
AC_SUBST([gsl_lib])
AC_SUBST([HDF5])
AC_SUBST([HDF5_LDFLAGS])
AC_SUBST([h5_inc])
AC_SUBST([h5_lib])

#AC_SUBST([LIBTOOL])



AC_PROG_MAKE_SET

AC_CONFIG_FILES([Makefile])
AC_OUTPUT