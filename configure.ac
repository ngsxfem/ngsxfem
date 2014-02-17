AC_INIT([ngsxfem],[0.0.0.1-alpha],[])
AM_INIT_AUTOMAKE([-Wall -Werror foreign])

AC_CONFIG_MACRO_DIR([m4])

m4_pattern_allow([AM_PROG_AR])
AM_PROG_AR

AC_PROG_CC
AC_PROG_CXX

AC_PREFIX_DEFAULT(["/opt/netgen"])
AC_DISABLE_STATIC
AC_PROG_LIBTOOL

AC_SUBST([SOLVER], ["pardiso"])

AC_ARG_WITH([netgen], 
        [  --with-netgen=dir       use Netgen installed in directory dir],
        [netgendir=$withval],
        [netgendir="/opt/netgen"]
        )



# guion=true
# AC_ARG_ENABLE([gui],
#         [AC_HELP_STRING([--disable-gui],[don't build ngsolve with GUI])],
#         [if test "$enableval" = yes; then guion=true; else guion=false; fi])

# AM_CONDITIONAL([NGSGUI], [test x$guion = xtrue])

# # Tcl configuration:
# if test a$guion = atrue ; then
# TEA_INIT([3.9])
# TEA_PATH_TCLCONFIG
# TEA_LOAD_TCLCONFIG 
# TEA_PUBLIC_TCL_HEADERS
# fi


CXX=$netgendir/bin/ngscxx

AC_ARG_WITH([lapack],
        [  --with-lapack=lapacklib   linker flags for linking dense linear algebra library LAPACK],          
        [AC_SUBST([LAPACK_INCLUDES], ["-DLAPACK"])]
        [AC_SUBST([LAPACK_LIBS], [$withval])]
        )

# AC_ARG_ENABLE([navierstokes], 
#         [  --enable-navierstokes   install Navier Stokes package],
#         [nsemakefile="navierstokes/Makefile"]
# 	[AC_SUBST([NSEDIR], ["navierstokes"])]
#         [nse="yes"],
#         [nse="no"]
#         )

## [GM] commented this out since it isn't in use anymore
#AC_ARG_ENABLE([dev], 
#        [  --enable-dev   install Development packages (experimental stuff)],
#        [nsedevmakefile="navierstokes/dev/Makefile"]
#	[AC_SUBST([NSEDEVDIR], ["navierstokes/dev"])]
#        [nse="yes"]
#        )

#AC_ARG_ENABLE([pardiso],
#        [  --enable-pardiso           enable sparse direct solver pardiso],          
#        [AC_SUBST([PARDISO_INCLUDES], ["-DUSE_PARDISO"])]
#        )

AC_ARG_ENABLE([mumps],
        [  --enable-mumps           enable sparse direct solver mumps],          
        [AC_SUBST([SOLVER_FLAGS], ["-DUSE_MUMPS"])]
        [AC_SUBST([SOLVER], ["mumps"])]
        )

AC_ARG_WITH([superlu],
        [  --with-superlu=superlulib   linker flags for linking SUPERLU],          
        [AC_SUBST([SOLVER_FLAGS], ["-DUSE_SUPERLU $withval"])]
        [AC_SUBST([SOLVER_LIBS], ["-lsuperlu $withval"])]
        [AC_SUBST([SOLVER], ["superlu"])]
        )

AC_ARG_ENABLE([parallel],
        [  --enable-parallel           enable mpi parallelisation],          
        [AC_SUBST([PARALLEL_FLAGS], ["-DPARALLEL"])]
        [parallel="yes"],
        [parallel="no"]
        )


#AM_CONDITIONAL([COND_NSE], [test "$nse" = "yes"])
#AM_CONDITIONAL([COND_HDF], [test "$hdf" = "yes"])

AC_SUBST([NETGEN_INCLUDES], [" "])
AC_SUBST([NETGEN_LIBS], ["-L$netgendir/lib -lngsolve"])

AC_CONFIG_HEADERS(config.h) AC_CONFIG_FILES(Makefile utils/Makefile spacetime/Makefile spacetime/pde_demos/Makefile cutint/Makefile cutint/pde_demos/Makefile parabolic/Makefile parabolic/pde_demos/Makefile xfem/Makefile xfem/pde_demos/Makefile)

#programming_demos/Makefile
#navierstokes/dev/Makefile 

#$nsemakefile $hdfmakefile)

AC_OUTPUT

AC_MSG_RESULT([
------------------------------------------------------------------------
  $PACKAGE $VERSION:  Automatic configuration OK.

  Enabled functionality:

  Netgen/NGSolve: ...... $netgendir

  Building:

    Type 'make' to compile $PACKAGE.

    Type 'make install' to install $PACKAGE.
------------------------------------------------------------------------
])
