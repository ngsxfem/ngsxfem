AC_INIT([ngsxfem],[0.6.0-alpha],[])
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

pyon=false
AC_ARG_ENABLE([python],
	[AC_HELP_STRING([--enable-python],[enable python stuff])],
 	[CXXFLAGS+=" -DNGSX_PYTHON"]
    [pyon=true]
	)

AC_SUBST([NETGEN_INCLUDES], ["-I$netgendir/include"])
AC_SUBST([NETGEN_LIBS], ["-L$netgendir/lib -lngsolve"])

AC_CONFIG_HEADERS(config.h) AC_CONFIG_FILES(Makefile utils/Makefile spacetime/Makefile cutint/Makefile parabolic/Makefile xfem/Makefile meancurvflow/Makefile stokes/Makefile pde_tutorial/Makefile pde_experimental/Makefile python/Makefile)

AC_OUTPUT

AC_MSG_RESULT([
------------------------------------------------------------------------
  $PACKAGE $VERSION:  Automatic configuration OK.

  Enabled functionality:

  Netgen/NGSolve: ........ $netgendir
  Python functionality: .. $pyon

  Building:

    Type 'make' to compile $PACKAGE.

    Type 'make install' to install $PACKAGE.
------------------------------------------------------------------------
])
