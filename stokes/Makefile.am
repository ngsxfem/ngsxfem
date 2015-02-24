#SUBDIRS = pde_demos

ACLOCAL_AMFLAGS = -I m4

AM_CPPFLAGS = -I$(netgendir)/include $(NETGEN_INCLUDES) $(LAPACK_INCLUDES) $(TCL_INCLUDES) -DLINUX $(SOLVER_FLAGS)

lib_LTLIBRARIES = libngsxfem_xstokes.la

libngsxfem_xstokes_la_SOURCES =  xstokesIntegrators.hpp xstokesIntegrators.cpp xmeancurv.hpp xmeancurv.cpp xstokesspace.hpp xstokesspace.cpp
libngsxfem_xstokes_la_LIBADD = ../spacetime/libngsxfem_spacetime.la ../utils/libngsxfem_utils.la  ../cutint/libngsxfem_cutint.la ../xfem/libngsxfem_xfem.la
libngsxfem_xstokes_la_LDFLAGS = -avoid-version -rdynamic $(SOLVER_LIBS) $(NETGEN_LIBS)



