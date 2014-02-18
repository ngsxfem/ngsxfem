SUBDIRS = pde_demos

ACLOCAL_AMFLAGS = -I m4

AM_CPPFLAGS = -I$(netgendir)/include $(NETGEN_INCLUDES) $(LAPACK_INCLUDES) $(TCL_INCLUDES) -DLINUX $(SOLVER_FLAGS)

lib_LTLIBRARIES = libngsxfem_xfem.la

libngsxfem_xfem_la_SOURCES = xfiniteelement.cpp xfiniteelement.hpp xFESpace.cpp xFESpace.hpp xfemVisInts.cpp xfemVisInts.hpp xfemIntegrators.cpp xfemIntegrators.hpp xfemNitsche.cpp xfemNitsche.hpp
libngsxfem_xfem_la_LIBADD = ../spacetime/libngsxfem_spacetime.la ../utils/libngsxfem_utils.la  ../cutint/libngsxfem_cutint.la
libngsxfem_xfem_la_LDFLAGS = -avoid-version -rdynamic $(SOLVER_LIBS) $(NETGEN_LIBS)



