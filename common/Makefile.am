SUBDIRS = pde_demos

ACLOCAL_AMFLAGS = -I m4

AM_CPPFLAGS = -I$(netgendir)/include $(NETGEN_INCLUDES) $(LAPACK_INCLUDES) $(TCL_INCLUDES) -DLINUX $(SOLVER_FLAGS)

lib_LTLIBRARIES = libngsxfem_common.la

libngsxfem_common_la_SOURCES = spacetimefespace.cpp spacetimefespace.hpp spacetimefe.cpp spacetimefe.hpp spacetimeintegrators.cpp
#libngsxfem_common_la_LIBADD = timeint/libngsflow_timeint.la
libngsxfem_common_la_LDFLAGS = -avoid-version -rdynamic $(SOLVER_LIBS) $(NETGEN_LIBS)



