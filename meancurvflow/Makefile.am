#SUBDIRS = pde_demos

ACLOCAL_AMFLAGS = -I m4

AM_CPPFLAGS = -I$(netgendir)/include $(NETGEN_INCLUDES) $(LAPACK_INCLUDES) $(TCL_INCLUDES) -DLINUX $(SOLVER_FLAGS)

lib_LTLIBRARIES = libngsxfem_meancurv.la

libngsxfem_meancurv_la_SOURCES = meancurv.cpp
libngsxfem_meancurv_la_LDFLAGS = -avoid-version -rdynamic $(SOLVER_LIBS) $(NETGEN_LIBS)



