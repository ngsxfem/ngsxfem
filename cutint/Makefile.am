SUBDIRS = pde_demos

ACLOCAL_AMFLAGS = -I m4

AM_CPPFLAGS = -I$(netgendir)/include $(NETGEN_INCLUDES) $(LAPACK_INCLUDES) $(TCL_INCLUDES) -DLINUX $(SOLVER_FLAGS)

lib_LTLIBRARIES = libngsxfem_test.la

libngsxfem_test_la_SOURCES = npxtest.cpp xintegration.hpp xintegration.cpp xdecompose.cpp xdecompose.hpp
libngsxfem_test_la_LIBADD = ../utils/libngsxfem_utils.la
libngsxfem_test_la_LDFLAGS = -avoid-version -rdynamic $(SOLVER_LIBS) $(NETGEN_LIBS)



