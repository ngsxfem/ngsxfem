#SUBDIRS = pde_demos

ACLOCAL_AMFLAGS = -I m4

AM_CPPFLAGS = -I$(netgendir)/include $(NETGEN_INCLUDES) $(LAPACK_INCLUDES) $(TCL_INCLUDES) -DLINUX $(SOLVER_FLAGS)

lib_LTLIBRARIES = libngsxfem_utils.la

libngsxfem_utils_la_SOURCES = fieldeval.cpp fieldeval.hpp error.cpp error.hpp output.cpp output.hpp
#libngsxfem_utils_la_LIBADD = ../spacetime/libngsxfem_spacetime.la
libngsxfem_utils_la_LDFLAGS = -avoid-version -rdynamic $(SOLVER_LIBS) $(NETGEN_LIBS)



