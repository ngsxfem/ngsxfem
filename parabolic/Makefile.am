SUBDIRS = pde_demos

ACLOCAL_AMFLAGS = -I m4

AM_CPPFLAGS = -I$(netgendir)/include $(NETGEN_INCLUDES) $(LAPACK_INCLUDES) $(TCL_INCLUDES) -DLINUX $(SOLVER_FLAGS)

lib_LTLIBRARIES = libngsxfem_parabolic.la

libngsxfem_parabolic_la_SOURCES = solveinstat.cpp solveinstatx.cpp
libngsxfem_parabolic_la_LIBADD = ../spacetime/libngsxfem_spacetime.la ../xfem/libngsxfem_xfem.la
libngsxfem_parabolic_la_LDFLAGS = -avoid-version -rdynamic $(SOLVER_LIBS) $(NETGEN_LIBS)



