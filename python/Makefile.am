#SUBDIRS = pde_demos

ACLOCAL_AMFLAGS = -I m4

AM_CPPFLAGS = $(NETGEN_INCLUDES)

lib_LTLIBRARIES = libngsxfem_py.la

libngsxfem_py_la_SOURCES = python_xfem.cpp
libngsxfem_py_la_LDFLAGS = -avoid-version -rdynamic $(SOLVER_LIBS) $(NETGEN_LIBS)



