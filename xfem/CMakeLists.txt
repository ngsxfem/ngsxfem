add_library(ngsxfem_xfem ${NGS_LIB_TYPE}
  ghostpenalty.cpp
  numprocx.cpp
  sdxfemIntegrators.cpp
  stxfemIntegrators.cpp
  xFESpace.cpp
  xfemIntegrators.cpp
  xfemNitsche.cpp
  xfemVisInts.cpp
  xfiniteelement.cpp
  python_xfem.cpp
  )

set_target_properties(ngsxfem_xfem PROPERTIES SUFFIX ".so")

if(NOT WIN32)
  target_link_libraries(ngsxfem_xfem ngsxfem_spacetime ngsxfem_utils ngsxfem_cutint ${PYTHON_LIBS} ${LAPACK_LIBRARIES} ${NETGEN_LIBS})
  install( TARGETS ngsxfem_xfem ${ngs_install_dir})
endif(NOT WIN32)

install( FILES
  ghostpenalty.hpp
  setvaluesx.hpp
  stxfemIntegrators.hpp
  xFESpace.hpp
  xfemIntegrators.hpp
  xfemNitsche.hpp
  xfemVisInts.hpp
  xfiniteelement.hpp
  DESTINATION include
  )

# install (FILES ngfemcf.py DESTINATION python)


##############################################################################
# Unhandled Makefile.am settings
##############################################################################

# libngfem_la_LIBADD = $(top_builddir)/ngstd/libngstd.la \
# 	$(top_builddir)/basiclinalg/libngbla.la $(LAPACK_LIBS)


