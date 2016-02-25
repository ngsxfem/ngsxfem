add_library(ngsxfem_utils ${NGS_LIB_TYPE}
  fieldeval.cpp fieldeval.hpp error.cpp error.hpp output.cpp output.hpp
  )

set_target_properties(ngsxfem_utils PROPERTIES SUFFIX ".so")

if(NOT WIN32)
    target_link_libraries(ngsxfem_utils ${PYTHON_LIBS} ${LAPACK_LIBRARIES} ${NETGEN_LIBS})
  install( TARGETS ngsxfem_utils ${ngs_install_dir})
endif(NOT WIN32)

# install( FILES
#         DESTINATION include
#        )

# install (FILES ngfemcf.py DESTINATION python)


##############################################################################
# Unhandled Makefile.am settings
##############################################################################

# libngfem_la_LIBADD = $(top_builddir)/ngstd/libngstd.la \
# 	$(top_builddir)/basiclinalg/libngbla.la $(LAPACK_LIBS)


