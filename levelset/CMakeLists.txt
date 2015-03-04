add_library(ngsxfem_levelset ${NGS_LIB_TYPE}
  corrsource.cpp
  python_levelset.cpp
  )

set_target_properties(ngsxfem_levelset PROPERTIES SUFFIX ".so")

if(NOT WIN32)
    target_link_libraries(ngsxfem_levelset ${PYTHON_LIBS} ${LAPACK_LIBRARIES} ${NETGEN_LIBS})
    install( TARGETS ngsxfem_levelset ${ngs_install_dir})
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


