add_library(ngsxfem_meancurv ${NGS_LIB_TYPE}
          meancurv.cpp
        )

set_target_properties(ngsxfem_meancurv PROPERTIES SUFFIX ".so")

if(NOT WIN32)
    target_link_libraries(ngsxfem_meancurv ${PYTHON_LIBS} ${LAPACK_LIBRARIES} ${NETGEN_LIBS})
    install( TARGETS ngsxfem_meancurv ${ngs_install_dir})
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


