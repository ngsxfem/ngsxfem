add_library(ngsxfem_osmosis ${NGS_LIB_TYPE}
          meancurv.cpp
          fictext.cpp
	  python_osmosis.cpp
        )

set_target_properties(ngsxfem_osmosis PROPERTIES SUFFIX ".so")

if(NOT WIN32)
    target_link_libraries(ngsxfem_osmosis ngsxfem_spacetime ngsxfem_utils ngsxfem_cutint ngsxfem_xfem ${PYTHON_LIBS} ${LAPACK_LIBRARIES} ${NETGEN_LIBS})
    install( TARGETS ngsxfem_osmosis ${ngs_install_dir})
endif(NOT WIN32)

#file(COPY  "${CMAKE_CURRENT_SOURCE_DIR}/project" DESTINATION "${CMAKE_CURRENT_BINARY_DIR}")
