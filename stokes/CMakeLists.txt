add_library(ngsxfem_xstokes ${NGS_LIB_TYPE}
  xstokesIntegrators.cpp
  xmeancurv.cpp
  xstokesspace.cpp
  )

set_target_properties(ngsxfem_xstokes PROPERTIES SUFFIX ".so")

if(NOT WIN32)
    target_link_libraries(ngsxfem_xstokes ngsxfem_spacetime ngsxfem_utils ngsxfem_cutint ngsxfem_xfem ${PYTHON_LIBS} ${LAPACK_LIBRARIES} ${NETGEN_LIBS})
    install( TARGETS ngsxfem_xstokes ${ngs_install_dir})
endif(NOT WIN32)

install( FILES
  xstokesIntegrators.hpp
  xmeancurv.hpp
  xstokesspace.hpp
  DESTINATION include
  )


