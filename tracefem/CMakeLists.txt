add_library(ngsxfem_tracefem ${NGS_LIB_TYPE}
  traceintegrators.cpp
  )

set_target_properties(ngsxfem_tracefem PROPERTIES SUFFIX ".so")

if(NOT WIN32)
  target_link_libraries(ngsxfem_tracefem ngsxfem_xfem ngsxfem_spacetime ngsxfem_utils ngsxfem_cutint ${PYTHON_LIBS} ${LAPACK_LIBRARIES} ${NETGEN_LIBS})
  install( TARGETS ngsxfem_tracefem ${ngs_install_dir})
endif(NOT WIN32)

install( FILES
  traceintegrators.hpp
  DESTINATION include
  )
