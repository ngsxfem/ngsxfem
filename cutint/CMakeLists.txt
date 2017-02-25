add_library(ngsxfem_cutint ${NGS_LIB_TYPE}
  fieldeval.cpp fieldeval.hpp
  xintegration.hpp xintegration.cpp
  xdecompose.cpp xdecompose.hpp  
  straightcutrule.hpp straightcutrule.cpp
        )

set_target_properties(ngsxfem_cutint PROPERTIES SUFFIX ".so")

if(NOT WIN32)
    target_link_libraries(ngsxfem_cutint ${PYTHON_LIBS} ${LAPACK_LIBRARIES} ${NETGEN_LIBS})
    install( TARGETS ngsxfem_cutint ${ngs_install_dir})
endif(NOT WIN32)

