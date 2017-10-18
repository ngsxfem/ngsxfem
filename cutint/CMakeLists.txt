add_ngsolve_python_module(ngsxfem_cutint_py python_cutint.cpp)
set_target_properties(ngsxfem_cutint_py PROPERTIES INSTALL_RPATH "${NETGEN_RPATH_TOKEN}/../${NETGEN_PYTHON_RPATH}")
install(TARGETS ngsxfem_cutint_py DESTINATION ${NGSOLVE_INSTALL_DIR_PYTHON}/xfem)
target_link_libraries(ngsxfem_cutint_py ngsxfem_cutint)

add_library(ngsxfem_cutint ${NGS_LIB_TYPE}
  fieldeval.cpp fieldeval.hpp
  xintegration.hpp xintegration.cpp
  xdecompose.cpp xdecompose.hpp  
  straightcutrule.hpp straightcutrule.cpp
  spacetimecutrule.hpp spacetimecutrule.cpp
        )

set_target_properties(ngsxfem_cutint PROPERTIES SUFFIX ".so")

target_link_libraries(ngsxfem_cutint ngsxfem_utils ngsxfem_spacetime ${PYTHON_LIBS} ${LAPACK_LIBRARIES} ${NETGEN_LIBS})
install( TARGETS ngsxfem_cutint ${ngs_install_dir})

