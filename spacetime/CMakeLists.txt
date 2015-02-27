add_library(ngsxfem_spacetime ${NGS_LIB_TYPE}
  spacetimefespace.cpp spacetimefespace.hpp spacetimefe.cpp spacetimefe.hpp spacetimeintegrators.cpp spacetimeintegrators.hpp
  )

set_target_properties(ngsxfem_spacetime PROPERTIES SUFFIX ".so")

if(NOT WIN32)
    target_link_libraries(ngsxfem_spacetime ${PYTHON_LIBS} ${LAPACK_LIBRARIES} ${NETGEN_LIBS})
    install( TARGETS ngsxfem_spacetime ${ngs_install_dir})
endif(NOT WIN32)
