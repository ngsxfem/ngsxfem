add_ngsolve_python_module(ngsxfem_utils_py python_utils.cpp)
set_target_properties(ngsxfem_utils_py PROPERTIES INSTALL_RPATH "${NETGEN_RPATH_TOKEN}/../${NETGEN_PYTHON_RPATH}")
install(TARGETS ngsxfem_utils_py DESTINATION ${NGSOLVE_INSTALL_DIR_PYTHON}/xfem)
target_link_libraries(ngsxfem_utils_py ngsxfem_utils)

#ngsxfem_lsetcurving ngsxfem_utils ngsxfem_cutint ngsxfem_xfem)



add_library(ngsxfem_utils ${NGS_LIB_TYPE}
  restrictedblf.cpp
  bitarraycf.cpp
  ngsxstd.cpp
)

set_target_properties(ngsxfem_utils PROPERTIES SUFFIX ".so")

target_link_libraries(ngsxfem_utils ${PYTHON_LIBS} ${LAPACK_LIBRARIES} ${NETGEN_LIBS})
install( TARGETS ngsxfem_utils ${ngs_install_dir})


install( FILES
  restrictedblf.hpp
  ngsxstd.hpp
  bitarraycf.hpp
  DESTINATION include
)


if(USE_PYTHON)
    install (FILES
            utils.py 
            DESTINATION ${PYTHON_PACKAGES_INSTALL_DIR}/xfem
            COMPONENT xfem
            )
endif(USE_PYTHON)




