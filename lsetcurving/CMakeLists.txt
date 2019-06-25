# add_ngsolve_python_module(ngsxfem_lsetcurving_py python_lsetcurving.cpp)
# set_target_properties(ngsxfem_lsetcurving_py PROPERTIES INSTALL_RPATH "${NETGEN_RPATH_TOKEN}/../${NETGEN_PYTHON_RPATH}")
# install(TARGETS ngsxfem_lsetcurving_py DESTINATION ${NGSOLVE_INSTALL_DIR_PYTHON}/xfem)
# target_link_libraries(ngsxfem_lsetcurving_py ngsxfem_lsetcurving ngsxfem_utils)

add_library(ngsxfem_lsetcurving ${NGS_LIB_TYPE}
  shiftintegrators.cpp calcpointshift.cpp calcgeomerrors.cpp lsetrefine.cpp
  projshift.cpp shiftedevaluate.cpp
  )


install( FILES
  shiftintegrators.hpp calcpointshift.hpp calcgeomerrors.hpp lsetrefine.hpp
  projshift.hpp shiftedevaluate.hpp
  DESTINATION include
  )

target_link_libraries(ngsxfem_lsetcurving ngsxfem_xfem ${PYTHON_LIBS} ${LAPACK_LIBRARIES} ${NETGEN_LIBS})
target_link_libraries(ngsxfem_lsetcurving ngsxfem_xfem ${PYTHON_LIBS} ${LAPACK_LIBRARIES} ${NETGEN_LIBS})
install( TARGETS ngsxfem_lsetcurving ${ngs_install_dir})

# install ( FILES lsetcurv.py DESTINATION python)

install (FILES
        lsetcurv.py 
        DESTINATION ${NGSOLVE_INSTALL_DIR_PYTHON}/xfem
        COMPONENT xfem
        )

add_subdirectory(py_demos)

