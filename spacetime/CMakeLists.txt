add_ngsolve_python_module(ngsxfem_spacetime_py python_spacetime.cpp)
set_target_properties(ngsxfem_spacetime_py PROPERTIES INSTALL_RPATH "${NETGEN_RPATH_TOKEN}/../${NETGEN_PYTHON_RPATH}")
install(TARGETS ngsxfem_spacetime_py DESTINATION ${NGSOLVE_INSTALL_DIR_PYTHON}/xfem)
target_link_libraries(ngsxfem_spacetime_py ngsxfem_spacetime ngsxfem_cutint ngsxfem_lsetcurving)

add_library(ngsxfem_spacetime ${NGS_LIB_TYPE}
  SpaceTimeFE.cpp
  SpaceTimeFESpace.cpp
  diffopDt.cpp
  timecf.cpp
)

set_target_properties(ngsxfem_spacetime PROPERTIES SUFFIX ".so")

if(NOT WIN32)
    target_link_libraries(ngsxfem_spacetime ${PYTHON_LIBS} ${LAPACK_LIBRARIES} ${NETGEN_LIBS})
    install( TARGETS ngsxfem_spacetime ${ngs_install_dir})
endif(NOT WIN32)


install( FILES
  SpaceTimeFE.hpp
  SpaceTimeFESpace.hpp
  diffopDt.hpp
  timecf.hpp
  DESTINATION include
)

install (FILES
        lset_spacetime.py
        DESTINATION ${NGSOLVE_INSTALL_DIR_PYTHON}/xfem
        COMPONENT xfem
        )


install( FILES
  py_demos/spacetime_P1P1.py
  DESTINATION share/ngsxfem
)
