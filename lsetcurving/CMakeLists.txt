install( FILES
  shiftintegrators.hpp calcpointshift.hpp calcgeomerrors.hpp lsetrefine.hpp
  projshift.hpp shiftedevaluate.hpp
  DESTINATION ${NGSOLVE_INSTALL_DIR_INCLUDE}
  )

add_subdirectory(py_demos)

