install( FILES
  shiftintegrators.hpp calcpointshift.hpp calcgeomerrors.hpp lsetrefine.hpp
  projshift.hpp shiftedevaluate.hpp
  DESTINATION ${NGSOLVE_INSTALL_DIR_INCLUDE}
  )

install (FILES
        lsetcurv.py 
        DESTINATION ${INST_DIR}
        COMPONENT xfem
        )

add_subdirectory(py_demos)

