install( FILES
  shiftintegrators.hpp calcpointshift.hpp calcgeomerrors.hpp lsetrefine.hpp
  projshift.hpp shiftedevaluate.hpp
  DESTINATION include
  )

install (FILES
        lsetcurv.py 
        DESTINATION ${NGSOLVE_INSTALL_DIR_PYTHON}/xfem
        COMPONENT xfem
        )

add_subdirectory(py_demos)

