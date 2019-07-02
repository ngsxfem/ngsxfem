install( FILES
  shiftintegrators.hpp calcpointshift.hpp calcgeomerrors.hpp lsetrefine.hpp
  projshift.hpp shiftedevaluate.hpp
  DESTINATION include
  )

install (FILES
        lsetcurv.py 
        DESTINATION xfem
        COMPONENT xfem
        )

add_subdirectory(py_demos)

