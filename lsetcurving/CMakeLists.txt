add_library(ngsxfem_lsetcurving ${NGS_LIB_TYPE}
  p1interpol.cpp shiftintegrators.cpp calcpointshift.cpp setdeformation.cpp calcgeomerrors.cpp lsetrefine.cpp
  projshift.cpp
  )

set_target_properties(ngsxfem_lsetcurving PROPERTIES SUFFIX ".so")


install( FILES
  p1interpol.hpp shiftintegrators.hpp calcpointshift.hpp calcgeomerrors.hpp lsetrefine.hpp
  projshift.hpp
  DESTINATION include
  )

if(NOT WIN32)
    target_link_libraries(ngsxfem_lsetcurving ${PYTHON_LIBS} ${LAPACK_LIBRARIES} ${NETGEN_LIBS})
  install( TARGETS ngsxfem_lsetcurving ${ngs_install_dir})
endif(NOT WIN32)

install ( FILES lsetcurv.py DESTINATION python)
