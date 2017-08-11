add_library(ngsxfem_lsetcurving ${NGS_LIB_TYPE}
  p1interpol.cpp shiftintegrators.cpp calcpointshift.cpp calcgeomerrors.cpp lsetrefine.cpp
  projshift.cpp shiftedevaluate.cpp
  )


install( FILES
  p1interpol.hpp shiftintegrators.hpp calcpointshift.hpp calcgeomerrors.hpp lsetrefine.hpp
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
