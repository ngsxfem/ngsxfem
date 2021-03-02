if(NETGEN_USE_PYTHON)
    install (FILES
        ngs_check.py
        DESTINATION ${NGSOLVE_INSTALL_DIR_PYTHON}/xfem
        COMPONENT xfem
    )
endif(NETGEN_USE_PYTHON)

install( FILES
  p1interpol.hpp 
  xprolongation.hpp
  restrictedblf.hpp
  restrictedfespace.hpp
  ngsxstd.hpp
  bitarraycf.hpp
  DESTINATION ${NGSOLVE_INSTALL_DIR_INCLUDE}
)
