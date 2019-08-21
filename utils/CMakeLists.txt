install( FILES
  p1interpol.hpp 
  xprolongation.hpp
  restrictedblf.hpp
  ngsxstd.hpp
  bitarraycf.hpp
  DESTINATION ${NGSOLVE_INSTALL_DIR_INCLUDE}
)


if(USE_PYTHON)
    install (FILES
            utils.py 
            DESTINATION ${INST_DIR}
            COMPONENT xfem
            )
endif(USE_PYTHON)




