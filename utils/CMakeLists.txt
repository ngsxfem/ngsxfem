install( FILES
  p1interpol.hpp 
  xprolongation.hpp
  restrictedblf.hpp
  ngsxstd.hpp
  bitarraycf.hpp
  DESTINATION include
)


if(USE_PYTHON)
    install (FILES
            utils.py 
            DESTINATION xfem
            COMPONENT xfem
            )
endif(USE_PYTHON)




