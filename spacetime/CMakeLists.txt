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
  py_demos/spaceP1_timeDGP1.py
  DESTINATION share/ngsxfem
)
