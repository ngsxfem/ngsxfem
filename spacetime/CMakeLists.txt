install( FILES
  SpaceTimeFE.hpp
  SpaceTimeFESpace.hpp
  diffopDt.hpp
  timecf.hpp
  spacetime_vtk.hpp
  DESTINATION ${NGSOLVE_INSTALL_DIR_INCLUDE}
)


install( FILES
  py_demos/spaceP1_timeDGP1.py
  DESTINATION ${NGSOLVE_INSTALL_DIR_RES}/ngsxfem/spacetime
)
