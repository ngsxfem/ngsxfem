install( FILES
  cutfem.py
  fictdom_dg_ghostpen.py
  fictdom_ghostpen.py
  lsetgeoms.py
  mlset_poisson.py
  moving_domain.py
  nxfem.py
  stokescutfem.py
  tracefem.py
  DESTINATION ${NGSOLVE_INSTALL_DIR_RES}/ngsxfem
)

install( FILES
  mpi/mpi_nxfem.py
  DESTINATION ${NGSOLVE_INSTALL_DIR_RES}/ngsxfem/mpi
)

install( FILES
  spacetime/spacetime_vtk.py
  spacetime/spacetimeCG_unfitted.py
  spacetime/spacetimeDG_unfitted.py
  spacetime/spacetimeDG_fitted.py
  DESTINATION ${NGSOLVE_INSTALL_DIR_RES}/ngsxfem/spacetime
)
