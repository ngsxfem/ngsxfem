install( FILES
  fictdom.py
  fictdom_dg.py
  fictdom_mlset.py
  lsetgeoms.py
  moving_domain.py
  stokescutfem.py
  tracefem.py
  unf_interf_prob.py
  DESTINATION ${XFEM_INSTALL_DIR_RES}/ngsxfem
)

install( FILES
  mpi/mpi_nxfem.py
  DESTINATION ${XFEM_INSTALL_DIR_RES}/ngsxfem/mpi
)

install( FILES
  spacetime/spacetime_vtk.py
  spacetime/spacetimeCG_unfitted.py
  spacetime/spacetimeDG_unfitted.py
  spacetime/spacetimeDG_fitted.py
  DESTINATION ${XFEM_INSTALL_DIR_RES}/ngsxfem/spacetime
)
