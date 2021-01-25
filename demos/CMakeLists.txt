install( FILES
  cutfem.py
  fictdom_dg_ghostpen.py
  fictdom_ghostpen.py
  fictdom_ghostpen_restrictedblf.py
  moving_domain.py
  nxfem.py
  nxfem_higher_order.py
  stokescutfem.py
  traceconverg.py
  tracefem.py
  tracefem3d.py
  xfespaces.py
  DESTINATION ${NGSOLVE_INSTALL_DIR_RES}/ngsxfem
)

install( FILES
  mpi/mpi_nxfem.py
  DESTINATION ${NGSOLVE_INSTALL_DIR_RES}/ngsxfem/mpi
)

install( FILES
  spacetime/test_RestrictToTime.py
  spacetime/spacetime_vtk.py
  spacetime/spaceP1_timeDGP1.py
  spacetime/spaceP1_timeCGP1.py
  spacetime/shiftedeval_debug.py
  spacetime/set_test.py
  spacetime/model_spacetime.py
  DESTINATION ${NGSOLVE_INSTALL_DIR_RES}/ngsxfem/spacetime
)

install( FILES
  mlset/poisson_on_square.py
  mlset/TensorFunctions.py
  DESTINATION ${NGSOLVE_INSTALL_DIR_RES}/ngsxfem/mlset
)
