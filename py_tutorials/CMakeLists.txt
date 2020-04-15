install( FILES
  cutfem.py
  fictdom_dg_ghostpen.py
  fictdom_ghostpen.py
  fictdom_ghostpen_restrictedblf.py
  moving_domain.py
  nxfem.py
  nxfem_higher_order.py
  stokesxfem.py
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
  mlset/poisson_on_square.py
  mlset/TensorFunctions.py
  DESTINATION ${NGSOLVE_INSTALL_DIR_RES}/ngsxfem/mlset
)
