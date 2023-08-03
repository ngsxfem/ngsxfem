# Demonstration files

The python files in this directory provide several examples that demonstrate the usage of `ngsxfem` features. For an introduction to these features we refer to the [ngsxfem-jupyter tutorials](https://github.com/ngsxfem/ngsxfem-jupyter). In the next table you can see the available demonstration files and a short description:

| filename | description | 
|-|-|
| [lsetgeoms.py](lsetgeoms.py) | Examples of smooth domains described by level set functions and their geometrical approximation. |
| [unf_interf_prob.py](unf_interf_prob.py) | Stationary scalar unfitted interface problem with an unfitted isoparametric FEM discretization and Nitsche. |
| [fictdom.py](fictdom.py) | Poisson problem on an unfitted (smooth) domain described by a level set function, discretized with CutFEM, Ghost Penalties and Nitsche. |
| [fictdom_dg.py](fictdom_dg.py) | Poisson problem on an unfitted (smooth) domain described by a level set function, discretized with Cut-DG (interior penalty), Ghost Penalties and Nitsche. |
| [fictdom_mlset.py](fictdom_mlset.py) | Poisson problem on an unfitted domain described by multiple level sets, discretized with CutFEM and Nitsche. |
| [stokescutfem.py](stokescutfem.py) | Stationary Stokes interface problem with an unfitted isoparametric Taylor-Hood-Nitsche discretization |
| [tracefem.py](tracefem.py) | stationary surface PDE problem with a TraceFEM discretization |
| [fictdom_dg.py](fictdom_dg.py) |   Poisson problem on an unfitted (smooth) domain described by a level set function, discretized with a Cut Discontinuous Galerkin method and Nitsche.|
| [moving_domain.py](moving_domain.py) | Parabolic PDE on a moving domain solved by an Eulerian time stepping method using isoparametric CutFEM and Ghost penalty stabilization. |
| [aggregates/aggfem_shapetester.py](aggregates/aggfem_shapetester.py) | Explanation of aggregated FEM basis functions as a means to deal with (small) cut elements. |
| [aggregates/fictdom_aggfem.py](aggregates/fictdom_aggfem.py) | FEM basis aggregation to solve a fictitious domain problem (without ghost penalties). |
| [aggregates/fictdom_dg_aggfem.py](aggregates/fictdom_dg_aggfem.py) | DG basis aggregation (cell merging) to solve a fictitious domain problem (without ghost penalties). |
| [spacetime/spacetimeDG_fitted.py](spacetime/spacetimeDG_fitted.py) | Parabolic PDE on a stationary fitted domain solved by a space-time DG-in-time method. |
| [spacetime/spacetimeDG_unfitted.py](spacetime/spacetimeDG_unfitted.py) | Parabolic PDE on a moving domain solved by a an unfitted isoparametric space-time DG-in-time CG-in-time method. |
| [spacetime/spacetimeCG_unfitted.py](spacetime/spacetimeCG_unfitted.py) | Parabolic PDE on a moving domain solved by a an unfitted isoparametric space-time continuous-in-time CG-in-time Petrov-Galerkin method. |
| [spacetime/spaceDGtimeDG_unfitted.py](spacetime/spacetimeDG_unfitted.py) | Parabolic PDE on a moving domain solved by a an unfitted isoparametric space-time DG-in-time DG-in-time method. |
| [spacetime/spaceDGtimeDG_unfitted_3D.py](spacetime/spacetimeDG_unfitted.py) | Parabolic PDE on a moving domain solved by a an unfitted isoparametric space-time DG-in-time DG-in-time method (3D). |
| [spacetime/spacetime_vtk.py](spacetime/spacetime_vtk.py) | Demonstration on how to export functions on a (2D+time) time slab to a vtk output. |
| [spdes/surfstokestracefem.py](spdes/surfstokestracefem.py) | TraceFEM discretization of the Stokes equations on a level set surface. |
| [mpi/mpi_nxfem.py](mpi/mpi_nxfem.py) | MPI version of a solution of a scalar unfitted interface problem in Nitsche-XFEM formulation. |
