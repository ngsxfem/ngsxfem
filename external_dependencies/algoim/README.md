# Algoim (third-party, bundled)

This directory contains a subset of **Algoim**, a C++ library by Robert I. Saye
(Lawrence Berkeley National Laboratory) providing *algorithms for implicitly
defined geometry, level set methods, and Voronoi implicit interface methods*.

* Project page: <https://algoim.github.io/>
* Repository: <https://github.com/algoim/algoim>
* Bundled version: commit `da1d81499608e1d499695d255f0233140b8c81e8`
  ("minor fixes", June 2023)
* License: modified BSD license of the Regents of the University of California,
  see [LICENSE](LICENSE) (copied verbatim from the Algoim repository)

ngsxfem uses Algoim's quadrature algorithms for domains implicitly defined by
multivariate polynomials in Bernstein form (`quadrature_multipoly.hpp`) as an
optional quadrature method on cut elements (see `cutint/algoimcutrule.hpp` and
`ngsxfemglobals.cut_quadrature`).

## Citation

If you use the Algoim-based quadrature of ngsxfem in your research, please cite
Algoim and the underlying publications, as requested on the Algoim project page:

* R. I. Saye, *High-order quadrature on multi-component domains implicitly
  defined by multivariate polynomials*, Journal of Computational Physics, 448,
  110720 (2022). doi: [10.1016/j.jcp.2021.110720](https://doi.org/10.1016/j.jcp.2021.110720)
* R. I. Saye, *High-order quadrature methods for implicitly defined surfaces and
  volumes in hyperrectangles*, SIAM Journal on Scientific Computing, 37(2),
  A993–A1019 (2015). doi: [10.1137/140966290](https://doi.org/10.1137/140966290)

## Contents

Only the headers required by `quadrature_multipoly.hpp` are bundled:

```
algoim/bernstein.hpp     algoim/polyset.hpp               algoim/tanhsinh.hpp
algoim/binomial.hpp      algoim/quadrature_multipoly.hpp  algoim/utility.hpp
algoim/booluarray.hpp    algoim/real.hpp                  algoim/uvector.hpp
algoim/gaussquad.hpp     algoim/sparkstack.hpp            algoim/xarray.hpp
algoim/multiloop.hpp
```

## Notes for ngsxfem developers

* Algoim is header-only, but some headers define non-inline functions. Hence,
  the headers must only be included in **one** translation unit of ngsxfem,
  namely `cutint/algoim_wrapper.cpp`. All other code uses the plain C++
  interface in `cutint/algoim_wrapper.hpp`.
* Algoim needs two LAPACKE routines (`LAPACKE_dgesvd`, `LAPACKE_dggevx`).
  ngsxfem does not link against LAPACK(E) directly. Instead,
  `cutint/algoim_lapacke/lapacke.h` provides self-contained replacements of
  exactly these two routines; the directory is put first on the include path of
  `cutint/algoim_wrapper.cpp`.
* To update Algoim, replace the headers listed above by the new versions, update
  the commit hash in this file and re-apply the modifications listed below.

## Modifications with respect to the upstream version

None.
