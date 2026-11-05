# Release Notes

## 2.1.2608 (unreleased)

Requires NGSolve 6.2.2608.

### Breaking changes

- **Level sets that are not piecewise linear GridFunctions** (higher order GridFunctions,
  general CoefficientFunctions) are no longer silently interpolated piecewise linearly
  on cut elements: by default (`ngsxfemglobals.cut_quadrature = "auto"`) the
  Algoim-based quadrature w.r.t. the level set function itself is used (more accurate,
  different results). With `cut_quadrature = "standard"` such level sets raise an
  exception, unless `subdivlvl > 0` (subdivision approach) is used. For the previous
  behaviour, interpolate the level set function with `InterpolateToP1`.

- **Installed header location changed.** Headers are no longer installed flat
  into the NGSolve include directory but into the Python package, keeping the
  source layout: `<site-packages>/xfem/include/{cutint,xfem,spacetime,lsetcurving,utils}/`.
  Code that includes ngsxfem headers needs to be adapted, e.g.
  `#include <xfem/aggregates.hpp>`, `#include <cutint/lsetintdomain.hpp>`.
  The recommended way is to use the new CMake package (see below).
- `NoDeformation(mesh, ...)`: `mesh` is now a required argument.

### New features

- **Algoim-based cut quadrature.** The quadrature of the [Algoim](https://algoim.github.io/)
  library (R. I. Saye, bundled in `external_dependencies/algoim`) is used on cut
  elements for level set functions that are not piecewise linear GridFunctions. It works
  directly with higher order level set functions (no piecewise linear approximation, no
  mesh deformation), on segments, triangles, quadrilaterals, tetrahedra and hexahedra.
  Selection:
  ```python
  ngsxfemglobals.cut_quadrature = "auto"   # default; or "standard" / "algoim"
  dx = dCut(lset_ho, NEG, cut_quadrature="algoim")   # per integral
  ci = CutInfo(mesh, lset_ho, cut_quadrature="algoim")
  ```
  `"auto"` (default) uses the standard rules for piecewise linear level set GridFunctions
  (unchanged, also in the isoparametric approach) and Algoim for all other level sets; the
  use of Algoim is reported once per run in the terminal (with link and reference to cite),
  mentioning an automatic choice. Algoim is used for element,
  boundary (`vb=BND`, including interfaces on the boundary), facet (`skeleton=True`) and
  element boundary (`element_boundary=True`) integrals and for multiple level sets
  (`DomainTypeArray`, also on quadrilaterals and hexahedra, `MultiLevelsetCutInfo`).
  Algoim also supports prisms, which the standard rules do not support; on prisms (and for
  multiple level sets on quadrilaterals, hexahedra and prisms) `"auto"` uses Algoim also for
  piecewise linear level sets. See `demos/fictdom_algoim.py` and `CUT_QUADRATURE`.

- **Algoim-based space-time quadrature** (1D+t, 2D+t, 3D+t) for space-time level sets of
  higher order in space (`GridFunction` on `ScalarTimeFE(k_t) * H1(mesh, order=k_s)`) and
  for `CoefficientFunction`s of the reference time `tref`, e.g.
  ```python
  lset = GridFunction(ScalarTimeFE(2) * H1(mesh, order=2))
  lset.Set(levelset)   # levelset depends on t = told + delta_t * tref
  dQ = delta_t * dCut(lset, NEG, time_order=4)
  dOmnew = dCut(lset, NEG, tref=1)
  ```
  No space-time mesh deformation (`LevelSetMeshAdaptation_Spacetime`) is needed. Supported
  are element integrals, skeleton facet integrals (including interfaces on faces in 3D+t,
  which the standard space-time rules do not support), integrals at a fixed reference time
  (`tref`) and `CutInfo(..., time_order=...)`. The quadrature works in space-time (time
  preferably as the outermost integration variable); on tetrahedra, time slicing at the
  sign changes of the level set in the vertices is combined with the spatial rules at fixed
  times. Space-time level sets that are piecewise linear in space keep the standard rules
  with `cut_quadrature="auto"`. Multiple level sets (`DomainTypeArray`) are supported in
  space-time and at fixed reference times as well (which the standard rules do not support;
  3D+t by time slicing), also by `MultiLevelsetCutInfo(mesh, lsets, time_order=...)`. The Algoim-based space-time rules are considerably more
  expensive than the standard (isoparametric) ones (in the space-time demos: about 10-20x
  in 2D+t, 15x in 3D+t).

- `Integrate` on level set domains and `CutInfo` increase their local heap automatically
  on overflow.

- **Reusable C++ core library.** The C++ part is split into a shared library
  `xfem_core` and the Python extension module `xfem`, which links against it.
  Other projects (e.g. further NGSolve add-ons) can link against ngsxfem:
  ```cmake
  execute_process(COMMAND ${Python3_EXECUTABLE} -m xfem.config
                  OUTPUT_VARIABLE xfem_DIR OUTPUT_STRIP_TRAILING_WHITESPACE)
  find_package(xfem CONFIG REQUIRED)
  target_link_libraries(my_target PRIVATE xfem::core)
  ```
  `xfem.config` also provides `get_include_dir()` and `get_library_dir()`.
  See `INSTALLATION.md`.
- **Windows**: the C++ core exports its symbols as a DLL (`NGSXFEM_API`); Windows
  wheels are now built and tested like the Linux/macOS ones (also for CPython 3.10).

### Bug fixes

- Space-time level set GridFunctions of higher order in space were treated as piecewise
  linear by the standard space-time rules (`SpaceTimeFESpace` does not report the order of
  its spatial space), which gave wrong integrals without notice. They now use the
  Algoim-based rules (`cut_quadrature="auto"`) or raise an exception (`"standard"`).
  The same holds for level sets in `Discontinuous(H1(mesh, order=k))`, k >= 2.
- `XFESpace` with a higher order level set in 3D: local heap overflow in the update of its
  `CutInfo` (fixed heap size).

- Element boundary cut integrals (`element_boundary=True`) of NEG/POS parts on 3D
  elements were scaled wrongly; interfaces on boundary faces in 3D (`dCut(lset, IF,
  vb=BND)`) gave NaN.

- `LevelSetMeshAdaptation` with `boundary_tangential`: elements that touch
  several selected boundaries (corners) are excluded from the `qn` projection.
  This avoids an over-correction there.
- Space-time VTK output: the output file is now closed after writing.
- `RestrictedDiffOp`: the flux on elements without dofs was zeroed with the
  wrong matrix shape.
- Integration of SIMD coefficient functions (`Integrate` and the cut
  integrals): corrected matrix shape (worked before only because of the identical
  memory layout).
- Compatibility fixes for current NGSolve (`FV<double>`, `constexpr if`).

### Build, packaging and CI

- ngsolve is now pinned exactly (`ngsolve==6.2.2607`) for building and at runtime;
  the custom dependency provider `_ngsolve_deps_provider.py` was removed.
- `ngsolve_addon.cmake` updated to the current NGSolve add-on template
  (`add_ngsolve_addon_library`, relative RPATHs, exported CMake targets);
  the C++ standard remains C++20.
- PyPI workflow reworked: wheels are tested with `pytest` in CI, builds can be
  triggered manually (`ngsolve_version`, `dry_run`), and new NGSolve releases
  are picked up automatically via `ngsolve-release-watcher`.
- Test extras: `pip install xfem[test]` now also installs `scipy` and `psutil`.

### Documentation

- Literature page rebuilt from a single data file (`doc/literature.yaml`), with
  about 30 new entries, updated journal references and links to
  reproduction code.
- Compiled PDFs (`literature.pdf`, `feature-details.pdf`) are published next
  to the Sphinx pages and linked directly from the README.
