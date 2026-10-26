# Release Notes

## 2.1.2608 (unreleased)

Requires NGSolve 6.2.2608.

### Breaking changes

- **Installed header location changed.** Headers are no longer installed flat
  into the NGSolve include directory but into the Python package, keeping the
  source layout: `<site-packages>/xfem/include/{cutint,xfem,spacetime,lsetcurving,utils}/`.
  Code that includes ngsxfem headers needs to be adapted, e.g.
  `#include <xfem/aggregates.hpp>`, `#include <cutint/lsetintdomain.hpp>`.
  The recommended way is to use the new CMake package (see below).
- `NoDeformation(mesh, ...)`: `mesh` is now a required argument.

### New features

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
