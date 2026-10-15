add_ngsolve_addon_library(xfem_core
      PACKAGE xfem
      EXPORT_NAME core
      SOURCES
      ../cutint/lsetintdomain.cpp
      ../cutint/fieldeval.cpp
      ../cutint/spacetimecutrule.cpp
      ../cutint/straightcutrule.cpp
      ../cutint/xdecompose.cpp
      ../cutint/mlsetintegration.cpp
      ../cutint/xintegration.cpp
      ../cutint/cutintegral.cpp
      ../lsetcurving/calcpointshift.cpp
      ../lsetcurving/lsetrefine.cpp
      ../lsetcurving/projshift.cpp
      ../lsetcurving/shiftedevaluate.cpp
      ../lsetcurving/shiftintegrators.cpp
      ../spacetime/diffopDt.cpp
      ../spacetime/diffopFixt.cpp
      ../spacetime/SpaceTimeFE.cpp
      ../spacetime/SpaceTimeFESpace.cpp
      ../spacetime/timecf.cpp
      ../spacetime/spacetime_vtk.cpp
      ../utils/bitarraycf.cpp
      ../utils/ngsxstd.cpp
      ../utils/p1interpol.cpp
      ../utils/restrictedblf.cpp
      ../utils/restrictedfespace.cpp
      ../utils/xprolongation.cpp
      ../xfem/cutinfo.cpp
      ../xfem/ghostpenalty.cpp
      ../xfem/sFESpace.cpp
      ../xfem/symboliccutbfi.cpp
      ../xfem/symboliccutlfi.cpp
      ../xfem/xfemdiffops.cpp
      ../xfem/xFESpace.cpp
      ../xfem/xfiniteelement.cpp
      ../xfem/aggregates.cpp
)
target_include_directories(xfem_core PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>)
target_compile_definitions(xfem_core PRIVATE NGSXFEM_EXPORTS
      PUBLIC $<$<BOOL:${SPACETIME_SANITY_CHECKS}>:SPACETIME_SANITY_CHECKS>)
set_target_properties(xfem_core PROPERTIES WINDOWS_EXPORT_ALL_SYMBOLS ON)

target_sources(xfem PRIVATE python_ngsxfem.cpp)
target_link_libraries(xfem PRIVATE xfem_core)
ngsolve_addon_set_relative_rpath(xfem ".")

install (FILES
        __init__.py config.py cutmg.py mlset.py lset_spacetime.py lset_smoothblend.py lsetcurv.py utils.py ngs_check.py _version.py
        DESTINATION ${XFEM_INSTALL_DIR_PYTHON}/xfem
        COMPONENT xfem
)
