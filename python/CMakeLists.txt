target_sources(xfem PRIVATE
      python_ngsxfem.cpp
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
      ../python/python_ngsxfem.cpp
      ../spacetime/diffopDt.cpp
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
install (FILES
        __init__.py cutmg.py mlset.py lset_spacetime.py lset_smoothblend.py lsetcurv.py utils.py ngs_check.py _version.py
        DESTINATION ${XFEM_INSTALL_DIR_PYTHON}/xfem
        COMPONENT xfem
)
