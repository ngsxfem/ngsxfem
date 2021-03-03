if(NETGEN_USE_PYTHON)
    add_ngsolve_python_module(ngsxfem_py 
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
    )
    set_target_properties(ngsxfem_py PROPERTIES INSTALL_RPATH "${NETGEN_RPATH_TOKEN}/../${NETGEN_PYTHON_RPATH}")
    install(TARGETS ngsxfem_py DESTINATION ${NGSOLVE_INSTALL_DIR_PYTHON}/xfem COMPONENT xfem)

    target_include_directories(ngsxfem_py PUBLIC ${NGSOLVE_INCLUDE_DIRS})
    # target_link_libraries(ngsxfem_py ngsxfem_lsetcurving ngsxfem_utils ngsxfem_cutint ngsxfem_xfem ngsxfem_spacetime)

    install (FILES
        __init__.py cutmg.py mlset.py lset_spacetime.py lsetcurv.py utils.py ngs_check.py
        DESTINATION ${NGSOLVE_INSTALL_DIR_PYTHON}/xfem
        COMPONENT xfem
    )

# build stub files for pybind11 packages
if(BUILD_STUB_FILES)
  execute_process(COMMAND ${NETGEN_PYTHON_EXECUTABLE} -c "import pybind11_stubgen; print(pybind11_stubgen.__file__)" OUTPUT_VARIABLE stubgen_path RESULT_VARIABLE pybind11_stubgen)
if(pybind11_stubgen AND NOT ${pybind11_stubgen} EQUAL 0)
  message(WARNING "pybind11-stubgen not found, if you want to create stub files
for better autocompletion support install it with pip.")
else()
  message("-- Found pybind11-stubgen: ${stubgen_path}")
  install(CODE "execute_process(COMMAND ${NETGEN_PYTHON_EXECUTABLE} -m pybind11_stubgen --no-setup-py xfem)") 
  install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../stubs/xfem-stubs/ DESTINATION ${NGSOLVE_INSTALL_DIR_PYTHON}/xfem/ COMPONENT xfem)
endif()
endif(BUILD_STUB_FILES)

endif(NETGEN_USE_PYTHON)
