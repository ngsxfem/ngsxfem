include(CTest)

if(APPLE)
    add_test(NAME cpu_usage COMMAND top "-l 1")
else(APPLE)
    add_test(NAME cpu_usage COMMAND top "-n1" "-b" )
endif(APPLE)

add_test(NAME py_demo_nxfem_lo COMMAND ${NETGEN_PYTHON_EXECUTABLE} 
  "${PROJECT_SOURCE_DIR}/demos/nxfem.py" "testmode")

add_test(NAME py_demo_cutfem_lo COMMAND ${NETGEN_PYTHON_EXECUTABLE} 
  "${PROJECT_SOURCE_DIR}/demos/cutfem.py" "testmode")

add_test(NAME py_demo_nxfem_ho COMMAND ${NETGEN_PYTHON_EXECUTABLE} 
  "${PROJECT_SOURCE_DIR}/demos/nxfem_higher_order.py" "testmode")

add_test(NAME py_demo_tracefem COMMAND ${NETGEN_PYTHON_EXECUTABLE} 
  "${PROJECT_SOURCE_DIR}/demos/tracefem.py" "testmode")

add_test(NAME py_demo_stokescutfem COMMAND ${NETGEN_PYTHON_EXECUTABLE} 
  "${PROJECT_SOURCE_DIR}/demos/stokescutfem.py" "testmode")

add_test(NAME py_demo_fictdom COMMAND ${NETGEN_PYTHON_EXECUTABLE} 
  "${PROJECT_SOURCE_DIR}/demos/fictdom_ghostpen.py" "testmode")

add_test(NAME py_demo_fictdom_dg COMMAND ${NETGEN_PYTHON_EXECUTABLE} 
  "${PROJECT_SOURCE_DIR}/demos/fictdom_dg_ghostpen.py" "testmode")

add_test(NAME pytests_xfes_ndof COMMAND ${NETGEN_PYTHON_EXECUTABLE} -m pytest
  "${PROJECT_SOURCE_DIR}/tests/pytests/test_xfes_ndof.py" WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}/tests")

add_test(NAME pytests_nonlin COMMAND ${NETGEN_PYTHON_EXECUTABLE} -m pytest
  "${PROJECT_SOURCE_DIR}/tests/pytests/test_nonlin.py" WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}/tests")

add_test(NAME pytests_calc_linearized COMMAND ${NETGEN_PYTHON_EXECUTABLE} -m pytest
  "${PROJECT_SOURCE_DIR}/tests/pytests/test_calc_linearized.py" WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}/tests")

add_test(NAME pytests_cutint COMMAND ${NETGEN_PYTHON_EXECUTABLE} -m pytest
  "${PROJECT_SOURCE_DIR}/tests/pytests/test_localgeom.py" WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}/tests")
  
add_test(NAME pytests_straightcutrule COMMAND ${NETGEN_PYTHON_EXECUTABLE} -m pytest
  "${PROJECT_SOURCE_DIR}/tests/pytests/test_straightcutrule.py" WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}/tests")
  
add_test(NAME pytests_backupcutrule COMMAND ${NETGEN_PYTHON_EXECUTABLE} -m pytest
  "${PROJECT_SOURCE_DIR}/tests/pytests/test_backuprule.py" WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}/tests")

add_test(NAME pytests_intcurved COMMAND ${NETGEN_PYTHON_EXECUTABLE} -m pytest
  "${PROJECT_SOURCE_DIR}/tests/pytests/test_intcurved.py" WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}/tests")

add_test(NAME pytests_nxfem COMMAND ${NETGEN_PYTHON_EXECUTABLE} -m pytest
  "${PROJECT_SOURCE_DIR}/tests/pytests/test_nxfem.py" WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}/tests")

add_test(NAME pytests_shifteval COMMAND ${NETGEN_PYTHON_EXECUTABLE} -m pytest
  "${PROJECT_SOURCE_DIR}/tests/pytests/test_shifteval.py" WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}/tests")

add_test(NAME pytests_quadchallenges COMMAND ${NETGEN_PYTHON_EXECUTABLE} -m pytest
  "${PROJECT_SOURCE_DIR}/tests/pytests/test_quadchallenges.py" WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}/tests")

add_test(NAME pytests_num_nxfem COMMAND ${NETGEN_PYTHON_EXECUTABLE} -m pytest
  "${PROJECT_SOURCE_DIR}/tests/pytests/test_num_nxfem.py" WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}/tests")
  
add_test(NAME pytests_spacetimecutrule COMMAND ${NETGEN_PYTHON_EXECUTABLE} -m pytest
  "${PROJECT_SOURCE_DIR}/tests/pytests/test_spacetimecutrule.py" WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}/tests")

add_test(NAME pytests_spacetime_lsetcurv COMMAND ${NETGEN_PYTHON_EXECUTABLE} -m pytest
  "${PROJECT_SOURCE_DIR}/tests/pytests/test_spacetime_lsetcurv.py" WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}/tests")
  
add_test(NAME pytests_spacetime_set COMMAND ${NETGEN_PYTHON_EXECUTABLE} -m pytest
  "${PROJECT_SOURCE_DIR}/tests/pytests/test_spacetime_set.py" WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}/tests")
    
add_test(NAME pytests_spacetime_vtk COMMAND ${NETGEN_PYTHON_EXECUTABLE} -m pytest
  "${PROJECT_SOURCE_DIR}/tests/pytests/test_spacetime_vtk.py" WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}/tests")

add_test(NAME pytests_navier_stokes COMMAND ${NETGEN_PYTHON_EXECUTABLE} -m pytest
  "${PROJECT_SOURCE_DIR}/tests/pytests/test_navier_stokes_2d1.py" WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}/tests")

add_test(NAME pytests_mlset_module COMMAND ${NETGEN_PYTHON_EXECUTABLE} -m pytest
  "${PROJECT_SOURCE_DIR}/tests/pytests/test_mlset_module_functionality.py" WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}/tests")

add_test(NAME pytests_mlset_integrate COMMAND ${NETGEN_PYTHON_EXECUTABLE} -m pytest
  "${PROJECT_SOURCE_DIR}/tests/pytests/test_mlset_integration.py" WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}/tests")

add_test(NAME pytests_mlset_pde COMMAND ${NETGEN_PYTHON_EXECUTABLE} -m pytest
  "${PROJECT_SOURCE_DIR}/tests/pytests/test_mlset_poisson.py" WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}/tests")

add_test(NAME pytests_memory COMMAND ${NETGEN_PYTHON_EXECUTABLE} -m pytest
  "${PROJECT_SOURCE_DIR}/tests/pytests/test_memory_usage.py" WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}/tests")

install( FILES
  ngsxfem_report.py
  DESTINATION ${NGSOLVE_INSTALL_DIR_RES}/ngsxfem/report
)
