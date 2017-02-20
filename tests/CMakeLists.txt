include(CTest)

add_test(NAME cpu_usage COMMAND top "-n1" "-b" )

add_test(NAME py_tutorial_nxfem_lo COMMAND ${PYTHON_EXECUTABLE} 
  "${PROJECT_SOURCE_DIR}/py_tutorials/nxfem.py")

add_test(NAME py_tutorial_nxfem_ho COMMAND ${PYTHON_EXECUTABLE} 
  "${PROJECT_SOURCE_DIR}/py_tutorials/nxfem_higher_order.py")

add_test(NAME py_tutorial_symbcutints1 COMMAND ${PYTHON_EXECUTABLE} 
  "${PROJECT_SOURCE_DIR}/py_tutorials/symboliccut.py")

add_test(NAME py_tutorial_symbcutints2 COMMAND ${PYTHON_EXECUTABLE} 
  "${PROJECT_SOURCE_DIR}/py_tutorials/symboliccut2.py")

add_test(NAME py_tutorial_symbcutints3 COMMAND ${PYTHON_EXECUTABLE} 
  "${PROJECT_SOURCE_DIR}/py_tutorials/symbolic_tracefem.py")

add_test(NAME py_tutorial_nxfem_old COMMAND ${PYTHON_EXECUTABLE} 
  "${PROJECT_SOURCE_DIR}/tests/nxfem_old_style.py")

add_test(NAME pytests_xfes_ndof COMMAND ${PYTHON_EXECUTABLE} -m pytest
  "${PROJECT_SOURCE_DIR}/tests/pytests/test_xfes_ndof.py")

add_test(NAME pytests_cutint COMMAND ${PYTHON_EXECUTABLE} -m pytest
  "${PROJECT_SOURCE_DIR}/tests/pytests/test_localgeom.py")

install( FILES
  ngsxfem_report.py
  DESTINATION share/ngsxfem/report
)
