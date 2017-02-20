include(CTest)

add_test(NAME cpu_usage COMMAND top "-n1" "-b" )

add_test(NAME py_tutorial_nxfem COMMAND ${PYTHON_EXECUTABLE} 
  "${PROJECT_SOURCE_DIR}/py_tutorials/nxfem.py")

add_test(NAME py_tutorial_nxfem_ho COMMAND ${PYTHON_EXECUTABLE} 
  "${PROJECT_SOURCE_DIR}/py_tutorials/nxfem_higher_order.py")

add_test(NAME py_tutorial_symbcutints1 COMMAND ${PYTHON_EXECUTABLE} 
  "${PROJECT_SOURCE_DIR}/py_tutorials/symboliccut.py")

add_test(NAME py_tutorial_symbcutints2 COMMAND ${PYTHON_EXECUTABLE} 
  "${PROJECT_SOURCE_DIR}/py_tutorials/symboliccut2.py")

add_test(NAME py_tutorial_symbcutints3 COMMAND ${PYTHON_EXECUTABLE} 
  "${PROJECT_SOURCE_DIR}/py_tutorials/symbolic_tracefem.py")

add_test(NAME py_tests_nxfem_old COMMAND ${PYTHON_EXECUTABLE} 
  "${PROJECT_SOURCE_DIR}/tests/nxfem_old_style.py")

