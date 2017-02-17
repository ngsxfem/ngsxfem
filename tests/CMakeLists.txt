include(CTest)

add_test(NAME cpu_usage COMMAND top "-n1" "-b" )

add_test(NAME py_nxfem COMMAND ${PYTHON_EXECUTABLE} 
  "${PROJECT_SOURCE_DIR}/py_tutorials/nitschexfem.py")

add_test(NAME py_symbcutints1 COMMAND ${PYTHON_EXECUTABLE} 
  "${PROJECT_SOURCE_DIR}/py_tutorials/symboliccut.py")

add_test(NAME py_symbcutints2 COMMAND ${PYTHON_EXECUTABLE} 
  "${PROJECT_SOURCE_DIR}/py_tutorials/symboliccut2.py")

add_test(NAME py_symbcutints3 COMMAND ${PYTHON_EXECUTABLE} 
  "${PROJECT_SOURCE_DIR}/py_tutorials/symbolic_tracefem.py")

