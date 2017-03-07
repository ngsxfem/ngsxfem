include(CTest)

add_test(NAME cpu_usage COMMAND top "-n1" "-b" )

add_test(NAME py_tutorial_nxfem_lo COMMAND ${PYTHON_EXECUTABLE} 
  "${PROJECT_SOURCE_DIR}/py_tutorials/nxfem.py")

add_test(NAME py_tutorial_cutfem_lo COMMAND ${PYTHON_EXECUTABLE} 
  "${PROJECT_SOURCE_DIR}/py_tutorials/cutfem.py")

add_test(NAME py_tutorial_nxfem_ho COMMAND ${PYTHON_EXECUTABLE} 
  "${PROJECT_SOURCE_DIR}/py_tutorials/nxfem_higher_order.py")

add_test(NAME py_tutorial_tracefem COMMAND ${PYTHON_EXECUTABLE} 
  "${PROJECT_SOURCE_DIR}/py_tutorials/tracefem.py")

add_test(NAME pytests_xfes_ndof COMMAND ${PYTHON_EXECUTABLE} -m pytest
  "${PROJECT_SOURCE_DIR}/tests/pytests/test_xfes_ndof.py" WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}/tests")

add_test(NAME pytests_cutint COMMAND ${PYTHON_EXECUTABLE} -m pytest
  "${PROJECT_SOURCE_DIR}/tests/pytests/test_localgeom.py" WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}/tests")
  
add_test(NAME pytests_straightcutrule COMMAND ${PYTHON_EXECUTABLE} -m pytest
  "${PROJECT_SOURCE_DIR}/tests/pytests/test_straightcutrule.py" WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}/tests")
  
add_test(NAME pytests_backupcutrule COMMAND ${PYTHON_EXECUTABLE} -m pytest
  "${PROJECT_SOURCE_DIR}/tests/pytests/test_backuprule.py" WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}/tests")

add_test(NAME pytests_nxfem COMMAND ${PYTHON_EXECUTABLE} -m pytest
  "${PROJECT_SOURCE_DIR}/tests/pytests/test_nxfem.py" WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}/tests")

install( FILES
  ngsxfem_report.py
  DESTINATION share/ngsxfem/report
)
