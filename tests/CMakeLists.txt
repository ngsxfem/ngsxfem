file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/../ngs_solved "echo CTEST_FULL_OUTPUT && ngs \$1 | tee test_out \n FOUND=`grep \"Solved\" test_out` \n if [ \"\${FOUND}\" != \"\" ] \n then \n cat ng.prof \n exit 0 \n else \n exit 1 \n fi\n")
file(COPY ${CMAKE_CURRENT_BINARY_DIR}/../ngs_solved
  DESTINATION ${CMAKE_CURRENT_BINARY_DIR}
  FILE_PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ
  GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)
file(REMOVE ${CMAKE_CURRENT_BINARY_DIR}/../ngs_solved)

install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/ngs_solved DESTINATION ${CMAKE_CURRENT_BINARY_DIR})


include(CTest)

add_test("pde_d1_xnitsche" "ngs_solved" 
  "${PROJECT_SOURCE_DIR}/pde_tutorial/d1_approx.pde")

add_test("pde_d2_xnitsche" "ngs_solved" 
  "${PROJECT_SOURCE_DIR}/pde_tutorial/d2_xnitsche.pde")

add_test("pde_d7_stokes" "ngs_solved" 
  "${PROJECT_SOURCE_DIR}/stokes/pde_demos/d7_stokes.pde")
