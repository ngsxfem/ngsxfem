file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/../ngs_solved "echo CTEST_FULL_OUTPUT\nngs \$1 &> test_out \nret_code=\$? \ncat test_out\necho \"ret_code=\${ret_code}\"\nFOUND=`grep \"Solved\" test_out` \nif [ \"\${FOUND}\" != \"\" ] \nthen \n cat ng.prof \n exit \${ret_code} \nelse \n exit 1 \nfi\n")
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

add_test("pde_d8_stokes3D" "ngs_solved" 
  "${PROJECT_SOURCE_DIR}/stokes/pde_demos/d8_stokes3D.pde")

# cmake_minimum_required(VERSION 2.8)

find_package(Boost COMPONENTS unit_test_framework REQUIRED)
add_executable(test_xlocalgeom test_xlocalgeom.cpp)
target_link_libraries(test_xlocalgeom ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY} ${NETGEN_LIBS} ngsxfem_spacetime ngsxfem_utils ngsxfem_cutint ngsxfem_xfem)

add_test("xlocalgeom" "test_xlocalgeom")

