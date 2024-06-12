cmake_minimum_required(VERSION 3.16)

set(addon_name xfem)
project(${addon_name})

set(CMAKE_BUILD_TYPE RELWITHDEBINFO CACHE STRING "release type")

# Download ngsolve_addon.cmake from template project
set(ngsolve_addon_version 2264321bac8fb565eb4d38e5232369178646cd5f)
set(ngsolve_addon_md5 1deafa0b768789ce07de3e149aeaeb51)
file(DOWNLOAD https://raw.githubusercontent.com/NGSolve/ngsolve-addon-template/${ngsolve_addon_version}/ngsolve_addon.cmake
  ${CMAKE_CURRENT_SOURCE_DIR}/ngsolve_addon.cmake
  EXPECTED_MD5 ${ngsolve_addon_md5}
)
include(ngsolve_addon.cmake)

set(SPACETIME_SANITY_CHECKS ON CACHE BOOL "use space-time sanity checks (abuses ngsolve's intrule nr)")
if(SPACETIME_SANITY_CHECKS)
  add_definitions(-DSPACETIME_SANITY_CHECKS)
endif(SPACETIME_SANITY_CHECKS)

set(CMAKE_INCLUDE_CURRENT_DIR ON)
include_directories ("${PROJECT_SOURCE_DIR}/include")
include_directories ("${PROJECT_BINARY_DIR}")

if(SKBUILD_PROJECT_VERSION)
  set(NGSXFEM_VERSION ${SKBUILD_PROJECT_VERSION})
  set(XFEM_INSTALL_DIR_PYTHON .)
  set(XFEM_INSTALL_DIR_INCLUDE ${SKBUILD_DATA_DIR}/include/xfem)
  set(XFEM_INSTALL_DIR_RES ${SKBUILD_DATA_DIR}/share/xfem)
else()
  set(NGSXFEM_VERSION "2.1.2303.dev0")
  set(XFEM_INSTALL_DIR_PYTHON ${NETGEN_INSTALL_DIR_PYTHON})
  set(XFEM_INSTALL_DIR_INCLUDE ${NETGEN_INSTALL_DIR_INCLUDE})
  set(XFEM_INSTALL_DIR_RES ${NETGEN_INSTALL_DIR_RES})
  file(WRITE ${CMAKE_CURRENT_SOURCE_DIR}/python/_version.py "__version__ = \"${NGSXFEM_VERSION}\"\n")
endif()

add_definitions(-DNGSXFEM_VERSION="${NGSXFEM_VERSION}")
message(STATUS "NGSXFEM_VERSION: ${NGSXFEM_VERSION}")

add_ngsolve_addon(${addon_name})
#######################################################################
# CTest
enable_testing()
include(CTest)

add_subdirectory(python)
add_subdirectory(utils)
add_subdirectory(cutint)
add_subdirectory(xfem)
add_subdirectory(demos)
add_subdirectory(lsetcurving)
add_subdirectory(tests)
add_subdirectory(spacetime)

install(TARGETS ${addon_name} DESTINATION ${XFEM_INSTALL_DIR_PYTHON}/${addon_name})

# generate stub files for autocomplete in IDEs
# this must be done at the very end (such that the stubgen generation happens after the python modules are installed)
#install(CODE ${stubgen_generation_code})
#install(DIRECTORY ${stubgen_directory} DESTINATION ${XFEM_INSTALL_DIR_PYTHON}/${addon_name})

# uninstall target
configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake_modules/cmake_uninstall.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
    IMMEDIATE @ONLY)
add_custom_target(uninstall_xfem
    COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake)
