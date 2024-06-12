cmake_minimum_required(VERSION 3.16)

project(xfem)

# Set default build type to Release if not specified
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE RELWITHDEBINFO CACHE STRING "Choose the type of build." FORCE)
endif()

# Download ngsolve_addon.cmake from template project
set(ngsolve_addon_version bbfb9648f2ff7a952a3a8edca58f460c198ca462)
set(ngsolve_addon_md5 d1b208d6aa63edce64a264e5e0370d3a)
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
else()
  set(NGSXFEM_VERSION "2.1.2303.dev0")
  file(WRITE ${CMAKE_CURRENT_SOURCE_DIR}/python/_version.py "__version__ = \"${NGSXFEM_VERSION}\"\n")
endif()

set(XFEM_INSTALL_DIR_PYTHON ${ADDON_INSTALL_DIR_PYTHON})
set(XFEM_INSTALL_DIR_INCLUDE ${ADDON_INSTALL_DIR_INCLUDE})
set(XFEM_INSTALL_DIR_RES ${ADDON_INSTALL_DIR_RES})

add_definitions(-DNGSXFEM_VERSION="${NGSXFEM_VERSION}")
message(STATUS "NGSXFEM_VERSION: ${NGSXFEM_VERSION}")

add_ngsolve_addon(xfem)
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

install(TARGETS xfem DESTINATION ${XFEM_INSTALL_DIR_PYTHON}/xfem)

# ngsolve_generate_stub_files(xfem)

# uninstall target
configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake_modules/cmake_uninstall.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
    IMMEDIATE @ONLY)
add_custom_target(uninstall_xfem
    COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake)
