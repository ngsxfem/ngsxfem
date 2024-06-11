cmake_minimum_required(VERSION 3.16)

set(addon_name xfem)
project(${addon_name})

file(DOWNLOAD https://raw.githubusercontent.com/NGSolve/ngsolve-addon-template/317ecdf599286663cb64672874e43ed672b149a3/ngsolve_addon.cmake ${CMAKE_CURRENT_SOURCE_DIR}/ngsolve_addon.cmake) #EXPECTED_MD5 450caf3fa2b4bdee67624cfe448353a2)
include(ngsolve_addon.cmake)

set(CMAKE_INCLUDE_CURRENT_DIR ON)
include_directories ("${PROJECT_SOURCE_DIR}/include")
include_directories ("${PROJECT_BINARY_DIR}")

if(SKBUILD_PROJECT_VERSION)
  set(NGSXFEM_VERSION ${SKBUILD_PROJECT_VERSION})
  set(XFEM_INSTALL_DIR_INCLUDE ${SKBUILD_DATA_DIR}/include/xfem)
  set(XFEM_INSTALL_DIR_RES ${SKBUILD_DATA_DIR}/share/xfem)
else()
  set(NGSXFEM_VERSION "0.0.1")
  set(XFEM_INSTALL_DIR_INCLUDE ${NGSOLVE_INSTALL_DIR_INCLUDE})
  set(XFEM_INSTALL_DIR_RES ${NGSOLVE_INSTALL_DIR_RES})
endif()

add_definitions(-DNGSXFEM_VERSION="${NGSXFEM_VERSION}")
message(STATUS "NGSXFEM_VERSION: ${NGSXFEM_VERSION}")

add_ngsolve_addon(${addon_name})

add_subdirectory(python)
add_subdirectory(utils)
add_subdirectory(cutint)
add_subdirectory(xfem)
add_subdirectory(demos)
add_subdirectory(lsetcurving)
add_subdirectory(tests)
add_subdirectory(spacetime)

install(TARGETS ${addon_name} DESTINATION ${addon_name})

# generate stub files for autocomplete in IDEs
# this must be done at the very end (such that the stubgen generation happens after the python modules are installed)
install(CODE ${stubgen_generation_code})
install(DIRECTORY ${stubgen_directory} DESTINATION ${addon_name})

