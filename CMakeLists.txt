cmake_minimum_required(VERSION 2.8)

if(NOT CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "${INSTALL_DIR}" CACHE INTERNAL "Prefix prepended to install directories" FORCE)
  if(NOT INSTALL_DIR)
    message( FATAL_ERROR "Please provide an CMAKE_INSTALL_PREFIX, CMake will exit." )
  endif(NOT INSTALL_DIR)
endif(NOT CMAKE_INSTALL_PREFIX)

set(BUILD_NGSOLVE ON CACHE BOOL "build NGSolve from scratch")
set(CMAKE_BUILD_TYPE RELWITHDEBINFO CACHE STRING "release type")
set(USE_CCACHE OFF CACHE BOOL "use ccache")
set(USE_GUI ON CACHE BOOL "use Netgen GUI")

if(CMAKE_GENERATOR STREQUAL "Ninja")
  set(ADDNINJA "-GNinja")
  set(BUILDEXE "ninja")
else(CMAKE_GENERATOR STREQUAL "Ninja")         
  set(ADDNINJA "")
  set(BUILDEXE "make")
endif(CMAKE_GENERATOR STREQUAL "Ninja")

if(BUILD_NGSOLVE)
  message(STATUS "Use NGSolve from submodule, updating modules...")
  execute_process(COMMAND git submodule update --init --recursive WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

  set(BUILD_NGSOLVE_THREADS "1" CACHE STRING "build NGSolve with n threads")
  set(USE_UMFPACK ON CACHE BOOL "use umfpack library (only if ngsolve is build)")
  set(USE_MKL ON CACHE BOOL "use mkl library (only if ngsolve is build)")
  set(MKL_ROOT "/opt/intel/mkl" CACHE STRING "mkl root dir (only if ngsolve is build)")
  set(TCL_INCLUDE_PATH "/usr/include" CACHE STRING "tcl include path")

  execute_process(COMMAND mkdir -p ngsolve)
  execute_process(COMMAND cmake
    -DINSTALL_DIR=${CMAKE_INSTALL_PREFIX}
    -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}
    -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
    -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
    -DCMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS}
    -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
    -DUSE_GUI=${USE_GUI}
    -DUSE_CCACHE=${USE_CCACHE}
    -DUSE_UMFPACK=${USE_UMFPACK}
    -DUSE_MKL=${USE_MKL}
    -DMKL_ROOT=${MKL_ROOT}
    -DTCL_INCLUDE_PATH=${TCL_INCLUDE_PATH} ${ADDNINJA}
    ${CMAKE_CURRENT_SOURCE_DIR}/external_dependencies/ngsolve/ WORKING_DIRECTORY ngsolve)

  # LINUX ONLY:
  if(CMAKE_GENERATOR STREQUAL "Ninja")
    execute_process(COMMAND ninja -j${BUILD_NGSOLVE_THREADS} WORKING_DIRECTORY ngsolve)
    execute_process(COMMAND ninja install WORKING_DIRECTORY ngsolve)
  else(CMAKE_GENERATOR STREQUAL "Ninja")         
    execute_process(COMMAND make -j${BUILD_NGSOLVE_THREADS} WORKING_DIRECTORY ngsolve)
    execute_process(COMMAND make install WORKING_DIRECTORY ngsolve)
  endif(CMAKE_GENERATOR STREQUAL "Ninja")


endif(BUILD_NGSOLVE)

if(USE_CCACHE)
  find_program(CCACHE_FOUND ccache)
  if(CCACHE_FOUND)
    set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE ccache)
    set_property(GLOBAL PROPERTY RULE_LAUNCH_LINK ccache)
  endif(CCACHE_FOUND)
endif(USE_CCACHE)


project(xfem)
find_package(NGSolve REQUIRED CONFIG)
set(NGSXFEM_VERSION_MAJOR 1)
set(NGSXFEM_VERSION_MINOR 0.0)
set(NGSXFEM_VERSION "${NGSXFEM_VERSION_MAJOR}.${NGSXFEM_VERSION_MINOR}-dev")
set(PACKAGE_VERSION ${NGSXFEM_VERSION})


add_definitions(-DNGSXFEM_VERSION="${NGSXFEM_VERSION}")

if(WIN32)
    set(NGS_LIB_TYPE STATIC)
else(WIN32)
    set(NGS_LIB_TYPE SHARED)
endif(WIN32)

set(NETGEN_LIBS solve)


if (ADDITIONAL_PATHS)
    list(REVERSE CMAKE_SYSTEM_PREFIX_PATH)
    list(APPEND CMAKE_SYSTEM_PREFIX_PATH ${ADDITIONAL_PATHS})
    list(REVERSE CMAKE_SYSTEM_PREFIX_PATH)
endif (ADDITIONAL_PATHS)

# use rpath
set(CMAKE_MACOSX_RPATH TRUE)
set(CMAKE_INSTALL_RPATH "${NETGEN_RPATH_TOKEN};${NETGEN_RPATH_TOKEN}/${NETGEN_RPATH}")

set(ADDITIONAL_PATHS "" CACHE PATH "List of paths to additional libraries in non-standard locations, seperated by ';'")

#######################################################################
# CTest
enable_testing()
include(CTest)

set(ngs_install_dir RUNTIME DESTINATION ${NGSOLVE_INSTALL_DIR_BIN} LIBRARY DESTINATION ${NGSOLVE_INSTALL_DIR_LIB} ARCHIVE DESTINATION ${NGSOLVE_INSTALL_DIR_LIB})


include_directories ("${PROJECT_SOURCE_DIR}/include")
include_directories ("${PROJECT_BINARY_DIR}")

set(CMAKE_INCLUDE_CURRENT_DIR ON)

add_subdirectory(cutint)
add_subdirectory(python)
add_subdirectory(utils)
add_subdirectory(xfem)
add_subdirectory(py_tutorials)
add_subdirectory(lsetcurving)
add_subdirectory(tests)
add_subdirectory(spacetime)

# get_directory_property( DirDefs COMPILE_DEFINITIONS )

set(CPACK_GENERATOR "DEB")
set(CPACK_PACKAGE_VERSION ${NGSXFEM_VERSION_MAJOR}.${NGSXFEM_VERSION_MINOR}-dev)

set(CPACK_PACKAGE_NAME xfem)
set(CPACK_DEBIAN_PACKAGE_ARCHITECTURE amd64) # TODO should be set via: dpkg --print-architecture

set(CPACK_DEBIAN_PACKAGE_MAINTAINER Christoph Lehrenfeld <christoph.lehrenfeld@gmail.com>)
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY Short description of xfem)
set(CPACK_PACKAGE_DESCRIPTION Long description of xfem)
set(CPACK_DEBIAN_PACKAGE_DEPENDS "ngsolve (>= 6.2)")
# TODO (use dpkg -s <packagename> to retrieve version)

set(CPACK_DEBIAN_PACKAGE_SECTION Science)
set(CPACK_DEBIAN_PACKAGE_VERSION ${CPACK_PACKAGE_VERSION})
include(CPack)

# uninstall target
configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake_modules/cmake_uninstall.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
    IMMEDIATE @ONLY)
add_custom_target(uninstall_xfem
    COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake)

if(NOT WIN32)
    string(ASCII 27 Esc)
    set(ColorReset "${Esc}[m")
    set(BoldBlue  "${Esc}[1;34m")
endif(NOT WIN32)

function(message)
    _MESSAGE("${BoldBlue}${ARGV}${ColorReset}")
endfunction()

string(TOUPPER "CMAKE_CXX_FLAGS_${CMAKE_BUILD_TYPE}" name)
set(flags "${${name}} ${CMAKE_CXX_FLAGS}")

message("
------------------------------------------------------------------------
  ${PROJECT_NAME} ${PACKAGE_VERSION}:  Automatic configuration OK.

  Install directory:
    ${CMAKE_INSTALL_PREFIX}

  Compiler:
    CMAKE_CXX_COMPILER = ${CMAKE_CXX_COMPILER}
    
  Compiler Flags:
    flags = ${flags}

  Generator:
    ${CMAKE_GENERATOR}

  Build Type:
    type = ${CMAKE_BUILD_TYPE}

  Build NGSolve from scratch: ${BUILD_NGSOLVE}
  Build xfem (and NGSolve) with ccache: ${USE_CCACHE}

  Building:

    Type '${BUILDEXE}' to compile ${PROJECT_NAME}.

    Type '${BUILDEXE} install' to install ${PROJECT_NAME}.

------------------------------------------------------------------------
")


