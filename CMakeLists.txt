cmake_minimum_required(VERSION 2.8)

message(STATUS "Use Ngsolve from submodule, updating modules...")
execute_process(COMMAND git submodule update --init --recursive WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
#set(USE_UMFPACK ON)

#add_custom_target(check_submodules_start_xfem ALL cmake -P cmake_modules/check_submodules.cmake WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} )

#include_directories("/usr/include/tcl8.5/tk-private/generic/")
#include_directories("/usr/include/tcl8.5/tk-private/unix/")

set(CXX_STANDARD 14)

#add_subdirectory(external_dependencies/ngsolve)
execute_process(COMMAND mkdir -p ngsolve)
execute_process(COMMAND cmake -DINSTALL_DIR=${INSTALL_DIR} -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE} -DCXX_STANDARD=${CXX_STANDARD} -DUSE_CCACHE=${USE_CCACHE} -DTCL_INCLUDE_PATH=/usr/include/tcl8.5/ ${CMAKE_CURRENT_SOURCE_DIR}/external_dependencies/ngsolve/ WORKING_DIRECTORY ngsolve)
#execute_process(COMMAND make WORKING_DIRECTORY ngsolve)
execute_process(COMMAND make install -j2 WORKING_DIRECTORY ngsolve)

set(BUILD_NETGEN ON)
set(NETGENDIR "${INSTALL_DIR}/bin")

set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG -std=c++14 -fPIC  -DHAVE_CONFIG_H -DUSE_TIMEOFDAY -DTCL -DUSE_MKL -DUSE_PARDISO -DLAPACK -DNGS_PYTHON -DUSE_UMFPACK")
set(CMAKE_CXX_FLAGS "-O3 -DNDEBUG -std=c++14 -fPIC  -DHAVE_CONFIG_H -DUSE_TIMEOFDAY -DTCL -DUSE_MKL -DUSE_PARDISO -DLAPACK -DNGS_PYTHON -DUSE_UMFPACK")

link_directories("${INSTALL_DIR}/lib/")
#link_directories("${CMAKE_CURRENT_SOURCE_DIR}/external_dependencies/ngsolve/solve")
#link_directories("${CMAKE_CURRENT_SOURCE_DIR}/external_dependencies/ngsolve/comp")
#link_directories("${CMAKE_CURRENT_SOURCE_DIR}/external_dependencies/ngsolve/fem")
#link_directories("${CMAKE_CURRENT_SOURCE_DIR}/external_dependencies/ngsolve/ngstd")
#link_directories("${CMAKE_CURRENT_SOURCE_DIR}/external_dependencies/ngsolve/basiclinalg")
#link_directories("${CMAKE_CURRENT_SOURCE_DIR}/external_dependencies/ngsolve/linalg")

include_directories("${INSTALL_DIR}/include/")
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/external_dependencies/ngsolve/external_dependencies/netgen/libsrc/visualization")
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/external_dependencies/ngsolve/external_dependencies/netgen/libsrc/general")
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/external_dependencies/ngsolve/external_dependencies/netgen/libsrc/include")
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/external_dependencies/ngsolve/include")
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/external_dependencies/ngsolve/ngstd")
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/external_dependencies/ngsolve/fem")
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/external_dependencies/ngsolve/external_dependencies/netgen/external_dependencies/pybind11/include/")

#######################################################################
# find netgen
#if(NOT NETGENDIR)
#    set(NETGENDIR $ENV{NETGENDIR} CACHE PATH "Path to netgen binary")
#else(NOT NETGENDIR)
#    set(NETGENDIR ${INSTALL_DIR}/bin CACHE PATH "Path to netgen binary")
#endif(NOT NETGENDIR)
get_filename_component(NETGEN_BASEDIR ${NETGENDIR} PATH)
if(NOT INSTALL_DIR)
  set(INSTALL_DIR ${NETGEN_BASEDIR} CACHE PATH "Path to install dir")
endif(NOT INSTALL_DIR)

set(CMAKE_INSTALL_PREFIX "${INSTALL_DIR}" CACHE INTERNAL "Prefix prepended to install directories" FORCE)

#SET(CMAKE_CXX_COMPILER "ccache g++")
if(USE_CCACHE)
  find_program(CCACHE_FOUND ccache)
  if(CCACHE_FOUND)
    set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE ccache)
    set_property(GLOBAL PROPERTY RULE_LAUNCH_LINK ccache)
  endif(CCACHE_FOUND)
endif(USE_CCACHE)

project(xfem)
set(NGSXFEM_VERSION_MAJOR 0)
set(NGSXFEM_VERSION_MINOR 6.0)
set(NGSXFEM_VERSION "${NGSXFEM_VERSION_MAJOR}.${NGSXFEM_VERSION_MINOR}-dev")
set(PACKAGE_VERSION ${NGSXFEM_VERSION})


add_definitions(-DNGSX_PYTHON)

# SET(CXX_COMPILER "${INSTALL_DIR}/bin/ngscxx" CACHE STRING "" FORCE)
# SET(CXX "${INSTALL_DIR}/bin/ngscxx" CACHE STRING "" FORCE)


set(NGS_LIB_TYPE SHARED)
set(NETGEN_LIBS "-L${NETGEN_BASEDIR}/lib -lngsolve -lngcomp -lngfem -lngstd -lngbla")
#set(NETGEN_LIBS "-L${CMAKE_CURRENT_SOURCE_DIR}/external_dependencies/ngsolve/solve -lngsolve")
#set(NETGEN_LIBS "${NETGEN_LIBS} -L${CMAKE_CURRENT_SOURCE_DIR}/external_dependencies/ngsolve/comp -lngcomp")
#set(NETGEN_LIBS "${NETGEN_LIBS} -L${CMAKE_CURRENT_SOURCE_DIR}/external_dependencies/ngsolve/fem -lngfem")
#set(NETGEN_LIBS "${NETGEN_LIBS} -L${CMAKE_CURRENT_SOURCE_DIR}/external_dependencies/ngsolve/ngstd -lngstd")
#set(NETGEN_LIBS "${NETGEN_LIBS} -L${CMAKE_CURRENT_SOURCE_DIR}/external_dependencies/ngsolve/basiclinalg -lngbla")

if (ADDITIONAL_PATHS)
    list(REVERSE CMAKE_SYSTEM_PREFIX_PATH)
    list(APPEND CMAKE_SYSTEM_PREFIX_PATH ${ADDITIONAL_PATHS})
    list(REVERSE CMAKE_SYSTEM_PREFIX_PATH)
endif (ADDITIONAL_PATHS)

# use rpath
SET(CMAKE_SKIP_BUILD_RPATH  FALSE)
SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE) 
SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
LIST(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
IF("${isSystemDir}" STREQUAL "-1")
   SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
ENDIF("${isSystemDir}" STREQUAL "-1")

set(ADDITIONAL_PATHS "" CACHE PATH "List of paths to additional libraries in non-standard locations, seperated by ';'")

SET(USE_PYTHON TRUE)

if (USE_PYTHON)
    set(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH}" "${PROJECT_SOURCE_DIR}/cmake_modules/python")
    set(PYTHON_VERSION "3" CACHE STRING "Python version (only Python >= 3.0 supported)")
    set(Python_ADDITIONAL_VERSIONS 3.5)
    if( PYTHON_VERSION VERSION_LESS 3 )
        message(FATAL_ERROR "NGSolve/ngsxfem supports only Python 3")
    endif( PYTHON_VERSION VERSION_LESS 3 )
    find_package(PythonInterp ${PYTHON_VERSION} REQUIRED)
    find_package(PythonLibs ${PYTHON_VERSION}  REQUIRED)
    add_definitions(-DNGS_PYTHON)

    include_directories(${PYTHON_INCLUDE_DIRS})
    set(PYTHON_LIBS "${PYTHON_LIBRARIES}")
    execute_process(COMMAND ${PYTHON_EXECUTABLE} -c "from distutils.sysconfig import get_python_lib; print(get_python_lib(1,0,''))" OUTPUT_VARIABLE PYTHON_PACKAGES_INSTALL_DIR OUTPUT_STRIP_TRAILING_WHITESPACE)
     message("python install path: ${PYTHON_PACKAGES_INSTALL_DIR}")
endif (USE_PYTHON)

#######################################################################
# CTest
enable_testing()
include(CTest)

set(ngs_install_dir RUNTIME DESTINATION bin LIBRARY DESTINATION lib ARCHIVE DESTINATION lib)

include_directories ("${PROJECT_SOURCE_DIR}/include")
include_directories ("${PROJECT_BINARY_DIR}")

set(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH}" "${PROJECT_SOURCE_DIR}/cmake_modules")

set(CMAKE_INCLUDE_CURRENT_DIR ON)

include (CheckIncludeFiles)
check_include_files (dlfcn.h HAVE_DLFCN_H)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/config.h.cmake ${CMAKE_CURRENT_BINARY_DIR}/config.h)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/config.h DESTINATION include)
# add_definitions(-DHAVE_CONFIG_H)

add_subdirectory(cutint)
#add_subdirectory(meancurvflow)
add_subdirectory(levelset)
add_subdirectory(python)
add_subdirectory(stokes)
add_subdirectory(utils)
add_subdirectory(xfem)
add_subdirectory(pde_tutorial)
add_subdirectory(tracefem)
add_subdirectory(lsetcurving)
add_subdirectory(osmosis)
add_subdirectory(tests)

# get_directory_property( DirDefs COMPILE_DEFINITIONS )

install (PROGRAMS DESTINATION bin)

set(CPACK_GENERATOR "DEB")
set(CPACK_PACKAGE_VERSION 0.6.1-dev)

set(CPACK_PACKAGE_NAME xfem)
set(CPACK_DEBIAN_PACKAGE_ARCHITECTURE amd64) # TODO should be set via: dpkg --print-architecture

set(CPACK_DEBIAN_PACKAGE_MAINTAINER Christoph Lehrenfeld <christoph.lehrenfeld@gmail.com>)
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY Short description of xfem)
set(CPACK_PACKAGE_DESCRIPTION Long description of xfem)
set(CPACK_DEBIAN_PACKAGE_DEPENDS "ngsolve (>= 6.0)")
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

##############################################################################
# Unhandled Makefile.am settings
##############################################################################

# echo "creating compiler wrapper ngscxx"
# file(WRITE ngscxx "${CMAKE_CXX_COMPILER}  -I. -I${netgendir}/../include ${CMAKE_CXX_FLAGS} -L${netgendir}/../lib -fpic  \$*  -Wl,--rpath=${NETGENDIR}/../lib\n") # $LAPACK_INCLUDES $TCL_INCLUDES $MUMPS_INCLUDES $PARDISO_INCLUDES, $LAPACK_LIBS 
# install(PROGRAMS "ngscxx" DESTINATION bin)
# echo "$CXX -I. $GOLD_INCLUDE $CXXFLAGS -L$netgendir/lib -fpic $LAPACK_INCLUDES $TCL_INCLUDES $MUMPS_INCLUDES $PARDISO_INCLUDES \$*  $LAPACK_LIBS -Wl,--rpath=$netgendir/lib" > ngscxx
# chmod a+x ngscxx

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

  Building:

    Type 'make' to compile ${PROJECT_NAME}.

    Type 'make install' to install ${PROJECT_NAME}.

  Install directory:
    ${CMAKE_INSTALL_PREFIX}

  Compiler:
    CMAKE_CXX_COMPILER = ${CMAKE_CXX_COMPILER}
    
  Compiler Flags:
    flags = ${flags}
  
  Build Type:
    type = ${CMAKE_BUILD_TYPE}

------------------------------------------------------------------------
")
#TODO:   Example programs will be built but not installed.


