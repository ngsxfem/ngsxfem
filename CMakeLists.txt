cmake_minimum_required(VERSION 2.8)

#######################################################################
# find netgen
if(NOT NETGENDIR)
    set(NETGENDIR $ENV{NETGENDIR} CACHE PATH "Path to netgen binary")
else(NOT NETGENDIR)
    set(NETGENDIR ${INSTALL_DIR}/bin CACHE PATH "Path to netgen binary")
endif(NOT NETGENDIR)
get_filename_component(NETGEN_BASEDIR ${NETGENDIR} PATH)
if(NOT INSTALL_DIR)
  set(INSTALL_DIR ${NETGEN_BASEDIR} CACHE PATH "Path to install dir")
endif(NOT INSTALL_DIR)

set(CMAKE_INSTALL_PREFIX "${INSTALL_DIR}" CACHE INTERNAL "Prefix prepended to install directories" FORCE)

SET(CMAKE_CXX_COMPILER "${NETGENDIR}/ngscxx" CACHE STRING "" FORCE)

project(xfem)
set(NGSXFEM_VERSION_MAJOR 0)
set(NGSXFEM_VERSION_MINOR 6.0)
set(NGSXFEM_VERSION "${NGSXFEM_VERSION_MAJOR}.${NGSXFEM_VERSION_MINOR}-dev")
set(PACKAGE_VERSION ${NGSXFEM_VERSION})


add_definitions(-DNGSX_PYTHON)

# SET(CXX_COMPILER "${INSTALL_DIR}/bin/ngscxx" CACHE STRING "" FORCE)
# SET(CXX "${INSTALL_DIR}/bin/ngscxx" CACHE STRING "" FORCE)


set(NGS_LIB_TYPE SHARED)
set(NETGEN_LIBS "-L${NETGENDIR}/lib -lngsolve -lngcomp -lngfem -lngstd -lngbla")

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
add_subdirectory(meancurvflow)
add_subdirectory(levelset)
add_subdirectory(parabolic)
add_subdirectory(python)
add_subdirectory(spacetime)
add_subdirectory(stokes)
add_subdirectory(utils)
add_subdirectory(xfem)
add_subdirectory(pde_tutorial)
add_subdirectory(tracefem)
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
add_custom_target(uninstall
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


# SET(CMAKE_CXX_FLAGS "" CACHE STRING "" FORCE) 
# SET(CMAKE_CXX_FLAGS "")

function(message)
    _MESSAGE("${BoldBlue}${ARGV}${ColorReset}")
endfunction()

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

------------------------------------------------------------------------
")
#TODO:   Example programs will be built but not installed.


