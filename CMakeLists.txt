cmake_minimum_required(VERSION 2.8)

project(xfem)

set(BUILD_NGSOLVE ON CACHE BOOL "build NGSolve from scratch")
set(CHECK_NGSOLVE_VERSION ON CACHE BOOL "Check if installed NGSolve (if exist) is compatible")
set(CMAKE_BUILD_TYPE RELWITHDEBINFO CACHE STRING "release type")
set(USE_CCACHE OFF CACHE BOOL "use ccache")
set(USE_GUI ON CACHE BOOL "use Netgen GUI")
set(BUILD_STUB_FILES OFF CACHE BOOL "Build stub files for better autocompletion")
set(SPACETIME_SANITY_CHECKS ON CACHE BOOL "use space-time sanity checks (abuses ngsolve's intrule nr)")

set(INST_DIR xfem CACHE PATH "Install dir appended to prefix")

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

  if(NOT CMAKE_INSTALL_PREFIX)
    set(CMAKE_INSTALL_PREFIX "${INSTALL_DIR}" CACHE INTERNAL "Prefix prepended to install directories" FORCE)
    if(NOT INSTALL_DIR)
      message( FATAL_ERROR "Please provide an CMAKE_INSTALL_PREFIX, CMake will exit." )
    endif(NOT INSTALL_DIR)
  endif(NOT CMAKE_INSTALL_PREFIX)

  set(BUILD_NGSOLVE_THREADS "1" CACHE STRING "build NGSolve with n threads")
  set(USE_UMFPACK ON CACHE BOOL "use umfpack library (only if ngsolve is build)")
  set(USE_MKL OFF CACHE BOOL "use mkl library (only if ngsolve is build)")
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

else(BUILD_NGSOLVE)
  find_package(PythonInterp 3 REQUIRED)
  find_package(PythonLibs 3 REQUIRED)
  # check if CMAKE_INSTALL_PREFIX is set by user, otherwise install to user-local python modules directory
  if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    execute_process(COMMAND ${PYTHON_EXECUTABLE} -c "import site; print(site.USER_SITE)" OUTPUT_VARIABLE install_dir OUTPUT_STRIP_TRAILING_WHITESPACE)
    set(CMAKE_INSTALL_PREFIX ${install_dir} CACHE PATH "Install dir" FORCE)
    message(STATUS "install dir ${install_dir}")
  endif(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  find_package(NGSolve CONFIG REQUIRED
    HINTS /usr $ENV{NETGENDIR}/.. /opt/netgen/ /Applications/Netgen.app/Contents/Resources/CMake C:/netgen $ENV{NETGENDIR}/../lib/cmake/ngsolve/ $ENV{NETGENDIR}/../Resources/CMake
  )


  if (CHECK_NGSOLVE_VERSION)
    execute_process(COMMAND ${NETGEN_PYTHON_EXECUTABLE} -c "import sys; sys.path.insert(1,\"${PROJECT_SOURCE_DIR}/python/\"); import ngs_check; sys.exit(ngs_check.check_if_ngsolve_newer_than(ngs_check.__ngsolve_required__))" RESULT_VARIABLE ngsfits)
    if (NOT ${ngsfits} EQUAL 1)
      message(FATAL_ERROR "
NGSolve version is smaller than required.
Update NGSolve or downgrade ngsxfem!
If you know what you are doing you can turn off this check by setting the cmake option CHECK_NGSOLVE_VERSION off: \"cmake -DCHECK_NGSOLVE_VERSION=OFF ...\"
")
    endif (NOT ${ngsfits} EQUAL 1)
  endif (CHECK_NGSOLVE_VERSION)
  set(CMAKE_CXX_COMPILER "ngscxx" CACHE STRING "Default CXX compiler" FORCE)
  set(CMAKE_LINKER "ngsld" CACHE STRING "Default linker" FORCE)

endif(BUILD_NGSOLVE)

if(USE_CCACHE)
  find_program(CCACHE_FOUND ccache)
  if(CCACHE_FOUND)
    set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE ccache)
    set_property(GLOBAL PROPERTY RULE_LAUNCH_LINK ccache)
  endif(CCACHE_FOUND)
endif(USE_CCACHE)

if(SPACETIME_SANITY_CHECKS)
  add_definitions(-DSPACETIME_SANITY_CHECKS)
endif(SPACETIME_SANITY_CHECKS)


find_package(NGSolve REQUIRED CONFIG)
set(NGSXFEM_VERSION_MAJOR 2)
set(NGSXFEM_VERSION_MINOR 0)
set(NGSXFEM_VERSION_CYCLE 2104)
set(NGSXFEM_VERSION "${NGSXFEM_VERSION_MAJOR}.${NGSXFEM_VERSION_MINOR}.${NGSXFEM_VERSION_CYCLE}")
set(PACKAGE_VERSION ${NGSXFEM_VERSION})


add_definitions(-DNGSXFEM_VERSION="${NGSXFEM_VERSION}")

# use rpath
set(CMAKE_MACOSX_RPATH TRUE)
set(CMAKE_INSTALL_RPATH "${NETGEN_RPATH_TOKEN};${NETGEN_RPATH_TOKEN}/${NETGEN_RPATH}")

#######################################################################
# CTest
enable_testing()
include(CTest)

set(ngs_install_dir RUNTIME DESTINATION ${NGSOLVE_INSTALL_DIR_BIN} LIBRARY DESTINATION ${NGSOLVE_INSTALL_DIR_LIB} ARCHIVE DESTINATION ${NGSOLVE_INSTALL_DIR_LIB})


include_directories ("${PROJECT_SOURCE_DIR}/include")
include_directories ("${PROJECT_BINARY_DIR}")

set(CMAKE_INCLUDE_CURRENT_DIR ON)

add_subdirectory(python)
add_subdirectory(utils)
add_subdirectory(cutint)
add_subdirectory(xfem)
add_subdirectory(demos)
add_subdirectory(lsetcurving)
add_subdirectory(tests)
add_subdirectory(spacetime)

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
    
  Linker:
    CMAKE_LINKER = ${CMAKE_LINKER}
    
  Compiler Flags:
    flags = ${flags}

    * NGSolve compile flags:
      NGSOLVE_COMPILE_FLAGS = ${NGSOLVE_COMPILE_FLAGS}

  Generator:
    ${CMAKE_GENERATOR}

  Build Type:
    type = ${CMAKE_BUILD_TYPE}

  Build with space-time quadrature sanity checks (abuses ngsolve's intrule nr): ${SPACETIME_SANITY_CHECKS}

  Build NGSolve from scratch: ${BUILD_NGSOLVE}
  Build xfem (and NGSolve) with ccache: ${USE_CCACHE}
  Build xfem with python stub files: ${BUILD_STUB_FILES}

  Building:

    Type '${BUILDEXE}' to compile ${PROJECT_NAME}.

    Type '${BUILDEXE} install' to install ${PROJECT_NAME}.

------------------------------------------------------------------------
")


