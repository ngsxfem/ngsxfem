###############################################################################
# This file was taken from https://github.com/NGSolve/ngsolve-addon-template
# Make sure to check for updates regularly.
# Don't change anything here (unless you know what you are doing!)
###############################################################################
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Find NGSolve and Netgen using python
if(CMAKE_VERSION VERSION_LESS "3.18")
  find_package(Python3 REQUIRED COMPONENTS Interpreter Development)
else()
  find_package(Python3 REQUIRED COMPONENTS Interpreter Development.Module)
endif()


set(Netgen_DIR "" CACHE PATH "Path to directory containing NetgenConfig.cmake")
set(NGSolve_DIR "" CACHE PATH "Path to directory containing NGSolveConfig.cmake")

execute_process(COMMAND ${Python3_EXECUTABLE} -m netgen.config OUTPUT_STRIP_TRAILING_WHITESPACE OUTPUT_VARIABLE Netgen_DIR)
execute_process(COMMAND ${Python3_EXECUTABLE} -m ngsolve.config OUTPUT_STRIP_TRAILING_WHITESPACE OUTPUT_VARIABLE NGSolve_DIR)

find_package(NGSolve CONFIG REQUIRED)

macro(add_ngsolve_addon module_name)
  # Create the module
  add_library(${module_name} SHARED ${ARGN})
  target_link_libraries(${module_name} PUBLIC ngsolve Python3::Module)
  set_target_properties(${module_name} PROPERTIES PREFIX "" CXX_STANDARD 17)

  # Python does not recognize .dll (Windows) and .dylib (MacOS) file endings as modules
  if(WIN32)
    set_target_properties(${module_name} PROPERTIES SUFFIX ".pyd" )
  else(WIN32)
    set_target_properties(${module_name} PROPERTIES SUFFIX ".so")
  endif(WIN32)
endmacro()

execute_process(COMMAND ${Python3_EXECUTABLE} -c "import sys,sysconfig,os.path; print(os.path.relpath(sysconfig.get_path('platlib'), sys.prefix))"
  OUTPUT_STRIP_TRAILING_WHITESPACE OUTPUT_VARIABLE python3_library_dir
)

if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  if(NETGEN_BULID_FOR_CONDA)
    # Netgen is assumed to be installed as python package
    # Set install prefix to user-base if a user site is available, sys.prefix otherwise
    execute_process(COMMAND ${Python3_EXECUTABLE} -c "import sys; print(sys.prefix)"
      OUTPUT_STRIP_TRAILING_WHITESPACE OUTPUT_VARIABLE install_prefix
    )
    execute_process(COMMAND ${Python3_EXECUTABLE} -m site --user-base
      OUTPUT_STRIP_TRAILING_WHITESPACE OUTPUT_VARIABLE user_base RESULT_VARIABLE ret
    )
    if (ret EQUAL 0)
      set(install_prefix ${user_base})
    endif()
    set(CMAKE_INSTALL_PREFIX ${install_prefix}/${python3_library_dir} CACHE PATH "Install dir" FORCE)
  else()
    # Netgen is self-compiled -> install into Netgen directory
    set(CMAKE_INSTALL_PREFIX ${NETGEN_DIR} CACHE PATH "Install prefix" FORCE)
  endif()
  set(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT OFF)
endif(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)

if(SKBUILD_PLATLIB_DIR)
  set(ADDON_INSTALL_DIR_PYTHON ${SKBUILD_PLATLIB_DIR})
  set(ADDON_INSTALL_DIR_BIN ${SKBUILD_DATA_DIR}/bin)
  set(ADDON_INSTALL_DIR_INCLUDE ${SKBUILD_DATA_DIR}/include)
  set(ADDON_INSTALL_DIR_LIB ${SKBUILD_DATA_DIR}/lib)
  set(ADDON_INSTALL_DIR_RES ${SKBUILD_DATA_DIR}/share)
  set(ADDON_INSTALL_DIR_CMAKE ${SKBUILD_DATA_DIR}/lib/cmake)
  set(stubgen_working_dir ${SKBUILD_PLATLIB_DIR})
elseif(NETGEN_BUILD_FOR_CONDA)
  set(stubgen_working_dir ${CMAKE_INSTALL_PREFIX})
  set(ADDON_INSTALL_DIR_PYTHON ${python3_library_dir})
  set(ADDON_INSTALL_DIR_BIN bin)
  set(ADDON_INSTALL_DIR_INCLUDE include)
  set(ADDON_INSTALL_DIR_LIB lib)
  set(ADDON_INSTALL_DIR_RES share)
  set(ADDON_INSTALL_DIR_CMAKE lib/cmake)
else()
  set(stubgen_working_dir ${CMAKE_INSTALL_PREFIX}/${NETGEN_INSTALL_DIR_PYTHON})
  set(ADDON_INSTALL_DIR_PYTHON ${NETGEN_INSTALL_DIR_PYTHON})
  set(ADDON_INSTALL_DIR_BIN ${NETGEN_INSTALL_DIR_BIN})
  set(ADDON_INSTALL_DIR_INCLUDE ${NETGEN_INSTALL_DIR_INCLUDE})
  set(ADDON_INSTALL_DIR_LIB ${NETGEN_INSTALL_DIR_LIB})
  set(ADDON_INSTALL_DIR_RES ${NETGEN_INSTALL_DIR_RES})
  set(ADDON_INSTALL_DIR_CMAKE ${NETGEN_INSTALL_DIR_CMAKE})
endif()

macro(ngsolve_generate_stub_files module_name)
  set(stubgen_generation_code "execute_process(WORKING_DIRECTORY ${stubgen_working_dir} COMMAND ${Python3_EXECUTABLE} -m pybind11_stubgen --ignore-all-errors -o ${CMAKE_CURRENT_BINARY_DIR}/stubs ${module_name})")
  set(stubgen_directory "${CMAKE_CURRENT_BINARY_DIR}/stubs/${module_name}/")

  install(CODE ${stubgen_generation_code})
  install(DIRECTORY ${stubgen_directory} DESTINATION ${ADDON_INSTALL_DIR_PYTHON}/${module_name})
endmacro()

message(STATUS "Install dir: ${CMAKE_INSTALL_PREFIX}")
