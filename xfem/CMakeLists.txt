add_library(ngsxfem_xfem ${NGS_LIB_TYPE}
  ghostpenalty.cpp
  xFESpace.cpp
  xfemIntegrators.cpp
  xfemNitsche.cpp
  xfemVisInts.cpp
  xfiniteelement.cpp
  symboliccutbfi.cpp
  symboliccutlfi.cpp
  )

# stxfemIntegrators.cpp

set_target_properties(ngsxfem_xfem PROPERTIES SUFFIX ".so")

if(NOT WIN32)
  target_link_libraries(ngsxfem_xfem ngsxfem_cutint ${PYTHON_LIBS} ${LAPACK_LIBRARIES} ${NETGEN_LIBS})
  install( TARGETS ngsxfem_xfem ${ngs_install_dir})
endif(NOT WIN32)

install( FILES
  ghostpenalty.hpp
  xFESpace.hpp
  xfemIntegrators.hpp
  xfemNitsche.hpp
  xfemVisInts.hpp
  xfiniteelement.hpp
  symboliccutbfi.hpp
  symboliccutlfi.hpp
  DESTINATION include
  )

# stxfemIntegrators.hpp

if(USE_PYTHON)
    install (FILES
            __init__.py 
            DESTINATION ${PYTHON_PACKAGES_INSTALL_DIR}/xfem
            COMPONENT xfem
            )
endif(USE_PYTHON)



