add_library(ngsxfem_py ${NGS_LIB_TYPE}
          python_xfem.cpp
        )

set_target_properties(ngsxfem_py PROPERTIES SUFFIX ".so")

if(NOT WIN32)
    target_link_libraries(ngsxfem_py ngsxfem_lsetcurving  ${PYTHON_LIBS} ${LAPACK_LIBRARIES} ${NETGEN_LIBS})
    install( TARGETS ngsxfem_py ${ngs_install_dir})
endif(NOT WIN32)

if(USE_PYTHON)
    install (FILES
            __init__.py 
            DESTINATION ${PYTHON_PACKAGES_INSTALL_DIR}/xfem
            COMPONENT xfem
            )
endif(USE_PYTHON)
