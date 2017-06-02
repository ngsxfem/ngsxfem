if(NETGEN_USE_PYTHON)
    add_ngsolve_python_module(ngsxfem_py python_xfem.cpp)
    set_target_properties(ngsxfem_py PROPERTIES INSTALL_RPATH "${NETGEN_RPATH_TOKEN}/../${NETGEN_PYTHON_RPATH}")
    install(TARGETS ngsxfem_py DESTINATION ${NGSOLVE_INSTALL_DIR_PYTHON}/xfem)

    if(NOT WIN32)
        target_link_libraries(ngsxfem_py ngsxfem_lsetcurving ngsxfem_utils)
        install( TARGETS ngsxfem_py DESTINATION ${NGSOLVE_INSTALL_DIR_PYTHON}/xfem )
    endif(NOT WIN32)

    install (FILES
        __init__.py 
        DESTINATION ${NGSOLVE_INSTALL_DIR_PYTHON}/xfem
        COMPONENT xfem
    )
endif(NETGEN_USE_PYTHON)
