add_library(ngsxfem_utils ${NGS_LIB_TYPE}
  bitarraycf.cpp
  ngsxstd.cpp
)

set_target_properties(ngsxfem_utils PROPERTIES SUFFIX ".so")

target_link_libraries(ngsxfem_utils ${PYTHON_LIBS} ${LAPACK_LIBRARIES} ${NETGEN_LIBS})
install( TARGETS ngsxfem_utils ${ngs_install_dir})


install( FILES
  ngsxstd.hpp
  bitarraycf.hpp
  DESTINATION include
)


if(USE_PYTHON)
    install (FILES
            utils.py 
            DESTINATION ${PYTHON_PACKAGES_INSTALL_DIR}/xfem
            COMPONENT xfem
            )
endif(USE_PYTHON)




