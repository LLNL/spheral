set(PYTHON_PREFIX ${CMAKE_CURRENT_BINARY_DIR}/${lib_name})
set(PYTHON_CACHE "${CACHE_DIR}/Python-2.7.15.tgz")
set(PYTHON_SRC_DIR ${PYTHON_PREFIX}/src/python)
set(PYTHON_INSTALL_DIR ${${lib_name}_DIR})
set(PYTHON_EXE ${PYTHON_INSTALL_DIR}/bin/python)
set(PYTHON_URL "http://www.python.org/ftp/python/2.7.15/Python-2.7.15.tgz")
set(PYTHON_SITE_PACKAGE_DIR ${PYTHON_INSTALL_DIR}/lib/python2.7/site-packages)

set(${lib_name}_INCLUDES ${PYTHON_INSTALL_DIR}/include/python2.7/)

if(${lib_name}_BUILD)

  if (EXISTS ${PYTHON_CACHE})
    set(PYTHON_URL ${PYTHON_CACHE})
  endif()

  ExternalProject_add(${lib_name}
    PREFIX ${PYTHON_PREFIX}
    URL ${PYTHON_URL} 
    DOWNLOAD_DIR ${CACHE_DIR}
    CONFIGURE_COMMAND env CC=${CMAKE_C_COMPILER} CXX=${CMAKE_CXX_COMPILER} ${PYTHON_SRC_DIR}/configure
                      --with-cxx-main='${CMAKE_CXX_COMPILER}'
                      --disable-ipv6
                      --prefix=${PYTHON_INSTALL_DIR}
    BUILD_COMMAND make -j${TPL_PARALLEL_BUILD}
    INSTALL_COMMAND make install

    LOG_DOWNLOAD ${OUT_PROTOCOL_EP}
    LOG_CONFIGURE ${OUT_PROTOCOL_EP}
    LOG_BUILD ${OUT_PROTOCOL_EP}
    LOG_INSTALL ${OUT_PROTOCOL_EP}
  )
  set(${lib_name}-install-dep ${lib_name})

else()
  set(${lib_name}_ADD_BLT_TARGET OFF)
endif()

add_custom_target(
  ${lib_name}-install
  COMMAND ${PYTHON_EXE} -V &> python-version.log
  DEPENDS ${${lib_name}-install-dep}
)

list(APPEND spheral_py_depends 
  python-install)


set(PYTHON_EXE ${PYTHON_EXE} PARENT_SCOPE)
set(PYTHON_INSTALL_DIR ${PYTHON_INSTALL_DIR} PARENT_SCOPE)
set(PYTHON_SITE_PACKAGE_DIR ${PYTHON_SITE_PACKAGE_DIR} PARENT_SCOPE)
