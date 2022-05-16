set(PYTHON_PREFIX ${CMAKE_CURRENT_BINARY_DIR}/${lib_name})
set(PYTHON_VERSION "2.7.18")
set(PYTHON_CACHE "${CACHE_DIR}/Python-${PYTHON_VERSION}.tgz")
set(PYTHON_SRC_DIR ${PYTHON_PREFIX}/src/python)
set(PYTHON_INSTALL_DIR ${${lib_name}_DIR})
set(PYTHON_EXE ${PYTHON_INSTALL_DIR}/bin/python)
set(PYTHON_URL "http://www.python.org/ftp/python/${PYTHON_VERSION}/Python-${PYTHON_VERSION}.tgz")
set(PYTHON_SITE_PACKAGE_DIR ${PYTHON_INSTALL_DIR}/lib/python2.7/site-packages)

set(${lib_name}_libs libpython2.7.so)
set(${lib_name}_INCLUDES $<BUILD_INTERFACE:${PYTHON_INSTALL_DIR}/include/python2.7/>)

set(PYTHON_C_COMPILER ${CMAKE_C_COMPILER})
set(PYTHON_CXX_COMPILER ${CMAKE_CXX_COMPILER})

if(PYTHON_C_COMPILER MATCHES "mpicc$")
  execute_process(COMMAND ${CMAKE_C_COMPILER} -show
                  OUTPUT_VARIABLE MPICC_SHOW )
  string (REPLACE " " ";" MPICC_SHOW_LIST "${MPICC_SHOW}")
  list(GET MPICC_SHOW_LIST 0 PYTHON_C_COMPILER)
  message("${PYTHON_C_COMPILER}")
endif()
if(PYTHON_CXX_COMPILER MATCHES "mpicxx$")
  execute_process(COMMAND ${CMAKE_CXX_COMPILER} -show
                  OUTPUT_VARIABLE MPICXX_SHOW )
  string (REPLACE " " ";" MPICXX_SHOW_LIST "${MPICXX_SHOW}")
  list(GET MPICXX_SHOW_LIST 0 PYTHON_CXX_COMPILER)
  message("${PYTHON_CXX_COMPILER}")
endif()

add_custom_target(
  ${lib_name}-install
  COMMAND ${PYTHON_EXE} -V &> python-version.log
  DEPENDS ${${lib_name}-install-dep}
)

list(APPEND spheral_py_depends 
  python-install)


set(${lib_name}_BUILD ${${lib_name}_BUILD} PARENT_SCOPE)
set(PYTHON_EXE ${PYTHON_EXE} PARENT_SCOPE)
set(PYTHON_INSTALL_DIR ${PYTHON_INSTALL_DIR} PARENT_SCOPE)
set(PYTHON_SITE_PACKAGE_DIR ${PYTHON_SITE_PACKAGE_DIR} PARENT_SCOPE)
