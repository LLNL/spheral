set(PYTHON_INSTALL_DIR ${${lib_name}_DIR})
set(PYTHON_EXE ${PYTHON_INSTALL_DIR}/bin/python)
set(PYTHON_SITE_PACKAGE_DIR ${PYTHON_INSTALL_DIR}/lib/python2.7/site-packages)

set(${lib_name}_libs libpython2.7.so)
set(${lib_name}_INCLUDES $<BUILD_INTERFACE:${PYTHON_INSTALL_DIR}/include/python2.7/>)

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
