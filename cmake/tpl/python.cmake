set(PYTHON_INSTALL_DIR ${${lib_name}_DIR})
set(PYTHON_EXE ${PYTHON_INSTALL_DIR}/bin/python)
set(PYTHON_SITE_PACKAGE_DIR ${PYTHON_INSTALL_DIR}/lib/python3.9/site-packages)

list(APPEND ${lib_name}_libs ${Python3_LIBRARIES})
list(APPEND ${lib_name}_INCLUDES ${Python3_INCLUDE_DIRS})
message("** ${lib_name}_libs : ${${lib_name}_libs}")
message("** ${lib_name}_INCLUDES : ${${lib_name}_INCLUDES}")

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
