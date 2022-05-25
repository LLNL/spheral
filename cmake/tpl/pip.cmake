if (NOT DEFINED PIP_EXE)
  set(PIP_EXE ${PYTHON_INSTALL_DIR}/bin/pip2.7)
endif()

set(${lib_name}_ADD_BLT_TARGET OFF)

add_custom_target(
  ${lib_name}-install
  COMMAND ${CMAKE_COMMAND} -E env PYTHONPATH="${setuptools_DIR}:${pip_DIR}" ${PYTHON_EXE} ${PIP_EXE} -V &> pip-version.log
  DEPENDS ${${lib_name}-install-dep}
)

list(APPEND ${dep_list}
  pip-install)

set(PIP_EXE ${PIP_EXE} PARENT_SCOPE)
