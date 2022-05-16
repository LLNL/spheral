set(PIP_PREFIX ${CMAKE_CURRENT_BINARY_DIR}/${lib_name})
set(PIP_DIST pip-9.0.1-py2.py3-none-any.whl)
set(PIP_URL "https://pypi.python.org/packages/b6/ac/7015eb97dc749283ffdec1c3a88ddb8ae03b8fad0f0e611408f196358da3/pip-9.0.1-py2.py3-none-any.whl")
set(PIP_CACHE ${CACHE_DIR}/${PIP_DIST})
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
