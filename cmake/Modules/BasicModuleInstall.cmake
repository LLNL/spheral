add_custom_target(
  ${lib_name}
  COMMAND ${PIP_EXE} ${OUT_PROTOCOL_PIP} download --no-binary :all -d ${CACHE_DIR} ${lib_name}
  COMMAND ${PIP_EXE} ${OUT_PROTOCOL_PIP} install --upgrade ${lib_name} --no-index --find-links ${CACHE_DIR}
  DEPENDS pip-install setuptools wheel ${${lib_name}_EXTRA_DEPENDS}
)
