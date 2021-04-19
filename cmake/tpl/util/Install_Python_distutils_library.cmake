macro(Install_Python_distutils_library lib_name lib_file lib_url lib_dir)
  add_custom_command(
    OUTPUT ${lib_file}
    COMMAND echo "-- downloading ${lib_file}"
    COMMAND wget ${lib_url} --no-check-certificate -O ${lib_file} > ${lib_file}-download.log
    )
  add_custom_target(${lib_name}
    COMMAND echo "-- pip installing ${lib_list}"
    COMMAND ${PYTHON_EXE} ${PIP_EXE} ${OUT_PROTOCOL_PIP} install ${${lib_list}} --no-index --find-links ${CACHE_DIR}
    DEPENDS pip-install ${${lib_list}_DEPENDS} ${lib_file}
    )
  list(APPEND spheral_py_depends ${lib_name})
endmacro()
