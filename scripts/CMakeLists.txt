if (NOT ENABLE_CXXONLY)

  # We Use the spack generated build time python path to
  # pick up all of the python modules needed for our runtime.
  # Some modules need to install from directories other that
  # site-packages so we strip away site-packages pack and 
  # copy the full contents of the python library prefix to
  # our virtual env in spheral-setup-venv.sh
  string(REGEX REPLACE "lib\/python3.9\/site-packages\/?[A-Za-z]*:" "* " VIRTUALENV_PYTHONPATH_COPY "${SPACK_PYTHONPATH}:")

  set(SPHERAL_ATS_BUILD_CONFIG_ARGS )

  if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    list(APPEND SPHERAL_ATS_BUILD_CONFIG_ARGS "--filter='\"level<100\"'")
  endif()

  if (NOT ENABLE_MPI)
    list(APPEND SPHERAL_ATS_BUILD_CONFIG_ARGS "--filter='\"np<2\"'")
  endif()

  if (NOT SPHERAL_ENABLE_GSPH)
     list(APPEND SPHERAL_ATS_BUILD_CONFIG_ARGS "--filter='\"not gsph\"'")
  endif()

  if (NOT SPHERAL_ENABLE_SVPH)
     list(APPEND SPHERAL_ATS_BUILD_CONFIG_ARGS "--filter='\"not svph\"'")
  endif()

  if ($ENV{SYS_TYPE} MATCHES ".*blueos.*")
    list(APPEND SPHERAL_ATS_BUILD_CONFIG_ARGS "--addOp --smpi_off")
  endif()

  string(REPLACE ";" " " SPHERAL_ATS_BUILD_CONFIG_ARGS_STRING "${SPHERAL_ATS_BUILD_CONFIG_ARGS}")
  
  configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/spheral-setup-venv.in"
    "${CMAKE_CURRENT_BINARY_DIR}/spheral-setup-venv.sh"
    )

  configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/spheral-env.in"
    "${CMAKE_CURRENT_BINARY_DIR}/spheral-env.sh"
    )

  configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/atstest.in"
    "${CMAKE_CURRENT_BINARY_DIR}/atstest.sh"
    )

  configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/lcatstest.in"
    "${CMAKE_CURRENT_BINARY_DIR}/lcatstest.sh"
    )

  install(FILES 
      "${CMAKE_CURRENT_BINARY_DIR}/spheral-setup-venv.sh"
      "${CMAKE_CURRENT_BINARY_DIR}/spheral-env.sh"
      "${CMAKE_CURRENT_BINARY_DIR}/atstest.sh"
      "${CMAKE_CURRENT_BINARY_DIR}/lcatstest.sh"
      "${CMAKE_CURRENT_SOURCE_DIR}/lc/lcats"
    DESTINATION "${CMAKE_INSTALL_PREFIX}/scripts"
    )

  install(CODE "execute_process( \
    COMMAND bash ${CMAKE_CURRENT_BINARY_DIR}/spheral-setup-venv.sh \ 
    )")
endif()
