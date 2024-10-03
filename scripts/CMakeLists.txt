if (NOT ENABLE_CXXONLY)

  # We use the spack generated build time python path to
  # pick up all of the python modules needed for our runtime.
  # Since some modules need to install from directories other than
  # site-packages, we strip away site-packages and pack and
  # copy the full contents of the python library prefix to
  # our virtual env in spheral-setup-venv.sh
  string(REGEX REPLACE "lib\/python3.9\/site-packages\/?[A-Za-z]*:" ";" VIRTUALENV_PYTHONPATH_COPY "${SPACK_PYTHONPATH}:")

  set(SPHERAL_ATS_BUILD_CONFIG_ARGS )

  if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    list(APPEND SPHERAL_ATS_BUILD_CONFIG_ARGS "--filter='\"level<100\"'")
  endif()

  if (NOT ENABLE_MPI)
    list(APPEND SPHERAL_ATS_BUILD_CONFIG_ARGS "--filter='\"np<2\"'")
  endif()

  if (NOT SPHERAL_ENABLE_FSISPH)
     list(APPEND SPHERAL_ATS_BUILD_CONFIG_ARGS "--filter='\"not fsisph\"'")
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

  configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/performance/performance.py.in"
    "${CMAKE_CURRENT_BINARY_DIR}/performance/performance.py"
  )

  install(FILES
      "${CMAKE_CURRENT_BINARY_DIR}/spheral-setup-venv.sh"
      "${CMAKE_CURRENT_BINARY_DIR}/spheral-env.sh"
      "${CMAKE_CURRENT_BINARY_DIR}/atstest.sh"
      "${CMAKE_CURRENT_BINARY_DIR}/lcatstest.sh"
      "${CMAKE_CURRENT_SOURCE_DIR}/lc/lcats"
    DESTINATION "${CMAKE_INSTALL_PREFIX}/scripts"
  )

  install(FILES
      "${CMAKE_CURRENT_BINARY_DIR}/performance/performance.py"
    DESTINATION "${CMAKE_INSTALL_PREFIX}/tests"
  )

  install(CODE "execute_process( \
    COMMAND env PYTHONPATH=${SPACK_PYTHONPATH} ${PYTHON_EXE} -m venv .venv --without-pip --prompt \
    'Spheral>')"
  )

  foreach(_venv_dir ${VIRTUALENV_PYTHONPATH_COPY})
    if(NOT ${_venv_dir} MATCHES "sphinx")
      install(DIRECTORY ${_venv_dir}
        USE_SOURCE_PERMISSIONS
        MESSAGE_NEVER
        DESTINATION "${CMAKE_INSTALL_PREFIX}/.venv"
        PATTERN "*\/.spack*" EXCLUDE
        PATTERN "*\/tests\/*" EXCLUDE
        PATTERN "*.pyc" EXCLUDE
      )
    endif()
  endforeach()

  install(CODE "execute_process( \
    COMMAND bash ${CMAKE_CURRENT_BINARY_DIR}/spheral-setup-venv.sh)"
  )
endif()
