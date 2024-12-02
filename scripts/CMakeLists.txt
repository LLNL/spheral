if (NOT ENABLE_CXXONLY)

<<<<<<< HEAD
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
=======
  # We use the spack generated build time python path to
  # pick up all of the python modules needed for our runtime.
  # Since some modules need to install from directories other than
  # site-packages, we strip away site-packages and pack and
  # copy the full contents of the python library prefix to
  # our virtual env in spheral-setup-venv.sh
  string(REGEX REPLACE "lib\/python3.9\/site-packages\/?[A-Za-z]*:" ";" VIRTUALENV_PYTHONPATH_COPY "${SPACK_PYTHONPATH}:")
>>>>>>> develop

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

  install(FILES 
      "${CMAKE_CURRENT_BINARY_DIR}/spheral-setup-venv.sh"
      "${CMAKE_CURRENT_BINARY_DIR}/spheral-env.sh"
      "${CMAKE_CURRENT_BINARY_DIR}/atstest.sh"
      "${CMAKE_CURRENT_SOURCE_DIR}/spheral_ats.py"
      "${CMAKE_CURRENT_SOURCE_DIR}/gitlab/performance_analysis.py"
      "${CMAKE_CURRENT_SOURCE_DIR}/spheralutils.py"
    DESTINATION "${CMAKE_INSTALL_PREFIX}/scripts"
  )

<<<<<<< HEAD
  install(FILES
      "${CMAKE_CURRENT_BINARY_DIR}/performance/performance.py"
    DESTINATION "${CMAKE_INSTALL_PREFIX}/tests"
  )

  configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/runtime-requirements.txt.in"
    "${CMAKE_CURRENT_BINARY_DIR}/runtime-requirements.txt"
  )

  # Install Runtime python dependencies into the Install directory.
  Spheral_Python_Env(python_runtime_env 
    REQUIREMENTS runtime-requirements.txt 
    PREFIX ${CMAKE_INSTALL_PREFIX}
  )

  install(FILES ${polytope_DIR}/${SPHERAL_SITE_PACKAGES_PATH}/polytope/polytope.so
    DESTINATION ${CMAKE_INSTALL_PREFIX}/.venv/${SPHERAL_SITE_PACKAGES_PATH}/polytope/
  )
=======
  install(CODE "execute_process( \
    COMMAND env PYTHONPATH=${SPACK_PYTHONPATH} ${PYTHON_EXE} -m venv .venv --without-pip --prompt \
    'Spheral>')"
  )

  # Copy over all of the python TPL files, with a few exceptions
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
>>>>>>> develop

  install(CODE "execute_process( \
    COMMAND bash ${CMAKE_CURRENT_BINARY_DIR}/spheral-setup-venv.sh)"
  )
endif()
